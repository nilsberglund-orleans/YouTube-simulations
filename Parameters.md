### Parameter values for YouTube simulations ###

Created by **Nils Berglund** and optimized by **Marco Mancini**

C code for videos on YouTube Channel https://www.youtube.com/c/NilsBerglund

Below are parameter values used for different simulations, as well as initial conditions used in 
function animation. Some simulations use variants of the published code. The list is going to be 
updated gradually. 


### 12 April 22 - Interacting dipoles on a surface of genus 2 ###

**Program:** `lennardjones.c` 

```
#define MOVIE 1         /* set to 1 to generate movie */
#define DOUBLE_MOVIE 1  /* set to 1 to produce movies for wave height and energy simultaneously */

#define TIME_LAPSE 0     /* set to 1 to add a time-lapse movie at the end */
                         /* so far incompatible with double movie */
#define TIME_LAPSE_FACTOR 4    /* factor of time-lapse movie */


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

#define ADD_FIXED_OBSTACLES 1   /* set to 1 do add fixed circular obstacles */
#define OBSTACLE_PATTERN 2  /* pattern of obstacles, see list in global_ljones.c */

#define TWO_TYPES 0         /* set to 1 to have two types of particles */
#define TPYE_PROPORTION 0.8 /* proportion of particles of first type */
#define SYMMETRIZE_FORCE 0  /* set to 1 to symmetrize two-particle interaction, only needed if particles are not all the same */
#define CENTER_PX 0         /* set to 1 to center horizontal momentum */
#define CENTER_PY 0         /* set to 1 to center vertical momentum */
#define CENTER_PANGLE 0     /* set to 1 to center angular momentum */

#define INTERACTION 5       /* particle interaction, see list in global_ljones.c */
#define INTERACTION_B 1     /* particle interaction for second type of particle, see list in global_ljones.c */
#define SPIN_INTER_FREQUENCY 2.0 /* angular frequency of spin-spin interaction */
#define SPIN_INTER_FREQUENCY_B 2.0 /* angular frequency of spin-spin interaction for second particle type */

#define P_PERCOL 0.25       /* probability of having a circle in C_RAND_PERCOL arrangement */
#define NPOISSON 100        /* number of points for Poisson C_RAND_POISSON arrangement */
#define PDISC_DISTANCE 1.75  /* minimal distance in Poisson disc process, controls density of particles */
#define PDISC_CANDIDATES 100 /* number of candidates in construction of Poisson disc process */
#define RANDOM_POLY_ANGLE 0 /* set to 1 to randomize angle of polygons */

#define LAMBDA 2.0	    /* parameter controlling the dimensions of domain */
#define MU 0.045 	    /* parameter controlling radius of particles */
#define MU_B 0.0254         /* parameter controlling radius of particles of second type */
#define NPOLY 3             /* number of sides of polygon */
#define APOLY 0.125         /* angle by which to turn polygon, in units of Pi/2 */ 
#define MDEPTH 4            /* depth of computation of Menger gasket */
#define MRATIO 3            /* ratio defining Menger gasket */
#define MANDELLEVEL 1000    /* iteration level for Mandelbrot set */
#define MANDELLIMIT 10.0    /* limit value for approximation of Mandelbrot set */
#define FOCI 1              /* set to 1 to draw focal points of ellipse */
#define NGRIDX 30           /* number of grid point for grid of disks */
#define NGRIDY 20           /* number of grid point for grid of disks */
#define EHRENFEST_RADIUS 0.9    /* radius of container for Ehrenfest urn configuration */
#define EHRENFEST_WIDTH 0.035     /* width of tube for Ehrenfest urn configuration */

#define X_SHOOTER -0.2
#define Y_SHOOTER -0.6
#define X_TARGET 0.4
#define Y_TARGET 0.7        /* shooter and target positions in laser fight */

/* Parameters for length and speed of simulation */

#define NSTEPS 3400      /* number of frames of movie */
#define NVID 200          /* number of iterations between images displayed on screen */
#define NSEG 250         /* number of segments of boundary */
#define INITIAL_TIME 10    /* time after which to start saving frames */
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

#define BOUNDARY_COND 14

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

#define ENERGY_HUE_MIN 330.0      /* color of original particle */
#define ENERGY_HUE_MAX 50.0        /* color of saturated particle */
#define PARTICLE_HUE_MIN 359.0      /* color of original particle */
#define PARTICLE_HUE_MAX 0.0        /* color of saturated particle */
#define PARTICLE_EMAX 1.0e3           /* energy of particle with hottest color */
#define HUE_TYPE0 280.0     /* hue of particles of type 0 */
#define HUE_TYPE1 135.0      /* hue of particles of type 1 */

#define RANDOM_RADIUS 0     /* set to 1 for random circle radius */
#define DT_PARTICLE 3.0e-7    /* time step for particle displacement */
#define KREPEL 12.0          /* constant in repelling force between particles */
#define EQUILIBRIUM_DIST 5.0    /* Lennard-Jones equilibrium distance */
#define EQUILIBRIUM_DIST_B 5.0  /* Lennard-Jones equilibrium distance for second type of particle */
#define REPEL_RADIUS 20.0    /* radius in which repelling force acts (in units of particle radius) */
#define DAMPING 0.0     /* damping coefficient of particles */
#define PARTICLE_MASS 1.0   /* mass of particle of radius MU */
#define PARTICLE_MASS_B 1.0   /* mass of particle of radius MU */
#define PARTICLE_INERTIA_MOMENT 0.2     /* moment of inertia of particle */
#define PARTICLE_INERTIA_MOMENT_B 0.02     /* moment of inertia of second type of particle */
#define V_INITIAL 10.0        /* initial velocity range */
#define OMEGA_INITIAL 10.0        /* initial angular velocity range */

#define THERMOSTAT 1        /* set to 1 to switch on thermostat */
#define SIGMA 5.0           /* noise intensity in thermostat */
#define BETA 0.0002          /* initial inverse temperature */
#define MU_XI 0.01           /* friction constant in thermostat */
#define KSPRING_BOUNDARY 5.0e9    /* confining harmonic potential outside simulation region */
#define KSPRING_OBSTACLE 5.0e8    /* harmonic potential of obstacles */
#define NBH_DIST_FACTOR 3.0       /* radius in which to count neighbours */
#define GRAVITY 0.0            /* gravity acting on all particles */
#define INCREASE_GRAVITY 0     /* set to 1 to increase gravity during the simulation */
#define GRAVITY_FACTOR 100.0   /* factor by which to increase gravity */
#define GRAVITY_RESTORE_TIME 750    /* time at end of simulation with gravity restored to initial value */

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
#define BETA_FACTOR 50.0   /* factor by which to change BETA during simulation */
#define N_TOSCILLATIONS 2.25   /* number of temperature oscillations in BETA schedule */
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
#define N_T_AVERAGE 50       /* size of temperature averaging window */
#define MAX_PRESSURE 3.0e10  /* pressure shown in "hottest" color */
#define PARTIAL_THERMO_COUPLING 0   /* set to 1 to couple only particles to the right of obstacle to thermostat */
#define PARTIAL_THERMO_SHIFT 0.5    /* distance from obstacle at the right of which particles are coupled to thermostat */

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
#define POSITION_Y_DEPENDENCE 0     /* set to 1 for the separation between particles to be vertical */
#define PRINT_ENTROPY 0     /* set to 1 to compute entropy */

#define PRINT_PARTICLE_NUMBER 0     /* set to 1 to print total number of particles */

#define EHRENFEST_COPY 0    /* set to 1 to add equal number of larger particles (for Ehrenfest model) */

#define LID_MASS 1000.0     /* mass of lid for BC_RECTANGLE_LID b.c. */
#define LID_WIDTH 0.1       /* width of lid for BC_RECTANGLE_LID b.c. */
#define WALL_MASS 2000.0    /* mass of wall for BC_RECTANGLE_WALL b.c. */
#define WALL_FRICTION 0.0   /* friction on wall for BC_RECTANGLE_WALL b.c. */
#define WALL_WIDTH 0.1      /* width of wall for BC_RECTANGLE_WALL b.c. */
#define WALL_VMAX 100.0     /* max speed of wall */
#define WALL_TIME 500       /* time during which to keep wall */

#define FLOOR_FORCE 1      /* set to 1 to limit force on particle to FMAX */
#define FMAX 1.0e9         /* maximal force */
#define FLOOR_OMEGA 1      /* set to 1 to limit particle momentum to PMAX */
#define PMAX 1000.0        /* maximal force */

#define HASHX 36   /* size of hashgrid in x direction */
#define HASHY 20   /* size of hashgrid in y direction */
#define HASHMAX 100  /* maximal number of particles per hashgrid cell */
#define HASHGRID_PADDING 0.1    /* padding of hashgrid outside simulation window */

#define DRAW_COLOR_SCHEME 0     /* set to 1 to plot the color scheme */
#define COLORBAR_RANGE 8.0    /* scale of color scheme bar */
#define COLORBAR_RANGE_B 12.0    /* scale of color scheme bar for 2nd part */
#define ROTATE_COLOR_SCHEME 0   /* set to 1 to draw color scheme horizontally */

```

### 11 April 22 - 3D waves in a 2D Julia set ###

**Program:** `wave_3d.c` 

**Initial condition in function `animation()`:** 

```
    init_circular_wave(0.2, 0.4, phi, psi, xy_in);
    add_circular_wave(-1.0, -0.2, -0.4, phi, psi, xy_in);
```

```
#define MOVIE 1         /* set to 1 to generate movie */
#define DOUBLE_MOVIE 1  /* set to 1 to produce movies for wave height and energy simultaneously */

/* General geometrical parameters */

#define WINWIDTH 	1920  /* window width */
#define WINHEIGHT 	1000  /* window height */
#define NX 1920          /* number of grid points on x axis */
#define NY 1000          /* number of grid points on y axis */

#define XMIN -2.0
#define XMAX 2.0	/* x interval  */
#define YMIN -1.041666667
#define YMAX 1.041666667	/* y interval for 9/16 aspect ratio */

#define HIGHRES 0        /* set to 1 if resolution of grid is double that of displayed image */

#define JULIA_SCALE 0.8 /* scaling for Julia sets */

/* Choice of the billiard table */

#define B_DOMAIN 16        /* choice of domain shape, see list in global_pdes.c */

#define CIRCLE_PATTERN 201   /* pattern of circles or polygons, see list in global_pdes.c */

#define P_PERCOL 0.25       /* probability of having a circle in C_RAND_PERCOL arrangement */
#define NPOISSON 300        /* number of points for Poisson C_RAND_POISSON arrangement */
#define RANDOM_POLY_ANGLE 1 /* set to 1 to randomize angle of polygons */

#define LAMBDA 0.6	    /* parameter controlling the dimensions of domain */
#define MU 0.6              /* parameter controlling the dimensions of domain */
#define NPOLY 6             /* number of sides of polygon */
#define APOLY 0.0           /* angle by which to turn polygon, in units of Pi/2 */ 
#define MDEPTH 3            /* depth of computation of Menger gasket */
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

#define ISO_XSHIFT_LEFT -1.0 
#define ISO_XSHIFT_RIGHT 0.9
#define ISO_YSHIFT_LEFT 0.05
#define ISO_YSHIFT_RIGHT 0.05 
#define ISO_SCALE 0.85           /* coordinates for isospectral billiards */

/* You can add more billiard tables by adapting the functions */
/* xy_in_billiard and draw_billiard below */

/* Physical parameters of wave equation */

#define TWOSPEEDS 0          /* set to 1 to replace hardcore boundary by medium with different speed */
#define OSCILLATE_LEFT 0     /* set to 1 to add oscilating boundary condition on the left */
#define OSCILLATE_TOPBOT 0   /* set to 1 to enforce a planar wave on top and bottom boundary */

#define OMEGA 0.005        /* frequency of periodic excitation */
#define AMPLITUDE 0.8      /* amplitude of periodic excitation */ 
#define COURANT 0.06       /* Courant number */
#define COURANTB 0.03      /* Courant number in medium B */
#define GAMMA 0.0          /* damping factor in wave equation */
#define GAMMAB 1.0e-7        /* damping factor in wave equation */
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
#define OSCILLATING_SOURCE_PERIOD 20    /* period of oscillating source */

/* Boundary conditions, see list in global_pdes.c  */

#define B_COND 2

/* Parameters for length and speed of simulation */

#define NSTEPS 2500        /* number of frames of movie */
#define NVID 10           /* number of iterations between images displayed on screen */
#define NSEG 1000         /* number of segments of boundary */
#define INITIAL_TIME 0      /* time after which to start saving frames */
#define BOUNDARY_WIDTH 3    /* width of billiard boundary */

#define PAUSE 200       /* number of frames after which to pause */
#define PSLEEP 2         /* sleep time during pause */
#define SLEEP1  1        /* initial sleeping time */
#define SLEEP2  1        /* final sleeping time */
#define MID_FRAMES 200    /* number of still frames between parts of two-part movie */
#define END_FRAMES 100   /* number of still frames at end of movie */
#define FADE 1           /* set to 1 to fade at end of movie */

/* Parameters of initial condition */

#define INITIAL_AMP 0.5         /* amplitude of initial condition */
#define INITIAL_VARIANCE 0.0005  /* variance of initial condition */
#define INITIAL_WAVELENGTH  0.1  /* wavelength of initial condition */

/* Plot type, see list in global_pdes.c  */

#define PLOT 103

#define PLOT_B 104        /* plot type for second movie */

#define P_3D_AMPLITUDE  101     /* color depends on amplitude */
#define P_3D_ANGLE 102          /* color depends on angle with fixed direction */
#define P_3D_AMP_ANGLE 103      /* color depends on amplitude, luminosity depends on angle */
#define P_3D_ENERGY 104         /* color depends on energy, luminosity depends on angle */

#define AMPLITUDE_HIGH_RES 1    /* set to 1 to increase resolution of P_3D_AMPLITUDE plot */
#define FADE_IN_OBSTACLE 1      /* set to 1 to fade color inside obstacles */

#define PLOT_SCALE_ENERGY 0.05      /* vertical scaling in energy plot */

/* 3D representation */

#define REPRESENTATION_3D 1     /* choice of 3D representation */ 

#define REP_AXO_3D 0        /* linear projection (axonometry) */
#define REP_PROJ_3D 1       /* projection on plane orthogonal to observer line of sight */

/* Color schemes */

#define COLOR_PALETTE 14    /* Color palette, see list in global_pdes.c  */
#define COLOR_PALETTE_B 11     /* Color palette, see list in global_pdes.c  */

#define BLACK 1          /* background */

#define COLOR_SCHEME 3   /* choice of color scheme, see list in global_pdes.c  */

#define SCALE 0          /* set to 1 to adjust color scheme to variance of field */
#define SLOPE 1.0       /* sensitivity of color on wave amplitude */
#define VSCALE_AMPLITUDE 0.2     /* additional scaling factor for color scheme P_3D_AMPLITUDE */
#define VSCALE_ENERGY 0.35       /* additional scaling factor for color scheme P_3D_ENERGY */
#define ATTENUATION 0.0  /* exponential attenuation coefficient of contrast with time */
#define E_SCALE 200.0     /* scaling factor for energy representation */
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
#define COLORBAR_RANGE 3.0     /* scale of color scheme bar */
#define COLORBAR_RANGE_B 5.0    /* scale of color scheme bar for 2nd part */
#define ROTATE_COLOR_SCHEME 0   /* set to 1 to draw color scheme horizontally */

#define SAVE_TIME_SERIES 0      /* set to 1 to save wave time series at a point */

/* For debugging purposes only */
#define FLOOR 0         /* set to 1 to limit wave amplitude to VMAX */
#define VMAX 10.0       /* max value of wave amplitude */

/* Parameters controlling 3D projection */

double u_3d[2] = {0.75, -0.45};     /* projections of basis vectors for REP_AXO_3D representation */
double v_3d[2] = {-0.75, -0.45};
double w_3d[2] = {0.0, 0.015};
double light[3] = {0.816496581, -0.40824829, 0.40824829};      /* vector of "light" direction for P_3D_ANGLE color scheme */
double observer[3] = {10.0, 6.0, 8.5};    /* location of observer for REP_PROJ_3D representation */ 

#define Z_SCALING_FACTOR 0.018     /* overall scaling factor of z axis for REP_PROJ_3D representation */
#define XY_SCALING_FACTOR 3.75     /* overall scaling factor for on-screen (x,y) coordinates after projection */
#define ZMAX_FACTOR 1.0           /* max value of z coordinate for REP_PROJ_3D representation */
#define XSHIFT_3D 0.0             /* overall x shift for REP_PROJ_3D representation */
#define YSHIFT_3D 0.0             /* overall y shift for REP_PROJ_3D representation */

```

```
double julia_x = 0.37468, julia_y = 0.21115; 

```

### 10 April 22 - Brownian particles on a surface of genus 2 ###

**Program:** `lennardjones.c` 

```
#define MOVIE 1         /* set to 1 to generate movie */
#define DOUBLE_MOVIE 0  /* set to 1 to produce movies for wave height and energy simultaneously */

#define TIME_LAPSE 1     /* set to 1 to add a time-lapse movie at the end */
                         /* so far incompatible with double movie */
#define TIME_LAPSE_FACTOR 3    /* factor of time-lapse movie */

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

#define ADD_FIXED_OBSTACLES 1   /* set to 1 do add fixed circular obstacles */
#define OBSTACLE_PATTERN 2  /* pattern of obstacles, see list in global_ljones.c */

#define TWO_TYPES 0         /* set to 1 to have two types of particles */
#define TPYE_PROPORTION 0.8 /* proportion of particles of first type */
#define SYMMETRIZE_FORCE 1  /* set to 1 to symmetrize two-particle interaction, only needed if particles are not all the same */
#define CENTER_PX 0         /* set to 1 to center horizontal momentum */
#define CENTER_PY 0         /* set to 1 to center vertical momentum */
#define CENTER_PANGLE 0     /* set to 1 to center angular momentum */

#define INTERACTION 1       /* particle interaction, see list in global_ljones.c */
#define INTERACTION_B 1     /* particle interaction for second type of particle, see list in global_ljones.c */
#define SPIN_INTER_FREQUENCY 5.0 /* angular frequency of spin-spin interaction */
#define SPIN_INTER_FREQUENCY_B 2.0 /* angular frequency of spin-spin interaction for second particle type */

#define P_PERCOL 0.25       /* probability of having a circle in C_RAND_PERCOL arrangement */
#define NPOISSON 100        /* number of points for Poisson C_RAND_POISSON arrangement */
#define PDISC_DISTANCE 5.75  /* minimal distance in Poisson disc process, controls density of particles */
#define PDISC_CANDIDATES 100 /* number of candidates in construction of Poisson disc process */
#define RANDOM_POLY_ANGLE 0 /* set to 1 to randomize angle of polygons */

#define LAMBDA 2.0	    /* parameter controlling the dimensions of domain */
#define MU 0.015  	    /* parameter controlling radius of particles */
#define MU_B 0.0254         /* parameter controlling radius of particles of second type */
#define NPOLY 3             /* number of sides of polygon */
#define APOLY 0.125         /* angle by which to turn polygon, in units of Pi/2 */ 
#define MDEPTH 4            /* depth of computation of Menger gasket */
#define MRATIO 3            /* ratio defining Menger gasket */
#define MANDELLEVEL 1000    /* iteration level for Mandelbrot set */
#define MANDELLIMIT 10.0    /* limit value for approximation of Mandelbrot set */
#define FOCI 1              /* set to 1 to draw focal points of ellipse */
#define NGRIDX 30           /* number of grid point for grid of disks */
#define NGRIDY 20           /* number of grid point for grid of disks */
#define EHRENFEST_RADIUS 0.9    /* radius of container for Ehrenfest urn configuration */
#define EHRENFEST_WIDTH 0.035     /* width of tube for Ehrenfest urn configuration */

#define X_SHOOTER -0.2
#define Y_SHOOTER -0.6
#define X_TARGET 0.4
#define Y_TARGET 0.7        /* shooter and target positions in laser fight */

/* Parameters for length and speed of simulation */

#define NSTEPS 5100      /* number of frames of movie */
#define NVID 200          /* number of iterations between images displayed on screen */
#define NSEG 250         /* number of segments of boundary */
#define INITIAL_TIME 20     /* time after which to start saving frames */
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

#define BOUNDARY_COND 14

/* Plot type, see list in global_ljones.c  */

#define PLOT 5
#define PLOT_B 4        /* plot type for second movie */

#define COLOR_BONDS 1   /* set to 1 to color bonds according to length */

/* Color schemes */

#define COLOR_PALETTE 10     /* Color palette, see list in global_ljones.c  */

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

#define ENERGY_HUE_MIN 330.0      /* color of original particle */
#define ENERGY_HUE_MAX 50.0        /* color of saturated particle */
#define PARTICLE_HUE_MIN 359.0      /* color of original particle */
#define PARTICLE_HUE_MAX 0.0        /* color of saturated particle */
#define PARTICLE_EMAX 1.0e3           /* energy of particle with hottest color */
#define HUE_TYPE0 280.0     /* hue of particles of type 0 */
#define HUE_TYPE1 135.0      /* hue of particles of type 1 */
#define HUE_TYPE2 70.0      /* hue of particles of type 1 */
#define HUE_TYPE3 210.0      /* hue of particles of type 1 */

#define RANDOM_RADIUS 0     /* set to 1 for random circle radius */
#define DT_PARTICLE 5.0e-7    /* time step for particle displacement */
#define KREPEL 12.0          /* constant in repelling force between particles */
#define EQUILIBRIUM_DIST 5.0    /* Lennard-Jones equilibrium distance */
#define EQUILIBRIUM_DIST_B 5.0  /* Lennard-Jones equilibrium distance for second type of particle */
#define REPEL_RADIUS 20.0    /* radius in which repelling force acts (in units of particle radius) */
#define DAMPING 0.0     /* damping coefficient of particles */
#define PARTICLE_MASS 1.0   /* mass of particle of radius MU */
#define PARTICLE_MASS_B 1.0   /* mass of particle of radius MU */
#define PARTICLE_INERTIA_MOMENT 0.2     /* moment of inertia of particle */
#define PARTICLE_INERTIA_MOMENT_B 0.02     /* moment of inertia of second type of particle */
#define V_INITIAL 10.0        /* initial velocity range */
#define OMEGA_INITIAL 10.0        /* initial angular velocity range */

#define THERMOSTAT 1        /* set to 1 to switch on thermostat */
#define SIGMA 5.0           /* noise intensity in thermostat */
#define BETA 0.0001          /* initial inverse temperature */
#define MU_XI 0.01           /* friction constant in thermostat */
#define KSPRING_BOUNDARY 5.0e9    /* confining harmonic potential outside simulation region */
#define KSPRING_OBSTACLE 5.0e8    /* harmonic potential of obstacles */
#define NBH_DIST_FACTOR 4.0       /* radius in which to count neighbours */
#define GRAVITY 0.0            /* gravity acting on all particles */
#define INCREASE_GRAVITY 0     /* set to 1 to increase gravity during the simulation */
#define GRAVITY_FACTOR 100.0   /* factor by which to increase gravity */
#define GRAVITY_RESTORE_TIME 750    /* time at end of simulation with gravity restored to initial value */

#define ROTATION 0           /* set to 1 to include rotation of particles */
#define COUPLE_ANGLE_TO_THERMOSTAT 0    /* set to 1 to couple angular degrees of freedom to thermostat */
#define DIMENSION_FACTOR 1.0  /* scaling factor taking into account number of degrees of freedom */  
#define KTORQUE 50.0          /* force constant in angular dynamics */
#define KTORQUE_B 10.0        /* force constant in angular dynamics */
#define KTORQUE_DIFF 150.0    /* force constant in angular dynamics for different particles */
#define DRAW_SPIN 0           /* set to 1 to draw spin vectors of particles */
#define DRAW_SPIN_B 0         /* set to 1 to draw spin vectors of particles */
#define DRAW_CROSS 1          /* set to 1 to draw cross on particles of second type */
#define SPIN_RANGE 7.0       /* range of spin-spin interaction */
#define SPIN_RANGE_B 5.0     /* range of spin-spin interaction for second type of particle */
#define QUADRUPOLE_RATIO 0.6  /* anisotropy in quadrupole potential */ 

#define INCREASE_BETA 0  /* set to 1 to increase BETA during simulation */
#define BETA_FACTOR 20.0   /* factor by which to change BETA during simulation */
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
#define OBSTACLE_RADIUS 0.12  /* radius of obstacle for circle boundary conditions */
#define FUNNEL_WIDTH  0.25   /* funnel width for funnel boundary conditions */
#define OBSTACLE_XMIN 0.0    /* initial position of obstacle */
#define OBSTACLE_XMAX 3.0    /* final position of obstacle */
#define RECORD_PRESSURES 0   /* set to 1 to record pressures on obstacle */
#define N_PRESSURES 100      /* number of intervals to record pressure */
#define N_P_AVERAGE 100      /* size of pressure averaging window */
#define N_T_AVERAGE 50       /* size of temperature averaging window */
#define MAX_PRESSURE 3.0e10  /* pressure shown in "hottest" color */
#define PARTIAL_THERMO_COUPLING 0   /* set to 1 to couple only particles to the right of obstacle to thermostat */
#define PARTIAL_THERMO_SHIFT 0.5    /* distance from obstacle at the right of which particles are coupled to thermostat */

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
#define N_TRACER_PARTICLES 3    /* number of tracer particles */
#define TRAJECTORY_LENGTH 8000   /* length of recorded trajectory */
#define TRACER_PARTICLE_MASS 4.0    /* relative mass of tracer particle */
#define TRAJECTORY_WIDTH 3      /* width of tracer particle trajectory */

#define POSITION_DEPENDENT_TYPE 0   /* set to 1 to make particle type depend on initial position */
#define POSITION_Y_DEPENDENCE 0     /* set to 1 for the separation between particles to be vertical */
#define PRINT_ENTROPY 0     /* set to 1 to compute entropy */

#define PRINT_PARTICLE_NUMBER 0     /* set to 1 to print total number of particles */

#define EHRENFEST_COPY 0    /* set to 1 to add equal number of larger particles (for Ehrenfest model) */

#define LID_MASS 1000.0     /* mass of lid for BC_RECTANGLE_LID b.c. */
#define LID_WIDTH 0.1       /* width of lid for BC_RECTANGLE_LID b.c. */
#define WALL_MASS 2000.0    /* mass of wall for BC_RECTANGLE_WALL b.c. */
#define WALL_FRICTION 0.0   /* friction on wall for BC_RECTANGLE_WALL b.c. */
#define WALL_WIDTH 0.1      /* width of wall for BC_RECTANGLE_WALL b.c. */
#define WALL_VMAX 100.0     /* max speed of wall */
#define WALL_TIME 500       /* time during which to keep wall */

#define FLOOR_FORCE 1      /* set to 1 to limit force on particle to FMAX */
#define FMAX 1.0e9         /* maximal force */
#define FLOOR_OMEGA 1      /* set to 1 to limit particle momentum to PMAX */
#define PMAX 1000.0        /* maximal force */

#define HASHX 34   /* size of hashgrid in x direction */
#define HASHY 18   /* size of hashgrid in y direction */
#define HASHMAX 100  /* maximal number of particles per hashgrid cell */
#define HASHGRID_PADDING 0.1    /* padding of hashgrid outside simulation window */

#define DRAW_COLOR_SCHEME 0     /* set to 1 to plot the color scheme */
#define COLORBAR_RANGE 8.0    /* scale of color scheme bar */
#define COLORBAR_RANGE_B 12.0    /* scale of color scheme bar for 2nd part */
#define ROTATE_COLOR_SCHEME 0   /* set to 1 to draw color scheme horizontally */

```

### 9 April 22 - A hexagonal parabolic resonator, in 3D ###

**Program:** `wave_3d.c` 

**Initial condition in function `animation()`:** `init_circular_wave(0.0, 0.0, phi, psi, xy_in);`

```
#define MOVIE 1         /* set to 1 to generate movie */
#define DOUBLE_MOVIE 1  /* set to 1 to produce movies for wave height and energy simultaneously */

/* General geometrical parameters */

#define WINWIDTH 	1920  /* window width */
#define WINHEIGHT 	1000  /* window height */
#define NX 1920          /* number of grid points on x axis */
#define NY 1000          /* number of grid points on y axis */

#define XMIN -1.8
#define XMAX 2.2	/* x interval  */
#define YMIN -1.041666667
#define YMAX 1.041666667	/* y interval for 9/16 aspect ratio */

#define HIGHRES 0        /* set to 1 if resolution of grid is double that of displayed image */

#define JULIA_SCALE 1.0 /* scaling for Julia sets */

/* Choice of the billiard table */

#define B_DOMAIN 32        /* choice of domain shape, see list in global_pdes.c */

#define CIRCLE_PATTERN 201   /* pattern of circles or polygons, see list in global_pdes.c */

#define P_PERCOL 0.25       /* probability of having a circle in C_RAND_PERCOL arrangement */
#define NPOISSON 300        /* number of points for Poisson C_RAND_POISSON arrangement */
#define RANDOM_POLY_ANGLE 1 /* set to 1 to randomize angle of polygons */

#define LAMBDA 0.0	    /* parameter controlling the dimensions of domain */
#define MU 1.0              /* parameter controlling the dimensions of domain */
#define NPOLY 6             /* number of sides of polygon */
#define APOLY 0.0           /* angle by which to turn polygon, in units of Pi/2 */ 
#define MDEPTH 3            /* depth of computation of Menger gasket */
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

#define TWOSPEEDS 0          /* set to 1 to replace hardcore boundary by medium with different speed */
#define OSCILLATE_LEFT 0     /* set to 1 to add oscilating boundary condition on the left */
#define OSCILLATE_TOPBOT 0   /* set to 1 to enforce a planar wave on top and bottom boundary */

#define OMEGA 0.005        /* frequency of periodic excitation */
#define AMPLITUDE 0.8      /* amplitude of periodic excitation */ 
#define COURANT 0.06       /* Courant number */
#define COURANTB 0.03      /* Courant number in medium B */
#define GAMMA 0.0          /* damping factor in wave equation */
#define GAMMAB 1.0e-7        /* damping factor in wave equation */
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
#define OSCILLATING_SOURCE_PERIOD 20    /* period of oscillating source */

/* Boundary conditions, see list in global_pdes.c  */

#define B_COND 2

/* Parameters for length and speed of simulation */

#define NSTEPS 2350        /* number of frames of movie */
#define NVID 18          /* number of iterations between images displayed on screen */
#define NSEG 1000         /* number of segments of boundary */
#define INITIAL_TIME 0      /* time after which to start saving frames */
#define BOUNDARY_WIDTH 3    /* width of billiard boundary */

#define PAUSE 200       /* number of frames after which to pause */
#define PSLEEP 2         /* sleep time during pause */
#define SLEEP1  1        /* initial sleeping time */
#define SLEEP2  1        /* final sleeping time */
#define MID_FRAMES 20    /* number of still frames between parts of two-part movie */
#define END_FRAMES 100    /* number of still frames at end of movie */

/* Parameters of initial condition */

#define INITIAL_AMP 0.5         /* amplitude of initial condition */
#define INITIAL_VARIANCE 0.0025  /* variance of initial condition */
#define INITIAL_WAVELENGTH  0.5  /* wavelength of initial condition */

/* Plot type, see list in global_pdes.c  */

#define PLOT 103

#define PLOT_B 104        /* plot type for second movie */

#define P_3D_AMPLITUDE  101     /* color depends on amplitude */
#define P_3D_ANGLE 102          /* color depends on angle with fixed direction */
#define P_3D_AMP_ANGLE 103      /* color depends on amplitude, luminosity depends on angle */
#define P_3D_ENERGY 104         /* color depends on energy, luminosity depends on angle */

#define AMPLITUDE_HIGH_RES 1    /* set to 1 to increase resolution of P_3D_AMPLITUDE plot */
#define FADE_IN_OBSTACLE 1      /* set to 1 to fade color inside obstacles */

#define PLOT_SCALE_ENERGY 0.15      /* vertical scaling in energy plot */

/* 3D representation */

#define REPRESENTATION_3D 1     /* choice of 3D representation */ 

#define REP_AXO_3D 0        /* linear projection (axonometry) */
#define REP_PROJ_3D 1       /* projection on plane orthogonal to observer line of sight */

/* Color schemes */

#define COLOR_PALETTE 17     /* Color palette, see list in global_pdes.c  */
#define COLOR_PALETTE_B 14     /* Color palette, see list in global_pdes.c  */

#define BLACK 1          /* background */

#define COLOR_SCHEME 3   /* choice of color scheme, see list in global_pdes.c  */

#define SCALE 0          /* set to 1 to adjust color scheme to variance of field */
#define SLOPE 1.0       /* sensitivity of color on wave amplitude */
#define VSCALE_AMPLITUDE 0.1    /* additional scaling factor for color scheme P_3D_AMPLITUDE */
#define VSCALE_ENERGY 0.25       /* additional scaling factor for color scheme P_3D_ENERGY */
#define ATTENUATION 0.0  /* exponential attenuation coefficient of contrast with time */
#define E_SCALE 200.0     /* scaling factor for energy representation */
#define LOG_SCALE 1.0     /* scaling factor for energy log representation */
#define LOG_SHIFT 1.0     /* shift of colors on log scale */
#define RESCALE_COLOR_IN_CENTER 0   /* set to 1 to decrease color intentiy in the center (for wave escaping ring) */

#define COLORHUE 260     /* initial hue of water color for scheme C_LUM */
#define COLORDRIFT 0.0   /* how much the color hue drifts during the whole simulation */
#define LUMMEAN 0.5      /* amplitude of luminosity variation for scheme C_LUM */
#define LUMAMP 0.3       /* amplitude of luminosity variation for scheme C_LUM */
#define HUEMEAN 180.0    /* mean value of hue for color scheme C_HUE */
#define HUEAMP -180.0      /* amplitude of variation of hue for color scheme C_HUE */

#define DRAW_COLOR_SCHEME 1     /* set to 1 to plot the color scheme */
#define COLORBAR_RANGE 3.0     /* scale of color scheme bar */
#define COLORBAR_RANGE_B 5.0    /* scale of color scheme bar for 2nd part */
#define ROTATE_COLOR_SCHEME 0   /* set to 1 to draw color scheme horizontally */

#define SAVE_TIME_SERIES 0      /* set to 1 to save wave time series at a point */

/* For debugging purposes only */
#define FLOOR 0         /* set to 1 to limit wave amplitude to VMAX */
#define VMAX 10.0       /* max value of wave amplitude */

/* Parameters controlling 3D projection */

double u_3d[2] = {0.75, -0.45};     /* projections of basis vectors for REP_AXO_3D representation */
double v_3d[2] = {-0.75, -0.45};
double w_3d[2] = {0.0, 0.015};
double light[3] = {0.40824829, 0.816496581, 0.40824829};      /* vector of "light" direction for P_3D_ANGLE color scheme */
double observer[3] = {8.0, 8.0, 8.0};    /* location of observer for REP_PROJ_3D representation */ 

#define Z_SCALING_FACTOR 0.0075   /* overall scaling factor of z axis for REP_PROJ_3D representation */
#define XY_SCALING_FACTOR 3.0     /* overall scaling factor for on-screen (x,y) coordinates after projection */
#define ZMAX_FACTOR 1.0           /* max value of z coordinate for REP_PROJ_3D representation */

```

### 8 April 22 - Crystal formation on a surface of genus 2 ###

**Program:** `lennardjones.c` 

```
#define MOVIE 1         /* set to 1 to generate movie */
#define DOUBLE_MOVIE 1  /* set to 1 to produce movies for wave height and energy simultaneously */

#define TIME_LAPSE 0     /* set to 1 to add a time-lapse movie at the end */
                         /* so far incompatible with double movie */
#define TIME_LAPSE_FACTOR 4    /* factor of time-lapse movie */


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

#define ADD_FIXED_OBSTACLES 1   /* set to 1 do add fixed circular obstacles */
#define OBSTACLE_PATTERN 2  /* pattern of obstacles, see list in global_ljones.c */

#define TWO_TYPES 0         /* set to 1 to have two types of particles */
#define TPYE_PROPORTION 0.8 /* proportion of particles of first type */
#define SYMMETRIZE_FORCE 0  /* set to 1 to symmetrize two-particle interaction, only needed if particles are not all the same */
#define CENTER_PX 0         /* set to 1 to center horizontal momentum */
#define CENTER_PY 0         /* set to 1 to center vertical momentum */
#define CENTER_PANGLE 0     /* set to 1 to center angular momentum */

#define INTERACTION 1       /* particle interaction, see list in global_ljones.c */
#define INTERACTION_B 1     /* particle interaction for second type of particle, see list in global_ljones.c */
#define SPIN_INTER_FREQUENCY 5.0 /* angular frequency of spin-spin interaction */
#define SPIN_INTER_FREQUENCY_B 2.0 /* angular frequency of spin-spin interaction for second particle type */

#define P_PERCOL 0.25       /* probability of having a circle in C_RAND_PERCOL arrangement */
#define NPOISSON 100        /* number of points for Poisson C_RAND_POISSON arrangement */
#define PDISC_DISTANCE 3.25 /* minimal distance in Poisson disc process, controls density of particles */
#define PDISC_CANDIDATES 100 /* number of candidates in construction of Poisson disc process */
#define RANDOM_POLY_ANGLE 0 /* set to 1 to randomize angle of polygons */

#define LAMBDA 2.0	    /* parameter controlling the dimensions of domain */
#define MU 0.021 	    /* parameter controlling radius of particles */
#define MU_B 0.0254         /* parameter controlling radius of particles of second type */
#define NPOLY 3             /* number of sides of polygon */
#define APOLY 0.125         /* angle by which to turn polygon, in units of Pi/2 */ 
#define MDEPTH 4            /* depth of computation of Menger gasket */
#define MRATIO 3            /* ratio defining Menger gasket */
#define MANDELLEVEL 1000    /* iteration level for Mandelbrot set */
#define MANDELLIMIT 10.0    /* limit value for approximation of Mandelbrot set */
#define FOCI 1              /* set to 1 to draw focal points of ellipse */
#define NGRIDX 30           /* number of grid point for grid of disks */
#define NGRIDY 20           /* number of grid point for grid of disks */
#define EHRENFEST_RADIUS 0.9    /* radius of container for Ehrenfest urn configuration */
#define EHRENFEST_WIDTH 0.035     /* width of tube for Ehrenfest urn configuration */

#define X_SHOOTER -0.2
#define Y_SHOOTER -0.6
#define X_TARGET 0.4
#define Y_TARGET 0.7        /* shooter and target positions in laser fight */

/* Parameters for length and speed of simulation */

#define NSTEPS 3400      /* number of frames of movie */
#define NVID 200          /* number of iterations between images displayed on screen */
#define NSEG 250         /* number of segments of boundary */
#define INITIAL_TIME 10    /* time after which to start saving frames */
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

#define BOUNDARY_COND 14

/* Plot type, see list in global_ljones.c  */

#define PLOT 0
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

#define ENERGY_HUE_MIN 330.0      /* color of original particle */
#define ENERGY_HUE_MAX 50.0        /* color of saturated particle */
#define PARTICLE_HUE_MIN 330.0      /* color of original particle */
#define PARTICLE_HUE_MAX 50.0        /* color of saturated particle */
#define PARTICLE_EMAX 1.0e3           /* energy of particle with hottest color */
#define HUE_TYPE0 280.0     /* hue of particles of type 0 */
#define HUE_TYPE1 135.0      /* hue of particles of type 1 */

#define RANDOM_RADIUS 0     /* set to 1 for random circle radius */
#define DT_PARTICLE 3.0e-7    /* time step for particle displacement */
#define KREPEL 12.0          /* constant in repelling force between particles */
#define EQUILIBRIUM_DIST 7.5    /* Lennard-Jones equilibrium distance */
#define EQUILIBRIUM_DIST_B 5.0  /* Lennard-Jones equilibrium distance for second type of particle */
#define REPEL_RADIUS 20.0    /* radius in which repelling force acts (in units of particle radius) */
#define DAMPING 0.0     /* damping coefficient of particles */
#define PARTICLE_MASS 1.0   /* mass of particle of radius MU */
#define PARTICLE_MASS_B 1.0   /* mass of particle of radius MU */
#define PARTICLE_INERTIA_MOMENT 0.2     /* moment of inertia of particle */
#define PARTICLE_INERTIA_MOMENT_B 0.02     /* moment of inertia of second type of particle */
#define V_INITIAL 10.0        /* initial velocity range */
#define OMEGA_INITIAL 10.0        /* initial angular velocity range */

#define THERMOSTAT 1        /* set to 1 to switch on thermostat */
#define SIGMA 5.0           /* noise intensity in thermostat */
#define BETA 0.001            /* initial inverse temperature */
#define MU_XI 0.01           /* friction constant in thermostat */
#define KSPRING_BOUNDARY 5.0e9    /* confining harmonic potential outside simulation region */
#define KSPRING_OBSTACLE 5.0e8    /* harmonic potential of obstacles */
#define NBH_DIST_FACTOR 5.5       /* radius in which to count neighbours */
#define GRAVITY 0.0            /* gravity acting on all particles */
#define INCREASE_GRAVITY 0     /* set to 1 to increase gravity during the simulation */
#define GRAVITY_FACTOR 100.0   /* factor by which to increase gravity */
#define GRAVITY_RESTORE_TIME 750    /* time at end of simulation with gravity restored to initial value */

#define ROTATION 0           /* set to 1 to include rotation of particles */
#define COUPLE_ANGLE_TO_THERMOSTAT 0    /* set to 1 to couple angular degrees of freedom to thermostat */
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
#define BETA_FACTOR 100.0   /* factor by which to change BETA during simulation */
#define N_TOSCILLATIONS 2.25   /* number of temperature oscillations in BETA schedule */
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
#define N_T_AVERAGE 50       /* size of temperature averaging window */
#define MAX_PRESSURE 3.0e10  /* pressure shown in "hottest" color */
#define PARTIAL_THERMO_COUPLING 0   /* set to 1 to couple only particles to the right of obstacle to thermostat */
#define PARTIAL_THERMO_SHIFT 0.5    /* distance from obstacle at the right of which particles are coupled to thermostat */

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
#define POSITION_Y_DEPENDENCE 0     /* set to 1 for the separation between particles to be vertical */
#define PRINT_ENTROPY 0     /* set to 1 to compute entropy */

#define PRINT_PARTICLE_NUMBER 0     /* set to 1 to print total number of particles */

#define EHRENFEST_COPY 0    /* set to 1 to add equal number of larger particles (for Ehrenfest model) */

#define LID_MASS 1000.0     /* mass of lid for BC_RECTANGLE_LID b.c. */
#define LID_WIDTH 0.1       /* width of lid for BC_RECTANGLE_LID b.c. */
#define WALL_MASS 2000.0    /* mass of wall for BC_RECTANGLE_WALL b.c. */
#define WALL_FRICTION 0.0   /* friction on wall for BC_RECTANGLE_WALL b.c. */
#define WALL_WIDTH 0.1      /* width of wall for BC_RECTANGLE_WALL b.c. */
#define WALL_VMAX 100.0     /* max speed of wall */
#define WALL_TIME 500       /* time during which to keep wall */

#define FLOOR_FORCE 1      /* set to 1 to limit force on particle to FMAX */
#define FMAX 1.0e9         /* maximal force */
#define FLOOR_OMEGA 1      /* set to 1 to limit particle momentum to PMAX */
#define PMAX 1000.0        /* maximal force */

#define HASHX 40   /* size of hashgrid in x direction */
#define HASHY 22   /* size of hashgrid in y direction */
#define HASHMAX 100  /* maximal number of particles per hashgrid cell */
#define HASHGRID_PADDING 0.1    /* padding of hashgrid outside simulation window */

#define DRAW_COLOR_SCHEME 0     /* set to 1 to plot the color scheme */
#define COLORBAR_RANGE 8.0    /* scale of color scheme bar */
#define COLORBAR_RANGE_B 12.0    /* scale of color scheme bar for 2nd part */
#define ROTATE_COLOR_SCHEME 0   /* set to 1 to draw color scheme horizontally */

```

### 7 April 22 - A diffraction grating in 3D ###

**Program:** `wave_3d.c` 

**Initial condition in function `animation()`:** `init_circular_wave(-1.0, 0.0, phi, psi, xy_in);`

```
#define MOVIE 1         /* set to 1 to generate movie */
#define DOUBLE_MOVIE 1  /* set to 1 to produce movies for wave height and energy simultaneously */

/* General geometrical parameters */

#define WINWIDTH 	1920  /* window width */
#define WINHEIGHT 	1000  /* window height */
#define NX 1920          /* number of grid points on x axis */
#define NY 1000          /* number of grid points on y axis */

#define XMIN -1.8
#define XMAX 2.2	/* x interval  */
#define YMIN -1.041666667
#define YMAX 1.041666667	/* y interval for 9/16 aspect ratio */

#define HIGHRES 0        /* set to 1 if resolution of grid is double that of displayed image */

#define JULIA_SCALE 1.0 /* scaling for Julia sets */

/* Choice of the billiard table */

#define B_DOMAIN 10        /* choice of domain shape, see list in global_pdes.c */

#define CIRCLE_PATTERN 201   /* pattern of circles or polygons, see list in global_pdes.c */

#define P_PERCOL 0.25       /* probability of having a circle in C_RAND_PERCOL arrangement */
#define NPOISSON 300        /* number of points for Poisson C_RAND_POISSON arrangement */
#define RANDOM_POLY_ANGLE 1 /* set to 1 to randomize angle of polygons */

#define LAMBDA 0.18	    /* parameter controlling the dimensions of domain */
#define MU 0.045              /* parameter controlling the dimensions of domain */
#define NPOLY 4             /* number of sides of polygon */
#define APOLY 0.5           /* angle by which to turn polygon, in units of Pi/2 */ 
#define MDEPTH 3            /* depth of computation of Menger gasket */
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

#define TWOSPEEDS 0          /* set to 1 to replace hardcore boundary by medium with different speed */
#define OSCILLATE_LEFT 0     /* set to 1 to add oscilating boundary condition on the left */
#define OSCILLATE_TOPBOT 0   /* set to 1 to enforce a planar wave on top and bottom boundary */

#define OMEGA 0.005        /* frequency of periodic excitation */
#define AMPLITUDE 0.8      /* amplitude of periodic excitation */ 
#define COURANT 0.06       /* Courant number */
#define COURANTB 0.03      /* Courant number in medium B */
#define GAMMA 0.0          /* damping factor in wave equation */
#define GAMMAB 1.0e-7        /* damping factor in wave equation */
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

/* Boundary conditions, see list in global_pdes.c  */

#define B_COND 2

/* Parameters for length and speed of simulation */

#define NSTEPS 1400        /* number of frames of movie */
#define NVID 12          /* number of iterations between images displayed on screen */
#define NSEG 1000         /* number of segments of boundary */
#define INITIAL_TIME 0      /* time after which to start saving frames */
#define BOUNDARY_WIDTH 2    /* width of billiard boundary */

#define PAUSE 200       /* number of frames after which to pause */
#define PSLEEP 2         /* sleep time during pause */
#define SLEEP1  1        /* initial sleeping time */
#define SLEEP2  1        /* final sleeping time */
#define MID_FRAMES 20    /* number of still frames between parts of two-part movie */
#define END_FRAMES 100    /* number of still frames at end of movie */

/* Parameters of initial condition */

#define INITIAL_AMP 0.75         /* amplitude of initial condition */
#define INITIAL_VARIANCE 0.0003  /* variance of initial condition */
#define INITIAL_WAVELENGTH  0.05  /* wavelength of initial condition */

/* Plot type, see list in global_pdes.c  */

#define PLOT 103

#define PLOT_B 104        /* plot type for second movie */

#define P_3D_AMPLITUDE  101     /* color depends on amplitude */
#define P_3D_ANGLE 102          /* color depends on angle with fixed direction */
#define P_3D_AMP_ANGLE 103      /* color depends on amplitude, luminosity depends on angle */
#define P_3D_ENERGY 104         /* color depends on energy, luminosity depends on angle */

#define AMPLITUDE_HIGH_RES 1    /* set to 1 to increase resolution of P_3D_AMPLITUDE plot */
#define FADE_IN_OBSTACLE 1      /* set to 1 to fade color inside obstacles */

#define PLOT_SCALE_ENERGY 0.1      /* vertical scaling in energy plot */

/* 3D representation */

#define REPRESENTATION_3D 1     /* choice of 3D representation */ 

#define REP_AXO_3D 0        /* linear projection (axonometry) */
#define REP_PROJ_3D 1       /* projection on plane orthogonal to observer line of sight */

/* Color schemes */

#define COLOR_PALETTE 16     /* Color palette, see list in global_pdes.c  */
#define COLOR_PALETTE_B 12     /* Color palette, see list in global_pdes.c  */

#define BLACK 1          /* background */

#define COLOR_SCHEME 3   /* choice of color scheme, see list in global_pdes.c  */

#define SCALE 0          /* set to 1 to adjust color scheme to variance of field */
#define SLOPE 1.0       /* sensitivity of color on wave amplitude */
#define VSCALE_AMPLITUDE 1.5    /* additional scaling factor for color scheme P_3D_AMPLITUDE */
#define VSCALE_ENERGY 3.0       /* additional scaling factor for color scheme P_3D_ENERGY */
#define ATTENUATION 0.0  /* exponential attenuation coefficient of contrast with time */
#define E_SCALE 200.0     /* scaling factor for energy representation */
#define LOG_SCALE 1.0     /* scaling factor for energy log representation */
#define LOG_SHIFT 1.0     /* shift of colors on log scale */
#define RESCALE_COLOR_IN_CENTER 0   /* set to 1 to decrease color intentiy in the center (for wave escaping ring) */

#define COLORHUE 260     /* initial hue of water color for scheme C_LUM */
#define COLORDRIFT 0.0   /* how much the color hue drifts during the whole simulation */
#define LUMMEAN 0.5      /* amplitude of luminosity variation for scheme C_LUM */
#define LUMAMP 0.3       /* amplitude of luminosity variation for scheme C_LUM */
#define HUEMEAN 180.0    /* mean value of hue for color scheme C_HUE */
#define HUEAMP -180.0      /* amplitude of variation of hue for color scheme C_HUE */

#define DRAW_COLOR_SCHEME 1     /* set to 1 to plot the color scheme */
#define COLORBAR_RANGE 4.0     /* scale of color scheme bar */
#define COLORBAR_RANGE_B 5.0    /* scale of color scheme bar for 2nd part */
#define ROTATE_COLOR_SCHEME 0   /* set to 1 to draw color scheme horizontally */

#define SAVE_TIME_SERIES 0      /* set to 1 to save wave time series at a point */

/* For debugging purposes only */
#define FLOOR 0         /* set to 1 to limit wave amplitude to VMAX */
#define VMAX 10.0       /* max value of wave amplitude */

/* Parameters controlling 3D projection */

double u_3d[2] = {0.75, -0.45};     /* projections of basis vectors for REP_AXO_3D representation */
double v_3d[2] = {-0.75, -0.45};
double w_3d[2] = {0.0, 0.015};
double light[3] = {0.40824829, 0.816496581, 0.40824829};      /* vector of "light" direction for P_3D_ANGLE color scheme */
double observer[3] = {8.0, 8.0, 8.0};    /* location of observer for REP_PROJ_3D representation */ 

#define Z_SCALING_FACTOR 0.0075   /* overall scaling factor of z axis for REP_PROJ_3D representation */
#define XY_SCALING_FACTOR 2.0     /* overall scaling factor for on-screen (x,y) coordinates after projection */
#define ZMAX_FACTOR 0.5           /* max value of z coordinate for REP_PROJ_3D representation */

```

### 6 April 22 - Lennard-Jones particles in increasing gravity, long version with gravity reset at the end ###

**Program:** `lennardjones.c` 

```
#define MOVIE 1         /* set to 1 to generate movie */
#define DOUBLE_MOVIE 1  /* set to 1 to produce movies for wave height and energy simultaneously */

#define TIME_LAPSE 0     /* set to 1 to add a time-lapse movie at the end */
                         /* so far incompatible with double movie */
#define TIME_LAPSE_FACTOR 4    /* factor of time-lapse movie */


/* General geometrical parameters */

#define WINWIDTH 	1280  /* window width */
#define WINHEIGHT 	720   /* window height */

#define XMIN -2.0
#define XMAX 2.0	/* x interval */
#define YMIN -1.125
#define YMAX 1.125	/* y interval for 9/16 aspect ratio */

#define INITXMIN -1.6
#define INITXMAX 1.6	/* x interval for initial condition */
#define INITYMIN -0.95
#define INITYMAX 0.65	/* y interval for initial condition */

#define BCXMIN -1.95
#define BCXMAX 1.95	/* x interval for boundary condition */
#define BCYMIN -1.05
#define BCYMAX 0.9	/* y interval for boundary condition */

#define OBSXMIN -2.0
#define OBSXMAX 2.0     /* x interval for motion of obstacle */

#define CIRCLE_PATTERN 8   /* pattern of circles, see list in global_ljones.c */

#define ADD_FIXED_OBSTACLES 0   /* set to 1 do add fixed circular obstacles */
#define OBSTACLE_PATTERN 0  /* pattern of obstacles, see list in global_ljones.c */

#define TWO_TYPES 0         /* set to 1 to have two types of particles */
#define TPYE_PROPORTION 0.8 /* proportion of particles of first type */
#define SYMMETRIZE_FORCE 0  /* set to 1 to symmetrize two-particle interaction, only needed if particles are not all the same */
#define CENTER_PX 0         /* set to 1 to center horizontal momentum */
#define CENTER_PY 0         /* set to 1 to center vertical momentum */
#define CENTER_PANGLE 0     /* set to 1 to center angular momentum */

#define INTERACTION 1       /* particle interaction, see list in global_ljones.c */
#define INTERACTION_B 1     /* particle interaction for second type of particle, see list in global_ljones.c */
#define SPIN_INTER_FREQUENCY 5.0 /* angular frequency of spin-spin interaction */
#define SPIN_INTER_FREQUENCY_B 2.0 /* angular frequency of spin-spin interaction for second particle type */

#define P_PERCOL 0.25       /* probability of having a circle in C_RAND_PERCOL arrangement */
#define NPOISSON 100        /* number of points for Poisson C_RAND_POISSON arrangement */
#define PDISC_DISTANCE 3.0  /* minimal distance in Poisson disc process, controls density of particles */
#define PDISC_CANDIDATES 100 /* number of candidates in construction of Poisson disc process */
#define RANDOM_POLY_ANGLE 0 /* set to 1 to randomize angle of polygons */

#define LAMBDA 2.0	    /* parameter controlling the dimensions of domain */
#define MU 0.015 	    /* parameter controlling radius of particles */
#define MU_B 0.0254         /* parameter controlling radius of particles of second type */
#define NPOLY 3             /* number of sides of polygon */
#define APOLY 0.125         /* angle by which to turn polygon, in units of Pi/2 */ 
#define MDEPTH 4            /* depth of computation of Menger gasket */
#define MRATIO 3            /* ratio defining Menger gasket */
#define MANDELLEVEL 1000    /* iteration level for Mandelbrot set */
#define MANDELLIMIT 10.0    /* limit value for approximation of Mandelbrot set */
#define FOCI 1              /* set to 1 to draw focal points of ellipse */
#define NGRIDX 30           /* number of grid point for grid of disks */
#define NGRIDY 20           /* number of grid point for grid of disks */
#define EHRENFEST_RADIUS 0.9    /* radius of container for Ehrenfest urn configuration */
#define EHRENFEST_WIDTH 0.035     /* width of tube for Ehrenfest urn configuration */

#define X_SHOOTER -0.2
#define Y_SHOOTER -0.6
#define X_TARGET 0.4
#define Y_TARGET 0.7        /* shooter and target positions in laser fight */

/* Parameters for length and speed of simulation */

#define NSTEPS 2750      /* number of frames of movie */
#define NVID 200          /* number of iterations between images displayed on screen */
#define NSEG 250         /* number of segments of boundary */
#define INITIAL_TIME 10    /* time after which to start saving frames */
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

#define BOUNDARY_COND 1

/* Plot type, see list in global_ljones.c  */

#define PLOT 0
#define PLOT_B 6        /* plot type for second movie */

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

#define ENERGY_HUE_MIN 330.0      /* color of original particle */
#define ENERGY_HUE_MAX 50.0        /* color of saturated particle */
#define PARTICLE_HUE_MIN 330.0      /* color of original particle */
#define PARTICLE_HUE_MAX 50.0        /* color of saturated particle */
#define PARTICLE_EMAX 2.0e3           /* energy of particle with hottest color */

#define RANDOM_RADIUS 0     /* set to 1 for random circle radius */
#define DT_PARTICLE 3.0e-7    /* time step for particle displacement */
#define KREPEL 12.0          /* constant in repelling force between particles */
#define EQUILIBRIUM_DIST 7.5    /* Lennard-Jones equilibrium distance */
#define EQUILIBRIUM_DIST_B 5.0  /* Lennard-Jones equilibrium distance for second type of particle */
#define REPEL_RADIUS 20.0    /* radius in which repelling force acts (in units of particle radius) */
#define DAMPING 0.0     /* damping coefficient of particles */
#define PARTICLE_MASS 1.0   /* mass of particle of radius MU */
#define PARTICLE_MASS_B 1.0   /* mass of particle of radius MU */
#define PARTICLE_INERTIA_MOMENT 0.2     /* moment of inertia of particle */
#define PARTICLE_INERTIA_MOMENT_B 0.02     /* moment of inertia of second type of particle */
#define V_INITIAL 10.0        /* initial velocity range */
#define OMEGA_INITIAL 10.0        /* initial angular velocity range */

#define THERMOSTAT 1        /* set to 1 to switch on thermostat */
#define SIGMA 5.0           /* noise intensity in thermostat */
#define BETA 0.005            /* initial inverse temperature */
#define MU_XI 0.01           /* friction constant in thermostat */
#define KSPRING_BOUNDARY 5.0e9    /* confining harmonic potential outside simulation region */
#define KSPRING_OBSTACLE 5.0e7    /* harmonic potential of obstacles */
#define NBH_DIST_FACTOR 5.5       /* radius in which to count neighbours */
#define GRAVITY 4000.0         /* gravity acting on all particles */
#define INCREASE_GRAVITY 1     /* set to 1 to increase gravity during the simulation */
#define GRAVITY_FACTOR 100.0   /* factor by which to increase gravity */
#define GRAVITY_RESTORE_TIME 750    /* time at end of simulation with gravity restored to initial value */

#define ROTATION 0           /* set to 1 to include rotation of particles */
#define COUPLE_ANGLE_TO_THERMOSTAT 0    /* set to 1 to couple angular degrees of freedom to thermostat */
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

#define INCREASE_BETA 0  /* set to 1 to increase BETA during simulation */
#define BETA_FACTOR 500.0   /* factor by which to change BETA during simulation */
#define N_TOSCILLATIONS 1.5   /* number of temperature oscillations in BETA schedule */
#define NO_OSCILLATION 1        /* set to 1 to have exponential BETA change only */
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
#define OBSTACLE_RADIUS 0.3  /* radius of obstacle for circle boundary conditions */
#define FUNNEL_WIDTH  0.25   /* funnel width for funnel boundary conditions */
#define OBSTACLE_XMIN 0.0    /* initial position of obstacle */
#define OBSTACLE_XMAX 3.0    /* final position of obstacle */
#define RECORD_PRESSURES 0   /* set to 1 to record pressures on obstacle */
#define N_PRESSURES 100      /* number of intervals to record pressure */
#define N_P_AVERAGE 100      /* size of pressure averaging window */
#define N_T_AVERAGE 50       /* size of temperature averaging window */
#define MAX_PRESSURE 3.0e10  /* pressure shown in "hottest" color */
#define PARTIAL_THERMO_COUPLING 0   /* set to 1 to couple only particles to the right of obstacle to thermostat */
#define PARTIAL_THERMO_SHIFT 0.5    /* distance from obstacle at the right of which particles are coupled to thermostat */

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
#define POSITION_Y_DEPENDENCE 0     /* set to 1 for the separation between particles to be vertical */
#define PRINT_ENTROPY 0     /* set to 1 to compute entropy */

#define PRINT_PARTICLE_NUMBER 0     /* set to 1 to print total number of particles */

#define EHRENFEST_COPY 0    /* set to 1 to add equal number of larger particles (for Ehrenfest model) */

#define LID_MASS 1000.0     /* mass of lid for BC_RECTANGLE_LID b.c. */
#define LID_WIDTH 0.1       /* width of lid for BC_RECTANGLE_LID b.c. */
#define WALL_MASS 2000.0    /* mass of wall for BC_RECTANGLE_WALL b.c. */
#define WALL_FRICTION 0.0   /* friction on wall for BC_RECTANGLE_WALL b.c. */
#define WALL_WIDTH 0.1      /* width of wall for BC_RECTANGLE_WALL b.c. */
#define WALL_VMAX 100.0     /* max speed of wall */
#define WALL_TIME 500       /* time during which to keep wall */

#define FLOOR_FORCE 1      /* set to 1 to limit force on particle to FMAX */
#define FMAX 1.0e9         /* maximal force */
#define FLOOR_OMEGA 1      /* set to 1 to limit particle momentum to PMAX */
#define PMAX 1000.0        /* maximal force */

#define HASHX 40   /* size of hashgrid in x direction */
#define HASHY 22   /* size of hashgrid in y direction */
#define HASHMAX 100  /* maximal number of particles per hashgrid cell */
#define HASHGRID_PADDING 0.1    /* padding of hashgrid outside simulation window */

#define DRAW_COLOR_SCHEME 0     /* set to 1 to plot the color scheme */
#define COLORBAR_RANGE 8.0    /* scale of color scheme bar */
#define COLORBAR_RANGE_B 12.0    /* scale of color scheme bar for 2nd part */
#define ROTATE_COLOR_SCHEME 0   /* set to 1 to draw color scheme horizontally */

```

### 5 April 22 - Young’s double-slit experiment with waves, in 3D ###

**Program:** `wave_3d.c` 

**Initial condition in function `animation()`:** `init_circular_wave(-1.0, 0.0, phi, psi, xy_in);`

```
#define MOVIE 1         /* set to 1 to generate movie */
#define DOUBLE_MOVIE 1  /* set to 1 to produce movies for wave height and energy simultaneously */

/* General geometrical parameters */

#define WINWIDTH 	1920  /* window width */
#define WINHEIGHT 	1000  /* window height */
#define NX 1920          /* number of grid points on x axis */
#define NY 1000          /* number of grid points on y axis */

#define XMIN -1.8
#define XMAX 2.2	/* x interval  */
#define YMIN -1.041666667
#define YMAX 1.041666667	/* y interval for 9/16 aspect ratio */

#define HIGHRES 0        /* set to 1 if resolution of grid is double that of displayed image */

#define JULIA_SCALE 1.0 /* scaling for Julia sets */

/* Choice of the billiard table */

#define B_DOMAIN 9        /* choice of domain shape, see list in global_pdes.c */

#define CIRCLE_PATTERN 201   /* pattern of circles or polygons, see list in global_pdes.c */

#define P_PERCOL 0.25       /* probability of having a circle in C_RAND_PERCOL arrangement */
#define NPOISSON 300        /* number of points for Poisson C_RAND_POISSON arrangement */
#define RANDOM_POLY_ANGLE 1 /* set to 1 to randomize angle of polygons */

#define LAMBDA 0.5	    /* parameter controlling the dimensions of domain */
#define MU 0.04              /* parameter controlling the dimensions of domain */
#define NPOLY 4             /* number of sides of polygon */
#define APOLY 0.5           /* angle by which to turn polygon, in units of Pi/2 */ 
#define MDEPTH 3            /* depth of computation of Menger gasket */
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

#define TWOSPEEDS 0          /* set to 1 to replace hardcore boundary by medium with different speed */
#define OSCILLATE_LEFT 0     /* set to 1 to add oscilating boundary condition on the left */
#define OSCILLATE_TOPBOT 0   /* set to 1 to enforce a planar wave on top and bottom boundary */

#define OMEGA 0.005        /* frequency of periodic excitation */
#define AMPLITUDE 0.8      /* amplitude of periodic excitation */ 
#define COURANT 0.06       /* Courant number */
#define COURANTB 0.03      /* Courant number in medium B */
#define GAMMA 0.0          /* damping factor in wave equation */
#define GAMMAB 1.0e-7        /* damping factor in wave equation */
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

/* Boundary conditions, see list in global_pdes.c  */

#define B_COND 2

/* Parameters for length and speed of simulation */

#define NSTEPS 1600        /* number of frames of movie */
#define NVID 12          /* number of iterations between images displayed on screen */
#define NSEG 1000         /* number of segments of boundary */
#define INITIAL_TIME 0      /* time after which to start saving frames */
#define BOUNDARY_WIDTH 2    /* width of billiard boundary */

#define PAUSE 200       /* number of frames after which to pause */
#define PSLEEP 2         /* sleep time during pause */
#define SLEEP1  1        /* initial sleeping time */
#define SLEEP2  1        /* final sleeping time */
#define MID_FRAMES 20    /* number of still frames between parts of two-part movie */
#define END_FRAMES 100    /* number of still frames at end of movie */

/* Parameters of initial condition */

#define INITIAL_AMP 0.75         /* amplitude of initial condition */
#define INITIAL_VARIANCE 0.0003  /* variance of initial condition */
#define INITIAL_WAVELENGTH  0.05  /* wavelength of initial condition */

/* Plot type, see list in global_pdes.c  */

// #define PLOT 104
#define PLOT 103

#define PLOT_B 104        /* plot type for second movie */

#define P_3D_AMPLITUDE  101     /* color depends on amplitude */
#define P_3D_ANGLE 102          /* color depends on angle with fixed direction */
#define P_3D_AMP_ANGLE 103      /* color depends on amplitude, luminosity depends on angle */
#define P_3D_ENERGY 104         /* color depends on energy, luminosity depends on angle */

#define AMPLITUDE_HIGH_RES 1    /* set to 1 to increase resolution of P_3D_AMPLITUDE plot */
#define FADE_IN_OBSTACLE 1      /* set to 1 to fade color inside obstacles */

#define PLOT_SCALE_ENERGY 0.08      /* vertical scaling in energy plot */

/* Color schemes */

#define COLOR_PALETTE 11     /* Color palette, see list in global_pdes.c  */
#define COLOR_PALETTE_B 13     /* Color palette, see list in global_pdes.c  */

#define BLACK 1          /* background */

#define COLOR_SCHEME 3   /* choice of color scheme, see list in global_pdes.c  */

#define SCALE 0          /* set to 1 to adjust color scheme to variance of field */
#define SLOPE 1.0       /* sensitivity of color on wave amplitude */
#define VSCALE_AMPLITUDE 1.5    /* additional scaling factor for color scheme P_3D_AMPLITUDE */
#define VSCALE_ENERGY 3.0       /* additional scaling factor for color scheme P_3D_ENERGY */
#define ATTENUATION 0.0  /* exponential attenuation coefficient of contrast with time */
#define E_SCALE 200.0     /* scaling factor for energy representation */
#define LOG_SCALE 1.0     /* scaling factor for energy log representation */
#define LOG_SHIFT 1.0     /* shift of colors on log scale */
#define RESCALE_COLOR_IN_CENTER 0   /* set to 1 to decrease color intentiy in the center (for wave escaping ring) */

#define COLORHUE 260     /* initial hue of water color for scheme C_LUM */
#define COLORDRIFT 0.0   /* how much the color hue drifts during the whole simulation */
#define LUMMEAN 0.5      /* amplitude of luminosity variation for scheme C_LUM */
#define LUMAMP 0.3       /* amplitude of luminosity variation for scheme C_LUM */
#define HUEMEAN 180.0    /* mean value of hue for color scheme C_HUE */
#define HUEAMP -180.0      /* amplitude of variation of hue for color scheme C_HUE */

#define DRAW_COLOR_SCHEME 1     /* set to 1 to plot the color scheme */
#define COLORBAR_RANGE 4.0     /* scale of color scheme bar */
#define COLORBAR_RANGE_B 5.0    /* scale of color scheme bar for 2nd part */
#define ROTATE_COLOR_SCHEME 0   /* set to 1 to draw color scheme horizontally */

#define SAVE_TIME_SERIES 0      /* set to 1 to save wave time series at a point */

/* For debugging purposes only */
#define FLOOR 0         /* set to 1 to limit wave amplitude to VMAX */
#define VMAX 10.0       /* max value of wave amplitude */

double u_3d[2] = {0.75, -0.45};     /* projections of basis vectors */
double v_3d[2] = {-0.75, -0.45};
double w_3d[2] = {0.0, 0.015};

```

### 4 April 22 - A wave traveling between parabolic reflectors, in 3D ###

**Program:** `wave_3d.c` 

**Initial condition in function `animation()`:** `init_circular_wave(LAMBDA, 0.0, phi, psi, xy_in);`

```
#define MOVIE 1         /* set to 1 to generate movie */
#define DOUBLE_MOVIE 1  /* set to 1 to produce movies for wave height and energy simultaneously */

/* General geometrical parameters */

#define WINWIDTH 	1920  /* window width */
#define WINHEIGHT 	1000  /* window height */
#define NX 1920          /* number of grid points on x axis */
#define NY 1000          /* number of grid points on y axis */

#define XMIN -1.8
#define XMAX 2.2	/* x interval  */
#define YMIN -1.041666667
#define YMAX 1.041666667	/* y interval for 9/16 aspect ratio */

#define HIGHRES 0        /* set to 1 if resolution of grid is double that of displayed image */

#define JULIA_SCALE 1.0 /* scaling for Julia sets */

/* Choice of the billiard table */

#define B_DOMAIN 19        /* choice of domain shape, see list in global_pdes.c */

#define CIRCLE_PATTERN 201   /* pattern of circles or polygons, see list in global_pdes.c */

#define P_PERCOL 0.25       /* probability of having a circle in C_RAND_PERCOL arrangement */
#define NPOISSON 300        /* number of points for Poisson C_RAND_POISSON arrangement */
#define RANDOM_POLY_ANGLE 1 /* set to 1 to randomize angle of polygons */

#define LAMBDA 1.1	    /* parameter controlling the dimensions of domain */
#define MU 0.5              /* parameter controlling the dimensions of domain */
#define NPOLY 4             /* number of sides of polygon */
#define APOLY 0.5           /* angle by which to turn polygon, in units of Pi/2 */ 
#define MDEPTH 3            /* depth of computation of Menger gasket */
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

#define TWOSPEEDS 0          /* set to 1 to replace hardcore boundary by medium with different speed */
#define OSCILLATE_LEFT 0     /* set to 1 to add oscilating boundary condition on the left */
#define OSCILLATE_TOPBOT 0   /* set to 1 to enforce a planar wave on top and bottom boundary */

#define OMEGA 0.005        /* frequency of periodic excitation */
#define AMPLITUDE 0.8      /* amplitude of periodic excitation */ 
#define COURANT 0.05       /* Courant number */
#define COURANTB 0.03      /* Courant number in medium B */
#define GAMMA 0.0          /* damping factor in wave equation */
#define GAMMAB 1.0e-7        /* damping factor in wave equation */
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

#define B_COND 2

/* Parameters for length and speed of simulation */

#define NSTEPS 1700        /* number of frames of movie */
#define NVID 15          /* number of iterations between images displayed on screen */
#define NSEG 1000         /* number of segments of boundary */
#define INITIAL_TIME 0      /* time after which to start saving frames */
#define BOUNDARY_WIDTH 2    /* width of billiard boundary */

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

#define PLOT 103

#define PLOT_B 104        /* plot type for second movie */

#define P_3D_AMPLITUDE  101     /* color depends on amplitude */
#define P_3D_ANGLE 102          /* color depends on angle with fixed direction */
#define P_3D_AMP_ANGLE 103      /* color depends on amplitude, luminosity depends on angle */
#define P_3D_ENERGY 104         /* color depends on energy, luminosity depends on angle */

#define AMPLITUDE_HIGH_RES 1    /* set to 1 to increase resolution of P_3D_AMPLITUDE plot */
#define FADE_IN_OBSTACLE 1      /* set to 1 to fade color inside obstacles */

/* Color schemes */

#define COLOR_PALETTE 11     /* Color palette, see list in global_pdes.c  */
#define COLOR_PALETTE_B 13     /* Color palette, see list in global_pdes.c  */

#define BLACK 1          /* background */

#define COLOR_SCHEME 3   /* choice of color scheme, see list in global_pdes.c  */

#define SCALE 0          /* set to 1 to adjust color scheme to variance of field */
#define SLOPE 1.0       /* sensitivity of color on wave amplitude */
#define VSCALE_AMPLITUDE 0.7    /* additional scaling factor for color scheme P_3D_AMPLITUDE */
#define VSCALE_ENERGY 0.7       /* additional scaling factor for color scheme P_3D_ENERGY */
#define ATTENUATION 0.0  /* exponential attenuation coefficient of contrast with time */
#define E_SCALE 200.0     /* scaling factor for energy representation */
#define LOG_SCALE 1.0     /* scaling factor for energy log representation */
#define LOG_SHIFT 1.0     /* shift of colors on log scale */
#define RESCALE_COLOR_IN_CENTER 0   /* set to 1 to decrease color intentiy in the center (for wave escaping ring) */

#define COLORHUE 260     /* initial hue of water color for scheme C_LUM */
#define COLORDRIFT 0.0   /* how much the color hue drifts during the whole simulation */
#define LUMMEAN 0.5      /* amplitude of luminosity variation for scheme C_LUM */
#define LUMAMP 0.3       /* amplitude of luminosity variation for scheme C_LUM */
#define HUEMEAN 180.0    /* mean value of hue for color scheme C_HUE */
#define HUEAMP -180.0      /* amplitude of variation of hue for color scheme C_HUE */

#define DRAW_COLOR_SCHEME 1     /* set to 1 to plot the color scheme */
#define COLORBAR_RANGE 4.0     /* scale of color scheme bar */
#define COLORBAR_RANGE_B 5.0    /* scale of color scheme bar for 2nd part */
#define ROTATE_COLOR_SCHEME 0   /* set to 1 to draw color scheme horizontally */

#define SAVE_TIME_SERIES 0      /* set to 1 to save wave time series at a point */

/* For debugging purposes only */
#define FLOOR 0         /* set to 1 to limit wave amplitude to VMAX */
#define VMAX 10.0       /* max value of wave amplitude */

double u_3d[2] = {-0.75, -0.45};     /* projections of basis vectors */
double v_3d[2] = {0.75, -0.45};
double w_3d[2] = {0.0, 0.025};

```

### 3 April 22 - Shock wave of a narrow wedge moving through & Lennard-Jones gas ###

**Program:** `lennardjones.c` 

```
#define MOVIE 1         /* set to 1 to generate movie */
#define DOUBLE_MOVIE 1  /* set to 1 to produce movies for wave height and energy simultaneously */

#define TIME_LAPSE 0     /* set to 1 to add a time-lapse movie at the end */
                         /* so far incompatible with double movie */
#define TIME_LAPSE_FACTOR 4    /* factor of time-lapse movie */


/* General geometrical parameters */

#define WINWIDTH 	1280  /* window width */
#define WINHEIGHT 	720   /* window height */

#define XMIN -2.75
#define XMAX 1.25	/* x interval */
#define YMIN -1.125
#define YMAX 1.125	/* y interval for 9/16 aspect ratio */

#define INITXMIN -3.65
#define INITXMAX 2.15	/* x interval for initial condition */
#define INITYMIN -1.0
#define INITYMAX 1.0	/* y interval for initial condition */

#define BCXMIN -3.75
#define BCXMAX 2.25	/* x interval for boundary condition */
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
#define CENTER_PX 0         /* set to 1 to center horizontal momentum */
#define CENTER_PY 0         /* set to 1 to center vertical momentum */
#define CENTER_PANGLE 0     /* set to 1 to center angular momentum */

#define INTERACTION 1       /* particle interaction, see list in global_ljones.c */
#define INTERACTION_B 1     /* particle interaction for second type of particle, see list in global_ljones.c */
#define SPIN_INTER_FREQUENCY 5.0 /* angular frequency of spin-spin interaction */
#define SPIN_INTER_FREQUENCY_B 2.0 /* angular frequency of spin-spin interaction for second particle type */

#define P_PERCOL 0.25       /* probability of having a circle in C_RAND_PERCOL arrangement */
#define NPOISSON 100        /* number of points for Poisson C_RAND_POISSON arrangement */
#define PDISC_DISTANCE 5.0  /* minimal distance in Poisson disc process, controls density of particles */
#define PDISC_CANDIDATES 100 /* number of candidates in construction of Poisson disc process */
#define RANDOM_POLY_ANGLE 0 /* set to 1 to randomize angle of polygons */

#define LAMBDA 2.0	    /* parameter controlling the dimensions of domain */
#define MU 0.012 	    /* parameter controlling radius of particles */
#define MU_B 0.0254         /* parameter controlling radius of particles of second type */
#define NPOLY 3             /* number of sides of polygon */
#define APOLY 0.125         /* angle by which to turn polygon, in units of Pi/2 */ 
#define MDEPTH 4            /* depth of computation of Menger gasket */
#define MRATIO 3            /* ratio defining Menger gasket */
#define MANDELLEVEL 1000    /* iteration level for Mandelbrot set */
#define MANDELLIMIT 10.0    /* limit value for approximation of Mandelbrot set */
#define FOCI 1              /* set to 1 to draw focal points of ellipse */
#define NGRIDX 30           /* number of grid point for grid of disks */
#define NGRIDY 20           /* number of grid point for grid of disks */
#define EHRENFEST_RADIUS 0.9    /* radius of container for Ehrenfest urn configuration */
#define EHRENFEST_WIDTH 0.035     /* width of tube for Ehrenfest urn configuration */

#define X_SHOOTER -0.2
#define Y_SHOOTER -0.6
#define X_TARGET 0.4
#define Y_TARGET 0.7        /* shooter and target positions in laser fight */

/* Parameters for length and speed of simulation */

// #define NSTEPS 500      /* number of frames of movie */
#define NSTEPS 3200      /* number of frames of movie */
#define NVID 200          /* number of iterations between images displayed on screen */
#define NSEG 250         /* number of segments of boundary */
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

#define BOUNDARY_COND 8

/* Plot type, see list in global_ljones.c  */

#define PLOT 0
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

#define ENERGY_HUE_MIN 330.0      /* color of original particle */
#define ENERGY_HUE_MAX 50.0        /* color of saturated particle */
#define PARTICLE_HUE_MIN 330.0      /* color of original particle */
#define PARTICLE_HUE_MAX 50.0        /* color of saturated particle */
#define PARTICLE_EMAX 1.0e3           /* energy of particle with hottest color */

#define RANDOM_RADIUS 0     /* set to 1 for random circle radius */
#define DT_PARTICLE 3.0e-7    /* time step for particle displacement */
#define KREPEL 12.0          /* constant in repelling force between particles */
#define EQUILIBRIUM_DIST 5.0  /* Lennard-Jones equilibrium distance for second type of particle */
#define EQUILIBRIUM_DIST_B 5.0  /* Lennard-Jones equilibrium distance for second type of particle */
#define REPEL_RADIUS 20.0    /* radius in which repelling force acts (in units of particle radius) */
#define DAMPING 0.0     /* damping coefficient of particles */
#define PARTICLE_MASS 1.0   /* mass of particle of radius MU */
#define PARTICLE_MASS_B 1.0   /* mass of particle of radius MU */
#define PARTICLE_INERTIA_MOMENT 0.2     /* moment of inertia of particle */
#define PARTICLE_INERTIA_MOMENT_B 0.02     /* moment of inertia of second type of particle */
#define V_INITIAL 10.0        /* initial velocity range */
#define OMEGA_INITIAL 10.0        /* initial angular velocity range */

#define THERMOSTAT 1        /* set to 1 to switch on thermostat */
#define SIGMA 5.0           /* noise intensity in thermostat */
#define BETA 0.01           /* initial inverse temperature */
#define MU_XI 0.01            /* friction constant in thermostat */
#define KSPRING_BOUNDARY 5.0e7    /* confining harmonic potential outside simulation region */
#define KSPRING_OBSTACLE 5.0e7    /* harmonic potential of obstacles */
#define NBH_DIST_FACTOR 5.5       /* radius in which to count neighbours */
#define GRAVITY 0.0            /* gravity acting on all particles */
#define INCREASE_GRAVITY 0     /* set to 1 to increase gravity during the simulation */
#define GRAVITY_FACTOR 100.0   /* factor by which to increase gravity */

#define ROTATION 0           /* set to 1 to include rotation of particles */
#define COUPLE_ANGLE_TO_THERMOSTAT 0    /* set to 1 to couple angular degrees of freedom to thermostat */
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

#define INCREASE_BETA 0  /* set to 1 to increase BETA during simulation */
#define BETA_FACTOR 500.0   /* factor by which to change BETA during simulation */
#define N_TOSCILLATIONS 1.5   /* number of temperature oscillations in BETA schedule */
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
#define OBSTACLE_RADIUS 0.3  /* radius of obstacle for circle boundary conditions */
#define FUNNEL_WIDTH  0.25   /* funnel width for funnel boundary conditions */
#define OBSTACLE_XMIN 0.0    /* initial position of obstacle */
#define OBSTACLE_XMAX 3.0    /* final position of obstacle */
#define RECORD_PRESSURES 0   /* set to 1 to record pressures on obstacle */
#define N_PRESSURES 100      /* number of intervals to record pressure */
#define N_P_AVERAGE 100      /* size of pressure averaging window */
#define N_T_AVERAGE 50       /* size of temperature averaging window */
#define MAX_PRESSURE 3.0e10  /* pressure shown in "hottest" color */
#define PARTIAL_THERMO_COUPLING 1   /* set to 1 to couple only particles to the right of obstacle to thermostat */
#define PARTIAL_THERMO_SHIFT 0.5    /* distance from obstacle at the right of which particles are coupled to thermostat */

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
#define POSITION_Y_DEPENDENCE 0     /* set to 1 for the separation between particles to be vertical */
#define PRINT_ENTROPY 0     /* set to 1 to compute entropy */

#define PRINT_PARTICLE_NUMBER 0     /* set to 1 to print total number of particles */

#define EHRENFEST_COPY 0    /* set to 1 to add equal number of larger particles (for Ehrenfest model) */

#define LID_MASS 1000.0     /* mass of lid for BC_RECTANGLE_LID b.c. */
#define LID_WIDTH 0.1       /* width of lid for BC_RECTANGLE_LID b.c. */
#define WALL_MASS 2000.0    /* mass of wall for BC_RECTANGLE_WALL b.c. */
#define WALL_FRICTION 0.0   /* friction on wall for BC_RECTANGLE_WALL b.c. */
#define WALL_WIDTH 0.1      /* width of wall for BC_RECTANGLE_WALL b.c. */
#define WALL_VMAX 100.0     /* max speed of wall */
#define WALL_TIME 500       /* time during which to keep wall */

#define FLOOR_FORCE 1      /* set to 1 to limit force on particle to FMAX */
#define FMAX 1.0e9         /* maximal force */
#define FLOOR_OMEGA 1      /* set to 1 to limit particle momentum to PMAX */
#define PMAX 1000.0        /* maximal force */

#define HASHX 90   /* size of hashgrid in x direction */
#define HASHY 30   /* size of hashgrid in y direction */
#define HASHMAX 100  /* maximal number of particles per hashgrid cell */
#define HASHGRID_PADDING 0.1    /* padding of hashgrid outside simulation window */

#define DRAW_COLOR_SCHEME 0     /* set to 1 to plot the color scheme */
#define COLORBAR_RANGE 8.0    /* scale of color scheme bar */
#define COLORBAR_RANGE_B 12.0    /* scale of color scheme bar for 2nd part */
#define ROTATE_COLOR_SCHEME 0   /* set to 1 to draw color scheme horizontally */

```

### 2 April 22 - 3D waves flowing over Sierpinski carpets ###

**Program:** `wave_3d_sierpinski.c` (variant of `wave_3d.c`) 

**Initial condition in function `animation()`:** `init_circular_wave_sierpinski(1.5, 0.0, phi, psi, xy_in, mdepth);`

```
#define MOVIE 1         /* set to 1 to generate movie */
#define DOUBLE_MOVIE 0  /* set to 1 to produce movies for wave height and energy simultaneously */

/* General geometrical parameters */

#define WINWIDTH 	1920  /* window width */
#define WINHEIGHT 	1000  /* window height */
#define NX 1920          /* number of grid points on x axis */
#define NY 1000          /* number of grid points on y axis */

#define XMIN -1.8
#define XMAX 2.2	/* x interval  */
#define YMIN -1.041666667
#define YMAX 1.041666667	/* y interval for 9/16 aspect ratio */

#define HIGHRES 0        /* set to 1 if resolution of grid is double that of displayed image */

#define JULIA_SCALE 1.0 /* scaling for Julia sets */

/* Choice of the billiard table */

#define B_DOMAIN 15        /* choice of domain shape, see list in global_pdes.c */

#define CIRCLE_PATTERN 201   /* pattern of circles or polygons, see list in global_pdes.c */

#define P_PERCOL 0.25       /* probability of having a circle in C_RAND_PERCOL arrangement */
#define NPOISSON 300        /* number of points for Poisson C_RAND_POISSON arrangement */
#define RANDOM_POLY_ANGLE 1 /* set to 1 to randomize angle of polygons */

#define LAMBDA 1.0	    /* parameter controlling the dimensions of domain */
#define MU 0.3              /* parameter controlling the dimensions of domain */
#define NPOLY 4             /* number of sides of polygon */
#define APOLY 0.5           /* angle by which to turn polygon, in units of Pi/2 */ 
#define MDEPTH 3            /* depth of computation of Menger gasket */
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
#define OSCILLATE_LEFT 0     /* set to 1 to add oscilating boundary condition on the left */
#define OSCILLATE_TOPBOT 0   /* set to 1 to enforce a planar wave on top and bottom boundary */

#define OMEGA 0.005        /* frequency of periodic excitation */
#define AMPLITUDE 0.8      /* amplitude of periodic excitation */ 
#define COURANT 0.05       /* Courant number */
#define COURANTB 0.03      /* Courant number in medium B */
#define GAMMA 0.0          /* damping factor in wave equation */
#define GAMMAB 5.0e-7        /* damping factor in wave equation */
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

#define B_COND 2

/* Parameters for length and speed of simulation */

#define NSTEPS 1000        /* number of frames of movie */
#define NVID 15          /* number of iterations between images displayed on screen */
#define NSEG 1000         /* number of segments of boundary */
#define INITIAL_TIME 0      /* time after which to start saving frames */
#define BOUNDARY_WIDTH 2    /* width of billiard boundary */

#define PAUSE 200       /* number of frames after which to pause */
#define PSLEEP 2         /* sleep time during pause */
#define SLEEP1  1        /* initial sleeping time */
#define SLEEP2  1        /* final sleeping time */
#define MID_FRAMES 20    /* number of still frames between parts of two-part movie */
#define END_FRAMES 50    /* number of still frames at end of movie */

/* Parameters of initial condition */

#define INITIAL_AMP 0.75         /* amplitude of initial condition */
#define INITIAL_VARIANCE 0.00025  /* variance of initial condition */
#define INITIAL_WAVELENGTH  0.015  /* wavelength of initial condition */

/* Plot type, see list in global_pdes.c  */

#define PLOT 103

#define PLOT_B 102        /* plot type for second movie */

#define P_3D_AMPLITUDE  101     /* color depends on amplitude */
#define P_3D_ANGLE 102          /* color depends on angle with fixed direction */
#define P_3D_AMP_ANGLE 103      /* color depends on amplitude, luminosity depends on angle */

#define AMPLITUDE_HIGH_RES 1    /* set to 1 to increase resolution of P_3D_AMPLITUDE plot */
#define FADE_IN_OBSTACLE 1      /* set to 1 to fade color inside obstacles */

/* Color schemes */

#define COLOR_PALETTE 11     /* Color palette, see list in global_pdes.c  */
#define COLOR_PALETTE_B 18     /* Color palette, see list in global_pdes.c  */

#define BLACK 1          /* background */

#define COLOR_SCHEME 3   /* choice of color scheme, see list in global_pdes.c  */

#define SCALE 0          /* set to 1 to adjust color scheme to variance of field */
#define SLOPE 1.0       /* sensitivity of color on wave amplitude */
#define VSCALE_AMPLITUDE 0.7    /* additional scaling factor for color scheme P_3D_AMPLITUDE */
#define ATTENUATION 0.0  /* exponential attenuation coefficient of contrast with time */
#define E_SCALE 300.0     /* scaling factor for energy representation */
#define LOG_SCALE 1.0     /* scaling factor for energy log representation */
#define LOG_SHIFT 1.0     /* shift of colors on log scale */
#define RESCALE_COLOR_IN_CENTER 0   /* set to 1 to decrease color intentiy in the center (for wave escaping ring) */

#define COLORHUE 260     /* initial hue of water color for scheme C_LUM */
#define COLORDRIFT 0.0   /* how much the color hue drifts during the whole simulation */
#define LUMMEAN 0.5      /* amplitude of luminosity variation for scheme C_LUM */
#define LUMAMP 0.3       /* amplitude of luminosity variation for scheme C_LUM */
#define HUEMEAN 180.0    /* mean value of hue for color scheme C_HUE */
#define HUEAMP -180.0      /* amplitude of variation of hue for color scheme C_HUE */

#define DRAW_COLOR_SCHEME 1     /* set to 1 to plot the color scheme */
#define COLORBAR_RANGE 4.0     /* scale of color scheme bar */
#define COLORBAR_RANGE_B 1.0    /* scale of color scheme bar for 2nd part */
#define ROTATE_COLOR_SCHEME 0   /* set to 1 to draw color scheme horizontally */

#define SAVE_TIME_SERIES 0      /* set to 1 to save wave time series at a point */

/* For debugging purposes only */
#define FLOOR 0         /* set to 1 to limit wave amplitude to VMAX */
#define VMAX 10.0       /* max value of wave amplitude */

```

### 1 April 22 - A novel way of generating energy from water ###

**Program:** `lennardjones.c` 

```
#define MOVIE 1         /* set to 1 to generate movie */
#define DOUBLE_MOVIE 0  /* set to 1 to produce movies for wave height and energy simultaneously */

#define TIME_LAPSE 1     /* set to 1 to add a time-lapse movie at the end */
                         /* so far incompatible with double movie */
#define TIME_LAPSE_FACTOR 4    /* factor of time-lapse movie */


/* General geometrical parameters */

#define WINWIDTH 	1280  /* window width */
#define WINHEIGHT 	720   /* window height */

#define XMIN -2.0
#define XMAX 2.0	/* x interval */
#define YMIN -1.125
#define YMAX 1.125	/* y interval for 9/16 aspect ratio */

#define INITXMIN -1.95
#define INITXMAX 1.95	/* x interval for initial condition */
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
#define CENTER_PX 0         /* set to 1 to center horizontal momentum */
#define CENTER_PY 0         /* set to 1 to center vertical momentum */
#define CENTER_PANGLE 0     /* set to 1 to center angular momentum */

#define INTERACTION 7       /* particle interaction, see list in global_ljones.c */
#define INTERACTION_B 1     /* particle interaction for second type of particle, see list in global_ljones.c */
#define SPIN_INTER_FREQUENCY 5.0 /* angular frequency of spin-spin interaction */
#define SPIN_INTER_FREQUENCY_B 2.0 /* angular frequency of spin-spin interaction for second particle type */

#define P_PERCOL 0.25       /* probability of having a circle in C_RAND_PERCOL arrangement */
#define NPOISSON 100        /* number of points for Poisson C_RAND_POISSON arrangement */
#define PDISC_DISTANCE 3.0   /* minimal distance in Poisson disc process, controls density of particles */
#define PDISC_CANDIDATES 100 /* number of candidates in construction of Poisson disc process */
#define RANDOM_POLY_ANGLE 0 /* set to 1 to randomize angle of polygons */

#define LAMBDA 2.0	    /* parameter controlling the dimensions of domain */
#define MU 0.035 	    /* parameter controlling radius of particles */
#define MU_B 0.0254         /* parameter controlling radius of particles of second type */
#define NPOLY 3             /* number of sides of polygon */
#define APOLY 0.0           /* angle by which to turn polygon, in units of Pi/2 */ 
#define MDEPTH 4            /* depth of computation of Menger gasket */
#define MRATIO 3            /* ratio defining Menger gasket */
#define MANDELLEVEL 1000    /* iteration level for Mandelbrot set */
#define MANDELLIMIT 10.0    /* limit value for approximation of Mandelbrot set */
#define FOCI 1              /* set to 1 to draw focal points of ellipse */
#define NGRIDX 30           /* number of grid point for grid of disks */
#define NGRIDY 20           /* number of grid point for grid of disks */
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

#define BOUNDARY_COND 0

/* Plot type, see list in global_ljones.c  */

#define PLOT 0
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

#define ENERGY_HUE_MIN 330.0      /* color of original particle */
#define ENERGY_HUE_MAX 50.0        /* color of saturated particle */
#define PARTICLE_HUE_MIN 330.0      /* color of original particle */
#define PARTICLE_HUE_MAX 50.0        /* color of saturated particle */
#define PARTICLE_EMAX 1.5e3           /* energy of particle with hottest color */

#define RANDOM_RADIUS 0     /* set to 1 for random circle radius */
#define DT_PARTICLE 1.0e-6    /* time step for particle displacement */
#define KREPEL 10.0          /* constant in repelling force between particles */
#define EQUILIBRIUM_DIST 3.75  /* Lennard-Jones equilibrium distance for second type of particle */
#define EQUILIBRIUM_DIST_B 5.0  /* Lennard-Jones equilibrium distance for second type of particle */
#define REPEL_RADIUS 20.0    /* radius in which repelling force acts (in units of particle radius) */
#define DAMPING 0.0     /* damping coefficient of particles */
#define PARTICLE_MASS 1.0   /* mass of particle of radius MU */
#define PARTICLE_MASS_B 1.0   /* mass of particle of radius MU */
#define PARTICLE_INERTIA_MOMENT 0.2     /* moment of inertia of particle */
#define PARTICLE_INERTIA_MOMENT_B 0.02     /* moment of inertia of second type of particle */
#define V_INITIAL 10.0        /* initial velocity range */
#define OMEGA_INITIAL 10.0        /* initial angular velocity range */

#define THERMOSTAT 1        /* set to 1 to switch on thermostat */
#define SIGMA 5.0           /* noise intensity in thermostat */
#define BETA 0.0005            /* initial inverse temperature */
#define MU_XI 0.05           /* friction constant in thermostat */
#define KSPRING_BOUNDARY 5.0e7    /* confining harmonic potential outside simulation region */
#define KSPRING_OBSTACLE 5.0e7    /* harmonic potential of obstacles */
#define NBH_DIST_FACTOR 5.5       /* radius in which to count neighbours */
#define GRAVITY 2000.0          /* gravity acting on all particles */
#define INCREASE_GRAVITY 0     /* set to 1 to increase gravity during the simulation */
#define GRAVITY_FACTOR 100.0   /* factor by which to increase gravity */

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
#define NO_OSCILLATION 1        /* set to 1 to have exponential BETA change only */
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
#define N_T_AVERAGE 50       /* size of temperature averaging window */
#define MAX_PRESSURE 3.0e10  /* pressure shown in "hottest" color */
#define PARTIAL_THERMO_COUPLING 0   /* set to 1 to couple only particles to the right of obstacle to thermostat */

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
#define POSITION_Y_DEPENDENCE 0     /* set to 1 for the separation between particles to be vertical */
#define PRINT_ENTROPY 0     /* set to 1 to compute entropy */

#define PRINT_PARTICLE_NUMBER 0     /* set to 1 to print total number of particles */

#define EHRENFEST_COPY 0    /* set to 1 to add equal number of larger particles (for Ehrenfest model) */

#define LID_MASS 1000.0     /* mass of lid for BC_RECTANGLE_LID b.c. */
#define LID_WIDTH 0.1       /* width of lid for BC_RECTANGLE_LID b.c. */
#define WALL_MASS 2000.0    /* mass of wall for BC_RECTANGLE_WALL b.c. */
#define WALL_FRICTION 0.0   /* friction on wall for BC_RECTANGLE_WALL b.c. */
#define WALL_WIDTH 0.1      /* width of wall for BC_RECTANGLE_WALL b.c. */
#define WALL_VMAX 100.0     /* max speed of wall */
#define WALL_TIME 500       /* time during which to keep wall */

#define FLOOR_FORCE 1      /* set to 1 to limit force on particle to FMAX */
#define FMAX 1.0e8         /* maximal force */
#define FLOOR_OMEGA 1      /* set to 1 to limit particle momentum to PMAX */
#define PMAX 1000.0        /* maximal force */

#define HASHX 25   /* size of hashgrid in x direction */
#define HASHY 15   /* size of hashgrid in y direction */
#define HASHMAX 100  /* maximal number of particles per hashgrid cell */
#define HASHGRID_PADDING 0.1    /* padding of hashgrid outside simulation window */

#define DRAW_COLOR_SCHEME 0     /* set to 1 to plot the color scheme */
#define COLORBAR_RANGE 8.0    /* scale of color scheme bar */
#define COLORBAR_RANGE_B 12.0    /* scale of color scheme bar for 2nd part */
#define ROTATE_COLOR_SCHEME 0   /* set to 1 to draw color scheme horizontally */

```

