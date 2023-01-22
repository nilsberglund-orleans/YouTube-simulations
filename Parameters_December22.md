### Parameter values for YouTube simulations ###

Created by **Nils Berglund** and optimized by **Marco Mancini**

C code for videos on YouTube Channel https://www.youtube.com/c/NilsBerglund

Below are parameter values used for different simulations, as well as initial conditions used in 
function animation. Some simulations use variants of the published code. The list is going to be 
updated gradually. 


### 31 December 22 - Explosive lensing ###

**Program:** `wave_billiard.c` 

**Initial condition in function `animation()`:** 
```
    init_circular_wave(2.0*LAMBDA*cos(APOLY*PID), 2.0*LAMBDA*sin(APOLY*PID), phi, psi, xy_in);
    angle = DPI/(double)NPOLY;
    for (j=1; j<NPOLY; j++)
        add_circular_wave(1.0, 2.0*LAMBDA*cos((double)j*angle + APOLY*PID), 2.0*LAMBDA*sin((double)j*angle + APOLY*PID), phi, psi, xy_in);
        
    if (ALTERNATE_OSCILLATING_SOURCE) sign = -sign;
    for (j=0; j<NPOLY; j++)
        add_circular_wave(sign, 2.0*LAMBDA*cos((double)j*angle + APOLY*PID), 2.0*LAMBDA*sin((double)j*angle + APOLY*PID), phi, psi, xy_in);
```

```
#define MOVIE 1         /* set to 1 to generate movie */
#define DOUBLE_MOVIE 1  /* set to 1 to produce movies for wave height and energy simultaneously */
#define SAVE_MEMORY 1   /* set to 1 to save memory when writing tiff images */
#define NO_EXTRA_BUFFER_SWAP 1    /* some OS require one less buffer swap when recording images */

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

#define B_DOMAIN 59        /* choice of domain shape, see list in global_pdes.c */

#define CIRCLE_PATTERN 1   /* pattern of circles or polygons, see list in global_pdes.c */

#define COMPARISON 0        /* set to 1 to compare two different patterns (beta) */
#define B_DOMAIN_B 20       /* second domain shape, for comparisons */
#define CIRCLE_PATTERN_B 0  /* second pattern of circles or polygons */

#define P_PERCOL 0.25       /* probability of having a circle in C_RAND_PERCOL arrangement */
#define NPOISSON 300        /* number of points for Poisson C_RAND_POISSON arrangement */
#define RANDOM_POLY_ANGLE 1 /* set to 1 to randomize angle of polygons */

#define LAMBDA 0.5	    /* parameter controlling the dimensions of domain */
#define MU 0.5              /* parameter controlling the dimensions of domain */
#define NPOLY 6             /* number of sides of polygon */
#define APOLY 0.0           /* angle by which to turn polygon, in units of Pi/2 */ 
#define MDEPTH 6            /* depth of computation of Menger gasket */
#define MRATIO 3            /* ratio defining Menger gasket */
#define MANDELLEVEL 1000    /* iteration level for Mandelbrot set */
#define MANDELLIMIT 10.0    /* limit value for approximation of Mandelbrot set */
#define FOCI 1              /* set to 1 to draw focal points of ellipse */
#define NGRIDX 12           /* number of grid point for grid of disks */
#define NGRIDY 15           /* number of grid point for grid of disks */

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
#define OSCILLATION_SCHEDULE 1  /* oscillation schedule, see list in global_pdes.c */

#define OMEGA 0.0005       /* frequency of periodic excitation */
#define AMPLITUDE 0.8      /* amplitude of periodic excitation */ 
#define ACHIRP 0.25        /* acceleration coefficient in chirp */
#define DAMPING 0.0        /* damping of periodic excitation */
#define COURANT 0.1        /* Courant number */
#define COURANTB 0.04658753  /* Courant number in medium B */
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
#define OSCILLATING_SOURCE_PERIOD 50    /* period of oscillating source */
#define ALTERNATE_OSCILLATING_SOURCE 1  /* set to 1 to alternate sign of oscillating source */

/* Boundary conditions, see list in global_pdes.c  */

#define B_COND 2

/* Parameters for length and speed of simulation */

#define NSTEPS 1800       /* number of frames of movie */
#define NVID 9           /* number of iterations between images displayed on screen */
#define NSEG 1000         /* number of segments of boundary */
#define INITIAL_TIME 0      /* time after which to start saving frames */
#define BOUNDARY_WIDTH 2    /* width of billiard boundary */
#define PRINT_SPEED 0       /* print speed of moving source */

#define PAUSE 200       /* number of frames after which to pause */
#define PSLEEP 1         /* sleep time during pause */
#define SLEEP1  1        /* initial sleeping time */
#define SLEEP2  1        /* final sleeping time */
#define MID_FRAMES 20    /* number of still frames between parts of two-part movie */
#define END_FRAMES 100    /* number of still frames at end of movie */
#define FADE 1           /* set to 1 to fade at end of movie */

/* Parameters of initial condition */

#define INITIAL_AMP 0.35         /* amplitude of initial condition */
#define INITIAL_VARIANCE 0.0002  /* variance of initial condition */
#define INITIAL_WAVELENGTH  0.01  /* wavelength of initial condition */

/* Plot type, see list in global_pdes.c  */

#define PLOT 1

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
#define E_SCALE 40.0     /* scaling factor for energy representation */
#define LOG_SCALE 1.0     /* scaling factor for energy log representation */
#define LOG_SHIFT 3.5     /* shift of colors on log scale */
#define FLUX_SCALE 1.0e3    /* scaling factor for enegy flux represtnation */
#define RESCALE_COLOR_IN_CENTER 0   /* set to 1 to decrease color intentiy in the center (for wave escaping ring) */

#define COLORHUE 260     /* initial hue of water color for scheme C_LUM */
#define COLORDRIFT 0.0   /* how much the color hue drifts during the whole simulation */
#define LUMMEAN 0.5      /* amplitude of luminosity variation for scheme C_LUM */
#define LUMAMP 0.3       /* amplitude of luminosity variation for scheme C_LUM */
#define HUEMEAN 180.0    /* mean value of hue for color scheme C_HUE */
#define HUEAMP -180.0      /* amplitude of variation of hue for color scheme C_HUE */

#define DRAW_COLOR_SCHEME 1    /* set to 1 to plot the color scheme */
#define COLORBAR_RANGE 0.75     /* scale of color scheme bar */
#define COLORBAR_RANGE_B 0.05   /* scale of color scheme bar for 2nd part */
#define ROTATE_COLOR_SCHEME 0   /* set to 1 to draw color scheme horizontally */

#define SAVE_TIME_SERIES 0      /* set to 1 to save wave time series at a point */

#define NXMAZE 7      /* width of maze */
#define NYMAZE 7      /* height of maze */
#define MAZE_MAX_NGBH 4     /* max number of neighbours of maze cell */
#define RAND_SHIFT 24        /* seed of random number generator */
#define MAZE_XSHIFT 0.0     /* horizontal shift of maze */

/* for compatibility with sub_wave and sub_maze */
#define ADD_POTENTIAL 0
#define POT_MAZE 7
#define POTENTIAL 0
/* end of constants only used by sub_wave and sub_maze */


/* For debugging purposes only */
#define FLOOR 0         /* set to 1 to limit wave amplitude to VMAX */
#define VMAX 10.0       /* max value of wave amplitude */

```

### 30 December 22 - Yet another, only marginally successful attempt to solve a maze with an incompressible fluid ###

**Program:** `rde.c` 

**Initial condition in function `animation()`:** `init_laminar_flow(flow_speed_schedule(0), LAMINAR_FLOW_MODULATION, 0.5, 0.0, phi, xy_in);`

```
#define MOVIE 1         /* set to 1 to generate movie */
#define DOUBLE_MOVIE 1  /* set to 1 to produce movies for wave height and energy simultaneously */
#define SAVE_MEMORY 1           /* set to 1 to save memory when writing tiff images */
#define NO_EXTRA_BUFFER_SWAP 1    /* some OS require one less buffer swap when recording images */

/* General geometrical parameters */

#define WINWIDTH 	1920  /* window width */
#define WINHEIGHT 	1150  /* window height */
#define NX 960          /* number of grid points on x axis */
#define NY 575          /* number of grid points on y axis */

#define XMIN -2.0
#define XMAX 2.0	/* x interval */
#define YMIN -1.197916667
#define YMAX 1.197916667	/* y interval for 9/16 aspect ratio */

/* Choice of simulated equation */

#define RDE_EQUATION 6  /* choice of reaction term, see list in global_3d.c */
#define NFIELDS 3       /* number of fields in reaction-diffusion equation */
#define NLAPLACIANS 1   /* number of fields for which to compute Laplacian */

#define ADD_POTENTIAL 0 /* set to 1 to add a potential (for Schrodinger equation) */
#define ADD_MAGNETIC_FIELD 0    /* set to 1 to add a magnetic field (for Schrodinger equation) - then set POTENTIAL 1 */
#define POTENTIAL 7     /* type of potential or vector potential, see list in global_3d.c  */

#define ANTISYMMETRIZE_WAVE_FCT 0   /* set tot 1 to make wave function antisymmetric */

#define JULIA_SCALE 0.5 /* scaling for Julia sets */

/* Choice of the billiard table */

#define B_DOMAIN 541          /* choice of domain shape, see list in global_pdes.c  */

#define CIRCLE_PATTERN 99    /* pattern of circles, see list in global_pdes.c */

#define P_PERCOL 0.25       /* probability of having a circle in C_RAND_PERCOL arrangement */
#define NPOISSON 300        /* number of points for Poisson C_RAND_POISSON arrangement */
#define RANDOM_POLY_ANGLE 0 /* set to 1 to randomize angle of polygons */

#define LAMBDA 0.3	    /* parameter controlling the dimensions of domain */
#define MU 0.5	            /* parameter controlling the dimensions of domain */
#define NPOLY 5             /* number of sides of polygon */
#define APOLY 2.0          /* angle by which to turn polygon, in units of Pi/2 */
#define MDEPTH 7            /* depth of computation of Menger gasket */
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

#define DT 0.0000012

#define VISCOSITY 2.0

#define RPSA 0.75         /* parameter in Rock-Paper-Scissors-type interaction */
#define RPSLZB 0.75       /* second parameter in Rock-Paper-Scissors-Lizard-Spock type interaction */

#define EPSILON 0.8     /* time scale separation */
#define DELTA 0.1       /* time scale separation */
#define FHNA 1.0        /* parameter in FHN equation */
#define FHNC -0.01      /* parameter in FHN equation */
#define K_HARMONIC 1.0  /* spring constant of harmonic potential */
#define K_COULOMB 0.5   /* constant in Coulomb potential */
#define V_MAZE 0.4      /* potential in walls of maze */
#define BZQ 0.0008      /* parameter in BZ equation */
#define BZF 1.2         /* parameter in BZ equation */
#define B_FIELD 10.0    /* magnetic field */
#define AB_RADIUS 0.2   /* radius of region with magnetic field for Aharonov-Bohm effect */
#define K_EULER 50.0    /* constant in stream function integration of Euler equation */

#define SMOOTHEN_VORTICITY 1    /* set to 1 to smoothen vorticity field in Euler equation */
#define SMOOTHEN_PERIOD 10      /* period between smoothenings */
#define SMOOTH_FACTOR 0.05      /* factor by which to smoothen */

#define ADD_TRACERS 1    /* set to 1 to add tracer particles (for Euler equations) */
#define N_TRACERS 1000    /* number of tracer particles */

#define T_OUT 2.0       /* outside temperature */
#define T_IN 0.0        /* inside temperature */
#define SPEED 0.0       /* speed of drift to the right */

#define ADD_NOISE 0     /* set to 1 to add noise, set to 2 to add noise in right half */
#define NOISE_INTENSITY 0.005      /* noise intensity */
#define CHANGE_NOISE 1      /* set to 1 to increase noise intensity */
#define NOISE_FACTOR 40.0   /* factor by which to increase noise intensity */
#define NOISE_INITIAL_TIME 100  /* initial time during which noise remains constant */

#define CHANGE_VISCOSITY 0      /* set to 1 to change the viscosity in the course of the simulation */
#define ADJUST_INTSTEP 0       /* set to 1 to decrease integration step when viscosity increases */
#define VISCOSITY_INITIAL_TIME 10  /* initial time during which viscosity remains constant */
#define VISCOSITY_FACTOR 100.0   /* factor by which to change viscosity */
#define VISCOSITY_MAX 2.0        /* max value of viscosity beyond which NVID is increased and integration step is decrase, 
                                    for numerical stability */
                                        
#define CHANGE_RPSLZB 0         /* set to 1 to change second parameter in Rock-Paper-Scissors-Lizard-Spock equation */
#define RPSLZB_CHANGE 0.75      /* factor by which to rpslzb parameter */
#define RPSLZB_INITIAL_TIME 0   /* initial time during which rpslzb remains constant */
#define RPSLZB_FINAL_TIME 500   /* final time during which rpslzb remains constant */
                                      
#define CHANGE_FLOW_SPEED 0     /* set to 1 to change speed of laminar flow */

/* Boundary conditions, see list in global_pdes.c  */

#define B_COND 1

#define IN_OUT_FLOW_BC 1         /* set to 1 to add in-flow/out-flow boundary conditions (for Euler equation) */
#define IN_OUT_FLOW_MIN_AMP 0.35  /* amplitude of in-flow/out-flow boundary conditions (for Euler equation) */
#define IN_OUT_FLOW_AMP 0.35      /* amplitude of in-flow/out-flow boundary conditions (for Euler equation) */
#define LAMINAR_FLOW_MODULATION 0.3     /* asymmetry of laminar flow */

#define EULER_GRADIENT_YSHIFT 0.0    /* y-shift in computation of gradient in Euler equation */

/* Parameters for length and speed of simulation */

#define NSTEPS 1600           /* number of frames of movie */
#define NVID 120          /* number of iterations between images displayed on screen */
#define ACCELERATION_FACTOR 1.0 /* factor by which to increase NVID in course of simulation */
#define DT_ACCELERATION_FACTOR 1.0 /* factor by which to increase time step in course of simulation  */
#define MAX_DT 0.024     /* maximal value of integration step */
#define NSEG 100         /* number of segments of boundary */
#define BOUNDARY_WIDTH 2    /* width of billiard boundary */

#define PAUSE 100       /* number of frames after which to pause */
#define PSLEEP 2         /* sleep time during pause */
#define SLEEP1  2        /* initial sleeping time */
#define SLEEP2  1        /* final sleeping time */
#define INITIAL_TIME 0  /* initial still time */
#define MID_FRAMES 50    /* number of still frames between parts of two-part movie */
#define END_FRAMES 50    /* number of still frames at end of movie */
#define FADE 1           /* set to 1 to fade at end of movie */

/* Visualisation */

#define PLOT_3D 0    /* controls whether plot is 2D or 3D */

#define ROTATE_VIEW 0       /* set to 1 to rotate position of observer */
#define ROTATE_ANGLE 360.0  /* total angle of rotation during simulation */

#define DRAW_PERIODICISED 1     /* set to 1 to repeat wave periodically in x and y directions */

/* Plot type - color scheme */

#define CPLOT 50
#define CPLOT_B 54

/* Plot type - height of 3D plot */

#define ZPLOT 52     /* z coordinate in 3D plot */
#define ZPLOT_B 51    /* z coordinate in second 3D plot */

#define AMPLITUDE_HIGH_RES 1    /* set to 1 to increase resolution of P_3D_AMPLITUDE plot */
#define SHADE_3D 1              /* set to 1 to change luminosity according to normal vector */
#define NON_DIRICHLET_BC 0      /* set to 1 to draw only facets in domain, if field is not zero on boundary */
#define WRAP_ANGLE 1            /* experimental: wrap angle to [0, 2Pi) for interpolation in angle schemes */
#define FADE_IN_OBSTACLE 0      /* set to 1 to fade color inside obstacles */
#define DRAW_OUTSIDE_GRAY 0     /* experimental - draw outside of billiard in gray */
#define ADD_POTENTIAL_TO_Z 1    /* set to 1 to add the external potential to z-coordinate of plot */
#define ADD_POT_CONSTANT 0.35   /* constant in front of added potential */

#define PLOT_SCALE_ENERGY 0.05      /* vertical scaling in energy plot */

#define PRINT_TIME 0        /* set to 1 to print running time */
#define PRINT_VISCOSITY 0   /* set to 1 to print viscosity */
#define PRINT_RPSLZB 0      /* set to 1 to print rpslzb parameter */
#define PRINT_PROBABILITIES 0   /* set to 1 to print probabilities (for Ehrenfest urn configuration) */
#define PRINT_NOISE 0       /* set to 1 to print noise intensity */
#define PRINT_FLOW_SPEED 1      /* set to 1 to print speed of flow */

#define DRAW_FIELD_LINES 0  /* set to 1 to draw field lines */
#define FIELD_LINE_WIDTH 1  /* width of field lines */
#define N_FIELD_LINES 120   /* number of field lines */
#define FIELD_LINE_FACTOR 120 /* factor controlling precision when computing origin of field lines */
#define DRAW_BILLIARD 1     /* set to 1 to draw boundary */
#define DRAW_BILLIARD_FRONT 1     /* set to 1 to draw boundary */
#define FILL_BILLIARD_COMPLEMENT 1  /* set to 1 to fill complement of billiard (for certain shapes only) */

/* 3D representation */

#define REPRESENTATION_3D 1     /* choice of 3D representation */ 

#define REP_AXO_3D 0        /* linear projection (axonometry) */
#define REP_PROJ_3D 1       /* projection on plane orthogonal to observer line of sight */

/* Color schemes, see list in global_pdes.c  */

#define COLOR_PALETTE 10       /* Color palette, see list in global_pdes.c  */
#define COLOR_PALETTE_B 11     /* Color palette, see list in global_pdes.c  */

#define BLACK 1          /* black background */

#define COLOR_SCHEME 3   /* choice of color scheme */

#define COLOR_PHASE_SHIFT 0.0   /* phase shift of color scheme, in units of Pi */

#define SCALE 0          /* set to 1 to adjust color scheme to variance of field */
#define SLOPE 1.0        /* sensitivity of color on wave amplitude */
#define VSCALE_AMPLITUDE 0.7      /* additional scaling factor for color scheme P_3D_AMPLITUDE */
#define ATTENUATION 0.0  /* exponential attenuation coefficient of contrast with time */
#define CURL_SCALE 0.000015   /* scaling factor for curl representation */
#define RESCALE_COLOR_IN_CENTER 0   /* set to 1 to decrease color intentiy in the center (for wave escaping ring) */
#define SLOPE_SCHROD_LUM 50.0       /* sensitivity of luminosity on module, for color scheme Z_ARGUMENT */
#define MIN_SCHROD_LUM 0.2       /* minimal luminosity in color scheme Z_ARGUMENT*/
#define VSCALE_PRESSURE 0.5      /* additional scaling factor for color scheme Z_EULER_PRESSURE */
#define PRESSURE_SHIFT 50.0        /* shift for color scheme Z_EULER_PRESSURE */
#define PRESSURE_LOG_SHIFT -2.0     /* shift for color scheme Z_EULER_PRESSURE */

#define COLORHUE 260     /* initial hue of water color for scheme C_LUM */
#define COLORDRIFT 0.0   /* how much the color hue drifts during the whole simulation */
#define LUMMEAN 0.5      /* amplitude of luminosity variation for scheme C_LUM */
#define LUMAMP 0.3       /* amplitude of luminosity variation for scheme C_LUM */
#define HUEMEAN 359.0    /* mean value of hue for color scheme C_HUE */
#define HUEAMP -359.0    /* amplitude of variation of hue for color scheme C_HUE */
#define E_SCALE 100.0    /* scaling factor for energy representation */
#define FLUX_SCALE 100.0 /* scaling factor for energy representation */
#define LOG_SCALE 0.5    /* scaling factor for energy log representation */
#define LOG_SHIFT 1.0   
#define LOG_MIN 1.0e-3   /* floor value for log vorticity plot */

#define NXMAZE 7      /* width of maze */
#define NYMAZE 7      /* height of maze */
#define MAZE_MAX_NGBH 4     /* max number of neighbours of maze cell */
#define RAND_SHIFT 0        /* seed of random number generator */
#define MAZE_XSHIFT 0.0     /* horizontal shift of maze */

#define DRAW_COLOR_SCHEME 1     /* set to 1 to plot the color scheme */
#define COLORBAR_RANGE 2.0      /* scale of color scheme bar */
#define COLORBAR_RANGE_B 2.0   /* scale of color scheme bar for 2nd part */
#define ROTATE_COLOR_SCHEME 0   /* set to 1 to draw color scheme horizontally */

/* only for compatibility with wave_common.c */
#define TWOSPEEDS 0          /* set to 1 to replace hardcore boundary by medium with different speed */
#define OMEGA 0.005        /* frequency of periodic excitation */
#define COURANT 0.08       /* Courant number */
#define COURANTB 0.03      /* Courant number in medium B */
#define INITIAL_AMP 0.5         /* amplitude of initial condition */
#define INITIAL_VARIANCE 0.0002  /* variance of initial condition */
#define INITIAL_WAVELENGTH  0.1  /* wavelength of initial condition */
#define VSCALE_ENERGY 200.0       /* additional scaling factor for color scheme P_3D_ENERGY */
#define PHASE_FACTOR 20.0       /* factor in computation of phase in color scheme P_3D_PHASE */
#define PHASE_SHIFT 0.0      /* shift of phase in color scheme P_3D_PHASE */
#define OSCILLATION_SCHEDULE 0  /* oscillation schedule, see list in global_pdes.c */
#define AMPLITUDE 0.8      /* amplitude of periodic excitation */ 
#define ACHIRP 0.2        /* acceleration coefficient in chirp */
#define DAMPING 0.0        /* damping of periodic excitation */
#define COMPARISON 0        /* set to 1 to compare two different patterns (beta) */
#define B_DOMAIN_B 20       /* second domain shape, for comparisons */
#define CIRCLE_PATTERN_B 0  /* second pattern of circles or polygons */
/* end of constants added only for compatibility with wave_common.c */


double u_3d[2] = {0.75, -0.45};     /* projections of basis vectors for REP_AXO_3D representation */
double v_3d[2] = {-0.75, -0.45};
double w_3d[2] = {0.0, 0.015};
double light[3] = {0.816496581, -0.40824829, 0.40824829};      /* vector of "light" direction for P_3D_ANGLE color scheme */
double observer[3] = {8.0, 8.0, 8.0};    /* location of observer for REP_PROJ_3D representation */ 
int reset_view = 0;         /* switch to reset 3D view parameters (for option ROTATE_VIEW) */

#define Z_SCALING_FACTOR 0.08  /* overall scaling factor of z axis for REP_PROJ_3D representation */
#define XY_SCALING_FACTOR 1.7  /* overall scaling factor for on-screen (x,y) coordinates after projection */
#define ZMAX_FACTOR 1.0        /* max value of z coordinate for REP_PROJ_3D representation */
#define XSHIFT_3D -0.1         /* overall x shift for REP_PROJ_3D representation */
#define YSHIFT_3D 0.1          /* overall y shift for REP_PROJ_3D representation */
#define BORDER_PADDING 0       /* distance from boundary at which to plot points, to avoid boundary effects due to gradient */

/* For debugging purposes only */
#define FLOOR 1         /* set to 1 to limit wave amplitude to VMAX */
#define VMAX 1000.0        /* max value of wave amplitude */
#define TEST_GRADIENT 0 /* print norm squared of gradient */

```

### 29 December 22 - Are circular mazes unilluminable? ###

**Program:** `particle_billiard.c` 

**Initial condition in function `animation()`:** `init_drop_config(-0.1, 0.1, 0.0, DPI, configs);`

```
#define MOVIE 1         /* set to 1 to generate movie */
#define SAVE_MEMORY 1           /* set to 1 to save memory when writing tiff images */

#define WINWIDTH 	1280  /* window width */
#define WINHEIGHT 	720   /* window height */

#define XMIN -2.0
#define XMAX 2.0	/* x interval */
#define YMIN -1.125
#define YMAX 1.125	/* y interval for 9/16 aspect ratio */

#define SCALING_FACTOR 1.0       /* scaling factor of drawing, needed for flower billiards, otherwise set to 1.0 */

/* Choice of the billiard table, see global_particles.c */

#define B_DOMAIN 31     /* choice of domain shape */

#define CIRCLE_PATTERN 1    /* pattern of circles */
#define POLYLINE_PATTERN 13  /* pattern of polyline */

#define ABSORBING_CIRCLES 0 /* set to 1 for circular scatterers to be absorbing */

#define NMAXCIRCLES 100000     /* total number of circles (must be at least NCX*NCY for square grid) */
#define NMAXPOLY 100000        /* total number of sides of polygonal line */   
#define NCX 30            /* number of circles in x direction */
#define NCY 20            /* number of circles in y direction */
#define NPOISSON 500        /* number of points for Poisson C_RAND_POISSON arrangement */
#define NGOLDENSPIRAL 2000  /* max number of points for C_GOLDEN_SPIRAL arrandement */
#define SDEPTH 1            /* Sierpinski gastket depth */

#define LAMBDA 1.5	/* parameter controlling shape of domain */
#define MU 0.005          /* second parameter controlling shape of billiard */
#define FOCI 1          /* set to 1 to draw focal points of ellipse */
#define NPOLY 6             /* number of sides of polygon */
#define APOLY 0.0           /* angle by which to turn polygon, in units of Pi/2 */ 
#define PENROSE_RATIO 2.5    /* parameter controlling the shape of small ellipses in Penrose room */

#define DRAW_BILLIARD 1     /* set to 1 to draw billiard */
#define DRAW_CONSTRUCTION_LINES 0   /* set to 1 to draw additional construction lines for billiard */
#define PERIODIC_BC 0       /* set to 1 to enforce periodic boundary conditions when drawing particles */

#define RESAMPLE 0      /* set to 1 if particles should be added when dispersion too large */
#define DEBUG 0         /* draw trajectories, for debugging purposes */

/* Simulation parameters */

#define NPART 10000    /* number of particles */
#define NPARTMAX 100000	/* maximal number of particles after resampling */
#define LMAX 0.01       /* minimal segment length triggering resampling */ 
#define DMIN 0.02       /* minimal distance to boundary for triggering resampling */ 
#define CYCLE 1         /* set to 1 for closed curve (start in all directions) */
#define SHOWTRAILS 1    /* set to 1 to keep trails of the particles */
#define SHOWZOOM 0      /* set to 1 to show zoom on specific area */
#define PRINT_PARTICLE_NUMBER 0 /* set to 1 to print number of particles */
#define PRINT_LEFT_RIGHT_PARTICLE_NUMBER 0 /* set to 1 to print number of particles on left and right side */
#define PRINT_CIRCLE_PARTICLE_NUMBER 1 /* set to 1 to print number of particles outside circular maze */
#define PRINT_COLLISION_NUMBER 0 /* set to 1 to print number of collisions */
#define TEST_ACTIVE 1   /* set to 1 to test whether particle is in billiard */

#define TEST_INITIAL_COND 0     /* set to 1 to allow only initial conditions that pass a test */

#define NSTEPS 7200      /* number of frames of movie */
#define TIME 1500        /* time between movie frames, for fluidity of real-time simulation */ 
#define DPHI 0.00002     /* integration step */
#define NVID 25          /* number of iterations between images displayed on screen */
#define END_FRAMES 25    /* number of still frames at the end of the movie */

/* Decreasing TIME accelerates the animation and the movie                               */
/* For constant speed of movie, TIME*DPHI should be kept constant                        */
/* However, increasing DPHI too much deterioriates quality of simulation                 */
/* NVID tells how often a picture is drawn in the animation, increase it for faster anim */
/* For a good quality movie, take for instance TIME = 400, DPHI = 0.00005, NVID = 100    */

/* Colors and other graphical parameters */

#define COLOR_PALETTE 18     /* Color palette, see list in global_pdes.c  */

#define NCOLORS 250      /* number of colors */
#define COLORSHIFT 0     /* hue of initial color */ 
#define COLOR_HUEMIN 30  /* minimal color hue */
#define COLOR_HUEMAX 330 /* maximal color hue */
#define RAINBOW_COLOR 1  /* set to 1 to use different colors for all particles */
#define FLOWER_COLOR 0   /* set to 1 to adapt initial colors to flower billiard (tracks vs core) */
#define NSEG 100         /* number of segments of boundary */
#define LENGTH 0.01       /* length of velocity vectors */
#define BILLIARD_WIDTH 2    /* width of billiard */
#define PARTICLE_WIDTH 2    /* width of particles */
#define FRONT_WIDTH 3       /* width of wave front */

#define BLACK 1             /* set to 1 for black background */
#define COLOR_OUTSIDE 0     /* set to 1 for colored outside */ 
#define OUTER_COLOR 270.0   /* color outside billiard */
#define PAINT_INT 0         /* set to 1 to paint interior in other color (for polygon/Reuleaux) */
#define PAINT_EXT 1         /* set to 1 to paint exterior */

#define PAUSE 1000       /* number of frames after which to pause */
#define PSLEEP 2         /* sleep time during pause */
#define SLEEP1  1        /* initial sleeping time */
#define SLEEP2  1       /* final sleeping time */

#define NXMAZE 8      /* width of maze */
#define NYMAZE 32      /* height of maze */
#define MAZE_MAX_NGBH 5     /* max number of neighbours of maze cell */
#define RAND_SHIFT 0        /* seed of random number generator */
#define MAZE_XSHIFT 0.0     /* horizontal shift of maze */
#define MAZE_RANDOM_FACTOR 0.1     /* randomization factor for S_MAZE_RANDOM */
#define MAZE_CORNER_RADIUS 0.5     /* radius of tounded corners in maze */

```

### 28 December 22 - A rainbow sprinkler ###

**Program:** `lennardjones.c` postprocessed with `ljones_movie.c` 

```
#define MOVIE 0         /* set to 1 to generate movie */
#define DOUBLE_MOVIE 0  /* set to 1 to produce movies for wave height and energy simultaneously */
#define SAVE_MEMORY 1   /* set to 1 to save memory while saving frames */

#define TIME_LAPSE 1     /* set to 1 to add a time-lapse movie at the end */
                         /* so far incompatible with double movie */
#define TIME_LAPSE_FACTOR 3    /* factor of time-lapse movie */
#define TIME_LAPSE_FIRST 1  /* set to 1 to show time-lapse version first */

#define SAVE_TIME_SERIES 1  /* set to 1 to save time series of particle positions */

/* General geometrical parameters */

#define WINWIDTH 	1280  /* window width */
#define WINHEIGHT 	720   /* window height */

#define XMIN -2.0
#define XMAX 2.0	/* x interval */
#define YMIN -1.125
#define YMAX 1.125	/* y interval for 9/16 aspect ratio */

#define INITXMIN -0.1
#define INITXMAX 0.1	/* x interval for initial condition */
#define INITYMIN 1.1
#define INITYMAX 1.2	/* y interval for initial condition */

#define BCXMIN -2.0
#define BCXMAX 2.0	/* x interval for boundary condition */
#define BCYMIN -1.125
#define BCYMAX 1.3	/* y interval for boundary condition */

#define OBSXMIN -2.0
#define OBSXMAX 2.0     /* x interval for motion of obstacle */

#define CIRCLE_PATTERN 1   /* pattern of circles, see list in global_ljones.c */

#define ADD_FIXED_OBSTACLES 0   /* set to 1 do add fixed circular obstacles */
#define OBSTACLE_PATTERN 3  /* pattern of obstacles, see list in global_ljones.c */

#define ADD_FIXED_SEGMENTS 1    /* set to 1 to add fixed segments as obstacles */
#define SEGMENT_PATTERN 1     /* pattern of repelling segments, see list in global_ljones.c */
#define ROCKET_SHAPE 2        /* shape of rocket combustion chamber, see list in global_ljones.c */
#define ROCKET_SHAPE_B 2      /* shape of second rocket */
#define NOZZLE_SHAPE 1        /* shape of nozzle, see list in global_ljones.c */
#define NOZZLE_SHAPE_B 0      /* shape of nozzle for second rocket, see list in global_ljones.c */

#define TWO_TYPES 0         /* set to 1 to have two types of particles */
#define TPYE_PROPORTION 0.5 /* proportion of particles of first type */
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
#define PDISC_DISTANCE 3.0  /* minimal distance in Poisson disc process, controls density of particles */
#define PDISC_CANDIDATES 100 /* number of candidates in construction of Poisson disc process */
#define RANDOM_POLY_ANGLE 0 /* set to 1 to randomize angle of polygons */

#define LAMBDA 0.3	    /* parameter controlling the dimensions of domain */
#define MU 0.009 	    /* parameter controlling radius of particles */
#define MU_B 0.009          /* parameter controlling radius of particles of second type */
#define NPOLY 3             /* number of sides of polygon */
#define APOLY 0.66           /* angle by which to turn polygon, in units of Pi/2 */ 
#define MDEPTH 4            /* depth of computation of Menger gasket */
#define MRATIO 3            /* ratio defining Menger gasket */
#define MANDELLEVEL 1000    /* iteration level for Mandelbrot set */
#define MANDELLIMIT 10.0    /* limit value for approximation of Mandelbrot set */
#define FOCI 1              /* set to 1 to draw focal points of ellipse */
#define NGRIDX 1            /* number of grid point for grid of disks */
#define NGRIDY 1           /* number of grid point for grid of disks */
#define EHRENFEST_RADIUS 0.9    /* radius of container for Ehrenfest urn configuration */
#define EHRENFEST_WIDTH 0.035     /* width of tube for Ehrenfest urn configuration */
#define TWO_CIRCLES_RADIUS_RATIO 0.8    /* ratio of radii for S_TWO_CIRCLES_EXT segment configuration */
#define DAM_WIDTH 0.05       /* width of dam for S_DAM segment configuration */

#define X_SHOOTER -0.2
#define Y_SHOOTER -0.6
#define X_TARGET 0.4
#define Y_TARGET 0.7        /* shooter and target positions in laser fight */

/* Parameters for length and speed of simulation */

#define NSTEPS 4000       /* number of frames of movie */
#define NVID 100         /* number of iterations between images displayed on screen */
#define NSEG 250         /* number of segments of boundary */
#define INITIAL_TIME 0   /* time after which to start saving frames */
#define OBSTACLE_INITIAL_TIME 350     /* time after which to start moving obstacle */
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

#define PLOT 8
#define PLOT_B 8        /* plot type for second movie */

#define DRAW_BONDS 0    /* set to 1 to draw bonds between neighbours */
#define COLOR_BONDS 1   /* set to 1 to color bonds according to length */
#define FILL_TRIANGLES 0    /* set to 1 to fill triangles between neighbours */
#define ALTITUDE_LINES 0    /* set to 1 to add horizontal lines to show altitude */
#define COLOR_SEG_GROUPS 1  /* set to 1 to collor segment groups differently */

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

#define ENERGY_HUE_MIN 330.0        /* color of original particle */
#define ENERGY_HUE_MAX 50.0         /* color of saturated particle */
#define PARTICLE_HUE_MIN 359.0      /* color of original particle */
#define PARTICLE_HUE_MAX 0.0        /* color of saturated particle */
#define PARTICLE_EMAX 3.0e2         /* energy of particle with hottest color */
#define HUE_TYPE0 70.0     /* hue of particles of type 0 */
#define HUE_TYPE1 280.0      /* hue of particles of type 1 */
#define HUE_TYPE2 .0      /* hue of particles of type 2 */
#define HUE_TYPE3 210.0     /* hue of particles of type 3 */

#define RANDOM_RADIUS 0     /* set to 1 for random circle radius */
#define DT_PARTICLE 3.0e-6    /* time step for particle displacement */
#define KREPEL 12.0          /* constant in repelling force between particles */
#define EQUILIBRIUM_DIST 4.0    /* Lennard-Jones equilibrium distance */
#define EQUILIBRIUM_DIST_B 2.5  /* Lennard-Jones equilibrium distance for second type of particle */
#define REPEL_RADIUS 20.0    /* radius in which repelling force acts (in units of particle radius) */
#define DAMPING 20.0          /* damping coefficient of particles */
#define INITIAL_DAMPING 40.0  /* damping coefficient of particles during initial phase */
#define PARTICLE_MASS 20.0    /* mass of particle of radius MU */
#define PARTICLE_MASS_B 1.0   /* mass of particle of radius MU_B */
#define PARTICLE_INERTIA_MOMENT 0.2     /* moment of inertia of particle */
#define PARTICLE_INERTIA_MOMENT_B 0.02     /* moment of inertia of second type of particle */
#define V_INITIAL 0.0        /* initial velocity range */
#define OMEGA_INITIAL 10.0        /* initial angular velocity range */

#define THERMOSTAT 0        /* set to 1 to switch on thermostat */
#define VARY_THERMOSTAT 0   /* set to 1 for time-dependent thermostat schedule */
#define SIGMA 5.0           /* noise intensity in thermostat */
#define BETA 0.02          /* initial inverse temperature */
#define MU_XI 0.01           /* friction constant in thermostat */
#define KSPRING_BOUNDARY 1.0e7    /* confining harmonic potential outside simulation region */
#define KSPRING_OBSTACLE 1.0e11    /* harmonic potential of obstacles */
#define NBH_DIST_FACTOR 7.5       /* radius in which to count neighbours */
#define GRAVITY 1500.0            /* gravity acting on all particles */
#define GRAVITY_X 0.0        /* horizontal gravity acting on all particles */
#define INCREASE_GRAVITY 0     /* set to 1 to increase gravity during the simulation */
#define GRAVITY_SCHEDULE 2     /* type of gravity schedule, see list in global_ljones.c */
#define GRAVITY_FACTOR 100.0    /* factor by which to increase gravity */
#define GRAVITY_INITIAL_TIME 200    /* time at start of simulation with constant gravity */
#define GRAVITY_RESTORE_TIME 700    /* time at end of simulation with gravity restored to initial value */

#define ROTATION 0           /* set to 1 to include rotation of particles */
#define COUPLE_ANGLE_TO_THERMOSTAT 0    /* set to 1 to couple angular degrees of freedom to thermostat */
#define DIMENSION_FACTOR 1.0  /* scaling factor taking into account number of degrees of freedom */  
#define KTORQUE 100.0          /* force constant in angular dynamics */
#define KTORQUE_B 10.0        /* force constant in angular dynamics */
#define KTORQUE_DIFF 150.0    /* force constant in angular dynamics for different particles */
#define DRAW_SPIN 0           /* set to 1 to draw spin vectors of particles */
#define DRAW_SPIN_B 0         /* set to 1 to draw spin vectors of particles */
#define DRAW_CROSS 1          /* set to 1 to draw cross on particles of second type */
#define SPIN_RANGE 7.0       /* range of spin-spin interaction */
#define SPIN_RANGE_B 5.0     /* range of spin-spin interaction for second type of particle */
#define QUADRUPOLE_RATIO 0.6  /* anisotropy in quadrupole potential */ 

#define INCREASE_BETA 0  /* set to 1 to increase BETA during simulation */
#define BETA_FACTOR 0.025   /* factor by which to change BETA during simulation */
#define N_TOSCILLATIONS 1.5   /* number of temperature oscillations in BETA schedule */
#define NO_OSCILLATION 1        /* set to 1 to have exponential BETA change only */
#define MIDDLE_CONSTANT_PHASE 1600   /* final phase in which temperature is constant */
#define FINAL_DECREASE_PHASE 1500    /* final phase in which temperature decreases */ 
#define FINAL_CONSTANT_PHASE -1     /* final phase in which temperature is constant */

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
#define PARTIAL_THERMO_REGION 1     /* region for partial thermostat coupling (see list in global_ljones.c) */
#define PARTIAL_THERMO_SHIFT 0.2    /* distance from obstacle at the right of which particles are coupled to thermostat */
#define PARTIAL_THERMO_WIDTH 0.5    /* vertical size of partial thermostat coupling */
#define PARTIAL_THERMO_HEIGHT 0.2   /* vertical size of partial thermostat coupling */

#define INCREASE_KREPEL 0   /* set to 1 to increase KREPEL during simulation */
#define KREPEL_FACTOR 1000.0   /* factor by which to change KREPEL during simulation */

#define PART_AT_BOTTOM 0     /* set to 1 to include "seed" particles at bottom */
#define MASS_PART_BOTTOM 10000.0 /* mass of particles at bottom */
#define NPART_BOTTOM 100     /* number of particles at the bottom */

#define ADD_PARTICLES 1    /* set to 1 to add particles */
#define ADD_TIME 0       /* time at which to add first particle */
#define ADD_PERIOD 2      /* time interval between adding further particles */
#define N_ADD_PARTICLES 20   /* number of particles to add */
#define FINAL_NOADD_PERIOD 700  /* final period where no particles are added */
#define SAFETY_FACTOR 2.0  /* no particles are added at distance less than MU*SAFETY_FACTOR of other particles */

#define TRACER_PARTICLE 0   /* set to 1 to have a tracer particle */
#define N_TRACER_PARTICLES 3    /* number of tracer particles */
#define TRAJECTORY_LENGTH 8000   /* length of recorded trajectory */
#define TRACER_PARTICLE_MASS 4.0    /* relative mass of tracer particle */
#define TRAJECTORY_WIDTH 3      /* width of tracer particle trajectory */

#define ROTATE_BOUNDARY 0           /* set to 1 to rotate the repelling segments */
#define SMOOTH_ROTATION 1           /* set to 1 to update segments at each time step (rather than at each movie frame) */
#define PERIOD_ROTATE_BOUNDARY 1000  /* period of rotating boundary */
#define ROTATE_INITIAL_TIME 0       /* initial time without rotation */
#define ROTATE_FINAL_TIME 100       /* final time without rotation */
#define ROTATE_CHANGE_TIME 0.33     /* relative duration of acceleration/deceleration phases */
#define OMEGAMAX 100.0              /* maximal rotation speed */
#define PRINT_OMEGA 0               /* set to 1 to print angular speed */
#define PRINT_PARTICLE_SPEEDS 0     /* set to 1 to print average speeds/momenta of particles */
#define PRINT_SEGMENTS_SPEEDS 0     /* set to 1 to print velocity of moving segments */

#define MOVE_BOUNDARY 0        /* set to 1 to move repelling segments, due to force from particles */
#define SEGMENTS_MASS 40.0     /* mass of collection of segments */
#define DEACTIVATE_SEGMENT 1    /* set to 1 to deactivate last segment after a certain time */
#define SEGMENT_DEACTIVATION_TIME 200   /* time at which to deactivate last segment */
#define RELEASE_ROCKET_AT_DEACTIVATION 1    /* set to 1 to limit segments velocity before segment release */
#define SEGMENTS_X0 1.0        /* initial position of segments */
#define SEGMENTS_Y0 0.8        /* initial position of segments */
#define SEGMENTS_VX0 0.0       /* initial velocity of segments */
#define SEGMENTS_VY0 0.0      /* initial velocity of segments */
#define DAMP_SEGS_AT_NEGATIVE_Y 0   /* set to 1 to dampen segments when y coordinate is negative */

#define MOVE_SEGMENT_GROUPS 0       /* set to 1 to group segments into moving units */
#define SEGMENT_GROUP_MASS 750.0    /* mass of segment group */
#define SEGMENT_GROUP_I 500.0       /* moment of inertia of segment group */
#define SEGMENT_GROUP_DAMPING 0.0   /* damping of segment groups */
#define GROUP_REPULSION 1           /* set to 1 for groups of segments to repel each other */
#define KSPRING_GROUPS 1.0e11       /* harmonic potential between segment groups */
#define GROUP_WIDTH 0.05            /* interaction width of groups */
#define GROUP_G_REPEL 0             /* set to 1 to add repulsion between centers of mass of groups */
#define GROUP_G_REPEL_RADIUS 1.2    /* radius within which centers of mass of groups repel each other */
#define TRACK_SEGMENT_GROUPS 0      /* set to 1 for view to track group of segments */
#define TRACK_X_PADDING 2.0         /* distance from x boundary where tracking starts */

#define POSITION_DEPENDENT_TYPE 1   /* set to 1 to make particle type depend on initial position */
#define POSITION_Y_DEPENDENCE 1     /* set to 1 for the separation between particles to be horizontal */
#define PRINT_ENTROPY 0     /* set to 1 to compute entropy */

#define REACTION_DIFFUSION 0    /* set to 1 to simulate a chemical reaction (particles may change type) */
#define RD_TYPES 3              /* number of types in reaction-diffusion equation */
#define REACTION_PROB 0.0045      /* probability controlling reaction term */ 

#define PRINT_PARTICLE_NUMBER 0     /* set to 1 to print total number of particles */
#define PRINT_LEFT 0        /* set to 1 to print certain parameters at the top left instead of right */
#define PLOT_SPEEDS 0       /* set to 1 to add a plot of obstacle speeds (e.g. for rockets) */
#define PLOT_TRAJECTORIES 0     /* set to 1 to add a plot of obstacle trajectories (e.g. for rockets) */
#define VMAX_PLOT_SPEEDS 0.6    /* vertical scale of plot of obstacle speeds */

#define EHRENFEST_COPY 0    /* set to 1 to add equal number of larger particles (for Ehrenfest model) */

#define LID_MASS 1000.0     /* mass of lid for BC_RECTANGLE_LID b.c. */
#define LID_WIDTH 0.1       /* width of lid for BC_RECTANGLE_LID b.c. */
#define WALL_MASS 2000.0    /* mass of wall for BC_RECTANGLE_WALL b.c. */
#define WALL_FRICTION 0.0   /* friction on wall for BC_RECTANGLE_WALL b.c. */
#define WALL_WIDTH 0.1      /* width of wall for BC_RECTANGLE_WALL b.c. */
#define WALL_VMAX 100.0     /* max speed of wall */
#define WALL_TIME 0         /* time during which to keep wall */

#define NXMAZE 10      /* width of maze */
#define NYMAZE 10      /* height of maze */
#define MAZE_MAX_NGBH 4     /* max number of neighbours of maze cell */
#define RAND_SHIFT 200      /* seed of random number generator */
#define MAZE_XSHIFT 0.0     /* horizontal shift of maze */

#define FLOOR_FORCE 1      /* set to 1 to limit force on particle to FMAX */
#define FMAX 1.0e12         /* maximal force */
#define FLOOR_OMEGA 0      /* set to 1 to limit particle momentum to PMAX */
#define PMAX 1000.0        /* maximal force */

#define HASHX 150   /* size of hashgrid in x direction */
#define HASHY 75    /* size of hashgrid in y direction */
#define HASHMAX 100  /* maximal number of particles per hashgrid cell */
#define HASHGRID_PADDING 0.1    /* padding of hashgrid outside simulation window */

#define DRAW_COLOR_SCHEME 0     /* set to 1 to plot the color scheme */
#define COLORBAR_RANGE 8.0    /* scale of color scheme bar */
#define COLORBAR_RANGE_B 12.0    /* scale of color scheme bar for 2nd part */
#define ROTATE_COLOR_SCHEME 0   /* set to 1 to draw color scheme horizontally */

```

### 27 December 22 - Turbulence (possibly numerical) in a flow of increasing speed around an obstacle ###

**Program:** `rde.c` 

**Initial condition in function `animation()`:** `init_laminar_flow(flow_speed_schedule(0), LAMINAR_FLOW_MODULATION, 0.5, -0.1, phi, xy_in);`

```
#define MOVIE 1         /* set to 1 to generate movie */
#define DOUBLE_MOVIE 1  /* set to 1 to produce movies for wave height and energy simultaneously */
#define SAVE_MEMORY 0           /* set to 1 to save memory when writing tiff images */

/* General geometrical parameters */

#define WINWIDTH 	1920  /* window width */
#define WINHEIGHT 	1000  /* window height */
#define NX 960          /* number of grid points on x axis */
#define NY 500          /* number of grid points on y axis */

#define XMIN -1.0
#define XMAX 3.0	/* x interval */
#define YMIN -1.041666667
#define YMAX 1.041666667	/* y interval for 9/16 aspect ratio */

/* Choice of simulated equation */

#define RDE_EQUATION 6  /* choice of reaction term, see list in global_3d.c */
#define NFIELDS 3       /* number of fields in reaction-diffusion equation */
#define NLAPLACIANS 1   /* number of fields for which to compute Laplacian */

#define ADD_POTENTIAL 0 /* set to 1 to add a potential (for Schrodinger equation) */
#define ADD_MAGNETIC_FIELD 0    /* set to 1 to add a magnetic field (for Schrodinger equation) - then set POTENTIAL 1 */
#define POTENTIAL 7     /* type of potential or vector potential, see list in global_3d.c  */

#define ANTISYMMETRIZE_WAVE_FCT 0   /* set tot 1 to make wave function antisymmetric */

#define JULIA_SCALE 0.5 /* scaling for Julia sets */

/* Choice of the billiard table */

#define B_DOMAIN 199          /* choice of domain shape, see list in global_pdes.c  */

#define CIRCLE_PATTERN 99    /* pattern of circles, see list in global_pdes.c */

#define P_PERCOL 0.25       /* probability of having a circle in C_RAND_PERCOL arrangement */
#define NPOISSON 300        /* number of points for Poisson C_RAND_POISSON arrangement */
#define RANDOM_POLY_ANGLE 0 /* set to 1 to randomize angle of polygons */

#define LAMBDA 0.3	    /* parameter controlling the dimensions of domain */
#define MU 0.5	            /* parameter controlling the dimensions of domain */
#define NPOLY 5             /* number of sides of polygon */
#define APOLY 2.0          /* angle by which to turn polygon, in units of Pi/2 */
#define MDEPTH 7            /* depth of computation of Menger gasket */
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

#define DT 0.0000012

#define VISCOSITY 2.0

#define RPSA 0.75         /* parameter in Rock-Paper-Scissors-type interaction */
#define RPSLZB 0.75       /* second parameter in Rock-Paper-Scissors-Lizard-Spock type interaction */

#define EPSILON 0.8     /* time scale separation */
#define DELTA 0.1       /* time scale separation */
#define FHNA 1.0        /* parameter in FHN equation */
#define FHNC -0.01      /* parameter in FHN equation */
#define K_HARMONIC 1.0  /* spring constant of harmonic potential */
#define K_COULOMB 0.5   /* constant in Coulomb potential */
#define V_MAZE 0.4      /* potential in walls of maze */
#define BZQ 0.0008      /* parameter in BZ equation */
#define BZF 1.2         /* parameter in BZ equation */
#define B_FIELD 10.0    /* magnetic field */
#define AB_RADIUS 0.2   /* radius of region with magnetic field for Aharonov-Bohm effect */
#define K_EULER 50.0    /* constant in stream function integration of Euler equation */

#define SMOOTHEN_VORTICITY 1    /* set to 1 to smoothen vorticity field in Euler equation */
#define SMOOTHEN_PERIOD 10      /* period between smoothenings */
#define SMOOTH_FACTOR 0.05      /* factor by which to smoothen */

#define ADD_TRACERS 1    /* set to 1 to add tracer particles (for Euler equations) */
#define N_TRACERS 1000    /* number of tracer particles */

#define T_OUT 2.0       /* outside temperature */
#define T_IN 0.0        /* inside temperature */
#define SPEED 0.0       /* speed of drift to the right */

#define ADD_NOISE 0     /* set to 1 to add noise, set to 2 to add noise in right half */
#define NOISE_INTENSITY 0.005      /* noise intensity */
#define CHANGE_NOISE 1      /* set to 1 to increase noise intensity */
#define NOISE_FACTOR 40.0   /* factor by which to increase noise intensity */
#define NOISE_INITIAL_TIME 100  /* initial time during which noise remains constant */

#define CHANGE_VISCOSITY 0      /* set to 1 to change the viscosity in the course of the simulation */
#define ADJUST_INTSTEP 0       /* set to 1 to decrease integration step when viscosity increases */
#define VISCOSITY_INITIAL_TIME 10  /* initial time during which viscosity remains constant */
#define VISCOSITY_FACTOR 100.0   /* factor by which to change viscosity */
#define VISCOSITY_MAX 2.0        /* max value of viscosity beyond which NVID is increased and integration step is decrase, 
                                    for numerical stability */
                                        
#define CHANGE_RPSLZB 0         /* set to 1 to change second parameter in Rock-Paper-Scissors-Lizard-Spock equation */
#define RPSLZB_CHANGE 0.75      /* factor by which to rpslzb parameter */
#define RPSLZB_INITIAL_TIME 0   /* initial time during which rpslzb remains constant */
#define RPSLZB_FINAL_TIME 500   /* final time during which rpslzb remains constant */
                                      
#define CHANGE_FLOW_SPEED 1     /* set to 1 to change speed of laminar flow */

/* Boundary conditions, see list in global_pdes.c  */

#define B_COND 1

#define IN_OUT_FLOW_BC 1         /* set to 1 to add in-flow/out-flow boundary conditions (for Euler equation) */
#define IN_OUT_FLOW_MIN_AMP 0.3  /* amplitude of in-flow/out-flow boundary conditions (for Euler equation) */
#define IN_OUT_FLOW_AMP 1.4      /* amplitude of in-flow/out-flow boundary conditions (for Euler equation) */
#define LAMINAR_FLOW_MODULATION 0.5     /* asymmetry of laminar flow */

#define EULER_GRADIENT_YSHIFT 0.0    /* y-shift in computation of gradient in Euler equation */

/* Parameters for length and speed of simulation */

#define NSTEPS 1600           /* number of frames of movie */
#define NVID 120          /* number of iterations between images displayed on screen */
#define ACCELERATION_FACTOR 1.0 /* factor by which to increase NVID in course of simulation */
#define DT_ACCELERATION_FACTOR 1.0 /* factor by which to increase time step in course of simulation  */
#define MAX_DT 0.024     /* maximal value of integration step */
#define NSEG 100         /* number of segments of boundary */
#define BOUNDARY_WIDTH 2    /* width of billiard boundary */

#define PAUSE 100       /* number of frames after which to pause */
#define PSLEEP 2         /* sleep time during pause */
#define SLEEP1  2        /* initial sleeping time */
#define SLEEP2  1        /* final sleeping time */
#define INITIAL_TIME 0  /* initial still time */
#define MID_FRAMES 50    /* number of still frames between parts of two-part movie */
#define END_FRAMES 50    /* number of still frames at end of movie */
#define FADE 1           /* set to 1 to fade at end of movie */

/* Visualisation */

#define PLOT_3D 0    /* controls whether plot is 2D or 3D */

#define ROTATE_VIEW 0       /* set to 1 to rotate position of observer */
#define ROTATE_ANGLE 360.0  /* total angle of rotation during simulation */

#define DRAW_PERIODICISED 1     /* set to 1 to repeat wave periodically in x and y directions */

/* Plot type - color scheme */

#define CPLOT 50
#define CPLOT_B 51

/* Plot type - height of 3D plot */

#define ZPLOT 52     /* z coordinate in 3D plot */
#define ZPLOT_B 51    /* z coordinate in second 3D plot */

#define AMPLITUDE_HIGH_RES 1    /* set to 1 to increase resolution of P_3D_AMPLITUDE plot */
#define SHADE_3D 1              /* set to 1 to change luminosity according to normal vector */
#define NON_DIRICHLET_BC 0      /* set to 1 to draw only facets in domain, if field is not zero on boundary */
#define WRAP_ANGLE 1            /* experimental: wrap angle to [0, 2Pi) for interpolation in angle schemes */
#define FADE_IN_OBSTACLE 0      /* set to 1 to fade color inside obstacles */
#define DRAW_OUTSIDE_GRAY 0     /* experimental - draw outside of billiard in gray */
#define ADD_POTENTIAL_TO_Z 1    /* set to 1 to add the external potential to z-coordinate of plot */
#define ADD_POT_CONSTANT 0.35   /* constant in front of added potential */

#define PLOT_SCALE_ENERGY 0.05      /* vertical scaling in energy plot */

#define PRINT_TIME 0        /* set to 1 to print running time */
#define PRINT_VISCOSITY 0   /* set to 1 to print viscosity */
#define PRINT_RPSLZB 0      /* set to 1 to print rpslzb parameter */
#define PRINT_PROBABILITIES 0   /* set to 1 to print probabilities (for Ehrenfest urn configuration) */
#define PRINT_NOISE 0       /* set to 1 to print noise intensity */
#define PRINT_FLOW_SPEED 1      /* set to 1 to print speed of flow */

#define DRAW_FIELD_LINES 0  /* set to 1 to draw field lines */
#define FIELD_LINE_WIDTH 1  /* width of field lines */
#define N_FIELD_LINES 120   /* number of field lines */
#define FIELD_LINE_FACTOR 120 /* factor controlling precision when computing origin of field lines */
#define DRAW_BILLIARD 1     /* set to 1 to draw boundary */
#define DRAW_BILLIARD_FRONT 1     /* set to 1 to draw boundary */
#define FILL_BILLIARD_COMPLEMENT 1  /* set to 1 to fill complement of billiard (for certain shapes only) */

/* 3D representation */

#define REPRESENTATION_3D 1     /* choice of 3D representation */ 

#define REP_AXO_3D 0        /* linear projection (axonometry) */
#define REP_PROJ_3D 1       /* projection on plane orthogonal to observer line of sight */

/* Color schemes, see list in global_pdes.c  */

#define COLOR_PALETTE 10       /* Color palette, see list in global_pdes.c  */
#define COLOR_PALETTE_B 18     /* Color palette, see list in global_pdes.c  */

#define BLACK 1          /* black background */

#define COLOR_SCHEME 3   /* choice of color scheme */

#define COLOR_PHASE_SHIFT 0.0   /* phase shift of color scheme, in units of Pi */

#define SCALE 0          /* set to 1 to adjust color scheme to variance of field */
#define SLOPE 1.0        /* sensitivity of color on wave amplitude */
#define VSCALE_AMPLITUDE 1.0      /* additional scaling factor for color scheme P_3D_AMPLITUDE */
#define ATTENUATION 0.0  /* exponential attenuation coefficient of contrast with time */
#define CURL_SCALE 0.000015   /* scaling factor for curl representation */
#define RESCALE_COLOR_IN_CENTER 0   /* set to 1 to decrease color intentiy in the center (for wave escaping ring) */
#define SLOPE_SCHROD_LUM 50.0       /* sensitivity of luminosity on module, for color scheme Z_ARGUMENT */
#define MIN_SCHROD_LUM 0.2       /* minimal luminosity in color scheme Z_ARGUMENT*/
#define VSCALE_PRESSURE 0.5      /* additional scaling factor for color scheme Z_EULER_PRESSURE */
#define PRESSURE_SHIFT 10.0        /* shift for color scheme Z_EULER_PRESSURE */
#define PRESSURE_LOG_SHIFT -2.5     /* shift for color scheme Z_EULER_PRESSURE */

#define COLORHUE 260     /* initial hue of water color for scheme C_LUM */
#define COLORDRIFT 0.0   /* how much the color hue drifts during the whole simulation */
#define LUMMEAN 0.5      /* amplitude of luminosity variation for scheme C_LUM */
#define LUMAMP 0.3       /* amplitude of luminosity variation for scheme C_LUM */
#define HUEMEAN 359.0    /* mean value of hue for color scheme C_HUE */
#define HUEAMP -359.0    /* amplitude of variation of hue for color scheme C_HUE */
#define E_SCALE 100.0    /* scaling factor for energy representation */
#define FLUX_SCALE 100.0 /* scaling factor for energy representation */
#define LOG_SCALE 0.5    /* scaling factor for energy log representation */
#define LOG_SHIFT 1.0   
#define LOG_MIN 1.0e-3   /* floor value for log vorticity plot */

#define NXMAZE 7      /* width of maze */
#define NYMAZE 7      /* height of maze */
#define MAZE_MAX_NGBH 4     /* max number of neighbours of maze cell */
#define RAND_SHIFT 0        /* seed of random number generator */
#define MAZE_XSHIFT 0.0     /* horizontal shift of maze */

#define DRAW_COLOR_SCHEME 1     /* set to 1 to plot the color scheme */
#define COLORBAR_RANGE 2.0      /* scale of color scheme bar */
#define COLORBAR_RANGE_B 1.7   /* scale of color scheme bar for 2nd part */
#define ROTATE_COLOR_SCHEME 0   /* set to 1 to draw color scheme horizontally */

/* only for compatibility with wave_common.c */
#define TWOSPEEDS 0          /* set to 1 to replace hardcore boundary by medium with different speed */
#define OMEGA 0.005        /* frequency of periodic excitation */
#define COURANT 0.08       /* Courant number */
#define COURANTB 0.03      /* Courant number in medium B */
#define INITIAL_AMP 0.5         /* amplitude of initial condition */
#define INITIAL_VARIANCE 0.0002  /* variance of initial condition */
#define INITIAL_WAVELENGTH  0.1  /* wavelength of initial condition */
#define VSCALE_ENERGY 200.0       /* additional scaling factor for color scheme P_3D_ENERGY */
#define PHASE_FACTOR 20.0       /* factor in computation of phase in color scheme P_3D_PHASE */
#define PHASE_SHIFT 0.0      /* shift of phase in color scheme P_3D_PHASE */
#define OSCILLATION_SCHEDULE 0  /* oscillation schedule, see list in global_pdes.c */
#define AMPLITUDE 0.8      /* amplitude of periodic excitation */ 
#define ACHIRP 0.2        /* acceleration coefficient in chirp */
#define DAMPING 0.0        /* damping of periodic excitation */
#define COMPARISON 0        /* set to 1 to compare two different patterns (beta) */
#define B_DOMAIN_B 20       /* second domain shape, for comparisons */
#define CIRCLE_PATTERN_B 0  /* second pattern of circles or polygons */
/* end of constants added only for compatibility with wave_common.c */


double u_3d[2] = {0.75, -0.45};     /* projections of basis vectors for REP_AXO_3D representation */
double v_3d[2] = {-0.75, -0.45};
double w_3d[2] = {0.0, 0.015};
double light[3] = {0.816496581, -0.40824829, 0.40824829};      /* vector of "light" direction for P_3D_ANGLE color scheme */
double observer[3] = {8.0, 8.0, 8.0};    /* location of observer for REP_PROJ_3D representation */ 
int reset_view = 0;         /* switch to reset 3D view parameters (for option ROTATE_VIEW) */

#define Z_SCALING_FACTOR 0.08  /* overall scaling factor of z axis for REP_PROJ_3D representation */
#define XY_SCALING_FACTOR 1.7  /* overall scaling factor for on-screen (x,y) coordinates after projection */
#define ZMAX_FACTOR 1.0        /* max value of z coordinate for REP_PROJ_3D representation */
#define XSHIFT_3D -0.1         /* overall x shift for REP_PROJ_3D representation */
#define YSHIFT_3D 0.1          /* overall y shift for REP_PROJ_3D representation */
#define BORDER_PADDING 0       /* distance from boundary at which to plot points, to avoid boundary effects due to gradient */

/* For debugging purposes only */
#define FLOOR 1         /* set to 1 to limit wave amplitude to VMAX */
#define VMAX 1000.0        /* max value of wave amplitude */
#define TEST_GRADIENT 0 /* print norm squared of gradient */

```

### 26 December 22 - Shining a rotating laser beam into a 16 by 16 maze with rounded corners ###

**Program:** `particle_trajectory.c` 

**Initial condition in function `animation()`:** 
```
    x = -1.65;
    y = 0.03;
    alpha = angle_schedule(time);        
    period = compute_trajectories_xy(x, y, alpha, alpha + DPI, configs, trajectory, traj_length, &xmax);
```

```
#define MOVIE 1         /* set to 1 to generate movie */

#define WINWIDTH 	1280  /* window width */
#define WINHEIGHT 	720   /* window height */
#define SAVE_MEMORY 1   /* set to 1 to save memory when writing tiff images */

#define XMIN -2.0
#define XMAX 2.0	/* x interval */
#define YMIN -1.125
#define YMAX 1.125	/* y interval for 9/16 aspect ratio */

#define SCALING_FACTOR 1.0       /* scaling factor of drawing, needed for flower billiards, otherwise set to 1.0 */

/* Choice of the billiard table, see global_particles.c */

#define B_DOMAIN 31     /* choice of domain shape */

#define CIRCLE_PATTERN 1   /* pattern of circles */
#define POLYLINE_PATTERN 10  /* pattern of polyline */

#define ABSORBING_CIRCLES 1 /* set to 1 for circular scatterers to be absorbing */

#define NMAXCIRCLES 50000        /* total number of circles (must be at least NCX*NCY for square grid) */
#define NMAXPOLY 50000        /* total number of sides of polygonal line */   
#define NCX 9             /* number of circles in x direction */
#define NCY 20            /* number of circles in y direction */
#define NPOISSON 500        /* number of points for Poisson C_RAND_POISSON arrangement */
#define NGOLDENSPIRAL 2000  /* max number of points for C_GOLDEN_SPIRAL arrandement */
#define SDEPTH 2            /* Sierpinski gastket depth */

#define LAMBDA -1.1	/* parameter controlling shape of domain */
#define MU 0.1          /* second parameter controlling shape of billiard */
#define FOCI 1          /* set to 1 to draw focal points of ellipse */
#define NPOLY 3             /* number of sides of polygon */
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
#define TRAJ_LENGTH 10000 /* length of trajectory */
#define LMAX 0.01       /* minimal segment length triggering resampling */ 
#define DMIN 0.02       /* minimal distance to boundary for triggering resampling */ 
#define CYCLE 1         /* set to 1 for closed curve (start in all directions) */
#define SHOWTRAILS 0    /* set to 1 to keep trails of the particles */
#define SHOWZOOM 0      /* set to 1 to show zoom on specific area */
#define TEST_ACTIVE 1   /* set to 1 to test whether particle is in billiard */
#define PRINT_TRAJECTORY_LENGTH 1   /* set to 1 to print length of trajectory 0 */
#define PRINT_TRAJECTORY_ANGLE 1    /* set to 1 to print angle of trajectory 0 */
#define PRINT_TRAJECTORY_PERIOD 0   /* set to 1 to print period of trajectory 0 */
#define DRAW_LENGTHS_PLOT 1         /* set to 1 to plot trajectory lengths */
#define LENGTHS_LOG_SCALE 1         /* set to 1 to use log scale for plot of lengths */
#define MIN_ANGLE 10.0         /* range of angles of trajectory */
#define MAX_ANGLE 50.0         /* range of angles of trajectory */
#define TEST_INITIAL_COND 1     /* set to 1 to allow only initial conditions that pass a test */

#define SLOW_AT_LONG_TRAJ 0     /* set to 1 to slow down movie for long trajectories */
#define ADD_SUCCESS_GALLERY 1   /* set to 1 to add gallery of successful trajectories at end of movie */
#define SUCCESS_GALLERY_FRAMES 25   /* number of frames per success */

#define NSTEPS 2000      /* number of frames of movie */
#define TIME 2500         /* time between movie frames, for fluidity of real-time simulation */ 
#define DPHI 0.00001     /* integration step */
#define NVID 150         /* number of iterations between images displayed on screen */

/* Decreasing TIME accelerates the animation and the movie                               */
/* For constant speed of movie, TIME*DPHI should be kept constant                        */
/* However, increasing DPHI too much deterioriates quality of simulation                 */
/* NVID tells how often a picture is drawn in the animation, increase it for faster anim */
/* For a good quality movie, take for instance TIME = 400, DPHI = 0.00005, NVID = 100    */

/* Colors and other graphical parameters */

#define COLOR_PALETTE 13     /* Color palette, see list in global_pdes.c  */

#define NCOLORS 500       /* number of colors */
#define COLORSHIFT 0     /* hue of initial color */ 
#define RAINBOW_COLOR 0  /* set to 1 to use different colors for all particles */
#define FLOWER_COLOR 0   /* set to 1 to adapt initial colors to flower billiard (tracks vs core) */
#define NSEG 100         /* number of segments of boundary */
#define LENGTH 0.03       /* length of velocity vectors */
#define BILLIARD_WIDTH 2    /* width of billiard */
#define PARTICLE_WIDTH 2    /* width of particles */
#define FRONT_WIDTH 3       /* width of wave front */
#define GRAPH_HUE 30.0      /* color hue of graph */
#define SUCCESS_HUE 90.0   /* color hue of success */

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

#define NXMAZE 16      /* width of maze */
#define NYMAZE 16     /* height of maze */
#define MAZE_MAX_NGBH 4     /* max number of neighbours of maze cell */
#define RAND_SHIFT 1       /* seed of random number generator */
#define MAZE_XSHIFT -0.5     /* horizontal shift of maze */
#define MAZE_RANDOM_FACTOR 0.1     /* randomization factor for S_MAZE_RANDOM */
#define MAZE_CORNER_RADIUS 0.4     /* radius of tounded corners in maze */

```

### 25 December 22 - A refracting star, viewed in 3D ###

**Program:** `wave_3d.c` 

**Initial condition in function `animation()`:** `init_circular_wave_mod(0.0, 0.0, phi, psi, xy_in);`

```
#define MOVIE 1         /* set to 1 to generate movie */
#define DOUBLE_MOVIE 1  /* set to 1 to produce movies for wave height and energy simultaneously */
#define SAVE_MEMORY 0   /* set to 1 to save memory when writing tiff images */

/* General geometrical parameters */

#define WINWIDTH 	1920  /* window width */
#define WINHEIGHT 	1000  /* window height */

#define HIGHRES 1        /* set to 1 if resolution of grid is double that of displayed image */

#define NX 1440          /* number of grid points on x axis */
#define NY 1440         /* number of grid points on y axis */

#define XMIN -1.5
#define XMAX 1.5	/* x interval  */
#define YMIN -1.5
#define YMAX 1.5	/* y interval  */

#define JULIA_SCALE 0.8 /* scaling for Julia sets */

/* Choice of the billiard table */

#define B_DOMAIN 42        /* choice of domain shape, see list in global_pdes.c */

#define CIRCLE_PATTERN 2   /* pattern of circles or polygons, see list in global_pdes.c */

#define COMPARISON 0        /* set to 1 to compare two different patterns */
#define B_DOMAIN_B 20       /* second domain shape, for comparisons */
#define CIRCLE_PATTERN_B 0  /* second pattern of circles or polygons */

#define VARIABLE_IOR 0      /* set to 1 for a variable index of refraction */
#define IOR 2               /* choice of index of refraction, see list in global_pdes.c */
#define MANDEL_IOR_SCALE -0.05   /* parameter controlling dependence of IoR on Mandelbrot escape speed */

#define P_PERCOL 0.25       /* probability of having a circle in C_RAND_PERCOL arrangement */
#define NPOISSON 300        /* number of points for Poisson C_RAND_POISSON arrangement */
#define RANDOM_POLY_ANGLE 1 /* set to 1 to randomize angle of polygons */

#define LAMBDA 0.7	    /* parameter controlling the dimensions of domain */
#define MU 0.2              /* parameter controlling the dimensions of domain */
#define NPOLY 24            /* number of sides of polygon */
#define APOLY 2.0           /* angle by which to turn polygon, in units of Pi/2 */ 
#define MDEPTH 7            /* depth of computation of Menger gasket */
#define MRATIO 3            /* ratio defining Menger gasket */
#define MANDELLEVEL 1000    /* iteration level for Mandelbrot set */
#define MANDELLIMIT 10.0    /* limit value for approximation of Mandelbrot set */
#define FOCI 1              /* set to 1 to draw focal points of ellipse */
#define NGRIDX 12           /* number of grid point for grid of disks */
#define NGRIDY 16           /* number of grid point for grid of disks */

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
#define COURANT 0.05        /* Courant number */
#define COURANTB 0.1       /* Courant number in medium B */
#define GAMMA 0.0          /* damping factor in wave equation */
#define GAMMAB 0.0            /* damping factor in wave equation */
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
#define OSCILLATING_SOURCE_PERIOD 100   /* period of oscillating source */
#define ALTERNATE_OSCILLATING_SOURCE 1  /* set to 1 to alternate sign of oscillating source */

/* Boundary conditions, see list in global_pdes.c  */

#define B_COND 2

#define PRECOMPUTE_BC 0     /* set to 1 to compute neighbours for Laplacian in advance */

/* Parameters for length and speed of simulation */

#define NSTEPS 2150        /* number of frames of movie */
#define NVID 5            /* number of iterations between images displayed on screen */
#define NSEG 1000         /* number of segments of boundary */
#define INITIAL_TIME 0      /* time after which to start saving frames */
#define BOUNDARY_WIDTH 2    /* width of billiard boundary */
#define PRINT_SPEED 0       /* set to 1 to print speed of moving source */

#define PAUSE 100       /* number of frames after which to pause */
#define PSLEEP 2         /* sleep time during pause */
#define SLEEP1  1        /* initial sleeping time */
#define SLEEP2  1        /* final sleeping time */
#define MID_FRAMES 100    /* number of still frames between parts of two-part movie */
#define END_FRAMES 100   /* number of still frames at end of movie */
#define FADE 1           /* set to 1 to fade at end of movie */
#define ROTATE_VIEW_WHILE_FADE 1    /* set to 1 to keep rotating viewpoint during fade */

/* Parameters of initial condition */

#define INITIAL_AMP 0.25         /* amplitude of initial condition */
#define INITIAL_VARIANCE 0.0005  /* variance of initial condition */
#define INITIAL_WAVELENGTH  0.02  /* wavelength of initial condition */


/* Plot type, see list in global_pdes.c  */

#define ZPLOT 103     /* wave height */
#define CPLOT 103     /* color scheme */

#define ZPLOT_B 107 
#define CPLOT_B 107        /* plot type for second movie */

#define AMPLITUDE_HIGH_RES 1    /* set to 1 to increase resolution of plot */
#define SHADE_3D 1              /* set to 1 to change luminosity according to normal vector */
#define NON_DIRICHLET_BC 0      /* set to 1 to draw only facets in domain, if field is not zero on boundary */
#define FLOOR_ZCOORD 1          /* set to 1 to draw only facets with z not too negative */
#define DRAW_BILLIARD 0         /* set to 1 to draw boundary */
#define DRAW_BILLIARD_FRONT 0   /* set to 1 to draw front of boundary after drawing wave */
#define DRAW_CONSTRUCTION_LINES 0   /* set to 1 to draw construction lines of certain domains */
#define FADE_IN_OBSTACLE 1      /* set to 1 to fade color inside obstacles */
#define DRAW_OUTSIDE_GRAY 0     /* experimental, draw outside of billiard in gray */

#define PLOT_SCALE_ENERGY 0.01         /* vertical scaling in energy plot */
#define PLOT_SCALE_LOG_ENERGY 0.25       /* vertical scaling in log energy plot */

/* 3D representation */

#define REPRESENTATION_3D 1     /* choice of 3D representation */ 

#define REP_AXO_3D 0        /* linear projection (axonometry) */
#define REP_PROJ_3D 1       /* projection on plane orthogonal to observer line of sight */

#define ROTATE_VIEW 1       /* set to 1 to rotate position of observer */
#define ROTATE_ANGLE 360.0  /* total angle of rotation during simulation */

/* Color schemes */

#define COLOR_PALETTE 14      /* Color palette, see list in global_pdes.c  */
#define COLOR_PALETTE_B 13    /* Color palette, see list in global_pdes.c  */

#define BLACK 1          /* background */

#define COLOR_SCHEME 3   /* choice of color scheme, see list in global_pdes.c  */

#define SCALE 0          /* set to 1 to adjust color scheme to variance of field */
#define SLOPE 1.0        /* sensitivity of color on wave amplitude */
#define VSCALE_AMPLITUDE 1.0   /* additional scaling factor for color scheme P_3D_AMPLITUDE */
#define VSCALE_ENERGY 100.0     /* additional scaling factor for color scheme P_3D_ENERGY */
#define PHASE_FACTOR 20.0       /* factor in computation of phase in color scheme P_3D_PHASE */
#define PHASE_SHIFT 0.0      /* shift of phase in color scheme P_3D_PHASE */
#define ATTENUATION 0.0    /* exponential attenuation coefficient of contrast with time */
#define E_SCALE 100.0      /* scaling factor for energy representation */
#define LOG_SCALE 0.75     /* scaling factor for energy log representation */
#define LOG_SHIFT 0.5      /* shift of colors on log scale */
#define LOG_ENERGY_FLOOR -10.0    /* floor value for log of (total) energy */
#define LOG_MEAN_ENERGY_SHIFT 1.0   /* additional shift for log of mean energy */
#define FLUX_SCALE 1.0e4    /* scaling factor for enegy flux represtnation */
#define RESCALE_COLOR_IN_CENTER 0   /* set to 1 to decrease color intentiy in the center (for wave escaping ring) */

#define COLORHUE 260     /* initial hue of water color for scheme C_LUM */
#define COLORDRIFT 0.0   /* how much the color hue drifts during the whole simulation */
#define LUMMEAN 0.5      /* amplitude of luminosity variation for scheme C_LUM */
#define LUMAMP 0.3       /* amplitude of luminosity variation for scheme C_LUM */
#define HUEMEAN 240.0    /* mean value of hue for color scheme C_HUE */
#define HUEAMP -200.0    /* amplitude of variation of hue for color scheme C_HUE */

#define NXMAZE 7      /* width of maze */
#define NYMAZE 7      /* height of maze */
#define MAZE_MAX_NGBH 4     /* max number of neighbours of maze cell */
#define RAND_SHIFT 24        /* seed of random number generator */
#define MAZE_XSHIFT 0.0     /* horizontal shift of maze */

#define DRAW_COLOR_SCHEME 0       /* set to 1 to plot the color scheme */
#define COLORBAR_RANGE 3.5      /* scale of color scheme bar */
#define COLORBAR_RANGE_B 3.0    /* scale of color scheme bar for 2nd part */
#define ROTATE_COLOR_SCHEME 0     /* set to 1 to draw color scheme horizontally */

#define SAVE_TIME_SERIES 0      /* set to 1 to save wave time series at a point */

/* for compatibility with sub_wave and sub_maze */
#define ADD_POTENTIAL 0
#define POTENTIAL 0
/* end of constants only used by sub_wave and sub_maze */


/* For debugging purposes only */
#define FLOOR 1         /* set to 1 to limit wave amplitude to VMAX */
#define VMAX 10.0       /* max value of wave amplitude */

/* Parameters controlling 3D projection */

double u_3d[2] = {0.75, -0.45};     /* projections of basis vectors for REP_AXO_3D representation */
double v_3d[2] = {-0.75, -0.45};
double w_3d[2] = {0.0, 0.015};
double light[3] = {0.816496581, -0.40824829, 0.40824829};      /* vector of "light" direction for P_3D_ANGLE color scheme */
double observer[3] = {8.0, 8.0, 12.0};    /* location of observer for REP_PROJ_3D representation */ 
int reset_view = 0;         /* switch to reset 3D view parameters (for option ROTATE_VIEW) */

#define Z_SCALING_FACTOR 0.6    /* overall scaling factor of z axis for REP_PROJ_3D representation */
#define XY_SCALING_FACTOR 2.4   /* overall scaling factor for on-screen (x,y) coordinates after projection */
#define ZMAX_FACTOR 1.0         /* max value of z coordinate for REP_PROJ_3D representation */
#define XSHIFT_3D -0.1          /* overall x shift for REP_PROJ_3D representation */
#define YSHIFT_3D 0.1           /* overall y shift for REP_PROJ_3D representation */

```

### 24 December 22 - A Christmas sandpile ###

**Program:** `lennardjones.c`, postprocessed with `ljones_movie.c`

```
#define MOVIE 0         /* set to 1 to generate movie */
#define DOUBLE_MOVIE 0  /* set to 1 to produce movies for wave height and energy simultaneously */
#define SAVE_MEMORY 1   /* set to 1 to save memory while saving frames */

#define TIME_LAPSE 1     /* set to 1 to add a time-lapse movie at the end */
                         /* so far incompatible with double movie */
#define TIME_LAPSE_FACTOR 3    /* factor of time-lapse movie */
#define TIME_LAPSE_FIRST 1  /* set to 1 to show time-lapse version first */

#define SAVE_TIME_SERIES 1  /* set to 1 to save time series of particle positions */

/* General geometrical parameters */

#define WINWIDTH 	1280  /* window width */
#define WINHEIGHT 	720   /* window height */

#define XMIN -2.0
#define XMAX 2.0	/* x interval */
#define YMIN -1.125
#define YMAX 1.125	/* y interval for 9/16 aspect ratio */

#define INITXMIN -0.1
#define INITXMAX 0.1	/* x interval for initial condition */
#define INITYMIN 1.1
#define INITYMAX 1.2	/* y interval for initial condition */

#define BCXMIN -2.0
#define BCXMAX 2.0	/* x interval for boundary condition */
#define BCYMIN -1.125
#define BCYMAX 1.3	/* y interval for boundary condition */

#define OBSXMIN -2.0
#define OBSXMAX 2.0     /* x interval for motion of obstacle */

#define CIRCLE_PATTERN 1   /* pattern of circles, see list in global_ljones.c */

#define ADD_FIXED_OBSTACLES 0   /* set to 1 do add fixed circular obstacles */
#define OBSTACLE_PATTERN 3  /* pattern of obstacles, see list in global_ljones.c */

#define ADD_FIXED_SEGMENTS 1    /* set to 1 to add fixed segments as obstacles */
#define SEGMENT_PATTERN 1     /* pattern of repelling segments, see list in global_ljones.c */
#define ROCKET_SHAPE 2        /* shape of rocket combustion chamber, see list in global_ljones.c */
#define ROCKET_SHAPE_B 2      /* shape of second rocket */
#define NOZZLE_SHAPE 1        /* shape of nozzle, see list in global_ljones.c */
#define NOZZLE_SHAPE_B 0      /* shape of nozzle for second rocket, see list in global_ljones.c */

#define TWO_TYPES 0         /* set to 1 to have two types of particles */
#define TPYE_PROPORTION 0.5 /* proportion of particles of first type */
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
#define PDISC_DISTANCE 3.0  /* minimal distance in Poisson disc process, controls density of particles */
#define PDISC_CANDIDATES 100 /* number of candidates in construction of Poisson disc process */
#define RANDOM_POLY_ANGLE 0 /* set to 1 to randomize angle of polygons */

#define LAMBDA 0.3	    /* parameter controlling the dimensions of domain */
#define MU 0.009 	    /* parameter controlling radius of particles */
#define MU_B 0.009          /* parameter controlling radius of particles of second type */
#define NPOLY 3             /* number of sides of polygon */
#define APOLY 0.66           /* angle by which to turn polygon, in units of Pi/2 */ 
#define MDEPTH 4            /* depth of computation of Menger gasket */
#define MRATIO 3            /* ratio defining Menger gasket */
#define MANDELLEVEL 1000    /* iteration level for Mandelbrot set */
#define MANDELLIMIT 10.0    /* limit value for approximation of Mandelbrot set */
#define FOCI 1              /* set to 1 to draw focal points of ellipse */
#define NGRIDX 1            /* number of grid point for grid of disks */
#define NGRIDY 1           /* number of grid point for grid of disks */
#define EHRENFEST_RADIUS 0.9    /* radius of container for Ehrenfest urn configuration */
#define EHRENFEST_WIDTH 0.035     /* width of tube for Ehrenfest urn configuration */
#define TWO_CIRCLES_RADIUS_RATIO 0.8    /* ratio of radii for S_TWO_CIRCLES_EXT segment configuration */
#define DAM_WIDTH 0.05       /* width of dam for S_DAM segment configuration */

#define X_SHOOTER -0.2
#define Y_SHOOTER -0.6
#define X_TARGET 0.4
#define Y_TARGET 0.7        /* shooter and target positions in laser fight */

/* Parameters for length and speed of simulation */

#define NSTEPS 3300      /* number of frames of movie */
// #define NSTEPS 1000       /* number of frames of movie */
#define NVID 100         /* number of iterations between images displayed on screen */
#define NSEG 250         /* number of segments of boundary */
#define INITIAL_TIME 0   /* time after which to start saving frames */
#define OBSTACLE_INITIAL_TIME 350     /* time after which to start moving obstacle */
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
#define PLOT_B 8        /* plot type for second movie */

#define DRAW_BONDS 0    /* set to 1 to draw bonds between neighbours */
#define COLOR_BONDS 1   /* set to 1 to color bonds according to length */
#define FILL_TRIANGLES 0    /* set to 1 to fill triangles between neighbours */
#define ALTITUDE_LINES 0    /* set to 1 to add horizontal lines to show altitude */
#define COLOR_SEG_GROUPS 1  /* set to 1 to collor segment groups differently */

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

#define ENERGY_HUE_MIN 330.0        /* color of original particle */
#define ENERGY_HUE_MAX 50.0         /* color of saturated particle */
#define PARTICLE_HUE_MIN 359.0      /* color of original particle */
#define PARTICLE_HUE_MAX 0.0        /* color of saturated particle */
#define PARTICLE_EMAX 3.0e2         /* energy of particle with hottest color */
#define HUE_TYPE0 70.0     /* hue of particles of type 0 */
#define HUE_TYPE1 280.0      /* hue of particles of type 1 */
#define HUE_TYPE2 .0      /* hue of particles of type 2 */
#define HUE_TYPE3 210.0     /* hue of particles of type 3 */

#define RANDOM_RADIUS 0     /* set to 1 for random circle radius */
#define DT_PARTICLE 3.0e-6    /* time step for particle displacement */
#define KREPEL 12.0          /* constant in repelling force between particles */
#define EQUILIBRIUM_DIST 4.0    /* Lennard-Jones equilibrium distance */
#define EQUILIBRIUM_DIST_B 2.5  /* Lennard-Jones equilibrium distance for second type of particle */
#define REPEL_RADIUS 20.0    /* radius in which repelling force acts (in units of particle radius) */
#define DAMPING 140.0          /* damping coefficient of particles */
#define INITIAL_DAMPING 40.0  /* damping coefficient of particles during initial phase */
#define PARTICLE_MASS 20.0    /* mass of particle of radius MU */
#define PARTICLE_MASS_B 1.0   /* mass of particle of radius MU_B */
#define PARTICLE_INERTIA_MOMENT 0.2     /* moment of inertia of particle */
#define PARTICLE_INERTIA_MOMENT_B 0.02     /* moment of inertia of second type of particle */
#define V_INITIAL 0.0        /* initial velocity range */
#define OMEGA_INITIAL 10.0        /* initial angular velocity range */

#define THERMOSTAT 0        /* set to 1 to switch on thermostat */
#define VARY_THERMOSTAT 0   /* set to 1 for time-dependent thermostat schedule */
#define SIGMA 5.0           /* noise intensity in thermostat */
#define BETA 0.02          /* initial inverse temperature */
#define MU_XI 0.01           /* friction constant in thermostat */
#define KSPRING_BOUNDARY 1.0e7    /* confining harmonic potential outside simulation region */
#define KSPRING_OBSTACLE 1.0e11    /* harmonic potential of obstacles */
#define NBH_DIST_FACTOR 7.5       /* radius in which to count neighbours */
#define GRAVITY 1500.0            /* gravity acting on all particles */
#define GRAVITY_X 0.0        /* horizontal gravity acting on all particles */
#define INCREASE_GRAVITY 0     /* set to 1 to increase gravity during the simulation */
#define GRAVITY_SCHEDULE 2     /* type of gravity schedule, see list in global_ljones.c */
#define GRAVITY_FACTOR 100.0    /* factor by which to increase gravity */
#define GRAVITY_INITIAL_TIME 200    /* time at start of simulation with constant gravity */
#define GRAVITY_RESTORE_TIME 700    /* time at end of simulation with gravity restored to initial value */

#define ROTATION 0           /* set to 1 to include rotation of particles */
#define COUPLE_ANGLE_TO_THERMOSTAT 0    /* set to 1 to couple angular degrees of freedom to thermostat */
#define DIMENSION_FACTOR 1.0  /* scaling factor taking into account number of degrees of freedom */  
#define KTORQUE 100.0          /* force constant in angular dynamics */
#define KTORQUE_B 10.0        /* force constant in angular dynamics */
#define KTORQUE_DIFF 150.0    /* force constant in angular dynamics for different particles */
#define DRAW_SPIN 0           /* set to 1 to draw spin vectors of particles */
#define DRAW_SPIN_B 0         /* set to 1 to draw spin vectors of particles */
#define DRAW_CROSS 1          /* set to 1 to draw cross on particles of second type */
#define SPIN_RANGE 7.0       /* range of spin-spin interaction */
#define SPIN_RANGE_B 5.0     /* range of spin-spin interaction for second type of particle */
#define QUADRUPOLE_RATIO 0.6  /* anisotropy in quadrupole potential */ 

#define INCREASE_BETA 0  /* set to 1 to increase BETA during simulation */
#define BETA_FACTOR 0.025   /* factor by which to change BETA during simulation */
#define N_TOSCILLATIONS 1.5   /* number of temperature oscillations in BETA schedule */
#define NO_OSCILLATION 1        /* set to 1 to have exponential BETA change only */
#define MIDDLE_CONSTANT_PHASE 1600   /* final phase in which temperature is constant */
#define FINAL_DECREASE_PHASE 1500    /* final phase in which temperature decreases */ 
#define FINAL_CONSTANT_PHASE -1     /* final phase in which temperature is constant */

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
#define PARTIAL_THERMO_REGION 1     /* region for partial thermostat coupling (see list in global_ljones.c) */
#define PARTIAL_THERMO_SHIFT 0.2    /* distance from obstacle at the right of which particles are coupled to thermostat */
#define PARTIAL_THERMO_WIDTH 0.5    /* vertical size of partial thermostat coupling */
#define PARTIAL_THERMO_HEIGHT 0.2   /* vertical size of partial thermostat coupling */

#define INCREASE_KREPEL 0   /* set to 1 to increase KREPEL during simulation */
#define KREPEL_FACTOR 1000.0   /* factor by which to change KREPEL during simulation */

#define PART_AT_BOTTOM 0     /* set to 1 to include "seed" particles at bottom */
#define MASS_PART_BOTTOM 10000.0 /* mass of particles at bottom */
#define NPART_BOTTOM 100     /* number of particles at the bottom */

#define ADD_PARTICLES 1    /* set to 1 to add particles */
#define ADD_TIME 0       /* time at which to add first particle */
#define ADD_PERIOD 2      /* time interval between adding further particles */
#define N_ADD_PARTICLES 20   /* number of particles to add */
#define FINAL_NOADD_PERIOD 700  /* final period where no particles are added */
#define SAFETY_FACTOR 2.0  /* no particles are added at distance less than MU*SAFETY_FACTOR of other particles */

#define TRACER_PARTICLE 0   /* set to 1 to have a tracer particle */
#define N_TRACER_PARTICLES 3    /* number of tracer particles */
#define TRAJECTORY_LENGTH 8000   /* length of recorded trajectory */
#define TRACER_PARTICLE_MASS 4.0    /* relative mass of tracer particle */
#define TRAJECTORY_WIDTH 3      /* width of tracer particle trajectory */

#define ROTATE_BOUNDARY 0           /* set to 1 to rotate the repelling segments */
#define SMOOTH_ROTATION 1           /* set to 1 to update segments at each time step (rather than at each movie frame) */
#define PERIOD_ROTATE_BOUNDARY 1000  /* period of rotating boundary */
#define ROTATE_INITIAL_TIME 0       /* initial time without rotation */
#define ROTATE_FINAL_TIME 100       /* final time without rotation */
#define ROTATE_CHANGE_TIME 0.33     /* relative duration of acceleration/deceleration phases */
#define OMEGAMAX 100.0              /* maximal rotation speed */
#define PRINT_OMEGA 0               /* set to 1 to print angular speed */
#define PRINT_PARTICLE_SPEEDS 0     /* set to 1 to print average speeds/momenta of particles */
#define PRINT_SEGMENTS_SPEEDS 0     /* set to 1 to print velocity of moving segments */

#define MOVE_BOUNDARY 0        /* set to 1 to move repelling segments, due to force from particles */
#define SEGMENTS_MASS 40.0     /* mass of collection of segments */
#define DEACTIVATE_SEGMENT 1    /* set to 1 to deactivate last segment after a certain time */
#define SEGMENT_DEACTIVATION_TIME 200   /* time at which to deactivate last segment */
#define RELEASE_ROCKET_AT_DEACTIVATION 1    /* set to 1 to limit segments velocity before segment release */
#define SEGMENTS_X0 1.0        /* initial position of segments */
#define SEGMENTS_Y0 0.8        /* initial position of segments */
#define SEGMENTS_VX0 0.0       /* initial velocity of segments */
#define SEGMENTS_VY0 0.0      /* initial velocity of segments */
#define DAMP_SEGS_AT_NEGATIVE_Y 0   /* set to 1 to dampen segments when y coordinate is negative */

#define MOVE_SEGMENT_GROUPS 0       /* set to 1 to group segments into moving units */
#define SEGMENT_GROUP_MASS 750.0    /* mass of segment group */
#define SEGMENT_GROUP_I 500.0       /* moment of inertia of segment group */
#define SEGMENT_GROUP_DAMPING 0.0   /* damping of segment groups */
#define GROUP_REPULSION 1           /* set to 1 for groups of segments to repel each other */
#define KSPRING_GROUPS 1.0e11       /* harmonic potential between segment groups */
#define GROUP_WIDTH 0.05            /* interaction width of groups */
#define GROUP_G_REPEL 0             /* set to 1 to add repulsion between centers of mass of groups */
#define GROUP_G_REPEL_RADIUS 1.2    /* radius within which centers of mass of groups repel each other */
#define TRACK_SEGMENT_GROUPS 0      /* set to 1 for view to track group of segments */
#define TRACK_X_PADDING 2.0         /* distance from x boundary where tracking starts */

#define POSITION_DEPENDENT_TYPE 1   /* set to 1 to make particle type depend on initial position */
#define POSITION_Y_DEPENDENCE 1     /* set to 1 for the separation between particles to be horizontal */
#define PRINT_ENTROPY 0     /* set to 1 to compute entropy */

#define REACTION_DIFFUSION 0    /* set to 1 to simulate a chemical reaction (particles may change type) */
#define RD_TYPES 3              /* number of types in reaction-diffusion equation */
#define REACTION_PROB 0.0045      /* probability controlling reaction term */ 

#define PRINT_PARTICLE_NUMBER 0     /* set to 1 to print total number of particles */
#define PRINT_LEFT 0        /* set to 1 to print certain parameters at the top left instead of right */
#define PLOT_SPEEDS 0       /* set to 1 to add a plot of obstacle speeds (e.g. for rockets) */
#define PLOT_TRAJECTORIES 0     /* set to 1 to add a plot of obstacle trajectories (e.g. for rockets) */
#define VMAX_PLOT_SPEEDS 0.6    /* vertical scale of plot of obstacle speeds */

#define EHRENFEST_COPY 0    /* set to 1 to add equal number of larger particles (for Ehrenfest model) */

#define LID_MASS 1000.0     /* mass of lid for BC_RECTANGLE_LID b.c. */
#define LID_WIDTH 0.1       /* width of lid for BC_RECTANGLE_LID b.c. */
#define WALL_MASS 2000.0    /* mass of wall for BC_RECTANGLE_WALL b.c. */
#define WALL_FRICTION 0.0   /* friction on wall for BC_RECTANGLE_WALL b.c. */
#define WALL_WIDTH 0.1      /* width of wall for BC_RECTANGLE_WALL b.c. */
#define WALL_VMAX 100.0     /* max speed of wall */
#define WALL_TIME 0         /* time during which to keep wall */

#define NXMAZE 10      /* width of maze */
#define NYMAZE 10      /* height of maze */
#define MAZE_MAX_NGBH 4     /* max number of neighbours of maze cell */
#define RAND_SHIFT 200      /* seed of random number generator */
#define MAZE_XSHIFT 0.0     /* horizontal shift of maze */

#define FLOOR_FORCE 1      /* set to 1 to limit force on particle to FMAX */
#define FMAX 1.0e12         /* maximal force */
#define FLOOR_OMEGA 0      /* set to 1 to limit particle momentum to PMAX */
#define PMAX 1000.0        /* maximal force */

#define HASHX 150   /* size of hashgrid in x direction */
#define HASHY 75    /* size of hashgrid in y direction */
#define HASHMAX 100  /* maximal number of particles per hashgrid cell */
#define HASHGRID_PADDING 0.1    /* padding of hashgrid outside simulation window */

#define DRAW_COLOR_SCHEME 0     /* set to 1 to plot the color scheme */
#define COLORBAR_RANGE 8.0    /* scale of color scheme bar */
#define COLORBAR_RANGE_B 12.0    /* scale of color scheme bar for 2nd part */
#define ROTATE_COLOR_SCHEME 0   /* set to 1 to draw color scheme horizontally */

```

### 23 December 22 - An airfoil in a laminar flow ###

**Program:** `rde.c` 

**Initial condition in function `animation()`:** `init_laminar_flow(IN_OUT_FLOW_AMP, 0.2, 0.5, 0.0, phi, xy_in);`

```
#define MOVIE 1         /* set to 1 to generate movie */
#define DOUBLE_MOVIE 1  /* set to 1 to produce movies for wave height and energy simultaneously */
#define SAVE_MEMORY 0           /* set to 1 to save memory when writing tiff images */

/* General geometrical parameters */

#define WINWIDTH 	1920  /* window width */
#define WINHEIGHT 	1000  /* window height */
#define NX 960          /* number of grid points on x axis */
#define NY 500          /* number of grid points on y axis */

#define XMIN -2.0
#define XMAX 2.0	/* x interval */
#define YMIN -1.041666667
#define YMAX 1.041666667	/* y interval for 9/16 aspect ratio */

/* Choice of simulated equation */

#define RDE_EQUATION 6  /* choice of reaction term, see list in global_3d.c */
#define NFIELDS 3       /* number of fields in reaction-diffusion equation */
#define NLAPLACIANS 1   /* number of fields for which to compute Laplacian */

#define ADD_POTENTIAL 0 /* set to 1 to add a potential (for Schrodinger equation) */
#define ADD_MAGNETIC_FIELD 0    /* set to 1 to add a magnetic field (for Schrodinger equation) - then set POTENTIAL 1 */
#define POTENTIAL 7     /* type of potential or vector potential, see list in global_3d.c  */

#define ANTISYMMETRIZE_WAVE_FCT 0   /* set tot 1 to make wave function antisymmetric */

#define JULIA_SCALE 0.5 /* scaling for Julia sets */

/* Choice of the billiard table */

#define B_DOMAIN 198          /* choice of domain shape, see list in global_pdes.c  */

#define CIRCLE_PATTERN 99    /* pattern of circles, see list in global_pdes.c */

#define P_PERCOL 0.25       /* probability of having a circle in C_RAND_PERCOL arrangement */
#define NPOISSON 300        /* number of points for Poisson C_RAND_POISSON arrangement */
#define RANDOM_POLY_ANGLE 0 /* set to 1 to randomize angle of polygons */

#define LAMBDA 0.7	    /* parameter controlling the dimensions of domain */
#define MU 0.15	            /* parameter controlling the dimensions of domain */
#define NPOLY 5             /* number of sides of polygon */
#define APOLY 2.0          /* angle by which to turn polygon, in units of Pi/2 */
#define MDEPTH 7            /* depth of computation of Menger gasket */
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

#define DT 0.0000012

#define VISCOSITY 2.0

#define RPSA 0.75         /* parameter in Rock-Paper-Scissors-type interaction */
#define RPSLZB 0.75       /* second parameter in Rock-Paper-Scissors-Lizard-Spock type interaction */

#define EPSILON 0.8     /* time scale separation */
#define DELTA 0.1       /* time scale separation */
#define FHNA 1.0        /* parameter in FHN equation */
#define FHNC -0.01      /* parameter in FHN equation */
#define K_HARMONIC 1.0  /* spring constant of harmonic potential */
#define K_COULOMB 0.5   /* constant in Coulomb potential */
#define V_MAZE 0.4      /* potential in walls of maze */
#define BZQ 0.0008      /* parameter in BZ equation */
#define BZF 1.2         /* parameter in BZ equation */
#define B_FIELD 10.0    /* magnetic field */
#define AB_RADIUS 0.2   /* radius of region with magnetic field for Aharonov-Bohm effect */
#define K_EULER 50.0    /* constant in stream function integration of Euler equation */

#define SMOOTHEN_VORTICITY 1    /* set to 1 to smoothen vorticity field in Euler equation */
#define SMOOTHEN_PERIOD 10      /* period between smoothenings */
#define SMOOTH_FACTOR 0.05      /* factor by which to smoothen */

#define ADD_TRACERS 1    /* set to 1 to add tracer particles (for Euler equations) */
#define N_TRACERS 1000    /* number of tracer particles */

#define T_OUT 2.0       /* outside temperature */
#define T_IN 0.0        /* inside temperature */
#define SPEED 0.0       /* speed of drift to the right */

#define ADD_NOISE 0     /* set to 1 to add noise, set to 2 to add noise in right half */
#define NOISE_INTENSITY 0.005      /* noise intensity */
#define CHANGE_NOISE 1      /* set to 1 to increase noise intensity */
#define NOISE_FACTOR 40.0   /* factor by which to increase noise intensity */
#define NOISE_INITIAL_TIME 100  /* initial time during which noise remains constant */

#define CHANGE_VISCOSITY 0      /* set to 1 to change the viscosity in the course of the simulation */
#define ADJUST_INTSTEP 0       /* set to 1 to decrease integration step when viscosity increases */
#define VISCOSITY_INITIAL_TIME 10  /* initial time during which viscosity remains constant */
#define VISCOSITY_FACTOR 100.0   /* factor by which to change viscosity */
#define VISCOSITY_MAX 2.0        /* max value of viscosity beyond which NVID is increased and integration step is decrase, 
                                    for numerical stability */
                                        
#define CHANGE_RPSLZB 0         /* set to 1 to change second parameter in Rock-Paper-Scissors-Lizard-Spock equation */
#define RPSLZB_CHANGE 0.75      /* factor by which to rpslzb parameter */
#define RPSLZB_INITIAL_TIME 0   /* initial time during which rpslzb remains constant */
#define RPSLZB_FINAL_TIME 500   /* final time during which rpslzb remains constant */
                                      

/* Boundary conditions, see list in global_pdes.c  */

#define B_COND 1

#define IN_OUT_FLOW_BC 1    /* set to 1 to add in-flow/out-flow boundary conditions (for Euler equation) */
#define IN_OUT_FLOW_AMP 1.0 /* amplitude of in-flow/out-flow boundary conditions (for Euler equation) */

#define EULER_GRADIENT_YSHIFT 0.0    /* y-shift in computation of gradient in Euler equation */

/* Parameters for length and speed of simulation */

#define NSTEPS 1500       /* number of frames of movie */
#define NVID 120          /* number of iterations between images displayed on screen */
#define ACCELERATION_FACTOR 1.0 /* factor by which to increase NVID in course of simulation */
#define DT_ACCELERATION_FACTOR 1.0 /* factor by which to increase time step in course of simulation  */
#define MAX_DT 0.024     /* maximal value of integration step */
#define NSEG 100         /* number of segments of boundary */
#define BOUNDARY_WIDTH 2    /* width of billiard boundary */

#define PAUSE 100       /* number of frames after which to pause */
#define PSLEEP 2         /* sleep time during pause */
#define SLEEP1  2        /* initial sleeping time */
#define SLEEP2  1        /* final sleeping time */
#define INITIAL_TIME 0  /* initial still time */
#define MID_FRAMES 50    /* number of still frames between parts of two-part movie */
#define END_FRAMES 50    /* number of still frames at end of movie */
#define FADE 1           /* set to 1 to fade at end of movie */

/* Visualisation */

#define PLOT_3D 0    /* controls whether plot is 2D or 3D */

#define ROTATE_VIEW 0       /* set to 1 to rotate position of observer */
#define ROTATE_ANGLE 360.0  /* total angle of rotation during simulation */

#define DRAW_PERIODICISED 1     /* set to 1 to repeat wave periodically in x and y directions */

/* Plot type - color scheme */

#define CPLOT 51
#define CPLOT_B 54

/* Plot type - height of 3D plot */

#define ZPLOT 52     /* z coordinate in 3D plot */
#define ZPLOT_B 51    /* z coordinate in second 3D plot */

#define AMPLITUDE_HIGH_RES 1    /* set to 1 to increase resolution of P_3D_AMPLITUDE plot */
#define SHADE_3D 1              /* set to 1 to change luminosity according to normal vector */
#define NON_DIRICHLET_BC 0      /* set to 1 to draw only facets in domain, if field is not zero on boundary */
#define WRAP_ANGLE 1            /* experimental: wrap angle to [0, 2Pi) for interpolation in angle schemes */
#define FADE_IN_OBSTACLE 0      /* set to 1 to fade color inside obstacles */
#define DRAW_OUTSIDE_GRAY 0     /* experimental - draw outside of billiard in gray */
#define ADD_POTENTIAL_TO_Z 1    /* set to 1 to add the external potential to z-coordinate of plot */
#define ADD_POT_CONSTANT 0.35   /* constant in front of added potential */

#define PLOT_SCALE_ENERGY 0.05      /* vertical scaling in energy plot */

#define PRINT_TIME 0        /* set to 1 to print running time */
#define PRINT_VISCOSITY 0   /* set to 1 to print viscosity */
#define PRINT_RPSLZB 0      /* set to 1 to print rpslzb parameter */
#define PRINT_PROBABILITIES 0   /* set to 1 to print probabilities (for Ehrenfest urn configuration) */
#define PRINT_NOISE 0       /* set to 1 to print noise intensity */

#define DRAW_FIELD_LINES 0  /* set to 1 to draw field lines */
#define FIELD_LINE_WIDTH 1  /* width of field lines */
#define N_FIELD_LINES 120   /* number of field lines */
#define FIELD_LINE_FACTOR 120 /* factor controlling precision when computing origin of field lines */
#define DRAW_BILLIARD 1     /* set to 1 to draw boundary */
#define DRAW_BILLIARD_FRONT 1     /* set to 1 to draw boundary */
#define FILL_BILLIARD_COMPLEMENT 1  /* set to 1 to fill complement of billiard (for certain shapes only) */

/* 3D representation */

#define REPRESENTATION_3D 1     /* choice of 3D representation */ 

#define REP_AXO_3D 0        /* linear projection (axonometry) */
#define REP_PROJ_3D 1       /* projection on plane orthogonal to observer line of sight */

/* Color schemes, see list in global_pdes.c  */

#define COLOR_PALETTE 10       /* Color palette, see list in global_pdes.c  */
#define COLOR_PALETTE_B 13     /* Color palette, see list in global_pdes.c  */

#define BLACK 1          /* black background */

#define COLOR_SCHEME 3   /* choice of color scheme */

#define COLOR_PHASE_SHIFT 0.0   /* phase shift of color scheme, in units of Pi */

#define SCALE 0          /* set to 1 to adjust color scheme to variance of field */
#define SLOPE 1.0        /* sensitivity of color on wave amplitude */
#define VSCALE_AMPLITUDE 0.5      /* additional scaling factor for color scheme P_3D_AMPLITUDE */
#define ATTENUATION 0.0  /* exponential attenuation coefficient of contrast with time */
#define CURL_SCALE 0.000015   /* scaling factor for curl representation */
#define RESCALE_COLOR_IN_CENTER 0   /* set to 1 to decrease color intentiy in the center (for wave escaping ring) */
#define SLOPE_SCHROD_LUM 50.0       /* sensitivity of luminosity on module, for color scheme Z_ARGUMENT */
#define MIN_SCHROD_LUM 0.2       /* minimal luminosity in color scheme Z_ARGUMENT*/
#define VSCALE_PRESSURE 0.5      /* additional scaling factor for color scheme Z_EULER_PRESSURE */
#define PRESSURE_SHIFT 10.0        /* shift for color scheme Z_EULER_PRESSURE */
#define PRESSURE_LOG_SHIFT -2.5     /* shift for color scheme Z_EULER_PRESSURE */

#define COLORHUE 260     /* initial hue of water color for scheme C_LUM */
#define COLORDRIFT 0.0   /* how much the color hue drifts during the whole simulation */
#define LUMMEAN 0.5      /* amplitude of luminosity variation for scheme C_LUM */
#define LUMAMP 0.3       /* amplitude of luminosity variation for scheme C_LUM */
#define HUEMEAN 359.0    /* mean value of hue for color scheme C_HUE */
#define HUEAMP -359.0    /* amplitude of variation of hue for color scheme C_HUE */
#define E_SCALE 100.0    /* scaling factor for energy representation */
#define FLUX_SCALE 100.0 /* scaling factor for energy representation */
#define LOG_SCALE 0.5    /* scaling factor for energy log representation */
#define LOG_SHIFT 1.0 
#define LOG_MIN 1.0e-3   /* floor value for log vorticity plot */

#define NXMAZE 7      /* width of maze */
#define NYMAZE 7      /* height of maze */
#define MAZE_MAX_NGBH 4     /* max number of neighbours of maze cell */
#define RAND_SHIFT 0        /* seed of random number generator */
#define MAZE_XSHIFT 0.0     /* horizontal shift of maze */

#define DRAW_COLOR_SCHEME 1     /* set to 1 to plot the color scheme */
#define COLORBAR_RANGE 2.0      /* scale of color scheme bar */
#define COLORBAR_RANGE_B 2.0   /* scale of color scheme bar for 2nd part */
#define ROTATE_COLOR_SCHEME 0   /* set to 1 to draw color scheme horizontally */

/* only for compatibility with wave_common.c */
#define TWOSPEEDS 0          /* set to 1 to replace hardcore boundary by medium with different speed */
#define OMEGA 0.005        /* frequency of periodic excitation */
#define COURANT 0.08       /* Courant number */
#define COURANTB 0.03      /* Courant number in medium B */
#define INITIAL_AMP 0.5         /* amplitude of initial condition */
#define INITIAL_VARIANCE 0.0002  /* variance of initial condition */
#define INITIAL_WAVELENGTH  0.1  /* wavelength of initial condition */
#define VSCALE_ENERGY 200.0       /* additional scaling factor for color scheme P_3D_ENERGY */
#define PHASE_FACTOR 20.0       /* factor in computation of phase in color scheme P_3D_PHASE */
#define PHASE_SHIFT 0.0      /* shift of phase in color scheme P_3D_PHASE */
#define OSCILLATION_SCHEDULE 0  /* oscillation schedule, see list in global_pdes.c */
#define AMPLITUDE 0.8      /* amplitude of periodic excitation */ 
#define ACHIRP 0.2        /* acceleration coefficient in chirp */
#define DAMPING 0.0        /* damping of periodic excitation */
#define COMPARISON 0        /* set to 1 to compare two different patterns (beta) */
#define B_DOMAIN_B 20       /* second domain shape, for comparisons */
#define CIRCLE_PATTERN_B 0  /* second pattern of circles or polygons */
/* end of constants added only for compatibility with wave_common.c */


double u_3d[2] = {0.75, -0.45};     /* projections of basis vectors for REP_AXO_3D representation */
double v_3d[2] = {-0.75, -0.45};
double w_3d[2] = {0.0, 0.015};
double light[3] = {0.816496581, -0.40824829, 0.40824829};      /* vector of "light" direction for P_3D_ANGLE color scheme */
double observer[3] = {8.0, 8.0, 8.0};    /* location of observer for REP_PROJ_3D representation */ 
int reset_view = 0;         /* switch to reset 3D view parameters (for option ROTATE_VIEW) */

#define Z_SCALING_FACTOR 0.08  /* overall scaling factor of z axis for REP_PROJ_3D representation */
#define XY_SCALING_FACTOR 1.7  /* overall scaling factor for on-screen (x,y) coordinates after projection */
#define ZMAX_FACTOR 1.0        /* max value of z coordinate for REP_PROJ_3D representation */
#define XSHIFT_3D -0.1         /* overall x shift for REP_PROJ_3D representation */
#define YSHIFT_3D 0.1          /* overall y shift for REP_PROJ_3D representation */
#define BORDER_PADDING 0       /* distance from boundary at which to plot points, to avoid boundary effects due to gradient */

/* For debugging purposes only */
#define FLOOR 1         /* set to 1 to limit wave amplitude to VMAX */
#define VMAX 1000.0        /* max value of wave amplitude */
#define TEST_GRADIENT 0 /* print norm squared of gradient */

```

### 22 December 22 - Deep impact 8: Dropping heavy particles (mass ratio 20:1) into a Lennard-Jones fluid ###

**Program:** `lennardjones.c` 

```
#define MOVIE 1         /* set to 1 to generate movie */
#define DOUBLE_MOVIE 0  /* set to 1 to produce movies for wave height and energy simultaneously */

#define TIME_LAPSE 1     /* set to 1 to add a time-lapse movie at the end */
                         /* so far incompatible with double movie */
#define TIME_LAPSE_FACTOR 3    /* factor of time-lapse movie */
#define TIME_LAPSE_FIRST 1  /* set to 1 to show time-lapse version first */

/* General geometrical parameters */

#define WINWIDTH 	1280  /* window width */
#define WINHEIGHT 	720   /* window height */

#define XMIN -2.0
#define XMAX 2.0	/* x interval */
#define YMIN -1.125
#define YMAX 1.125	/* y interval for 9/16 aspect ratio */

#define INITXMIN -1.99
#define INITXMAX 1.99	/* x interval for initial condition */
#define INITYMIN -1.1
#define INITYMAX 1.1	/* y interval for initial condition */

#define BCXMIN -2.0
#define BCXMAX 2.0	/* x interval for boundary condition */
#define BCYMIN -1.125
#define BCYMAX 1.125	/* y interval for boundary condition */

#define OBSXMIN -2.0
#define OBSXMAX 2.0     /* x interval for motion of obstacle */

#define CIRCLE_PATTERN 1   /* pattern of circles, see list in global_ljones.c */

#define ADD_FIXED_OBSTACLES 0   /* set to 1 do add fixed circular obstacles */
#define OBSTACLE_PATTERN 3  /* pattern of obstacles, see list in global_ljones.c */

#define ADD_FIXED_SEGMENTS 1    /* set to 1 to add fixed segments as obstacles */
#define SEGMENT_PATTERN 20    /* pattern of repelling segments, see list in global_ljones.c */
#define ROCKET_SHAPE 2        /* shape of rocket combustion chamber, see list in global_ljones.c */
#define ROCKET_SHAPE_B 2      /* shape of second rocket */
#define NOZZLE_SHAPE 1        /* shape of nozzle, see list in global_ljones.c */
#define NOZZLE_SHAPE_B 0      /* shape of nozzle for second rocket, see list in global_ljones.c */

#define TWO_TYPES 0         /* set to 1 to have two types of particles */
#define TPYE_PROPORTION 0.5 /* proportion of particles of first type */
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
#define PDISC_DISTANCE 3.0  /* minimal distance in Poisson disc process, controls density of particles */
#define PDISC_CANDIDATES 100 /* number of candidates in construction of Poisson disc process */
#define RANDOM_POLY_ANGLE 0 /* set to 1 to randomize angle of polygons */

#define LAMBDA 0.3	    /* parameter controlling the dimensions of domain */
#define MU 0.009 	    /* parameter controlling radius of particles */
#define MU_B 0.009          /* parameter controlling radius of particles of second type */
#define NPOLY 3             /* number of sides of polygon */
#define APOLY 1.0           /* angle by which to turn polygon, in units of Pi/2 */ 
#define MDEPTH 4            /* depth of computation of Menger gasket */
#define MRATIO 3            /* ratio defining Menger gasket */
#define MANDELLEVEL 1000    /* iteration level for Mandelbrot set */
#define MANDELLIMIT 10.0    /* limit value for approximation of Mandelbrot set */
#define FOCI 1              /* set to 1 to draw focal points of ellipse */
#define NGRIDX 100           /* number of grid point for grid of disks */
#define NGRIDY 75           /* number of grid point for grid of disks */
#define EHRENFEST_RADIUS 0.9    /* radius of container for Ehrenfest urn configuration */
#define EHRENFEST_WIDTH 0.035     /* width of tube for Ehrenfest urn configuration */
#define TWO_CIRCLES_RADIUS_RATIO 0.8    /* ratio of radii for S_TWO_CIRCLES_EXT segment configuration */
#define DAM_WIDTH 0.05       /* width of dam for S_DAM segment configuration */

#define X_SHOOTER -0.2
#define Y_SHOOTER -0.6
#define X_TARGET 0.4
#define Y_TARGET 0.7        /* shooter and target positions in laser fight */

/* Parameters for length and speed of simulation */

#define NSTEPS 2000      /* number of frames of movie */
#define NVID 100         /* number of iterations between images displayed on screen */
#define NSEG 250         /* number of segments of boundary */
#define INITIAL_TIME 200     /* time after which to start saving frames */
#define OBSTACLE_INITIAL_TIME 350     /* time after which to start moving obstacle */
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
#define PLOT_B 8        /* plot type for second movie */

#define DRAW_BONDS 0    /* set to 1 to draw bonds between neighbours */
#define COLOR_BONDS 1   /* set to 1 to color bonds according to length */
#define FILL_TRIANGLES 0    /* set to 1 to fill triangles between neighbours */
#define ALTITUDE_LINES 0    /* set to 1 to add horizontal lines to show altitude */
#define COLOR_SEG_GROUPS 1  /* set to 1 to collor segment groups differently */

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

#define ENERGY_HUE_MIN 330.0        /* color of original particle */
#define ENERGY_HUE_MAX 50.0         /* color of saturated particle */
#define PARTICLE_HUE_MIN 359.0      /* color of original particle */
#define PARTICLE_HUE_MAX 0.0        /* color of saturated particle */
#define PARTICLE_EMAX 3.0e2         /* energy of particle with hottest color */
#define HUE_TYPE0 70.0     /* hue of particles of type 0 */
#define HUE_TYPE1 280.0      /* hue of particles of type 1 */
#define HUE_TYPE2 .0      /* hue of particles of type 2 */
#define HUE_TYPE3 210.0     /* hue of particles of type 3 */

#define RANDOM_RADIUS 0     /* set to 1 for random circle radius */
#define DT_PARTICLE 3.0e-6    /* time step for particle displacement */
#define KREPEL 12.0          /* constant in repelling force between particles */
#define EQUILIBRIUM_DIST 4.0    /* Lennard-Jones equilibrium distance */
#define EQUILIBRIUM_DIST_B 2.5  /* Lennard-Jones equilibrium distance for second type of particle */
#define REPEL_RADIUS 20.0    /* radius in which repelling force acts (in units of particle radius) */
#define DAMPING 15.0          /* damping coefficient of particles */
#define INITIAL_DAMPING 40.0  /* damping coefficient of particles during initial phase */
#define PARTICLE_MASS 20.0    /* mass of particle of radius MU */
#define PARTICLE_MASS_B 1.0   /* mass of particle of radius MU_B */
#define PARTICLE_INERTIA_MOMENT 0.2     /* moment of inertia of particle */
#define PARTICLE_INERTIA_MOMENT_B 0.02     /* moment of inertia of second type of particle */
#define V_INITIAL 0.0        /* initial velocity range */
#define OMEGA_INITIAL 10.0        /* initial angular velocity range */

#define THERMOSTAT 0        /* set to 1 to switch on thermostat */
#define VARY_THERMOSTAT 0   /* set to 1 for time-dependent thermostat schedule */
#define SIGMA 5.0           /* noise intensity in thermostat */
#define BETA 0.02          /* initial inverse temperature */
#define MU_XI 0.01           /* friction constant in thermostat */
#define KSPRING_BOUNDARY 1.0e7    /* confining harmonic potential outside simulation region */
#define KSPRING_OBSTACLE 1.0e11    /* harmonic potential of obstacles */
#define NBH_DIST_FACTOR 7.5       /* radius in which to count neighbours */
#define GRAVITY 1500.0            /* gravity acting on all particles */
#define GRAVITY_X 0.0        /* horizontal gravity acting on all particles */
#define INCREASE_GRAVITY 0     /* set to 1 to increase gravity during the simulation */
#define GRAVITY_SCHEDULE 2     /* type of gravity schedule, see list in global_ljones.c */
#define GRAVITY_FACTOR 100.0    /* factor by which to increase gravity */
#define GRAVITY_INITIAL_TIME 200    /* time at start of simulation with constant gravity */
#define GRAVITY_RESTORE_TIME 700    /* time at end of simulation with gravity restored to initial value */

#define ROTATION 0           /* set to 1 to include rotation of particles */
#define COUPLE_ANGLE_TO_THERMOSTAT 0    /* set to 1 to couple angular degrees of freedom to thermostat */
#define DIMENSION_FACTOR 1.0  /* scaling factor taking into account number of degrees of freedom */  
#define KTORQUE 100.0          /* force constant in angular dynamics */
#define KTORQUE_B 10.0        /* force constant in angular dynamics */
#define KTORQUE_DIFF 150.0    /* force constant in angular dynamics for different particles */
#define DRAW_SPIN 0           /* set to 1 to draw spin vectors of particles */
#define DRAW_SPIN_B 0         /* set to 1 to draw spin vectors of particles */
#define DRAW_CROSS 1          /* set to 1 to draw cross on particles of second type */
#define SPIN_RANGE 7.0       /* range of spin-spin interaction */
#define SPIN_RANGE_B 5.0     /* range of spin-spin interaction for second type of particle */
#define QUADRUPOLE_RATIO 0.6  /* anisotropy in quadrupole potential */ 

#define INCREASE_BETA 0  /* set to 1 to increase BETA during simulation */
#define BETA_FACTOR 0.025   /* factor by which to change BETA during simulation */
#define N_TOSCILLATIONS 1.5   /* number of temperature oscillations in BETA schedule */
#define NO_OSCILLATION 1        /* set to 1 to have exponential BETA change only */
#define MIDDLE_CONSTANT_PHASE 1600   /* final phase in which temperature is constant */
#define FINAL_DECREASE_PHASE 1500    /* final phase in which temperature decreases */ 
#define FINAL_CONSTANT_PHASE -1     /* final phase in which temperature is constant */

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
#define PARTIAL_THERMO_REGION 1     /* region for partial thermostat coupling (see list in global_ljones.c) */
#define PARTIAL_THERMO_SHIFT 0.2    /* distance from obstacle at the right of which particles are coupled to thermostat */
#define PARTIAL_THERMO_WIDTH 0.5    /* vertical size of partial thermostat coupling */
#define PARTIAL_THERMO_HEIGHT 0.2   /* vertical size of partial thermostat coupling */

#define INCREASE_KREPEL 0   /* set to 1 to increase KREPEL during simulation */
#define KREPEL_FACTOR 1000.0   /* factor by which to change KREPEL during simulation */

#define PART_AT_BOTTOM 0     /* set to 1 to include "seed" particles at bottom */
#define MASS_PART_BOTTOM 10000.0 /* mass of particles at bottom */
#define NPART_BOTTOM 100     /* number of particles at the bottom */

#define ADD_PARTICLES 0    /* set to 1 to add particles */
#define ADD_TIME 0       /* time at which to add first particle */
#define ADD_PERIOD 10000       /* time interval between adding further particles */
#define FINAL_NOADD_PERIOD 200  /* final period where no particles are added */
#define SAFETY_FACTOR 2.0  /* no particles are added at distance less than MU*SAFETY_FACTOR of other particles */

#define TRACER_PARTICLE 0   /* set to 1 to have a tracer particle */
#define N_TRACER_PARTICLES 3    /* number of tracer particles */
#define TRAJECTORY_LENGTH 8000   /* length of recorded trajectory */
#define TRACER_PARTICLE_MASS 4.0    /* relative mass of tracer particle */
#define TRAJECTORY_WIDTH 3      /* width of tracer particle trajectory */

#define ROTATE_BOUNDARY 0           /* set to 1 to rotate the repelling segments */
#define SMOOTH_ROTATION 1           /* set to 1 to update segments at each time step (rather than at each movie frame) */
#define PERIOD_ROTATE_BOUNDARY 1000  /* period of rotating boundary */
#define ROTATE_INITIAL_TIME 0       /* initial time without rotation */
#define ROTATE_FINAL_TIME 100       /* final time without rotation */
#define ROTATE_CHANGE_TIME 0.33     /* relative duration of acceleration/deceleration phases */
#define OMEGAMAX 100.0              /* maximal rotation speed */
#define PRINT_OMEGA 0               /* set to 1 to print angular speed */
#define PRINT_PARTICLE_SPEEDS 0     /* set to 1 to print average speeds/momenta of particles */
#define PRINT_SEGMENTS_SPEEDS 0     /* set to 1 to print velocity of moving segments */

#define MOVE_BOUNDARY 0        /* set to 1 to move repelling segments, due to force from particles */
#define SEGMENTS_MASS 40.0     /* mass of collection of segments */
#define DEACTIVATE_SEGMENT 1    /* set to 1 to deactivate last segment after a certain time */
#define SEGMENT_DEACTIVATION_TIME 200   /* time at which to deactivate last segment */
#define RELEASE_ROCKET_AT_DEACTIVATION 1    /* set to 1 to limit segments velocity before segment release */
#define SEGMENTS_X0 1.0        /* initial position of segments */
#define SEGMENTS_Y0 0.8        /* initial position of segments */
#define SEGMENTS_VX0 0.0       /* initial velocity of segments */
#define SEGMENTS_VY0 0.0      /* initial velocity of segments */
#define DAMP_SEGS_AT_NEGATIVE_Y 0   /* set to 1 to dampen segments when y coordinate is negative */

#define MOVE_SEGMENT_GROUPS 1       /* set to 1 to group segments into moving units */
#define SEGMENT_GROUP_MASS 750.0    /* mass of segment group */
#define SEGMENT_GROUP_I 500.0       /* moment of inertia of segment group */
#define SEGMENT_GROUP_DAMPING 0.0   /* damping of segment groups */
#define GROUP_REPULSION 1           /* set to 1 for groups of segments to repel each other */
#define KSPRING_GROUPS 1.0e11       /* harmonic potential between segment groups */
#define GROUP_WIDTH 0.05            /* interaction width of groups */
#define GROUP_G_REPEL 0             /* set to 1 to add repulsion between centers of mass of groups */
#define GROUP_G_REPEL_RADIUS 1.2    /* radius within which centers of mass of groups repel each other */
#define TRACK_SEGMENT_GROUPS 0      /* set to 1 for view to track group of segments */
#define TRACK_X_PADDING 2.0         /* distance from x boundary where tracking starts */

#define POSITION_DEPENDENT_TYPE 1   /* set to 1 to make particle type depend on initial position */
#define POSITION_Y_DEPENDENCE 1     /* set to 1 for the separation between particles to be horizontal */
#define PRINT_ENTROPY 0     /* set to 1 to compute entropy */

#define REACTION_DIFFUSION 0    /* set to 1 to simulate a chemical reaction (particles may change type) */
#define RD_TYPES 3              /* number of types in reaction-diffusion equation */
#define REACTION_PROB 0.0045      /* probability controlling reaction term */ 

#define PRINT_PARTICLE_NUMBER 0     /* set to 1 to print total number of particles */
#define PRINT_LEFT 0        /* set to 1 to print certain parameters at the top left instead of right */
#define PLOT_SPEEDS 0       /* set to 1 to add a plot of obstacle speeds (e.g. for rockets) */
#define PLOT_TRAJECTORIES 0     /* set to 1 to add a plot of obstacle trajectories (e.g. for rockets) */
#define VMAX_PLOT_SPEEDS 0.6    /* vertical scale of plot of obstacle speeds */

#define EHRENFEST_COPY 0    /* set to 1 to add equal number of larger particles (for Ehrenfest model) */

#define LID_MASS 1000.0     /* mass of lid for BC_RECTANGLE_LID b.c. */
#define LID_WIDTH 0.1       /* width of lid for BC_RECTANGLE_LID b.c. */
#define WALL_MASS 2000.0    /* mass of wall for BC_RECTANGLE_WALL b.c. */
#define WALL_FRICTION 0.0   /* friction on wall for BC_RECTANGLE_WALL b.c. */
#define WALL_WIDTH 0.1      /* width of wall for BC_RECTANGLE_WALL b.c. */
#define WALL_VMAX 100.0     /* max speed of wall */
#define WALL_TIME 0         /* time during which to keep wall */

#define NXMAZE 10      /* width of maze */
#define NYMAZE 10      /* height of maze */
#define MAZE_MAX_NGBH 4     /* max number of neighbours of maze cell */
#define RAND_SHIFT 200      /* seed of random number generator */
#define MAZE_XSHIFT 0.0     /* horizontal shift of maze */

#define FLOOR_FORCE 1      /* set to 1 to limit force on particle to FMAX */
#define FMAX 1.0e12         /* maximal force */
#define FLOOR_OMEGA 0      /* set to 1 to limit particle momentum to PMAX */
#define PMAX 1000.0        /* maximal force */

#define HASHX 75   /* size of hashgrid in x direction */
#define HASHY 37    /* size of hashgrid in y direction */
#define HASHMAX 100  /* maximal number of particles per hashgrid cell */
#define HASHGRID_PADDING 0.1    /* padding of hashgrid outside simulation window */

#define DRAW_COLOR_SCHEME 0     /* set to 1 to plot the color scheme */
#define COLORBAR_RANGE 8.0    /* scale of color scheme bar */
#define COLORBAR_RANGE_B 12.0    /* scale of color scheme bar for 2nd part */
#define ROTATE_COLOR_SCHEME 0   /* set to 1 to draw color scheme horizontally */

```

### 21 December 22 - Solving a maze with rounded corners using laser beams ###

**Program:** `particle_billiard.c` 

**Initial condition in function `animation()`:** `init_drop_config(-1.15, -0.01, 0.15*PI, 0.3*PI, configs);`

```
#define MOVIE 1         /* set to 1 to generate movie */

#define WINWIDTH 	1280  /* window width */
#define WINHEIGHT 	720   /* window height */

#define XMIN -2.0
#define XMAX 2.0	/* x interval */
#define YMIN -1.125
#define YMAX 1.125	/* y interval for 9/16 aspect ratio */

#define SCALING_FACTOR 1.0       /* scaling factor of drawing, needed for flower billiards, otherwise set to 1.0 */

/* Choice of the billiard table, see global_particles.c */

#define B_DOMAIN 31     /* choice of domain shape */

#define CIRCLE_PATTERN 1    /* pattern of circles */
#define POLYLINE_PATTERN 10  /* pattern of polyline */

#define ABSORBING_CIRCLES 0 /* set to 1 for circular scatterers to be absorbing */

#define NMAXCIRCLES 100000     /* total number of circles (must be at least NCX*NCY for square grid) */
#define NMAXPOLY 100000        /* total number of sides of polygonal line */   
#define NCX 30            /* number of circles in x direction */
#define NCY 20            /* number of circles in y direction */
#define NPOISSON 500        /* number of points for Poisson C_RAND_POISSON arrangement */
#define NGOLDENSPIRAL 2000  /* max number of points for C_GOLDEN_SPIRAL arrandement */
#define SDEPTH 1            /* Sierpinski gastket depth */

#define LAMBDA 1.5	/* parameter controlling shape of domain */
#define MU 0.005          /* second parameter controlling shape of billiard */
#define FOCI 1          /* set to 1 to draw focal points of ellipse */
#define NPOLY 6             /* number of sides of polygon */
#define APOLY 0.0           /* angle by which to turn polygon, in units of Pi/2 */ 
#define PENROSE_RATIO 2.5    /* parameter controlling the shape of small ellipses in Penrose room */

#define DRAW_BILLIARD 1     /* set to 1 to draw billiard */
#define DRAW_CONSTRUCTION_LINES 0   /* set to 1 to draw additional construction lines for billiard */
#define PERIODIC_BC 0       /* set to 1 to enforce periodic boundary conditions when drawing particles */

#define RESAMPLE 0      /* set to 1 if particles should be added when dispersion too large */
#define DEBUG 0         /* draw trajectories, for debugging purposes */

/* Simulation parameters */

#define NPART 10000    /* number of particles */
#define NPARTMAX 100000	/* maximal number of particles after resampling */
#define LMAX 0.01       /* minimal segment length triggering resampling */ 
#define DMIN 0.02       /* minimal distance to boundary for triggering resampling */ 
#define CYCLE 1         /* set to 1 for closed curve (start in all directions) */
#define SHOWTRAILS 1    /* set to 1 to keep trails of the particles */
#define SHOWZOOM 0      /* set to 1 to show zoom on specific area */
#define PRINT_PARTICLE_NUMBER 0 /* set to 1 to print number of particles */
#define PRINT_LEFT_RIGHT_PARTICLE_NUMBER 1 /* set to 1 to print number of particles on left and right side */
#define PRINT_COLLISION_NUMBER 0 /* set to 1 to print number of collisions */
#define TEST_ACTIVE 1   /* set to 1 to test whether particle is in billiard */

#define TEST_INITIAL_COND 0     /* set to 1 to allow only initial conditions that pass a test */

#define NSTEPS 6900      /* number of frames of movie */
#define TIME 1500        /* time between movie frames, for fluidity of real-time simulation */ 
#define DPHI 0.00002     /* integration step */
#define NVID 25          /* number of iterations between images displayed on screen */
#define END_FRAMES 25    /* number of still frames at the end of the movie */

/* Decreasing TIME accelerates the animation and the movie                               */
/* For constant speed of movie, TIME*DPHI should be kept constant                        */
/* However, increasing DPHI too much deterioriates quality of simulation                 */
/* NVID tells how often a picture is drawn in the animation, increase it for faster anim */
/* For a good quality movie, take for instance TIME = 400, DPHI = 0.00005, NVID = 100    */

/* Colors and other graphical parameters */

#define COLOR_PALETTE 10     /* Color palette, see list in global_pdes.c  */

#define NCOLORS 250      /* number of colors */
#define COLORSHIFT 0     /* hue of initial color */ 
#define RAINBOW_COLOR 0  /* set to 1 to use different colors for all particles */
#define FLOWER_COLOR 0   /* set to 1 to adapt initial colors to flower billiard (tracks vs core) */
#define NSEG 100         /* number of segments of boundary */
#define LENGTH 0.04       /* length of velocity vectors */
#define BILLIARD_WIDTH 2    /* width of billiard */
#define PARTICLE_WIDTH 2    /* width of particles */
#define FRONT_WIDTH 3       /* width of wave front */

#define BLACK 1             /* set to 1 for black background */
#define COLOR_OUTSIDE 0     /* set to 1 for colored outside */ 
#define OUTER_COLOR 270.0   /* color outside billiard */
#define PAINT_INT 0         /* set to 1 to paint interior in other color (for polygon/Reuleaux) */
#define PAINT_EXT 1         /* set to 1 to paint exterior */

#define PAUSE 1000       /* number of frames after which to pause */
#define PSLEEP 2         /* sleep time during pause */
#define SLEEP1  1        /* initial sleeping time */
#define SLEEP2  1       /* final sleeping time */

#define NXMAZE 32      /* width of maze */
#define NYMAZE 32      /* height of maze */
#define MAZE_MAX_NGBH 4     /* max number of neighbours of maze cell */
#define RAND_SHIFT 2        /* seed of random number generator */
#define MAZE_XSHIFT 0.0     /* horizontal shift of maze */
#define MAZE_RANDOM_FACTOR 0.1     /* randomization factor for S_MAZE_RANDOM */
#define MAZE_CORNER_RADIUS 0.4     /* radius of tounded corners in maze */

```

### 20 December 22 - A refracting von Koch snowflake ###

**Program:** `wave_billiard.c` 

**Initial condition in function `animation()`:** 
```
    init_circular_wave(0.0, 0.0, phi, psi, xy_in);
    
    if (ALTERNATE_OSCILLATING_SOURCE) sign = -sign;
    add_circular_wave(sign, 0.0, 0.0, phi, psi, xy_in);
``` 
    
```
#define MOVIE 1         /* set to 1 to generate movie */
#define DOUBLE_MOVIE 1  /* set to 1 to produce movies for wave height and energy simultaneously */
#define SAVE_MEMORY 0   /* set to 1 to save memory when writing tiff images */

/* General geometrical parameters */

#define WINWIDTH 	1920  /* window width */
#define WINHEIGHT 	1000  /* window height */
#define NX 3840          /* number of grid points on x axis */
#define NY 2000          /* number of grid points on y axis */

#define XMIN -2.0
#define XMAX 2.0	/* x interval  */
#define YMIN -1.041666667
#define YMAX 1.041666667	/* y interval for 9/16 aspect ratio */

#define HIGHRES 1       /* set to 1 if resolution of grid is double that of displayed image */

#define JULIA_SCALE 1.0 /* scaling for Julia sets */

/* Choice of the billiard table */

#define B_DOMAIN 41        /* choice of domain shape, see list in global_pdes.c */

#define CIRCLE_PATTERN 1   /* pattern of circles or polygons, see list in global_pdes.c */

#define COMPARISON 0        /* set to 1 to compare two different patterns (beta) */
#define B_DOMAIN_B 20       /* second domain shape, for comparisons */
#define CIRCLE_PATTERN_B 0  /* second pattern of circles or polygons */

#define P_PERCOL 0.25       /* probability of having a circle in C_RAND_PERCOL arrangement */
#define NPOISSON 300        /* number of points for Poisson C_RAND_POISSON arrangement */
#define RANDOM_POLY_ANGLE 1 /* set to 1 to randomize angle of polygons */

#define LAMBDA 0.85	    /* parameter controlling the dimensions of domain */
#define MU 0.85             /* parameter controlling the dimensions of domain */
#define NPOLY 24            /* number of sides of polygon */
#define APOLY 0.333333333   /* angle by which to turn polygon, in units of Pi/2 */ 
#define MDEPTH 6            /* depth of computation of Menger gasket */
#define MRATIO 3            /* ratio defining Menger gasket */
#define MANDELLEVEL 1000    /* iteration level for Mandelbrot set */
#define MANDELLIMIT 10.0    /* limit value for approximation of Mandelbrot set */
#define FOCI 1              /* set to 1 to draw focal points of ellipse */
#define NGRIDX 12           /* number of grid point for grid of disks */
#define NGRIDY 15           /* number of grid point for grid of disks */

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
#define OSCILLATION_SCHEDULE 1  /* oscillation schedule, see list in global_pdes.c */

#define OMEGA 0.0005       /* frequency of periodic excitation */
#define AMPLITUDE 0.8      /* amplitude of periodic excitation */ 
#define ACHIRP 0.25        /* acceleration coefficient in chirp */
#define DAMPING 0.0        /* damping of periodic excitation */
#define COURANT 0.05        /* Courant number */
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

#define ADD_OSCILLATING_SOURCE 1        /* set to 1 to add an oscillating wave source */
#define OSCILLATING_SOURCE_PERIOD 50    /* period of oscillating source */
#define ALTERNATE_OSCILLATING_SOURCE 1  /* set to 1 to alternate sign of oscillating source */

/* Boundary conditions, see list in global_pdes.c  */

#define B_COND 2

/* Parameters for length and speed of simulation */

#define NSTEPS 1800       /* number of frames of movie */
#define NVID 12          /* number of iterations between images displayed on screen */
#define NSEG 1000         /* number of segments of boundary */
#define INITIAL_TIME 0      /* time after which to start saving frames */
#define BOUNDARY_WIDTH 2    /* width of billiard boundary */
#define PRINT_SPEED 0       /* print speed of moving source */

#define PAUSE 200       /* number of frames after which to pause */
#define PSLEEP 1         /* sleep time during pause */
#define SLEEP1  1        /* initial sleeping time */
#define SLEEP2  1        /* final sleeping time */
#define MID_FRAMES 20    /* number of still frames between parts of two-part movie */
#define END_FRAMES 100    /* number of still frames at end of movie */
#define FADE 1           /* set to 1 to fade at end of movie */

/* Parameters of initial condition */

#define INITIAL_AMP 0.35         /* amplitude of initial condition */
#define INITIAL_VARIANCE 0.0002  /* variance of initial condition */
#define INITIAL_WAVELENGTH  0.01  /* wavelength of initial condition */

/* Plot type, see list in global_pdes.c  */

#define PLOT 6

#define PLOT_B 7        /* plot type for second movie */

/* Color schemes */

#define COLOR_PALETTE 12     /* Color palette, see list in global_pdes.c  */
#define COLOR_PALETTE_B 18     /* Color palette, see list in global_pdes.c  */

#define BLACK 1          /* background */

#define COLOR_SCHEME 3   /* choice of color scheme, see list in global_pdes.c  */

#define SCALE 0          /* set to 1 to adjust color scheme to variance of field */
#define SLOPE 1.0        /* sensitivity of color on wave amplitude */
#define PHASE_FACTOR 1.0       /* factor in computation of phase in color scheme P_3D_PHASE */
#define PHASE_SHIFT 0.0      /* shift of phase in color scheme P_3D_PHASE */
#define ATTENUATION 0.0  /* exponential attenuation coefficient of contrast with time */
#define E_SCALE 100.0     /* scaling factor for energy representation */
#define LOG_SCALE 0.25     /* scaling factor for energy log representation */
#define LOG_SHIFT 0.0     /* shift of colors on log scale */
#define FLUX_SCALE 1.0e3    /* scaling factor for enegy flux represtnation */
#define RESCALE_COLOR_IN_CENTER 0   /* set to 1 to decrease color intentiy in the center (for wave escaping ring) */

#define COLORHUE 260     /* initial hue of water color for scheme C_LUM */
#define COLORDRIFT 0.0   /* how much the color hue drifts during the whole simulation */
#define LUMMEAN 0.5      /* amplitude of luminosity variation for scheme C_LUM */
#define LUMAMP 0.3       /* amplitude of luminosity variation for scheme C_LUM */
#define HUEMEAN 180.0    /* mean value of hue for color scheme C_HUE */
#define HUEAMP -180.0      /* amplitude of variation of hue for color scheme C_HUE */

#define DRAW_COLOR_SCHEME 1    /* set to 1 to plot the color scheme */
#define COLORBAR_RANGE 1.5     /* scale of color scheme bar */
#define COLORBAR_RANGE_B 2.5    /* scale of color scheme bar for 2nd part */
#define ROTATE_COLOR_SCHEME 0   /* set to 1 to draw color scheme horizontally */

#define SAVE_TIME_SERIES 0      /* set to 1 to save wave time series at a point */

#define NXMAZE 7      /* width of maze */
#define NYMAZE 7      /* height of maze */
#define MAZE_MAX_NGBH 4     /* max number of neighbours of maze cell */
#define RAND_SHIFT 24        /* seed of random number generator */
#define MAZE_XSHIFT 0.0     /* horizontal shift of maze */

/* for compatibility with sub_wave and sub_maze */
#define ADD_POTENTIAL 0
#define POT_MAZE 7
#define POTENTIAL 0
/* end of constants only used by sub_wave and sub_maze */


/* For debugging purposes only */
#define FLOOR 0         /* set to 1 to limit wave amplitude to VMAX */
#define VMAX 10.0       /* max value of wave amplitude */

```

### 19 December 22 - Bloopers 4: Blow-up in an incompressible Euler flow ###

**Program:** `rde.c` 

**Initial condition in function `animation()`:** `init_shear_flow(1.0, 0.02, 0.15, 1, 1, 0.0, phi, xy_in);`

```
#define MOVIE 1         /* set to 1 to generate movie */
#define DOUBLE_MOVIE 1  /* set to 1 to produce movies for wave height and energy simultaneously */
#define SAVE_MEMORY 0           /* set to 1 to save memory when writing tiff images */

/* General geometrical parameters */

#define WINWIDTH 	1920  /* window width */
#define WINHEIGHT 	1000  /* window height */
#define NX 960          /* number of grid points on x axis */
#define NY 500          /* number of grid points on y axis */

#define XMIN -2.0
#define XMAX 2.0	/* x interval */
#define YMIN -1.041666667
#define YMAX 1.041666667	/* y interval for 9/16 aspect ratio */

/* Choice of simulated equation */

#define RDE_EQUATION 6  /* choice of reaction term, see list in global_3d.c */
#define NFIELDS 3       /* number of fields in reaction-diffusion equation */
#define NLAPLACIANS 1   /* number of fields for which to compute Laplacian */

#define ADD_POTENTIAL 0 /* set to 1 to add a potential (for Schrodinger equation) */
#define ADD_MAGNETIC_FIELD 0    /* set to 1 to add a magnetic field (for Schrodinger equation) - then set POTENTIAL 1 */
#define POTENTIAL 7     /* type of potential or vector potential, see list in global_3d.c  */

#define ANTISYMMETRIZE_WAVE_FCT 0   /* set tot 1 to make wave function antisymmetric */

#define JULIA_SCALE 0.5 /* scaling for Julia sets */

/* Choice of the billiard table */

#define B_DOMAIN 3          /* choice of domain shape, see list in global_pdes.c  */

#define CIRCLE_PATTERN 99    /* pattern of circles, see list in global_pdes.c */

#define P_PERCOL 0.25       /* probability of having a circle in C_RAND_PERCOL arrangement */
#define NPOISSON 300        /* number of points for Poisson C_RAND_POISSON arrangement */
#define RANDOM_POLY_ANGLE 0 /* set to 1 to randomize angle of polygons */

#define LAMBDA 0.25	    /* parameter controlling the dimensions of domain */
#define MU 0.15	            /* parameter controlling the dimensions of domain */
#define NPOLY 5             /* number of sides of polygon */
#define APOLY 2.0          /* angle by which to turn polygon, in units of Pi/2 */
#define MDEPTH 7            /* depth of computation of Menger gasket */
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

#define DT 0.0000012

#define VISCOSITY 2.0

#define RPSA 0.75         /* parameter in Rock-Paper-Scissors-type interaction */
#define RPSLZB 0.75       /* second parameter in Rock-Paper-Scissors-Lizard-Spock type interaction */

#define EPSILON 0.8     /* time scale separation */
#define DELTA 0.1       /* time scale separation */
#define FHNA 1.0        /* parameter in FHN equation */
#define FHNC -0.01      /* parameter in FHN equation */
#define K_HARMONIC 1.0  /* spring constant of harmonic potential */
#define K_COULOMB 0.5   /* constant in Coulomb potential */
#define V_MAZE 0.4      /* potential in walls of maze */
#define BZQ 0.0008      /* parameter in BZ equation */
#define BZF 1.2         /* parameter in BZ equation */
#define B_FIELD 10.0    /* magnetic field */
#define AB_RADIUS 0.2   /* radius of region with magnetic field for Aharonov-Bohm effect */
#define K_EULER 50.0    /* constant in stream function integration of Euler equation */

#define SMOOTHEN_VORTICITY 1    /* set to 1 to smoothen vorticity field in Euler equation */
#define SMOOTHEN_PERIOD 10      /* period between smoothenings */
#define SMOOTH_FACTOR 0.03      /* factor by which to smoothen */

#define ADD_TRACERS 1    /* set to 1 to add tracer particles (for Euler equations) */
#define N_TRACERS 1000    /* number of tracer particles */

#define T_OUT 2.0       /* outside temperature */
#define T_IN 0.0        /* inside temperature */
#define SPEED 0.0       /* speed of drift to the right */

#define ADD_NOISE 0     /* set to 1 to add noise, set to 2 to add noise in right half */
#define NOISE_INTENSITY 0.005      /* noise intensity */
#define CHANGE_NOISE 1      /* set to 1 to increase noise intensity */
#define NOISE_FACTOR 40.0   /* factor by which to increase noise intensity */
#define NOISE_INITIAL_TIME 100  /* initial time during which noise remains constant */

#define CHANGE_VISCOSITY 0      /* set to 1 to change the viscosity in the course of the simulation */
#define ADJUST_INTSTEP 0       /* set to 1 to decrease integration step when viscosity increases */
#define VISCOSITY_INITIAL_TIME 10  /* initial time during which viscosity remains constant */
#define VISCOSITY_FACTOR 100.0   /* factor by which to change viscosity */
#define VISCOSITY_MAX 2.0        /* max value of viscosity beyond which NVID is increased and integration step is decrase, 
                                    for numerical stability */
                                        
#define CHANGE_RPSLZB 0         /* set to 1 to change second parameter in Rock-Paper-Scissors-Lizard-Spock equation */
#define RPSLZB_CHANGE 0.75      /* factor by which to rpslzb parameter */
#define RPSLZB_INITIAL_TIME 0   /* initial time during which rpslzb remains constant */
#define RPSLZB_FINAL_TIME 500   /* final time during which rpslzb remains constant */
                                      

/* Boundary conditions, see list in global_pdes.c  */

#define B_COND 1

#define EULER_GRADIENT_YSHIFT 0.0    /* y-shift in computation of gradient in Euler equation */

/* Parameters for length and speed of simulation */

#define NSTEPS 1400       /* number of frames of movie */
#define NVID 200          /* number of iterations between images displayed on screen */
#define ACCELERATION_FACTOR 1.0 /* factor by which to increase NVID in course of simulation */
#define DT_ACCELERATION_FACTOR 1.0 /* factor by which to increase time step in course of simulation  */
#define MAX_DT 0.024     /* maximal value of integration step */
#define NSEG 100         /* number of segments of boundary */
#define BOUNDARY_WIDTH 2    /* width of billiard boundary */

#define PAUSE 100       /* number of frames after which to pause */
#define PSLEEP 2         /* sleep time during pause */
#define SLEEP1  2        /* initial sleeping time */
#define SLEEP2  1        /* final sleeping time */
#define INITIAL_TIME 0  /* initial still time */
#define MID_FRAMES 50    /* number of still frames between parts of two-part movie */
#define END_FRAMES 50    /* number of still frames at end of movie */
#define FADE 1           /* set to 1 to fade at end of movie */

/* Visualisation */

#define PLOT_3D 0    /* controls whether plot is 2D or 3D */

#define ROTATE_VIEW 0       /* set to 1 to rotate position of observer */
#define ROTATE_ANGLE 360.0  /* total angle of rotation during simulation */

#define DRAW_PERIODICISED 1     /* set to 1 to repeat wave periodically in x and y directions */

/* Plot type - color scheme */

#define CPLOT 54
#define CPLOT_B 53

/* Plot type - height of 3D plot */

#define ZPLOT 52     /* z coordinate in 3D plot */
#define ZPLOT_B 51    /* z coordinate in second 3D plot */

#define AMPLITUDE_HIGH_RES 1    /* set to 1 to increase resolution of P_3D_AMPLITUDE plot */
#define SHADE_3D 1              /* set to 1 to change luminosity according to normal vector */
#define NON_DIRICHLET_BC 0      /* set to 1 to draw only facets in domain, if field is not zero on boundary */
#define WRAP_ANGLE 1            /* experimental: wrap angle to [0, 2Pi) for interpolation in angle schemes */
#define FADE_IN_OBSTACLE 0      /* set to 1 to fade color inside obstacles */
#define DRAW_OUTSIDE_GRAY 0     /* experimental - draw outside of billiard in gray */
#define ADD_POTENTIAL_TO_Z 1    /* set to 1 to add the external potential to z-coordinate of plot */
#define ADD_POT_CONSTANT 0.35   /* constant in front of added potential */

#define PLOT_SCALE_ENERGY 0.05      /* vertical scaling in energy plot */

#define PRINT_TIME 0        /* set to 1 to print running time */
#define PRINT_VISCOSITY 0   /* set to 1 to print viscosity */
#define PRINT_RPSLZB 0      /* set to 1 to print rpslzb parameter */
#define PRINT_PROBABILITIES 0   /* set to 1 to print probabilities (for Ehrenfest urn configuration) */
#define PRINT_NOISE 0       /* set to 1 to print noise intensity */

#define DRAW_FIELD_LINES 0  /* set to 1 to draw field lines */
#define FIELD_LINE_WIDTH 1  /* width of field lines */
#define N_FIELD_LINES 120   /* number of field lines */
#define FIELD_LINE_FACTOR 120 /* factor controlling precision when computing origin of field lines */
#define DRAW_BILLIARD 1     /* set to 1 to draw boundary */
#define DRAW_BILLIARD_FRONT 1     /* set to 1 to draw boundary */
#define FILL_BILLIARD_COMPLEMENT 1  /* set to 1 to fill complement of billiard (for certain shapes only) */

/* 3D representation */

#define REPRESENTATION_3D 1     /* choice of 3D representation */ 

#define REP_AXO_3D 0        /* linear projection (axonometry) */
#define REP_PROJ_3D 1       /* projection on plane orthogonal to observer line of sight */

/* Color schemes, see list in global_pdes.c  */

#define COLOR_PALETTE 10       /* Color palette, see list in global_pdes.c  */
#define COLOR_PALETTE_B 11     /* Color palette, see list in global_pdes.c  */

#define BLACK 1          /* black background */

#define COLOR_SCHEME 3   /* choice of color scheme */

#define COLOR_PHASE_SHIFT 0.0   /* phase shift of color scheme, in units of Pi */

#define SCALE 0          /* set to 1 to adjust color scheme to variance of field */
#define SLOPE 1.0        /* sensitivity of color on wave amplitude */
#define VSCALE_AMPLITUDE 0.5      /* additional scaling factor for color scheme P_3D_AMPLITUDE */
#define ATTENUATION 0.0  /* exponential attenuation coefficient of contrast with time */
#define CURL_SCALE 0.000015   /* scaling factor for curl representation */
#define RESCALE_COLOR_IN_CENTER 0   /* set to 1 to decrease color intentiy in the center (for wave escaping ring) */
#define SLOPE_SCHROD_LUM 50.0       /* sensitivity of luminosity on module, for color scheme Z_ARGUMENT */
#define MIN_SCHROD_LUM 0.2       /* minimal luminosity in color scheme Z_ARGUMENT*/
#define VSCALE_PRESSURE 0.35      /* additional scaling factor for color scheme Z_EULER_PRESSURE */
#define PRESSURE_SHIFT 250.0        /* shift for color scheme Z_EULER_PRESSURE */
#define PRESSURE_LOG_SHIFT -1.8     /* shift for color scheme Z_EULER_PRESSURE */

#define COLORHUE 260     /* initial hue of water color for scheme C_LUM */
#define COLORDRIFT 0.0   /* how much the color hue drifts during the whole simulation */
#define LUMMEAN 0.5      /* amplitude of luminosity variation for scheme C_LUM */
#define LUMAMP 0.3       /* amplitude of luminosity variation for scheme C_LUM */
#define HUEMEAN 359.0    /* mean value of hue for color scheme C_HUE */
#define HUEAMP -359.0    /* amplitude of variation of hue for color scheme C_HUE */
#define E_SCALE 100.0    /* scaling factor for energy representation */
#define FLUX_SCALE 100.0 /* scaling factor for energy representation */
#define LOG_SCALE 0.5    /* scaling factor for energy log representation */
#define LOG_SHIFT 1.0 
#define LOG_MIN 1.0e-3   /* floor value for log vorticity plot */

#define NXMAZE 7      /* width of maze */
#define NYMAZE 7      /* height of maze */
#define MAZE_MAX_NGBH 4     /* max number of neighbours of maze cell */
#define RAND_SHIFT 0        /* seed of random number generator */
#define MAZE_XSHIFT 0.0     /* horizontal shift of maze */

#define DRAW_COLOR_SCHEME 1     /* set to 1 to plot the color scheme */
#define COLORBAR_RANGE 2.0      /* scale of color scheme bar */
#define COLORBAR_RANGE_B 3.0   /* scale of color scheme bar for 2nd part */
#define ROTATE_COLOR_SCHEME 0   /* set to 1 to draw color scheme horizontally */

/* only for compatibility with wave_common.c */
#define TWOSPEEDS 0          /* set to 1 to replace hardcore boundary by medium with different speed */
#define OMEGA 0.005        /* frequency of periodic excitation */
#define COURANT 0.08       /* Courant number */
#define COURANTB 0.03      /* Courant number in medium B */
#define INITIAL_AMP 0.5         /* amplitude of initial condition */
#define INITIAL_VARIANCE 0.0002  /* variance of initial condition */
#define INITIAL_WAVELENGTH  0.1  /* wavelength of initial condition */
#define VSCALE_ENERGY 200.0       /* additional scaling factor for color scheme P_3D_ENERGY */
#define PHASE_FACTOR 20.0       /* factor in computation of phase in color scheme P_3D_PHASE */
#define PHASE_SHIFT 0.0      /* shift of phase in color scheme P_3D_PHASE */
#define OSCILLATION_SCHEDULE 0  /* oscillation schedule, see list in global_pdes.c */
#define AMPLITUDE 0.8      /* amplitude of periodic excitation */ 
#define ACHIRP 0.2        /* acceleration coefficient in chirp */
#define DAMPING 0.0        /* damping of periodic excitation */
#define COMPARISON 0        /* set to 1 to compare two different patterns (beta) */
#define B_DOMAIN_B 20       /* second domain shape, for comparisons */
#define CIRCLE_PATTERN_B 0  /* second pattern of circles or polygons */
/* end of constants added only for compatibility with wave_common.c */


double u_3d[2] = {0.75, -0.45};     /* projections of basis vectors for REP_AXO_3D representation */
double v_3d[2] = {-0.75, -0.45};
double w_3d[2] = {0.0, 0.015};
double light[3] = {0.816496581, -0.40824829, 0.40824829};      /* vector of "light" direction for P_3D_ANGLE color scheme */
double observer[3] = {8.0, 8.0, 8.0};    /* location of observer for REP_PROJ_3D representation */ 
int reset_view = 0;         /* switch to reset 3D view parameters (for option ROTATE_VIEW) */

#define Z_SCALING_FACTOR 0.08  /* overall scaling factor of z axis for REP_PROJ_3D representation */
#define XY_SCALING_FACTOR 1.7  /* overall scaling factor for on-screen (x,y) coordinates after projection */
#define ZMAX_FACTOR 1.0        /* max value of z coordinate for REP_PROJ_3D representation */
#define XSHIFT_3D -0.1         /* overall x shift for REP_PROJ_3D representation */
#define YSHIFT_3D 0.1          /* overall y shift for REP_PROJ_3D representation */
#define BORDER_PADDING 0       /* distance from boundary at which to plot points, to avoid boundary effects due to gradient */

/* For debugging purposes only */
#define FLOOR 1         /* set to 1 to limit wave amplitude to VMAX */
#define VMAX 1000.0        /* max value of wave amplitude */
#define TEST_GRADIENT 0 /* print norm squared of gradient */

```

### 18 December 22 - Deep impact 7: Dropping heavy particles (mass ratio 5:1) into a Lennard-Jones fluid ###

**Program:** `lennardjones.c` 

```
#define MOVIE 1         /* set to 1 to generate movie */
#define DOUBLE_MOVIE 0  /* set to 1 to produce movies for wave height and energy simultaneously */

#define TIME_LAPSE 1     /* set to 1 to add a time-lapse movie at the end */
                         /* so far incompatible with double movie */
#define TIME_LAPSE_FACTOR 3    /* factor of time-lapse movie */
#define TIME_LAPSE_FIRST 1  /* set to 1 to show time-lapse version first */

/* General geometrical parameters */

#define WINWIDTH 	1280  /* window width */
#define WINHEIGHT 	720   /* window height */

#define XMIN -2.0
#define XMAX 2.0	/* x interval */
#define YMIN -1.125
#define YMAX 1.125	/* y interval for 9/16 aspect ratio */

#define INITXMIN -1.99
#define INITXMAX 1.99	/* x interval for initial condition */
#define INITYMIN -1.1
#define INITYMAX 1.1	/* y interval for initial condition */

#define BCXMIN -2.0
#define BCXMAX 2.0	/* x interval for boundary condition */
#define BCYMIN -1.125
#define BCYMAX 1.125	/* y interval for boundary condition */

#define OBSXMIN -2.0
#define OBSXMAX 2.0     /* x interval for motion of obstacle */

#define CIRCLE_PATTERN 1   /* pattern of circles, see list in global_ljones.c */

#define ADD_FIXED_OBSTACLES 0   /* set to 1 do add fixed circular obstacles */
#define OBSTACLE_PATTERN 3  /* pattern of obstacles, see list in global_ljones.c */

#define ADD_FIXED_SEGMENTS 1    /* set to 1 to add fixed segments as obstacles */
#define SEGMENT_PATTERN 20    /* pattern of repelling segments, see list in global_ljones.c */
#define ROCKET_SHAPE 2        /* shape of rocket combustion chamber, see list in global_ljones.c */
#define ROCKET_SHAPE_B 2      /* shape of second rocket */
#define NOZZLE_SHAPE 1        /* shape of nozzle, see list in global_ljones.c */
#define NOZZLE_SHAPE_B 0      /* shape of nozzle for second rocket, see list in global_ljones.c */

#define TWO_TYPES 0         /* set to 1 to have two types of particles */
#define TPYE_PROPORTION 0.5 /* proportion of particles of first type */
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
#define PDISC_DISTANCE 3.0  /* minimal distance in Poisson disc process, controls density of particles */
#define PDISC_CANDIDATES 100 /* number of candidates in construction of Poisson disc process */
#define RANDOM_POLY_ANGLE 0 /* set to 1 to randomize angle of polygons */

#define LAMBDA 0.3	    /* parameter controlling the dimensions of domain */
#define MU 0.009 	    /* parameter controlling radius of particles */
#define MU_B 0.009          /* parameter controlling radius of particles of second type */
#define NPOLY 3             /* number of sides of polygon */
#define APOLY 1.0           /* angle by which to turn polygon, in units of Pi/2 */ 
#define MDEPTH 4            /* depth of computation of Menger gasket */
#define MRATIO 3            /* ratio defining Menger gasket */
#define MANDELLEVEL 1000    /* iteration level for Mandelbrot set */
#define MANDELLIMIT 10.0    /* limit value for approximation of Mandelbrot set */
#define FOCI 1              /* set to 1 to draw focal points of ellipse */
#define NGRIDX 100           /* number of grid point for grid of disks */
#define NGRIDY 75           /* number of grid point for grid of disks */
#define EHRENFEST_RADIUS 0.9    /* radius of container for Ehrenfest urn configuration */
#define EHRENFEST_WIDTH 0.035     /* width of tube for Ehrenfest urn configuration */
#define TWO_CIRCLES_RADIUS_RATIO 0.8    /* ratio of radii for S_TWO_CIRCLES_EXT segment configuration */
#define DAM_WIDTH 0.05       /* width of dam for S_DAM segment configuration */

#define X_SHOOTER -0.2
#define Y_SHOOTER -0.6
#define X_TARGET 0.4
#define Y_TARGET 0.7        /* shooter and target positions in laser fight */

/* Parameters for length and speed of simulation */

#define NSTEPS 1800      /* number of frames of movie */
#define NVID 100         /* number of iterations between images displayed on screen */
#define NSEG 250         /* number of segments of boundary */
#define INITIAL_TIME 200     /* time after which to start saving frames */
#define OBSTACLE_INITIAL_TIME 350     /* time after which to start moving obstacle */
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
#define PLOT_B 8        /* plot type for second movie */

#define DRAW_BONDS 0    /* set to 1 to draw bonds between neighbours */
#define COLOR_BONDS 1   /* set to 1 to color bonds according to length */
#define FILL_TRIANGLES 0    /* set to 1 to fill triangles between neighbours */
#define ALTITUDE_LINES 0    /* set to 1 to add horizontal lines to show altitude */
#define COLOR_SEG_GROUPS 1  /* set to 1 to collor segment groups differently */

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

#define ENERGY_HUE_MIN 330.0        /* color of original particle */
#define ENERGY_HUE_MAX 50.0         /* color of saturated particle */
#define PARTICLE_HUE_MIN 359.0      /* color of original particle */
#define PARTICLE_HUE_MAX 0.0        /* color of saturated particle */
#define PARTICLE_EMAX 3.0e2         /* energy of particle with hottest color */
#define HUE_TYPE0 70.0     /* hue of particles of type 0 */
#define HUE_TYPE1 280.0      /* hue of particles of type 1 */
#define HUE_TYPE2 .0      /* hue of particles of type 2 */
#define HUE_TYPE3 210.0     /* hue of particles of type 3 */

#define RANDOM_RADIUS 0     /* set to 1 for random circle radius */
#define DT_PARTICLE 3.0e-6    /* time step for particle displacement */
#define KREPEL 12.0          /* constant in repelling force between particles */
#define EQUILIBRIUM_DIST 4.0    /* Lennard-Jones equilibrium distance */
#define EQUILIBRIUM_DIST_B 2.5  /* Lennard-Jones equilibrium distance for second type of particle */
#define REPEL_RADIUS 20.0    /* radius in which repelling force acts (in units of particle radius) */
#define DAMPING 15.0          /* damping coefficient of particles */
#define INITIAL_DAMPING 40.0  /* damping coefficient of particles during initial phase */
#define PARTICLE_MASS 5.0    /* mass of particle of radius MU */
#define PARTICLE_MASS_B 1.0   /* mass of particle of radius MU_B */
#define PARTICLE_INERTIA_MOMENT 0.2     /* moment of inertia of particle */
#define PARTICLE_INERTIA_MOMENT_B 0.02     /* moment of inertia of second type of particle */
#define V_INITIAL 0.0        /* initial velocity range */
#define OMEGA_INITIAL 10.0        /* initial angular velocity range */

#define THERMOSTAT 0        /* set to 1 to switch on thermostat */
#define VARY_THERMOSTAT 0   /* set to 1 for time-dependent thermostat schedule */
#define SIGMA 5.0           /* noise intensity in thermostat */
#define BETA 0.02          /* initial inverse temperature */
#define MU_XI 0.01           /* friction constant in thermostat */
#define KSPRING_BOUNDARY 1.0e7    /* confining harmonic potential outside simulation region */
#define KSPRING_OBSTACLE 1.0e11    /* harmonic potential of obstacles */
#define NBH_DIST_FACTOR 7.5       /* radius in which to count neighbours */
#define GRAVITY 1500.0            /* gravity acting on all particles */
#define GRAVITY_X 0.0        /* horizontal gravity acting on all particles */
#define INCREASE_GRAVITY 0     /* set to 1 to increase gravity during the simulation */
#define GRAVITY_SCHEDULE 2     /* type of gravity schedule, see list in global_ljones.c */
#define GRAVITY_FACTOR 100.0    /* factor by which to increase gravity */
#define GRAVITY_INITIAL_TIME 200    /* time at start of simulation with constant gravity */
#define GRAVITY_RESTORE_TIME 700    /* time at end of simulation with gravity restored to initial value */

#define ROTATION 0           /* set to 1 to include rotation of particles */
#define COUPLE_ANGLE_TO_THERMOSTAT 0    /* set to 1 to couple angular degrees of freedom to thermostat */
#define DIMENSION_FACTOR 1.0  /* scaling factor taking into account number of degrees of freedom */  
#define KTORQUE 100.0          /* force constant in angular dynamics */
#define KTORQUE_B 10.0        /* force constant in angular dynamics */
#define KTORQUE_DIFF 150.0    /* force constant in angular dynamics for different particles */
#define DRAW_SPIN 0           /* set to 1 to draw spin vectors of particles */
#define DRAW_SPIN_B 0         /* set to 1 to draw spin vectors of particles */
#define DRAW_CROSS 1          /* set to 1 to draw cross on particles of second type */
#define SPIN_RANGE 7.0       /* range of spin-spin interaction */
#define SPIN_RANGE_B 5.0     /* range of spin-spin interaction for second type of particle */
#define QUADRUPOLE_RATIO 0.6  /* anisotropy in quadrupole potential */ 

#define INCREASE_BETA 0  /* set to 1 to increase BETA during simulation */
#define BETA_FACTOR 0.025   /* factor by which to change BETA during simulation */
#define N_TOSCILLATIONS 1.5   /* number of temperature oscillations in BETA schedule */
#define NO_OSCILLATION 1        /* set to 1 to have exponential BETA change only */
#define MIDDLE_CONSTANT_PHASE 1600   /* final phase in which temperature is constant */
#define FINAL_DECREASE_PHASE 1500    /* final phase in which temperature decreases */ 
#define FINAL_CONSTANT_PHASE -1     /* final phase in which temperature is constant */

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
#define PARTIAL_THERMO_REGION 1     /* region for partial thermostat coupling (see list in global_ljones.c) */
#define PARTIAL_THERMO_SHIFT 0.2    /* distance from obstacle at the right of which particles are coupled to thermostat */
#define PARTIAL_THERMO_WIDTH 0.5    /* vertical size of partial thermostat coupling */
#define PARTIAL_THERMO_HEIGHT 0.2   /* vertical size of partial thermostat coupling */

#define INCREASE_KREPEL 0   /* set to 1 to increase KREPEL during simulation */
#define KREPEL_FACTOR 1000.0   /* factor by which to change KREPEL during simulation */

#define PART_AT_BOTTOM 0     /* set to 1 to include "seed" particles at bottom */
#define MASS_PART_BOTTOM 10000.0 /* mass of particles at bottom */
#define NPART_BOTTOM 100     /* number of particles at the bottom */

#define ADD_PARTICLES 0    /* set to 1 to add particles */
#define ADD_TIME 0       /* time at which to add first particle */
#define ADD_PERIOD 10000       /* time interval between adding further particles */
#define FINAL_NOADD_PERIOD 200  /* final period where no particles are added */
#define SAFETY_FACTOR 2.0  /* no particles are added at distance less than MU*SAFETY_FACTOR of other particles */

#define TRACER_PARTICLE 0   /* set to 1 to have a tracer particle */
#define N_TRACER_PARTICLES 3    /* number of tracer particles */
#define TRAJECTORY_LENGTH 8000   /* length of recorded trajectory */
#define TRACER_PARTICLE_MASS 4.0    /* relative mass of tracer particle */
#define TRAJECTORY_WIDTH 3      /* width of tracer particle trajectory */

#define ROTATE_BOUNDARY 0           /* set to 1 to rotate the repelling segments */
#define SMOOTH_ROTATION 1           /* set to 1 to update segments at each time step (rather than at each movie frame) */
#define PERIOD_ROTATE_BOUNDARY 1000  /* period of rotating boundary */
#define ROTATE_INITIAL_TIME 0       /* initial time without rotation */
#define ROTATE_FINAL_TIME 100       /* final time without rotation */
#define ROTATE_CHANGE_TIME 0.33     /* relative duration of acceleration/deceleration phases */
#define OMEGAMAX 100.0              /* maximal rotation speed */
#define PRINT_OMEGA 0               /* set to 1 to print angular speed */
#define PRINT_PARTICLE_SPEEDS 0     /* set to 1 to print average speeds/momenta of particles */
#define PRINT_SEGMENTS_SPEEDS 0     /* set to 1 to print velocity of moving segments */

#define MOVE_BOUNDARY 0        /* set to 1 to move repelling segments, due to force from particles */
#define SEGMENTS_MASS 40.0     /* mass of collection of segments */
#define DEACTIVATE_SEGMENT 1    /* set to 1 to deactivate last segment after a certain time */
#define SEGMENT_DEACTIVATION_TIME 200   /* time at which to deactivate last segment */
#define RELEASE_ROCKET_AT_DEACTIVATION 1    /* set to 1 to limit segments velocity before segment release */
#define SEGMENTS_X0 1.0        /* initial position of segments */
#define SEGMENTS_Y0 0.8        /* initial position of segments */
#define SEGMENTS_VX0 0.0       /* initial velocity of segments */
#define SEGMENTS_VY0 0.0      /* initial velocity of segments */
#define DAMP_SEGS_AT_NEGATIVE_Y 0   /* set to 1 to dampen segments when y coordinate is negative */

#define MOVE_SEGMENT_GROUPS 1       /* set to 1 to group segments into moving units */
#define SEGMENT_GROUP_MASS 750.0    /* mass of segment group */
#define SEGMENT_GROUP_I 500.0       /* moment of inertia of segment group */
#define SEGMENT_GROUP_DAMPING 0.0   /* damping of segment groups */
#define GROUP_REPULSION 1           /* set to 1 for groups of segments to repel each other */
#define KSPRING_GROUPS 1.0e11       /* harmonic potential between segment groups */
#define GROUP_WIDTH 0.05            /* interaction width of groups */
#define GROUP_G_REPEL 0             /* set to 1 to add repulsion between centers of mass of groups */
#define GROUP_G_REPEL_RADIUS 1.2    /* radius within which centers of mass of groups repel each other */
#define TRACK_SEGMENT_GROUPS 0      /* set to 1 for view to track group of segments */
#define TRACK_X_PADDING 2.0         /* distance from x boundary where tracking starts */

#define POSITION_DEPENDENT_TYPE 1   /* set to 1 to make particle type depend on initial position */
#define POSITION_Y_DEPENDENCE 1     /* set to 1 for the separation between particles to be horizontal */
#define PRINT_ENTROPY 0     /* set to 1 to compute entropy */

#define REACTION_DIFFUSION 0    /* set to 1 to simulate a chemical reaction (particles may change type) */
#define RD_TYPES 3              /* number of types in reaction-diffusion equation */
#define REACTION_PROB 0.0045      /* probability controlling reaction term */ 

#define PRINT_PARTICLE_NUMBER 0     /* set to 1 to print total number of particles */
#define PRINT_LEFT 0        /* set to 1 to print certain parameters at the top left instead of right */
#define PLOT_SPEEDS 0       /* set to 1 to add a plot of obstacle speeds (e.g. for rockets) */
#define PLOT_TRAJECTORIES 0     /* set to 1 to add a plot of obstacle trajectories (e.g. for rockets) */
#define VMAX_PLOT_SPEEDS 0.6    /* vertical scale of plot of obstacle speeds */

#define EHRENFEST_COPY 0    /* set to 1 to add equal number of larger particles (for Ehrenfest model) */

#define LID_MASS 1000.0     /* mass of lid for BC_RECTANGLE_LID b.c. */
#define LID_WIDTH 0.1       /* width of lid for BC_RECTANGLE_LID b.c. */
#define WALL_MASS 2000.0    /* mass of wall for BC_RECTANGLE_WALL b.c. */
#define WALL_FRICTION 0.0   /* friction on wall for BC_RECTANGLE_WALL b.c. */
#define WALL_WIDTH 0.1      /* width of wall for BC_RECTANGLE_WALL b.c. */
#define WALL_VMAX 100.0     /* max speed of wall */
#define WALL_TIME 0         /* time during which to keep wall */

#define NXMAZE 10      /* width of maze */
#define NYMAZE 10      /* height of maze */
#define MAZE_MAX_NGBH 4     /* max number of neighbours of maze cell */
#define RAND_SHIFT 200      /* seed of random number generator */
#define MAZE_XSHIFT 0.0     /* horizontal shift of maze */

#define FLOOR_FORCE 1      /* set to 1 to limit force on particle to FMAX */
#define FMAX 1.0e12         /* maximal force */
#define FLOOR_OMEGA 0      /* set to 1 to limit particle momentum to PMAX */
#define PMAX 1000.0        /* maximal force */

#define HASHX 75   /* size of hashgrid in x direction */
#define HASHY 37    /* size of hashgrid in y direction */
#define HASHMAX 100  /* maximal number of particles per hashgrid cell */
#define HASHGRID_PADDING 0.1    /* padding of hashgrid outside simulation window */

#define DRAW_COLOR_SCHEME 0     /* set to 1 to plot the color scheme */
#define COLORBAR_RANGE 8.0    /* scale of color scheme bar */
#define COLORBAR_RANGE_B 12.0    /* scale of color scheme bar for 2nd part */
#define ROTATE_COLOR_SCHEME 0   /* set to 1 to draw color scheme horizontally */

```

### 17 December 22 - Mazes with randomized angles are harder to solve with laser beams ###

**Program:** `particle_billiard.c` 

**Initial condition in function `animation()`:** `init_drop_config(-1.15, -0.01, 0.15*PI, 0.3*PI, configs);`

```
#define MOVIE 1         /* set to 1 to generate movie */

#define WINWIDTH 	1280  /* window width */
#define WINHEIGHT 	720   /* window height */

#define XMIN -2.0
#define XMAX 2.0	/* x interval */
#define YMIN -1.125
#define YMAX 1.125	/* y interval for 9/16 aspect ratio */

#define SCALING_FACTOR 1.0       /* scaling factor of drawing, needed for flower billiards, otherwise set to 1.0 */

/* Choice of the billiard table, see global_particles.c */

#define B_DOMAIN 30     /* choice of domain shape */

#define CIRCLE_PATTERN 1    /* pattern of circles */
#define POLYLINE_PATTERN 12  /* pattern of polyline */

#define ABSORBING_CIRCLES 0 /* set to 1 for circular scatterers to be absorbing */

#define NMAXCIRCLES 100000     /* total number of circles (must be at least NCX*NCY for square grid) */
#define NMAXPOLY 100000        /* total number of sides of polygonal line */   
#define NCX 30            /* number of circles in x direction */
#define NCY 20            /* number of circles in y direction */
#define NPOISSON 500        /* number of points for Poisson C_RAND_POISSON arrangement */
#define NGOLDENSPIRAL 2000  /* max number of points for C_GOLDEN_SPIRAL arrandement */
#define SDEPTH 1            /* Sierpinski gastket depth */

#define LAMBDA 1.5	/* parameter controlling shape of domain */
#define MU 0.005          /* second parameter controlling shape of billiard */
#define FOCI 1          /* set to 1 to draw focal points of ellipse */
#define NPOLY 6             /* number of sides of polygon */
#define APOLY 0.0           /* angle by which to turn polygon, in units of Pi/2 */ 
#define PENROSE_RATIO 2.5    /* parameter controlling the shape of small ellipses in Penrose room */

#define DRAW_BILLIARD 1     /* set to 1 to draw billiard */
#define DRAW_CONSTRUCTION_LINES 0   /* set to 1 to draw additional construction lines for billiard */
#define PERIODIC_BC 0       /* set to 1 to enforce periodic boundary conditions when drawing particles */

#define RESAMPLE 0      /* set to 1 if particles should be added when dispersion too large */
#define DEBUG 0         /* draw trajectories, for debugging purposes */

/* Simulation parameters */

#define NPART 20000     /* number of particles */
#define NPARTMAX 100000	/* maximal number of particles after resampling */
#define LMAX 0.01       /* minimal segment length triggering resampling */ 
#define DMIN 0.02       /* minimal distance to boundary for triggering resampling */ 
#define CYCLE 1         /* set to 1 for closed curve (start in all directions) */
#define SHOWTRAILS 1    /* set to 1 to keep trails of the particles */
#define SHOWZOOM 0      /* set to 1 to show zoom on specific area */
#define PRINT_PARTICLE_NUMBER 0 /* set to 1 to print number of particles */
#define PRINT_LEFT_RIGHT_PARTICLE_NUMBER 1 /* set to 1 to print number of particles on left and right side */
#define PRINT_COLLISION_NUMBER 0 /* set to 1 to print number of collisions */
#define TEST_ACTIVE 1   /* set to 1 to test whether particle is in billiard */

#define TEST_INITIAL_COND 0     /* set to 1 to allow only initial conditions that pass a test */

#define NSTEPS 8000      /* number of frames of movie */
#define TIME 1500        /* time between movie frames, for fluidity of real-time simulation */ 
#define DPHI 0.00002     /* integration step */
#define NVID 25          /* number of iterations between images displayed on screen */
#define END_FRAMES 25    /* number of still frames at the end of the movie */

/* Decreasing TIME accelerates the animation and the movie                               */
/* For constant speed of movie, TIME*DPHI should be kept constant                        */
/* However, increasing DPHI too much deterioriates quality of simulation                 */
/* NVID tells how often a picture is drawn in the animation, increase it for faster anim */
/* For a good quality movie, take for instance TIME = 400, DPHI = 0.00005, NVID = 100    */

/* Colors and other graphical parameters */

#define COLOR_PALETTE 0     /* Color palette, see list in global_pdes.c  */

#define NCOLORS 500      /* number of colors */
#define COLORSHIFT 0     /* hue of initial color */ 
#define RAINBOW_COLOR 0  /* set to 1 to use different colors for all particles */
#define FLOWER_COLOR 0   /* set to 1 to adapt initial colors to flower billiard (tracks vs core) */
#define NSEG 100         /* number of segments of boundary */
#define LENGTH 0.04       /* length of velocity vectors */
#define BILLIARD_WIDTH 2    /* width of billiard */
#define PARTICLE_WIDTH 2    /* width of particles */
#define FRONT_WIDTH 3       /* width of wave front */

#define BLACK 1             /* set to 1 for black background */
#define COLOR_OUTSIDE 0     /* set to 1 for colored outside */ 
#define OUTER_COLOR 270.0   /* color outside billiard */
#define PAINT_INT 0         /* set to 1 to paint interior in other color (for polygon/Reuleaux) */
#define PAINT_EXT 1         /* set to 1 to paint exterior */

#define PAUSE 1000       /* number of frames after which to pause */
#define PSLEEP 2         /* sleep time during pause */
#define SLEEP1  1        /* initial sleeping time */
#define SLEEP2  1       /* final sleeping time */

#define NXMAZE 30      /* width of maze */
#define NYMAZE 30      /* height of maze */
#define MAZE_MAX_NGBH 4     /* max number of neighbours of maze cell */
#define RAND_SHIFT 2        /* seed of random number generator */
#define MAZE_XSHIFT 0.0     /* horizontal shift of maze */
#define MAZE_RANDOM_FACTOR 0.1     /* randomization factor for S_MAZE_RANDOM */

```

### 16 December 22 - A refracting star ###

**Program:** `wave_billiard.c` 

**Initial condition in function `animation()`:** 
```
    init_circular_wave(0.0, 0.0, phi, psi, xy_in);
    
    if (ALTERNATE_OSCILLATING_SOURCE) sign = -sign;
    add_circular_wave(sign, 0.0, 0.0, phi, psi, xy_in);
```

```
#define MOVIE 1         /* set to 1 to generate movie */
#define DOUBLE_MOVIE 1  /* set to 1 to produce movies for wave height and energy simultaneously */
#define SAVE_MEMORY 0   /* set to 1 to save memory when writing tiff images */

/* General geometrical parameters */

#define WINWIDTH 	1920  /* window width */
#define WINHEIGHT 	1000  /* window height */
#define NX 3840          /* number of grid points on x axis */
#define NY 2000          /* number of grid points on y axis */

#define XMIN -2.0
#define XMAX 2.0	/* x interval  */
#define YMIN -1.041666667
#define YMAX 1.041666667	/* y interval for 9/16 aspect ratio */

#define HIGHRES 1       /* set to 1 if resolution of grid is double that of displayed image */

#define JULIA_SCALE 1.0 /* scaling for Julia sets */

/* Choice of the billiard table */

#define B_DOMAIN 42        /* choice of domain shape, see list in global_pdes.c */

#define CIRCLE_PATTERN 1   /* pattern of circles or polygons, see list in global_pdes.c */

#define COMPARISON 0        /* set to 1 to compare two different patterns (beta) */
#define B_DOMAIN_B 20       /* second domain shape, for comparisons */
#define CIRCLE_PATTERN_B 0  /* second pattern of circles or polygons */

#define P_PERCOL 0.25       /* probability of having a circle in C_RAND_PERCOL arrangement */
#define NPOISSON 300        /* number of points for Poisson C_RAND_POISSON arrangement */
#define RANDOM_POLY_ANGLE 1 /* set to 1 to randomize angle of polygons */

#define LAMBDA 0.7	    /* parameter controlling the dimensions of domain */
#define MU 0.2              /* parameter controlling the dimensions of domain */
#define NPOLY 24            /* number of sides of polygon */
#define APOLY 0.0           /* angle by which to turn polygon, in units of Pi/2 */ 
#define MDEPTH 6            /* depth of computation of Menger gasket */
#define MRATIO 3            /* ratio defining Menger gasket */
#define MANDELLEVEL 1000    /* iteration level for Mandelbrot set */
#define MANDELLIMIT 10.0    /* limit value for approximation of Mandelbrot set */
#define FOCI 1              /* set to 1 to draw focal points of ellipse */
#define NGRIDX 12           /* number of grid point for grid of disks */
#define NGRIDY 15           /* number of grid point for grid of disks */

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
#define OSCILLATION_SCHEDULE 1  /* oscillation schedule, see list in global_pdes.c */

#define OMEGA 0.0005       /* frequency of periodic excitation */
#define AMPLITUDE 0.8      /* amplitude of periodic excitation */ 
#define ACHIRP 0.25        /* acceleration coefficient in chirp */
#define DAMPING 0.0        /* damping of periodic excitation */
#define COURANT 0.05        /* Courant number */
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

#define ADD_OSCILLATING_SOURCE 1        /* set to 1 to add an oscillating wave source */
#define OSCILLATING_SOURCE_PERIOD 50    /* period of oscillating source */
#define ALTERNATE_OSCILLATING_SOURCE 1  /* set to 1 to alternate sign of oscillating source */

/* Boundary conditions, see list in global_pdes.c  */

#define B_COND 2

/* Parameters for length and speed of simulation */

#define NSTEPS 2400       /* number of frames of movie */
#define NVID 12          /* number of iterations between images displayed on screen */
#define NSEG 1000         /* number of segments of boundary */
#define INITIAL_TIME 0      /* time after which to start saving frames */
#define BOUNDARY_WIDTH 2    /* width of billiard boundary */
#define PRINT_SPEED 0       /* print speed of moving source */

#define PAUSE 200       /* number of frames after which to pause */
#define PSLEEP 1         /* sleep time during pause */
#define SLEEP1  1        /* initial sleeping time */
#define SLEEP2  1        /* final sleeping time */
#define MID_FRAMES 20    /* number of still frames between parts of two-part movie */
#define END_FRAMES 100    /* number of still frames at end of movie */
#define FADE 1           /* set to 1 to fade at end of movie */

/* Parameters of initial condition */

#define INITIAL_AMP 0.35         /* amplitude of initial condition */
#define INITIAL_VARIANCE 0.0002  /* variance of initial condition */
#define INITIAL_WAVELENGTH  0.01  /* wavelength of initial condition */

/* Plot type, see list in global_pdes.c  */

#define PLOT 6

#define PLOT_B 7        /* plot type for second movie */

/* Color schemes */

#define COLOR_PALETTE 13     /* Color palette, see list in global_pdes.c  */
#define COLOR_PALETTE_B 18     /* Color palette, see list in global_pdes.c  */

#define BLACK 1          /* background */

#define COLOR_SCHEME 3   /* choice of color scheme, see list in global_pdes.c  */

#define SCALE 0          /* set to 1 to adjust color scheme to variance of field */
#define SLOPE 1.0        /* sensitivity of color on wave amplitude */
#define PHASE_FACTOR 1.0       /* factor in computation of phase in color scheme P_3D_PHASE */
#define PHASE_SHIFT 0.0      /* shift of phase in color scheme P_3D_PHASE */
#define ATTENUATION 0.0  /* exponential attenuation coefficient of contrast with time */
#define E_SCALE 100.0     /* scaling factor for energy representation */
#define LOG_SCALE 0.25     /* scaling factor for energy log representation */
#define LOG_SHIFT 0.0     /* shift of colors on log scale */
#define FLUX_SCALE 1.0e3    /* scaling factor for enegy flux represtnation */
#define RESCALE_COLOR_IN_CENTER 0   /* set to 1 to decrease color intentiy in the center (for wave escaping ring) */

#define COLORHUE 260     /* initial hue of water color for scheme C_LUM */
#define COLORDRIFT 0.0   /* how much the color hue drifts during the whole simulation */
#define LUMMEAN 0.5      /* amplitude of luminosity variation for scheme C_LUM */
#define LUMAMP 0.3       /* amplitude of luminosity variation for scheme C_LUM */
#define HUEMEAN 180.0    /* mean value of hue for color scheme C_HUE */
#define HUEAMP -180.0      /* amplitude of variation of hue for color scheme C_HUE */

#define DRAW_COLOR_SCHEME 1    /* set to 1 to plot the color scheme */
#define COLORBAR_RANGE 1.5     /* scale of color scheme bar */
#define COLORBAR_RANGE_B 2.5    /* scale of color scheme bar for 2nd part */
#define ROTATE_COLOR_SCHEME 0   /* set to 1 to draw color scheme horizontally */

#define SAVE_TIME_SERIES 0      /* set to 1 to save wave time series at a point */

#define NXMAZE 7      /* width of maze */
#define NYMAZE 7      /* height of maze */
#define MAZE_MAX_NGBH 4     /* max number of neighbours of maze cell */
#define RAND_SHIFT 24        /* seed of random number generator */
#define MAZE_XSHIFT 0.0     /* horizontal shift of maze */

/* for compatibility with sub_wave and sub_maze */
#define ADD_POTENTIAL 0
#define POT_MAZE 7
#define POTENTIAL 0
/* end of constants only used by sub_wave and sub_maze */


/* For debugging purposes only */
#define FLOOR 0         /* set to 1 to limit wave amplitude to VMAX */
#define VMAX 10.0       /* max value of wave amplitude */

```

### 15 December 22 - Pressure of an incompressible Euler flow with a circular obstacle ###

**Program:** `rde.c` 

**Initial condition in function `animation()`:** `init_shear_flow(1.0, 0.02, 0.15, 1, 1, 0.0, phi, xy_in);`

```
#define MOVIE 1         /* set to 1 to generate movie */
#define DOUBLE_MOVIE 0  /* set to 1 to produce movies for wave height and energy simultaneously */
#define SAVE_MEMORY 0           /* set to 1 to save memory when writing tiff images */

/* General geometrical parameters */

#define WINWIDTH 	1000  /* window width */
#define WINHEIGHT 	1000  /* window height */
#define NX 500          /* number of grid points on x axis */
#define NY 500          /* number of grid points on y axis */

#define XMIN -1.0
#define XMAX 1.0	/* x interval */
#define YMIN -1.0
#define YMAX 1.0	/* y interval for 9/16 aspect ratio */

/* Choice of simulated equation */

#define RDE_EQUATION 6  /* choice of reaction term, see list in global_3d.c */
#define NFIELDS 3       /* number of fields in reaction-diffusion equation */
#define NLAPLACIANS 1   /* number of fields for which to compute Laplacian */

#define ADD_POTENTIAL 0 /* set to 1 to add a potential (for Schrodinger equation) */
#define ADD_MAGNETIC_FIELD 0    /* set to 1 to add a magnetic field (for Schrodinger equation) - then set POTENTIAL 1 */
#define POTENTIAL 7     /* type of potential or vector potential, see list in global_3d.c  */

#define ANTISYMMETRIZE_WAVE_FCT 0   /* set tot 1 to make wave function antisymmetric */

#define JULIA_SCALE 0.5 /* scaling for Julia sets */

/* Choice of the billiard table */

#define B_DOMAIN 3          /* choice of domain shape, see list in global_pdes.c  */

#define CIRCLE_PATTERN 99    /* pattern of circles, see list in global_pdes.c */

#define P_PERCOL 0.25       /* probability of having a circle in C_RAND_PERCOL arrangement */
#define NPOISSON 300        /* number of points for Poisson C_RAND_POISSON arrangement */
#define RANDOM_POLY_ANGLE 0 /* set to 1 to randomize angle of polygons */

#define LAMBDA 0.3	    /* parameter controlling the dimensions of domain */
#define MU 0.15	            /* parameter controlling the dimensions of domain */
#define NPOLY 5             /* number of sides of polygon */
#define APOLY 2.0          /* angle by which to turn polygon, in units of Pi/2 */
#define MDEPTH 7            /* depth of computation of Menger gasket */
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

#define DT 0.0000012

#define VISCOSITY 2.0

#define RPSA 0.75         /* parameter in Rock-Paper-Scissors-type interaction */
#define RPSLZB 0.75       /* second parameter in Rock-Paper-Scissors-Lizard-Spock type interaction */

#define EPSILON 0.8     /* time scale separation */
#define DELTA 0.1       /* time scale separation */
#define FHNA 1.0        /* parameter in FHN equation */
#define FHNC -0.01      /* parameter in FHN equation */
#define K_HARMONIC 1.0  /* spring constant of harmonic potential */
#define K_COULOMB 0.5   /* constant in Coulomb potential */
#define V_MAZE 0.4      /* potential in walls of maze */
#define BZQ 0.0008      /* parameter in BZ equation */
#define BZF 1.2         /* parameter in BZ equation */
#define B_FIELD 10.0    /* magnetic field */
#define AB_RADIUS 0.2   /* radius of region with magnetic field for Aharonov-Bohm effect */
#define K_EULER 50.0    /* constant in stream function integration of Euler equation */

#define SMOOTHEN_VORTICITY 1    /* set to 1 to smoothen vorticity field in Euler equation */
#define SMOOTHEN_PERIOD 10      /* period between smoothenings */
#define SMOOTH_FACTOR 0.03      /* factor by which to smoothen */

#define ADD_TRACERS 1    /* set to 1 to add tracer particles (for Euler equations) */
#define N_TRACERS 1000    /* number of tracer particles */

#define T_OUT 2.0       /* outside temperature */
#define T_IN 0.0        /* inside temperature */
#define SPEED 0.0       /* speed of drift to the right */

#define ADD_NOISE 0     /* set to 1 to add noise, set to 2 to add noise in right half */
#define NOISE_INTENSITY 0.005      /* noise intensity */
#define CHANGE_NOISE 1      /* set to 1 to increase noise intensity */
#define NOISE_FACTOR 40.0   /* factor by which to increase noise intensity */
#define NOISE_INITIAL_TIME 100  /* initial time during which noise remains constant */

#define CHANGE_VISCOSITY 0      /* set to 1 to change the viscosity in the course of the simulation */
#define ADJUST_INTSTEP 0       /* set to 1 to decrease integration step when viscosity increases */
#define VISCOSITY_INITIAL_TIME 10  /* initial time during which viscosity remains constant */
#define VISCOSITY_FACTOR 100.0   /* factor by which to change viscosity */
#define VISCOSITY_MAX 2.0        /* max value of viscosity beyond which NVID is increased and integration step is decrase, 
                                    for numerical stability */
                                        
#define CHANGE_RPSLZB 0         /* set to 1 to change second parameter in Rock-Paper-Scissors-Lizard-Spock equation */
#define RPSLZB_CHANGE 0.75      /* factor by which to rpslzb parameter */
#define RPSLZB_INITIAL_TIME 0   /* initial time during which rpslzb remains constant */
#define RPSLZB_FINAL_TIME 500   /* final time during which rpslzb remains constant */
                                      
/* Boundary conditions, see list in global_pdes.c  */

#define B_COND 1

#define EULER_GRADIENT_YSHIFT 0.0    /* y-shift in computation of gradient in Euler equation */

/* Parameters for length and speed of simulation */

#define NSTEPS 1000       /* number of frames of movie */
#define NVID 200          /* number of iterations between images displayed on screen */
#define ACCELERATION_FACTOR 1.0 /* factor by which to increase NVID in course of simulation */
#define DT_ACCELERATION_FACTOR 1.0 /* factor by which to increase time step in course of simulation  */
#define MAX_DT 0.024     /* maximal value of integration step */
#define NSEG 100         /* number of segments of boundary */
#define BOUNDARY_WIDTH 2    /* width of billiard boundary */

#define PAUSE 100       /* number of frames after which to pause */
#define PSLEEP 2         /* sleep time during pause */
#define SLEEP1  2        /* initial sleeping time */
#define SLEEP2  1        /* final sleeping time */
#define INITIAL_TIME 0  /* initial still time */
#define MID_FRAMES 50    /* number of still frames between parts of two-part movie */
#define END_FRAMES 50    /* number of still frames at end of movie */
#define FADE 1           /* set to 1 to fade at end of movie */

/* Visualisation */

#define PLOT_3D 0    /* controls whether plot is 2D or 3D */

#define ROTATE_VIEW 0       /* set to 1 to rotate position of observer */
#define ROTATE_ANGLE 360.0  /* total angle of rotation during simulation */

#define DRAW_PERIODICISED 1     /* set to 1 to repeat wave periodically in x and y directions */

/* Plot type - color scheme */

#define CPLOT 54
#define CPLOT_B 51

/* Plot type - height of 3D plot */

#define ZPLOT 52     /* z coordinate in 3D plot */
// #define ZPLOT 32     /* z coordinate in 3D plot */
#define ZPLOT_B 51    /* z coordinate in second 3D plot */

#define AMPLITUDE_HIGH_RES 1    /* set to 1 to increase resolution of P_3D_AMPLITUDE plot */
#define SHADE_3D 1              /* set to 1 to change luminosity according to normal vector */
#define NON_DIRICHLET_BC 0      /* set to 1 to draw only facets in domain, if field is not zero on boundary */
#define WRAP_ANGLE 1            /* experimental: wrap angle to [0, 2Pi) for interpolation in angle schemes */
#define FADE_IN_OBSTACLE 0      /* set to 1 to fade color inside obstacles */
#define DRAW_OUTSIDE_GRAY 0     /* experimental - draw outside of billiard in gray */
#define ADD_POTENTIAL_TO_Z 1    /* set to 1 to add the external potential to z-coordinate of plot */
#define ADD_POT_CONSTANT 0.35   /* constant in front of added potential */

#define PLOT_SCALE_ENERGY 0.05      /* vertical scaling in energy plot */

#define PRINT_TIME 0        /* set to 1 to print running time */
#define PRINT_VISCOSITY 0   /* set to 1 to print viscosity */
#define PRINT_RPSLZB 0      /* set to 1 to print rpslzb parameter */
#define PRINT_PROBABILITIES 0   /* set to 1 to print probabilities (for Ehrenfest urn configuration) */
#define PRINT_NOISE 0       /* set to 1 to print noise intensity */

#define DRAW_FIELD_LINES 0  /* set to 1 to draw field lines */
#define FIELD_LINE_WIDTH 1  /* width of field lines */
#define N_FIELD_LINES 120   /* number of field lines */
#define FIELD_LINE_FACTOR 120 /* factor controlling precision when computing origin of field lines */
#define DRAW_BILLIARD 1     /* set to 1 to draw boundary */
#define DRAW_BILLIARD_FRONT 1     /* set to 1 to draw boundary */
#define FILL_BILLIARD_COMPLEMENT 1  /* set to 1 to fill complement of billiard (for certain shapes only) */

/* 3D representation */

#define REPRESENTATION_3D 1     /* choice of 3D representation */ 

#define REP_AXO_3D 0        /* linear projection (axonometry) */
#define REP_PROJ_3D 1       /* projection on plane orthogonal to observer line of sight */

/* Color schemes, see list in global_pdes.c  */

#define COLOR_PALETTE 10       /* Color palette, see list in global_pdes.c  */
#define COLOR_PALETTE_B 13     /* Color palette, see list in global_pdes.c  */

#define BLACK 1          /* black background */

#define COLOR_SCHEME 3   /* choice of color scheme */

#define COLOR_PHASE_SHIFT 0.0   /* phase shift of color scheme, in units of Pi */

#define SCALE 0          /* set to 1 to adjust color scheme to variance of field */
#define SLOPE 1.0        /* sensitivity of color on wave amplitude */
#define VSCALE_AMPLITUDE 1.5      /* additional scaling factor for color scheme P_3D_AMPLITUDE */
#define ATTENUATION 0.0  /* exponential attenuation coefficient of contrast with time */
#define CURL_SCALE 0.000015   /* scaling factor for curl representation */
#define RESCALE_COLOR_IN_CENTER 0   /* set to 1 to decrease color intentiy in the center (for wave escaping ring) */
#define SLOPE_SCHROD_LUM 50.0       /* sensitivity of luminosity on module, for color scheme Z_ARGUMENT */
#define MIN_SCHROD_LUM 0.2       /* minimal luminosity in color scheme Z_ARGUMENT*/
#define VSCALE_PRESSURE 0.35      /* additional scaling factor for color scheme Z_EULER_PRESSURE */
#define PRESSURE_SHIFT 250.0        /* shift for color scheme Z_EULER_PRESSURE */
#define PRESSURE_LOG_SHIFT -1.8     /* shift for color scheme Z_EULER_PRESSURE */

#define COLORHUE 260     /* initial hue of water color for scheme C_LUM */
#define COLORDRIFT 0.0   /* how much the color hue drifts during the whole simulation */
#define LUMMEAN 0.5      /* amplitude of luminosity variation for scheme C_LUM */
#define LUMAMP 0.3       /* amplitude of luminosity variation for scheme C_LUM */
#define HUEMEAN 359.0    /* mean value of hue for color scheme C_HUE */
#define HUEAMP -359.0    /* amplitude of variation of hue for color scheme C_HUE */
#define E_SCALE 100.0    /* scaling factor for energy representation */
#define FLUX_SCALE 100.0 /* scaling factor for energy representation */
#define LOG_SCALE 0.5    /* scaling factor for energy log representation */
#define LOG_SHIFT 1.0 
#define LOG_MIN 1.0e-3   /* floor value for log vorticity plot */

#define NXMAZE 7      /* width of maze */
#define NYMAZE 7      /* height of maze */
#define MAZE_MAX_NGBH 4     /* max number of neighbours of maze cell */
#define RAND_SHIFT 0        /* seed of random number generator */
#define MAZE_XSHIFT 0.0     /* horizontal shift of maze */

#define DRAW_COLOR_SCHEME 0     /* set to 1 to plot the color scheme */
#define COLORBAR_RANGE 3.0      /* scale of color scheme bar */
#define COLORBAR_RANGE_B 3.0   /* scale of color scheme bar for 2nd part */
#define ROTATE_COLOR_SCHEME 0   /* set to 1 to draw color scheme horizontally */

/* only for compatibility with wave_common.c */
#define TWOSPEEDS 0          /* set to 1 to replace hardcore boundary by medium with different speed */
#define OMEGA 0.005        /* frequency of periodic excitation */
#define COURANT 0.08       /* Courant number */
#define COURANTB 0.03      /* Courant number in medium B */
#define INITIAL_AMP 0.5         /* amplitude of initial condition */
#define INITIAL_VARIANCE 0.0002  /* variance of initial condition */
#define INITIAL_WAVELENGTH  0.1  /* wavelength of initial condition */
#define VSCALE_ENERGY 200.0       /* additional scaling factor for color scheme P_3D_ENERGY */
#define PHASE_FACTOR 20.0       /* factor in computation of phase in color scheme P_3D_PHASE */
#define PHASE_SHIFT 0.0      /* shift of phase in color scheme P_3D_PHASE */
#define OSCILLATION_SCHEDULE 0  /* oscillation schedule, see list in global_pdes.c */
#define AMPLITUDE 0.8      /* amplitude of periodic excitation */ 
#define ACHIRP 0.2        /* acceleration coefficient in chirp */
#define DAMPING 0.0        /* damping of periodic excitation */
#define COMPARISON 0        /* set to 1 to compare two different patterns (beta) */
#define B_DOMAIN_B 20       /* second domain shape, for comparisons */
#define CIRCLE_PATTERN_B 0  /* second pattern of circles or polygons */
/* end of constants added only for compatibility with wave_common.c */

double u_3d[2] = {0.75, -0.45};     /* projections of basis vectors for REP_AXO_3D representation */
double v_3d[2] = {-0.75, -0.45};
double w_3d[2] = {0.0, 0.015};
double light[3] = {0.816496581, -0.40824829, 0.40824829};      /* vector of "light" direction for P_3D_ANGLE color scheme */
double observer[3] = {8.0, 8.0, 8.0};    /* location of observer for REP_PROJ_3D representation */ 
int reset_view = 0;         /* switch to reset 3D view parameters (for option ROTATE_VIEW) */

#define Z_SCALING_FACTOR 0.08  /* overall scaling factor of z axis for REP_PROJ_3D representation */
#define XY_SCALING_FACTOR 1.7  /* overall scaling factor for on-screen (x,y) coordinates after projection */
#define ZMAX_FACTOR 1.0        /* max value of z coordinate for REP_PROJ_3D representation */
#define XSHIFT_3D -0.1         /* overall x shift for REP_PROJ_3D representation */
#define YSHIFT_3D 0.1          /* overall y shift for REP_PROJ_3D representation */
#define BORDER_PADDING 0       /* distance from boundary at which to plot points, to avoid boundary effects due to gradient */

/* For debugging purposes only */
#define FLOOR 1         /* set to 1 to limit wave amplitude to VMAX */
#define VMAX 1000.0        /* max value of wave amplitude */
#define TEST_GRADIENT 0 /* print norm squared of gradient */

```

### 14 December 22 - Mixing of sticky square grains of sand during a cave-in ###

**Program:** `lennardjones.c` 

```
#define MOVIE 1         /* set to 1 to generate movie */
#define DOUBLE_MOVIE 0  /* set to 1 to produce movies for wave height and energy simultaneously */

#define TIME_LAPSE 1     /* set to 1 to add a time-lapse movie at the end */
                         /* so far incompatible with double movie */
#define TIME_LAPSE_FACTOR 3    /* factor of time-lapse movie */
#define TIME_LAPSE_FIRST 1  /* set to 1 to show time-lapse version first */

/* General geometrical parameters */

#define WINWIDTH 	1280  /* window width */
#define WINHEIGHT 	720   /* window height */

#define XMIN -2.0
#define XMAX 2.0	/* x interval */
#define YMIN -1.125
#define YMAX 1.125	/* y interval for 9/16 aspect ratio */

#define INITXMIN -1.95
#define INITXMAX 1.95	/* x interval for initial condition */
#define INITYMIN 0.05
#define INITYMAX 1.1	/* y interval for initial condition */

#define BCXMIN -2.0
#define BCXMAX 2.0	/* x interval for boundary condition */
#define BCYMIN -1.125
#define BCYMAX 1.125	/* y interval for boundary condition */

#define OBSXMIN -2.0
#define OBSXMAX 2.0     /* x interval for motion of obstacle */

#define CIRCLE_PATTERN 8  /* pattern of circles, see list in global_ljones.c */

#define ADD_FIXED_OBSTACLES 0   /* set to 1 do add fixed circular obstacles */
#define OBSTACLE_PATTERN 3  /* pattern of obstacles, see list in global_ljones.c */

#define ADD_FIXED_SEGMENTS 1    /* set to 1 to add fixed segments as obstacles */
#define SEGMENT_PATTERN 18     /* pattern of repelling segments, see list in global_ljones.c */
#define ROCKET_SHAPE 2        /* shape of rocket combustion chamber, see list in global_ljones.c */
#define ROCKET_SHAPE_B 2      /* shape of second rocket */
#define NOZZLE_SHAPE 2        /* shape of nozzle, see list in global_ljones.c */
#define NOZZLE_SHAPE_B 4      /* shape of nozzle for second rocket, see list in global_ljones.c */

#define TWO_TYPES 0         /* set to 1 to have two types of particles */
#define TPYE_PROPORTION 0.5 /* proportion of particles of first type */
#define SYMMETRIZE_FORCE 1  /* set to 1 to symmetrize two-particle interaction, only needed if particles are not all the same */
#define CENTER_PX 1         /* set to 1 to center horizontal momentum */
#define CENTER_PY 0         /* set to 1 to center vertical momentum */
#define CENTER_PANGLE 0     /* set to 1 to center angular momentum */

#define INTERACTION 2       /* particle interaction, see list in global_ljones.c */
#define INTERACTION_B 1     /* particle interaction for second type of particle, see list in global_ljones.c */
#define SPIN_INTER_FREQUENCY 4.0 /* angular frequency of spin-spin interaction */
#define SPIN_INTER_FREQUENCY_B 2.0 /* angular frequency of spin-spin interaction for second particle type */

#define P_PERCOL 0.25       /* probability of having a circle in C_RAND_PERCOL arrangement */
#define NPOISSON 100        /* number of points for Poisson C_RAND_POISSON arrangement */
#define PDISC_DISTANCE 2.5  /* minimal distance in Poisson disc process, controls density of particles */
#define PDISC_CANDIDATES 100 /* number of candidates in construction of Poisson disc process */
#define RANDOM_POLY_ANGLE 0 /* set to 1 to randomize angle of polygons */

#define LAMBDA 0.8	    /* parameter controlling the dimensions of domain */
#define MU 0.014 	    /* parameter controlling radius of particles */
#define MU_B 0.018          /* parameter controlling radius of particles of second type */
#define NPOLY 25            /* number of sides of polygon */
#define APOLY 0.666666666   /* angle by which to turn polygon, in units of Pi/2 */ 
#define MDEPTH 4            /* depth of computation of Menger gasket */
#define MRATIO 3            /* ratio defining Menger gasket */
#define MANDELLEVEL 1000    /* iteration level for Mandelbrot set */
#define MANDELLIMIT 10.0    /* limit value for approximation of Mandelbrot set */
#define FOCI 1              /* set to 1 to draw focal points of ellipse */
#define NGRIDX 90           /* number of grid point for grid of disks */
#define NGRIDY 100           /* number of grid point for grid of disks */
#define EHRENFEST_RADIUS 0.9    /* radius of container for Ehrenfest urn configuration */
#define EHRENFEST_WIDTH 0.035     /* width of tube for Ehrenfest urn configuration */
#define TWO_CIRCLES_RADIUS_RATIO 0.8    /* ratio of radii for S_TWO_CIRCLES_EXT segment configuration */
#define DAM_WIDTH 0.05       /* width of dam for S_DAM segment configuration */

#define X_SHOOTER -0.2
#define Y_SHOOTER -0.6
#define X_TARGET 0.4
#define Y_TARGET 0.7        /* shooter and target positions in laser fight */

/* Parameters for length and speed of simulation */

#define NSTEPS 6000      /* number of frames of movie */
#define NVID 150         /* number of iterations between images displayed on screen */
#define NSEG 250         /* number of segments of boundary */
#define INITIAL_TIME 400     /* time after which to start saving frames */
#define OBSTACLE_INITIAL_TIME 200     /* time after which to start moving obstacle */
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

#define PLOT 11
#define PLOT_B 8        /* plot type for second movie */

#define DRAW_BONDS 0    /* set to 1 to draw bonds between neighbours */
#define COLOR_BONDS 1   /* set to 1 to color bonds according to length */
#define FILL_TRIANGLES 0    /* set to 1 to fill triangles between neighbours */
#define ALTITUDE_LINES 0    /* set to 1 to add horizontal lines to show altitude */
#define COLOR_SEG_GROUPS 0  /* set to 1 to collor segment groups differently */

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

#define ENERGY_HUE_MIN 330.0        /* color of original particle */
#define ENERGY_HUE_MAX 50.0         /* color of saturated particle */
#define PARTICLE_HUE_MIN 359.0      /* color of original particle */
#define PARTICLE_HUE_MAX 0.0        /* color of saturated particle */
#define PARTICLE_EMAX 2.0e2         /* energy of particle with hottest color */
#define HUE_TYPE0 70.0     /* hue of particles of type 0 */
#define HUE_TYPE1 280.0      /* hue of particles of type 1 */
#define HUE_TYPE2 .0      /* hue of particles of type 2 */
#define HUE_TYPE3 210.0     /* hue of particles of type 3 */

#define RANDOM_RADIUS 0     /* set to 1 for random circle radius */
#define DT_PARTICLE 4.5e-6    /* time step for particle displacement */
#define KREPEL 0.5          /* constant in repelling force between particles */
#define EQUILIBRIUM_DIST 3.5    /* Lennard-Jones equilibrium distance */
#define EQUILIBRIUM_DIST_B 3.5  /* Lennard-Jones equilibrium distance for second type of particle */
#define REPEL_RADIUS 15.0    /* radius in which repelling force acts (in units of particle radius) */
#define DAMPING 400.0          /* damping coefficient of particles */
#define INITIAL_DAMPING 750.0  /* damping coefficient of particles during initial phase */
#define PARTICLE_MASS 1.0    /* mass of particle of radius MU */
#define PARTICLE_MASS_B 1.0   /* mass of particle of radius MU */
#define PARTICLE_INERTIA_MOMENT 0.02       /* moment of inertia of particle */
#define PARTICLE_INERTIA_MOMENT_B 0.02     /* moment of inertia of second type of particle */
#define V_INITIAL 0.0        /* initial velocity range */
#define OMEGA_INITIAL 10.0        /* initial angular velocity range */

#define THERMOSTAT 0        /* set to 1 to switch on thermostat */
#define VARY_THERMOSTAT 0   /* set to 1 for time-dependent thermostat schedule */
#define SIGMA 5.0           /* noise intensity in thermostat */
#define BETA 0.02          /* initial inverse temperature */
#define MU_XI 0.01           /* friction constant in thermostat */
#define KSPRING_BOUNDARY 1.0e7    /* confining harmonic potential outside simulation region */
#define KSPRING_OBSTACLE 1.0e11    /* harmonic potential of obstacles */
#define NBH_DIST_FACTOR 7.5        /* radius in which to count neighbours */
#define GRAVITY 4000.0             /* gravity acting on all particles */
#define GRAVITY_X 0.0        /* horizontal gravity acting on all particles */
#define INCREASE_GRAVITY 0     /* set to 1 to increase gravity during the simulation */
#define GRAVITY_SCHEDULE 2     /* type of gravity schedule, see list in global_ljones.c */
#define GRAVITY_FACTOR 100.0    /* factor by which to increase gravity */
#define GRAVITY_INITIAL_TIME 200    /* time at start of simulation with constant gravity */
#define GRAVITY_RESTORE_TIME 700    /* time at end of simulation with gravity restored to initial value */

#define ROTATION 1           /* set to 1 to include rotation of particles */
#define COUPLE_ANGLE_TO_THERMOSTAT 0    /* set to 1 to couple angular degrees of freedom to thermostat */
#define DIMENSION_FACTOR 1.0  /* scaling factor taking into account number of degrees of freedom */  
#define KTORQUE 20.0          /* force constant in angular dynamics */
#define KTORQUE_B 10.0        /* force constant in angular dynamics */
#define KTORQUE_DIFF 150.0    /* force constant in angular dynamics for different particles */
#define DRAW_SPIN 0           /* set to 1 to draw spin vectors of particles */
#define DRAW_SPIN_B 0         /* set to 1 to draw spin vectors of particles */
#define DRAW_CROSS 1          /* set to 1 to draw cross on particles of second type */
#define SPIN_RANGE 7.0       /* range of spin-spin interaction */
#define SPIN_RANGE_B 5.0     /* range of spin-spin interaction for second type of particle */
#define QUADRUPOLE_RATIO 0.6  /* anisotropy in quadrupole potential */ 

#define INCREASE_BETA 0  /* set to 1 to increase BETA during simulation */
#define BETA_FACTOR 0.025   /* factor by which to change BETA during simulation */
#define N_TOSCILLATIONS 1.5   /* number of temperature oscillations in BETA schedule */
#define NO_OSCILLATION 1        /* set to 1 to have exponential BETA change only */
#define MIDDLE_CONSTANT_PHASE 2000   /* final phase in which temperature is constant */
#define FINAL_DECREASE_PHASE 1300    /* final phase in which temperature decreases */ 
#define FINAL_CONSTANT_PHASE -1     /* final phase in which temperature is constant */

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
#define PARTIAL_THERMO_REGION 1     /* region for partial thermostat coupling (see list in global_ljones.c) */
#define PARTIAL_THERMO_SHIFT 0.2    /* distance from obstacle at the right of which particles are coupled to thermostat */
#define PARTIAL_THERMO_WIDTH 0.5    /* vertical size of partial thermostat coupling */
#define PARTIAL_THERMO_HEIGHT 0.2   /* vertical size of partial thermostat coupling */

#define INCREASE_KREPEL 0   /* set to 1 to increase KREPEL during simulation */
#define KREPEL_FACTOR 1000.0   /* factor by which to change KREPEL during simulation */

#define PART_AT_BOTTOM 0     /* set to 1 to include "seed" particles at bottom */
#define MASS_PART_BOTTOM 10000.0 /* mass of particles at bottom */
#define NPART_BOTTOM 100     /* number of particles at the bottom */

#define ADD_PARTICLES 0    /* set to 1 to add particles */
#define ADD_TIME 0       /* time at which to add first particle */
#define ADD_PERIOD 10000       /* time interval between adding further particles */
#define FINAL_NOADD_PERIOD 200  /* final period where no particles are added */
#define SAFETY_FACTOR 2.0  /* no particles are added at distance less than MU*SAFETY_FACTOR of other particles */

#define TRACER_PARTICLE 0   /* set to 1 to have a tracer particle */
#define N_TRACER_PARTICLES 3    /* number of tracer particles */
#define TRAJECTORY_LENGTH 8000   /* length of recorded trajectory */
#define TRACER_PARTICLE_MASS 4.0    /* relative mass of tracer particle */
#define TRAJECTORY_WIDTH 3      /* width of tracer particle trajectory */

#define ROTATE_BOUNDARY 0           /* set to 1 to rotate the repelling segments */
#define SMOOTH_ROTATION 1           /* set to 1 to update segments at each time step (rather than at each movie frame) */
#define PERIOD_ROTATE_BOUNDARY 1000  /* period of rotating boundary */
#define ROTATE_INITIAL_TIME 0       /* initial time without rotation */
#define ROTATE_FINAL_TIME 100       /* final time without rotation */
#define ROTATE_CHANGE_TIME 0.33     /* relative duration of acceleration/deceleration phases */
#define OMEGAMAX 100.0              /* maximal rotation speed */
#define PRINT_OMEGA 0               /* set to 1 to print angular speed */
#define PRINT_PARTICLE_SPEEDS 0     /* set to 1 to print average speeds/momenta of particles */
#define PRINT_SEGMENTS_SPEEDS 1     /* set to 1 to print velocity of moving segments */

#define MOVE_BOUNDARY 0        /* set to 1 to move repelling segments, due to force from particles */
#define SEGMENTS_MASS 40.0     /* mass of collection of segments */
#define DEACTIVATE_SEGMENT 1    /* set to 1 to deactivate last segment after a certain time */
#define SEGMENT_DEACTIVATION_TIME 200   /* time at which to deactivate last segment */
#define RELEASE_ROCKET_AT_DEACTIVATION 1    /* set to 1 to limit segments velocity before segment release */
#define SEGMENTS_X0 1.5        /* initial position of segments */
#define SEGMENTS_Y0 0.0        /* initial position of segments */
#define SEGMENTS_VX0 0.0       /* initial velocity of segments */
#define SEGMENTS_VY0 0.0      /* initial velocity of segments */
#define DAMP_SEGS_AT_NEGATIVE_Y 0   /* set to 1 to dampen segments when y coordinate is negative */

#define MOVE_SEGMENT_GROUPS 1       /* set to 1 to group segments into moving units */
#define SEGMENT_GROUP_MASS 1000.0   /* mass of segment group */
#define SEGMENT_GROUP_I 1000.0      /* moment of inertia of segment group */
#define SEGMENT_GROUP_DAMPING 0.0   /* damping of segment groups */
#define GROUP_REPULSION 1           /* set to 1 for groups of segments to repel each other */
#define KSPRING_GROUPS 1.0e11       /* harmonic potential between segment groups */
#define GROUP_WIDTH 0.05            /* interaction width of groups */
#define GROUP_G_REPEL 1             /* set to 1 to add repulsion between centers of mass of groups */
#define GROUP_G_REPEL_RADIUS 1.2    /* radius within which centers of mass of groups repel each other */
#define TRACK_SEGMENT_GROUPS 1      /* set to 1 for view to track group of segments */
#define TRACK_X_PADDING 2.0         /* distance from x boundary where tracking starts */

#define POSITION_DEPENDENT_TYPE 0   /* set to 1 to make particle type depend on initial position */
#define POSITION_Y_DEPENDENCE 0     /* set to 1 for the separation between particles to be horizontal */
#define PRINT_ENTROPY 0     /* set to 1 to compute entropy */

#define REACTION_DIFFUSION 0    /* set to 1 to simulate a chemical reaction (particles may change type) */
#define RD_TYPES 3              /* number of types in reaction-diffusion equation */
#define REACTION_PROB 0.0045      /* probability controlling reaction term */ 

#define PRINT_PARTICLE_NUMBER 0     /* set to 1 to print total number of particles */
#define PRINT_LEFT 0        /* set to 1 to print certain parameters at the top left instead of right */
#define PLOT_SPEEDS 0       /* set to 1 to add a plot of obstacle speeds (e.g. for rockets) */
#define PLOT_TRAJECTORIES 0     /* set to 1 to add a plot of obstacle trajectories (e.g. for rockets) */
#define VMAX_PLOT_SPEEDS 0.6    /* vertical scale of plot of obstacle speeds */

#define EHRENFEST_COPY 0    /* set to 1 to add equal number of larger particles (for Ehrenfest model) */

#define LID_MASS 1000.0     /* mass of lid for BC_RECTANGLE_LID b.c. */
#define LID_WIDTH 0.1       /* width of lid for BC_RECTANGLE_LID b.c. */
#define WALL_MASS 2000.0    /* mass of wall for BC_RECTANGLE_WALL b.c. */
#define WALL_FRICTION 0.0   /* friction on wall for BC_RECTANGLE_WALL b.c. */
#define WALL_WIDTH 0.1      /* width of wall for BC_RECTANGLE_WALL b.c. */
#define WALL_VMAX 100.0     /* max speed of wall */
#define WALL_TIME 0         /* time during which to keep wall */

#define NXMAZE 10      /* width of maze */
#define NYMAZE 10      /* height of maze */
#define MAZE_MAX_NGBH 4     /* max number of neighbours of maze cell */
#define RAND_SHIFT 200      /* seed of random number generator */
#define MAZE_XSHIFT 0.0     /* horizontal shift of maze */

#define FLOOR_FORCE 1      /* set to 1 to limit force on particle to FMAX */
#define FMAX 1.0e12         /* maximal force */
#define FLOOR_OMEGA 0      /* set to 1 to limit particle momentum to PMAX */
#define PMAX 1000.0        /* maximal force */

#define HASHX 120   /* size of hashgrid in x direction */
#define HASHY 60    /* size of hashgrid in y direction */
#define HASHMAX 100  /* maximal number of particles per hashgrid cell */
#define HASHGRID_PADDING 0.1    /* padding of hashgrid outside simulation window */

#define DRAW_COLOR_SCHEME 0     /* set to 1 to plot the color scheme */
#define COLORBAR_RANGE 8.0    /* scale of color scheme bar */
#define COLORBAR_RANGE_B 12.0    /* scale of color scheme bar for 2nd part */
#define ROTATE_COLOR_SCHEME 0   /* set to 1 to draw color scheme horizontally */

```

### 13 December 22 - Laser beams starting in the center of a 32 by 32 maze ###

**Program:** `particle_billiard.c` 

**Initial condition in function `animation()`:** `init_drop_config(0.0, 0.05, 0.0, DPI, configs);`

```
#define MOVIE 1         /* set to 1 to generate movie */

#define WINWIDTH 	1280  /* window width */
#define WINHEIGHT 	720   /* window height */

#define XMIN -2.0
#define XMAX 2.0	/* x interval */
#define YMIN -1.125
#define YMAX 1.125	/* y interval for 9/16 aspect ratio */

#define SCALING_FACTOR 1.0       /* scaling factor of drawing, needed for flower billiards, otherwise set to 1.0 */

/* Choice of the billiard table, see global_particles.c */

#define B_DOMAIN 30     /* choice of domain shape */

#define CIRCLE_PATTERN 1    /* pattern of circles */
#define POLYLINE_PATTERN 10  /* pattern of polyline */

#define ABSORBING_CIRCLES 0 /* set to 1 for circular scatterers to be absorbing */

#define NMAXCIRCLES 100000     /* total number of circles (must be at least NCX*NCY for square grid) */
#define NMAXPOLY 100000        /* total number of sides of polygonal line */   
#define NCX 30            /* number of circles in x direction */
#define NCY 20            /* number of circles in y direction */
#define NPOISSON 500        /* number of points for Poisson C_RAND_POISSON arrangement */
#define NGOLDENSPIRAL 2000  /* max number of points for C_GOLDEN_SPIRAL arrandement */
#define SDEPTH 1            /* Sierpinski gastket depth */

#define LAMBDA 1.5	/* parameter controlling shape of domain */
#define MU 0.5          /* second parameter controlling shape of billiard */
#define FOCI 1          /* set to 1 to draw focal points of ellipse */
#define NPOLY 6             /* number of sides of polygon */
#define APOLY 0.0           /* angle by which to turn polygon, in units of Pi/2 */ 
#define PENROSE_RATIO 2.5    /* parameter controlling the shape of small ellipses in Penrose room */

#define DRAW_BILLIARD 1     /* set to 1 to draw billiard */
#define DRAW_CONSTRUCTION_LINES 0   /* set to 1 to draw additional construction lines for billiard */
#define PERIODIC_BC 0       /* set to 1 to enforce periodic boundary conditions when drawing particles */

#define RESAMPLE 0      /* set to 1 if particles should be added when dispersion too large */
#define DEBUG 0         /* draw trajectories, for debugging purposes */

/* Simulation parameters */

#define NPART 40000     /* number of particles */
#define NPARTMAX 100000	/* maximal number of particles after resampling */
#define LMAX 0.01       /* minimal segment length triggering resampling */ 
#define DMIN 0.02       /* minimal distance to boundary for triggering resampling */ 
#define CYCLE 1         /* set to 1 for closed curve (start in all directions) */
#define SHOWTRAILS 1    /* set to 1 to keep trails of the particles */
#define SHOWZOOM 0      /* set to 1 to show zoom on specific area */
#define PRINT_PARTICLE_NUMBER 0 /* set to 1 to print number of particles */
#define PRINT_LEFT_RIGHT_PARTICLE_NUMBER 1 /* set to 1 to print number of particles on left and right side */
#define PRINT_COLLISION_NUMBER 0 /* set to 1 to print number of collisions */
#define TEST_ACTIVE 1   /* set to 1 to test whether particle is in billiard */

#define TEST_INITIAL_COND 0     /* set to 1 to allow only initial conditions that pass a test */

#define NSTEPS 7500     /* number of frames of movie */
#define TIME 1500        /* time between movie frames, for fluidity of real-time simulation */ 
#define DPHI 0.00002     /* integration step */
#define NVID 25          /* number of iterations between images displayed on screen */
#define END_FRAMES 25    /* number of still frames at the end of the movie */

/* Decreasing TIME accelerates the animation and the movie                               */
/* For constant speed of movie, TIME*DPHI should be kept constant                        */
/* However, increasing DPHI too much deterioriates quality of simulation                 */
/* NVID tells how often a picture is drawn in the animation, increase it for faster anim */
/* For a good quality movie, take for instance TIME = 400, DPHI = 0.00005, NVID = 100    */

/* Colors and other graphical parameters */

#define COLOR_PALETTE 10     /* Color palette, see list in global_pdes.c  */

#define NCOLORS 128      /* number of colors */
#define COLORSHIFT 0     /* hue of initial color */ 
#define RAINBOW_COLOR 1  /* set to 1 to use different colors for all particles */
#define FLOWER_COLOR 0   /* set to 1 to adapt initial colors to flower billiard (tracks vs core) */
#define NSEG 100         /* number of segments of boundary */
#define LENGTH 0.04       /* length of velocity vectors */
#define BILLIARD_WIDTH 2    /* width of billiard */
#define PARTICLE_WIDTH 2    /* width of particles */
#define FRONT_WIDTH 3       /* width of wave front */

#define BLACK 1             /* set to 1 for black background */
#define COLOR_OUTSIDE 0     /* set to 1 for colored outside */ 
#define OUTER_COLOR 270.0   /* color outside billiard */
#define PAINT_INT 0         /* set to 1 to paint interior in other color (for polygon/Reuleaux) */
#define PAINT_EXT 1         /* set to 1 to paint exterior */

#define PAUSE 1000       /* number of frames after which to pause */
#define PSLEEP 2         /* sleep time during pause */
#define SLEEP1  1        /* initial sleeping time */
#define SLEEP2  1       /* final sleeping time */

#define NXMAZE 32      /* width of maze */
#define NYMAZE 32      /* height of maze */
#define MAZE_MAX_NGBH 4     /* max number of neighbours of maze cell */
#define RAND_SHIFT 2        /* seed of random number generator */
#define MAZE_XSHIFT 0.0     /* horizontal shift of maze */

```

### 12 December 22 - Randomly started waves refracted on a honeycomb lattice ###

**Program:** `wave_billiard.c` 

**Initial condition in function `animation()`:** 
```

```
    x = XMIN + (XMAX - XMIN)*rand()/RAND_MAX;
    y = YMIN + (YMAX - YMIN)*rand()/RAND_MAX;
    init_circular_wave(x, y, phi, psi, xy_in);
    
    if (ALTERNATE_OSCILLATING_SOURCE) sign = -sign;
    x = XMIN + (XMAX - XMIN)*rand()/RAND_MAX;
    y = YMIN + (YMAX - YMIN)*rand()/RAND_MAX;
    add_circular_wave(sign, x, y, phi, psi, xy_in);
```
#define MOVIE 1         /* set to 1 to generate movie */
#define DOUBLE_MOVIE 1  /* set to 1 to produce movies for wave height and energy simultaneously */
#define SAVE_MEMORY 0   /* set to 1 to save memory when writing tiff images */

/* General geometrical parameters */

#define WINWIDTH 	1920  /* window width */
#define WINHEIGHT 	1000  /* window height */
#define NX 3840          /* number of grid points on x axis */
#define NY 2000          /* number of grid points on y axis */

#define XMIN -2.0
#define XMAX 2.0	/* x interval  */
#define YMIN -1.041666667
#define YMAX 1.041666667	/* y interval for 9/16 aspect ratio */

#define HIGHRES 1       /* set to 1 if resolution of grid is double that of displayed image */

#define JULIA_SCALE 1.0 /* scaling for Julia sets */

/* Choice of the billiard table */

#define B_DOMAIN 57        /* choice of domain shape, see list in global_pdes.c */

#define CIRCLE_PATTERN 1   /* pattern of circles or polygons, see list in global_pdes.c */

#define COMPARISON 0        /* set to 1 to compare two different patterns (beta) */
#define B_DOMAIN_B 20       /* second domain shape, for comparisons */
#define CIRCLE_PATTERN_B 0  /* second pattern of circles or polygons */

#define P_PERCOL 0.25       /* probability of having a circle in C_RAND_PERCOL arrangement */
#define NPOISSON 300        /* number of points for Poisson C_RAND_POISSON arrangement */
#define RANDOM_POLY_ANGLE 1 /* set to 1 to randomize angle of polygons */

#define LAMBDA 0.25	    /* parameter controlling the dimensions of domain */
#define MU 0.0              /* parameter controlling the dimensions of domain */
#define NPOLY 6             /* number of sides of polygon */
#define APOLY 0.0           /* angle by which to turn polygon, in units of Pi/2 */ 
#define MDEPTH 6            /* depth of computation of Menger gasket */
#define MRATIO 3            /* ratio defining Menger gasket */
#define MANDELLEVEL 1000    /* iteration level for Mandelbrot set */
#define MANDELLIMIT 10.0    /* limit value for approximation of Mandelbrot set */
#define FOCI 1              /* set to 1 to draw focal points of ellipse */
#define NGRIDX 12           /* number of grid point for grid of disks */
#define NGRIDY 15           /* number of grid point for grid of disks */

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
#define OSCILLATION_SCHEDULE 1  /* oscillation schedule, see list in global_pdes.c */

#define OMEGA 0.0005       /* frequency of periodic excitation */
#define AMPLITUDE 0.8      /* amplitude of periodic excitation */ 
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
/* The Courant number is given by c*DT/DX, where DT is the time step and DX the lattice spacing */
/* The physical damping coefficient is given by GAMMA/(DT)^2 */
/* Increasing COURANT speeds up the simulation, but decreases accuracy */
/* For similar wave forms, COURANT^2*GAMMA should be kept constant */

#define ADD_OSCILLATING_SOURCE 1        /* set to 1 to add an oscillating wave source */
#define OSCILLATING_SOURCE_PERIOD 50    /* period of oscillating source */
#define ALTERNATE_OSCILLATING_SOURCE 1  /* set to 1 to alternate sign of oscillating source */

/* Boundary conditions, see list in global_pdes.c  */

#define B_COND 2

/* Parameters for length and speed of simulation */

#define NSTEPS 1480       /* number of frames of movie */
#define NVID 12          /* number of iterations between images displayed on screen */
#define NSEG 1000         /* number of segments of boundary */
#define INITIAL_TIME 0      /* time after which to start saving frames */
#define BOUNDARY_WIDTH 2    /* width of billiard boundary */
#define PRINT_SPEED 0       /* print speed of moving source */

#define PAUSE 200       /* number of frames after which to pause */
#define PSLEEP 1         /* sleep time during pause */
#define SLEEP1  1        /* initial sleeping time */
#define SLEEP2  1        /* final sleeping time */
#define MID_FRAMES 20    /* number of still frames between parts of two-part movie */
#define END_FRAMES 100    /* number of still frames at end of movie */
#define FADE 1           /* set to 1 to fade at end of movie */

/* Parameters of initial condition */

#define INITIAL_AMP 0.35         /* amplitude of initial condition */
#define INITIAL_VARIANCE 0.0002  /* variance of initial condition */
#define INITIAL_WAVELENGTH  0.01  /* wavelength of initial condition */

/* Plot type, see list in global_pdes.c  */

#define PLOT 0
#define PLOT_B 5        /* plot type for second movie */

/* Color schemes */

#define COLOR_PALETTE 13     /* Color palette, see list in global_pdes.c  */
#define COLOR_PALETTE_B 18     /* Color palette, see list in global_pdes.c  */

#define BLACK 1          /* background */

#define COLOR_SCHEME 3   /* choice of color scheme, see list in global_pdes.c  */

#define SCALE 0          /* set to 1 to adjust color scheme to variance of field */
#define SLOPE 1.0        /* sensitivity of color on wave amplitude */
#define PHASE_FACTOR 1.0       /* factor in computation of phase in color scheme P_3D_PHASE */
#define PHASE_SHIFT 0.0      /* shift of phase in color scheme P_3D_PHASE */
#define ATTENUATION 0.0  /* exponential attenuation coefficient of contrast with time */
#define E_SCALE 100.0     /* scaling factor for energy representation */
#define LOG_SCALE 0.25     /* scaling factor for energy log representation */
#define LOG_SHIFT 0.0     /* shift of colors on log scale */
#define FLUX_SCALE 1.0e3    /* scaling factor for enegy flux represtnation */
#define RESCALE_COLOR_IN_CENTER 0   /* set to 1 to decrease color intentiy in the center (for wave escaping ring) */

#define COLORHUE 260     /* initial hue of water color for scheme C_LUM */
#define COLORDRIFT 0.0   /* how much the color hue drifts during the whole simulation */
#define LUMMEAN 0.5      /* amplitude of luminosity variation for scheme C_LUM */
#define LUMAMP 0.3       /* amplitude of luminosity variation for scheme C_LUM */
#define HUEMEAN 180.0    /* mean value of hue for color scheme C_HUE */
#define HUEAMP -180.0      /* amplitude of variation of hue for color scheme C_HUE */

#define DRAW_COLOR_SCHEME 1    /* set to 1 to plot the color scheme */
#define COLORBAR_RANGE 2.5     /* scale of color scheme bar */
#define COLORBAR_RANGE_B 2.5    /* scale of color scheme bar for 2nd part */
#define ROTATE_COLOR_SCHEME 0   /* set to 1 to draw color scheme horizontally */

#define SAVE_TIME_SERIES 0      /* set to 1 to save wave time series at a point */

#define NXMAZE 7      /* width of maze */
#define NYMAZE 7      /* height of maze */
#define MAZE_MAX_NGBH 4     /* max number of neighbours of maze cell */
#define RAND_SHIFT 24        /* seed of random number generator */
#define MAZE_XSHIFT 0.0     /* horizontal shift of maze */

/* for compatibility with sub_wave and sub_maze */
#define ADD_POTENTIAL 0
#define POT_MAZE 7
#define POTENTIAL 0
/* end of constants only used by sub_wave and sub_maze */


/* For debugging purposes only */
#define FLOOR 0         /* set to 1 to limit wave amplitude to VMAX */
#define VMAX 10.0       /* max value of wave amplitude */

```

### 11 December 22 - Two airfoils on a torus ###

**Program:** `rde.c` 

**Initial condition in function `animation()`:** `init_shear_flow(-1.0, 0.1, 0.2, 1, 1, 0.2, phi, xy_in);`

```
#define MOVIE 1         /* set to 1 to generate movie */
#define DOUBLE_MOVIE 1  /* set to 1 to produce movies for wave height and energy simultaneously */
#define SAVE_MEMORY 0           /* set to 1 to save memory when writing tiff images */

/* General geometrical parameters */

#define WINWIDTH 	1920  /* window width */
#define WINHEIGHT 	1000  /* window height */
#define NX 960          /* number of grid points on x axis */
#define NY 500          /* number of grid points on y axis */
 
#define XMIN -2.0
#define XMAX 2.0	/* x interval */
#define YMIN -1.041666667
#define YMAX 1.041666667	/* y interval for 9/16 aspect ratio */

/* Choice of simulated equation */

#define RDE_EQUATION 6  /* choice of reaction term, see list in global_3d.c */
#define NFIELDS 2       /* number of fields in reaction-diffusion equation */
#define NLAPLACIANS 1   /* number of fields for which to compute Laplacian */

#define ADD_POTENTIAL 0 /* set to 1 to add a potential (for Schrodinger equation) */
#define ADD_MAGNETIC_FIELD 0    /* set to 1 to add a magnetic field (for Schrodinger equation) - then set POTENTIAL 1 */
#define POTENTIAL 7     /* type of potential or vector potential, see list in global_3d.c  */

#define ANTISYMMETRIZE_WAVE_FCT 0   /* set tot 1 to make wave function antisymmetric */

#define JULIA_SCALE 0.5 /* scaling for Julia sets */

/* Choice of the billiard table */

#define B_DOMAIN 58         /* choice of domain shape, see list in global_pdes.c  */

#define CIRCLE_PATTERN 99    /* pattern of circles, see list in global_pdes.c */

#define P_PERCOL 0.25       /* probability of having a circle in C_RAND_PERCOL arrangement */
#define NPOISSON 300        /* number of points for Poisson C_RAND_POISSON arrangement */
#define RANDOM_POLY_ANGLE 0 /* set to 1 to randomize angle of polygons */

#define LAMBDA 0.7	    /* parameter controlling the dimensions of domain */
#define MU 0.15	            /* parameter controlling the dimensions of domain */
#define NPOLY 5             /* number of sides of polygon */
#define APOLY 2.0          /* angle by which to turn polygon, in units of Pi/2 */
#define MDEPTH 7            /* depth of computation of Menger gasket */
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

#define DT 0.0000012

#define VISCOSITY 2.0

#define RPSA 0.75         /* parameter in Rock-Paper-Scissors-type interaction */
#define RPSLZB 0.75       /* second parameter in Rock-Paper-Scissors-Lizard-Spock type interaction */

#define EPSILON 0.8     /* time scale separation */
#define DELTA 0.1       /* time scale separation */
#define FHNA 1.0        /* parameter in FHN equation */
#define FHNC -0.01      /* parameter in FHN equation */
#define K_HARMONIC 1.0  /* spring constant of harmonic potential */
#define K_COULOMB 0.5   /* constant in Coulomb potential */
#define V_MAZE 0.4      /* potential in walls of maze */
#define BZQ 0.0008      /* parameter in BZ equation */
#define BZF 1.2         /* parameter in BZ equation */
#define B_FIELD 10.0    /* magnetic field */
#define AB_RADIUS 0.2   /* radius of region with magnetic field for Aharonov-Bohm effect */
#define K_EULER 50.0    /* constant in stream function integration of Euler equation */

#define SMOOTHEN_VORTICITY 1    /* set to 1 to smoothen vorticity field in Euler equation */
#define SMOOTHEN_PERIOD 10      /* period between smoothenings */
#define SMOOTH_FACTOR 0.03      /* factor by which to smoothen */

#define ADD_TRACERS 1    /* set to 1 to add tracer particles (for Euler equations) */
#define N_TRACERS 1000    /* number of tracer particles */

#define T_OUT 2.0       /* outside temperature */
#define T_IN 0.0        /* inside temperature */
#define SPEED 0.0       /* speed of drift to the right */

#define ADD_NOISE 0     /* set to 1 to add noise, set to 2 to add noise in right half */
#define NOISE_INTENSITY 0.005      /* noise intensity */
#define CHANGE_NOISE 1      /* set to 1 to increase noise intensity */
#define NOISE_FACTOR 40.0   /* factor by which to increase noise intensity */
#define NOISE_INITIAL_TIME 100  /* initial time during which noise remains constant */

#define CHANGE_VISCOSITY 0      /* set to 1 to change the viscosity in the course of the simulation */
#define ADJUST_INTSTEP 0       /* set to 1 to decrease integration step when viscosity increases */
#define VISCOSITY_INITIAL_TIME 10  /* initial time during which viscosity remains constant */
#define VISCOSITY_FACTOR 100.0   /* factor by which to change viscosity */
#define VISCOSITY_MAX 2.0        /* max value of viscosity beyond which NVID is increased and integration step is decrase, 
                                    for numerical stability */
                                        
#define CHANGE_RPSLZB 0         /* set to 1 to change second parameter in Rock-Paper-Scissors-Lizard-Spock equation */
#define RPSLZB_CHANGE 0.75      /* factor by which to rpslzb parameter */
#define RPSLZB_INITIAL_TIME 0   /* initial time during which rpslzb remains constant */
#define RPSLZB_FINAL_TIME 500   /* final time during which rpslzb remains constant */
                                      

/* Boundary conditions, see list in global_pdes.c  */

#define B_COND 1

#define EULER_GRADIENT_YSHIFT 0.0    /* y-shift in computation of gradient in Euler equation */

/* Parameters for length and speed of simulation */

#define NSTEPS 2150       /* number of frames of movie */
#define NVID 120          /* number of iterations between images displayed on screen */
#define ACCELERATION_FACTOR 1.0 /* factor by which to increase NVID in course of simulation */
#define DT_ACCELERATION_FACTOR 1.0 /* factor by which to increase time step in course of simulation  */
#define MAX_DT 0.024     /* maximal value of integration step */
#define NSEG 100         /* number of segments of boundary */
#define BOUNDARY_WIDTH 2    /* width of billiard boundary */

#define PAUSE 100       /* number of frames after which to pause */
#define PSLEEP 2         /* sleep time during pause */
#define SLEEP1  2        /* initial sleeping time */
#define SLEEP2  1        /* final sleeping time */
#define INITIAL_TIME 0  /* initial still time */
#define MID_FRAMES 50    /* number of still frames between parts of two-part movie */
#define END_FRAMES 50    /* number of still frames at end of movie */
#define FADE 1           /* set to 1 to fade at end of movie */

/* Visualisation */

#define PLOT_3D 0    /* controls whether plot is 2D or 3D */

#define ROTATE_VIEW 0       /* set to 1 to rotate position of observer */
#define ROTATE_ANGLE 360.0  /* total angle of rotation during simulation */

#define DRAW_PERIODICISED 1     /* set to 1 to repeat wave periodically in x and y directions */

/* Plot type - color scheme */

#define CPLOT 52
#define CPLOT_B 51

/* Plot type - height of 3D plot */

#define ZPLOT 52     /* z coordinate in 3D plot */
#define ZPLOT_B 51    /* z coordinate in second 3D plot */

#define AMPLITUDE_HIGH_RES 1    /* set to 1 to increase resolution of P_3D_AMPLITUDE plot */
#define SHADE_3D 1              /* set to 1 to change luminosity according to normal vector */
#define NON_DIRICHLET_BC 0      /* set to 1 to draw only facets in domain, if field is not zero on boundary */
#define WRAP_ANGLE 1            /* experimental: wrap angle to [0, 2Pi) for interpolation in angle schemes */
#define FADE_IN_OBSTACLE 0      /* set to 1 to fade color inside obstacles */
#define DRAW_OUTSIDE_GRAY 0     /* experimental - draw outside of billiard in gray */
#define ADD_POTENTIAL_TO_Z 1    /* set to 1 to add the external potential to z-coordinate of plot */
#define ADD_POT_CONSTANT 0.35   /* constant in front of added potential */

#define PLOT_SCALE_ENERGY 0.05      /* vertical scaling in energy plot */

#define PRINT_TIME 0        /* set to 1 to print running time */
#define PRINT_VISCOSITY 0   /* set to 1 to print viscosity */
#define PRINT_RPSLZB 0      /* set to 1 to print rpslzb parameter */
#define PRINT_PROBABILITIES 0   /* set to 1 to print probabilities (for Ehrenfest urn configuration) */
#define PRINT_NOISE 0       /* set to 1 to print noise intensity */

#define DRAW_FIELD_LINES 0  /* set to 1 to draw field lines */
#define FIELD_LINE_WIDTH 1  /* width of field lines */
#define N_FIELD_LINES 120   /* number of field lines */
#define FIELD_LINE_FACTOR 120 /* factor controlling precision when computing origin of field lines */
#define DRAW_BILLIARD 1     /* set to 1 to draw boundary */
#define DRAW_BILLIARD_FRONT 1     /* set to 1 to draw boundary */
#define FILL_BILLIARD_COMPLEMENT 1  /* set to 1 to fill complement of billiard (for certain shapes only) */

/* 3D representation */

#define REPRESENTATION_3D 1     /* choice of 3D representation */ 

#define REP_AXO_3D 0        /* linear projection (axonometry) */
#define REP_PROJ_3D 1       /* projection on plane orthogonal to observer line of sight */

/* Color schemes, see list in global_pdes.c  */

#define COLOR_PALETTE 10       /* Color palette, see list in global_pdes.c  */
#define COLOR_PALETTE_B 13     /* Color palette, see list in global_pdes.c  */

#define BLACK 1          /* black background */

#define COLOR_SCHEME 3   /* choice of color scheme */

#define COLOR_PHASE_SHIFT 0.0   /* phase shift of color scheme, in units of Pi */

#define SCALE 0          /* set to 1 to adjust color scheme to variance of field */
#define SLOPE 1.0        /* sensitivity of color on wave amplitude */
#define VSCALE_AMPLITUDE 1.5      /* additional scaling factor for color scheme P_3D_AMPLITUDE */
#define ATTENUATION 0.0  /* exponential attenuation coefficient of contrast with time */
#define CURL_SCALE 0.000015   /* scaling factor for curl representation */
#define RESCALE_COLOR_IN_CENTER 0   /* set to 1 to decrease color intentiy in the center (for wave escaping ring) */
#define SLOPE_SCHROD_LUM 50.0       /* sensitivity of luminosity on module, for color scheme Z_ARGUMENT */
#define MIN_SCHROD_LUM 0.2       /* minimal luminosity in color scheme Z_ARGUMENT*/

#define COLORHUE 260     /* initial hue of water color for scheme C_LUM */
#define COLORDRIFT 0.0   /* how much the color hue drifts during the whole simulation */
#define LUMMEAN 0.5      /* amplitude of luminosity variation for scheme C_LUM */
#define LUMAMP 0.3       /* amplitude of luminosity variation for scheme C_LUM */
#define HUEMEAN 359.0    /* mean value of hue for color scheme C_HUE */
#define HUEAMP -359.0    /* amplitude of variation of hue for color scheme C_HUE */
#define E_SCALE 100.0    /* scaling factor for energy representation */
#define FLUX_SCALE 100.0 /* scaling factor for energy representation */
#define LOG_SCALE 0.5    /* scaling factor for energy log representation */
#define LOG_SHIFT 1.0 
#define LOG_MIN 1.0e-3   /* floor value for log vorticity plot */

#define NXMAZE 7      /* width of maze */
#define NYMAZE 7      /* height of maze */
#define MAZE_MAX_NGBH 4     /* max number of neighbours of maze cell */
#define RAND_SHIFT 0        /* seed of random number generator */
#define MAZE_XSHIFT 0.0     /* horizontal shift of maze */

#define DRAW_COLOR_SCHEME 1     /* set to 1 to plot the color scheme */
#define COLORBAR_RANGE 3.0      /* scale of color scheme bar */
#define COLORBAR_RANGE_B 3.0   /* scale of color scheme bar for 2nd part */
#define ROTATE_COLOR_SCHEME 0   /* set to 1 to draw color scheme horizontally */

/* only for compatibility with wave_common.c */
#define TWOSPEEDS 0          /* set to 1 to replace hardcore boundary by medium with different speed */
#define OMEGA 0.005        /* frequency of periodic excitation */
#define COURANT 0.08       /* Courant number */
#define COURANTB 0.03      /* Courant number in medium B */
#define INITIAL_AMP 0.5         /* amplitude of initial condition */
#define INITIAL_VARIANCE 0.0002  /* variance of initial condition */
#define INITIAL_WAVELENGTH  0.1  /* wavelength of initial condition */
#define VSCALE_ENERGY 200.0       /* additional scaling factor for color scheme P_3D_ENERGY */
#define PHASE_FACTOR 20.0       /* factor in computation of phase in color scheme P_3D_PHASE */
#define PHASE_SHIFT 0.0      /* shift of phase in color scheme P_3D_PHASE */
#define OSCILLATION_SCHEDULE 0  /* oscillation schedule, see list in global_pdes.c */
#define AMPLITUDE 0.8      /* amplitude of periodic excitation */ 
#define ACHIRP 0.2        /* acceleration coefficient in chirp */
#define DAMPING 0.0        /* damping of periodic excitation */
#define COMPARISON 0        /* set to 1 to compare two different patterns (beta) */
#define B_DOMAIN_B 20       /* second domain shape, for comparisons */
#define CIRCLE_PATTERN_B 0  /* second pattern of circles or polygons */
/* end of constants added only for compatibility with wave_common.c */


double u_3d[2] = {0.75, -0.45};     /* projections of basis vectors for REP_AXO_3D representation */
double v_3d[2] = {-0.75, -0.45};
double w_3d[2] = {0.0, 0.015};
double light[3] = {0.816496581, -0.40824829, 0.40824829};      /* vector of "light" direction for P_3D_ANGLE color scheme */
double observer[3] = {8.0, 8.0, 8.0};    /* location of observer for REP_PROJ_3D representation */ 
int reset_view = 0;         /* switch to reset 3D view parameters (for option ROTATE_VIEW) */

#define Z_SCALING_FACTOR 0.08  /* overall scaling factor of z axis for REP_PROJ_3D representation */
#define XY_SCALING_FACTOR 1.7  /* overall scaling factor for on-screen (x,y) coordinates after projection */
#define ZMAX_FACTOR 1.0        /* max value of z coordinate for REP_PROJ_3D representation */
#define XSHIFT_3D -0.1         /* overall x shift for REP_PROJ_3D representation */
#define YSHIFT_3D 0.1          /* overall y shift for REP_PROJ_3D representation */
#define BORDER_PADDING 0       /* distance from boundary at which to plot points, to avoid boundary effects due to gradient */

/* For debugging purposes only */
#define FLOOR 1         /* set to 1 to limit wave amplitude to VMAX */
#define VMAX 10.0        /* max value of wave amplitude */
#define TEST_GRADIENT 0 /* print norm squared of gradient */

```

### 10 December 22 - Mixing (or lack thereof) of Lennard-Jones particles during a cave-in ###

**Program:** `lennardjones.c` 

```
#define MOVIE 1         /* set to 1 to generate movie */
#define DOUBLE_MOVIE 0  /* set to 1 to produce movies for wave height and energy simultaneously */

#define TIME_LAPSE 1     /* set to 1 to add a time-lapse movie at the end */
                         /* so far incompatible with double movie */
#define TIME_LAPSE_FACTOR 3    /* factor of time-lapse movie */
#define TIME_LAPSE_FIRST 1  /* set to 1 to show time-lapse version first */

/* General geometrical parameters */

#define WINWIDTH 	1280  /* window width */
#define WINHEIGHT 	720   /* window height */

#define XMIN -2.0
#define XMAX 2.0	/* x interval */
#define YMIN -1.125
#define YMAX 1.125	/* y interval for 9/16 aspect ratio */

#define INITXMIN -1.95
#define INITXMAX 1.95	/* x interval for initial condition */
#define INITYMIN 0.05
#define INITYMAX 1.1	/* y interval for initial condition */

#define BCXMIN -2.0
#define BCXMAX 2.0	/* x interval for boundary condition */
#define BCYMIN -1.125
#define BCYMAX 1.125	/* y interval for boundary condition */

#define OBSXMIN -2.0
#define OBSXMAX 2.0     /* x interval for motion of obstacle */

#define CIRCLE_PATTERN 8  /* pattern of circles, see list in global_ljones.c */

#define ADD_FIXED_OBSTACLES 0   /* set to 1 do add fixed circular obstacles */
#define OBSTACLE_PATTERN 3  /* pattern of obstacles, see list in global_ljones.c */

#define ADD_FIXED_SEGMENTS 1    /* set to 1 to add fixed segments as obstacles */
#define SEGMENT_PATTERN 18     /* pattern of repelling segments, see list in global_ljones.c */
#define ROCKET_SHAPE 2        /* shape of rocket combustion chamber, see list in global_ljones.c */
#define ROCKET_SHAPE_B 2      /* shape of second rocket */
#define NOZZLE_SHAPE 2        /* shape of nozzle, see list in global_ljones.c */
#define NOZZLE_SHAPE_B 4      /* shape of nozzle for second rocket, see list in global_ljones.c */

#define TWO_TYPES 0         /* set to 1 to have two types of particles */
#define TPYE_PROPORTION 0.5 /* proportion of particles of first type */
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
#define PDISC_DISTANCE 2.7  /* minimal distance in Poisson disc process, controls density of particles */
#define PDISC_CANDIDATES 100 /* number of candidates in construction of Poisson disc process */
#define RANDOM_POLY_ANGLE 0 /* set to 1 to randomize angle of polygons */

#define LAMBDA 0.8	    /* parameter controlling the dimensions of domain */
#define MU 0.008 	    /* parameter controlling radius of particles */
#define MU_B 0.018          /* parameter controlling radius of particles of second type */
#define NPOLY 25            /* number of sides of polygon */
#define APOLY 0.666666666   /* angle by which to turn polygon, in units of Pi/2 */ 
#define MDEPTH 4            /* depth of computation of Menger gasket */
#define MRATIO 3            /* ratio defining Menger gasket */
#define MANDELLEVEL 1000    /* iteration level for Mandelbrot set */
#define MANDELLIMIT 10.0    /* limit value for approximation of Mandelbrot set */
#define FOCI 1              /* set to 1 to draw focal points of ellipse */
#define NGRIDX 90           /* number of grid point for grid of disks */
#define NGRIDY 100           /* number of grid point for grid of disks */
#define EHRENFEST_RADIUS 0.9    /* radius of container for Ehrenfest urn configuration */
#define EHRENFEST_WIDTH 0.035     /* width of tube for Ehrenfest urn configuration */
#define TWO_CIRCLES_RADIUS_RATIO 0.8    /* ratio of radii for S_TWO_CIRCLES_EXT segment configuration */
#define DAM_WIDTH 0.05       /* width of dam for S_DAM segment configuration */

#define X_SHOOTER -0.2
#define Y_SHOOTER -0.6
#define X_TARGET 0.4
#define Y_TARGET 0.7        /* shooter and target positions in laser fight */

/* Parameters for length and speed of simulation */

#define NSTEPS 7000      /* number of frames of movie */
#define NVID 170         /* number of iterations between images displayed on screen */
#define NSEG 250         /* number of segments of boundary */
#define INITIAL_TIME 200     /* time after which to start saving frames */
#define OBSTACLE_INITIAL_TIME 200     /* time after which to start moving obstacle */
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

#define PLOT 11
#define PLOT_B 8        /* plot type for second movie */

#define DRAW_BONDS 0    /* set to 1 to draw bonds between neighbours */
#define COLOR_BONDS 1   /* set to 1 to color bonds according to length */
#define FILL_TRIANGLES 0    /* set to 1 to fill triangles between neighbours */
#define ALTITUDE_LINES 0    /* set to 1 to add horizontal lines to show altitude */
#define COLOR_SEG_GROUPS 0  /* set to 1 to collor segment groups differently */

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

#define ENERGY_HUE_MIN 330.0        /* color of original particle */
#define ENERGY_HUE_MAX 50.0         /* color of saturated particle */
#define PARTICLE_HUE_MIN 359.0      /* color of original particle */
#define PARTICLE_HUE_MAX 0.0        /* color of saturated particle */
#define PARTICLE_EMAX 1.0e2         /* energy of particle with hottest color */
#define HUE_TYPE0 70.0     /* hue of particles of type 0 */
#define HUE_TYPE1 280.0      /* hue of particles of type 1 */
#define HUE_TYPE2 .0      /* hue of particles of type 2 */
#define HUE_TYPE3 210.0     /* hue of particles of type 3 */

#define RANDOM_RADIUS 0     /* set to 1 for random circle radius */
#define DT_PARTICLE 4.5e-6    /* time step for particle displacement */
#define KREPEL 12.0          /* constant in repelling force between particles */
#define EQUILIBRIUM_DIST 3.7    /* Lennard-Jones equilibrium distance */
#define EQUILIBRIUM_DIST_B 3.5  /* Lennard-Jones equilibrium distance for second type of particle */
#define REPEL_RADIUS 15.0    /* radius in which repelling force acts (in units of particle radius) */
#define DAMPING 250.0          /* damping coefficient of particles */
#define INITIAL_DAMPING 500.0  /* damping coefficient of particles during initial phase */
#define PARTICLE_MASS 1.0    /* mass of particle of radius MU */
#define PARTICLE_MASS_B 1.0   /* mass of particle of radius MU */
#define PARTICLE_INERTIA_MOMENT 0.2     /* moment of inertia of particle */
#define PARTICLE_INERTIA_MOMENT_B 0.02     /* moment of inertia of second type of particle */
#define V_INITIAL 0.0        /* initial velocity range */
#define OMEGA_INITIAL 10.0        /* initial angular velocity range */

#define THERMOSTAT 0        /* set to 1 to switch on thermostat */
#define VARY_THERMOSTAT 0   /* set to 1 for time-dependent thermostat schedule */
#define SIGMA 5.0           /* noise intensity in thermostat */
#define BETA 0.02          /* initial inverse temperature */
#define MU_XI 0.01           /* friction constant in thermostat */
#define KSPRING_BOUNDARY 1.0e7    /* confining harmonic potential outside simulation region */
#define KSPRING_OBSTACLE 1.0e11    /* harmonic potential of obstacles */
#define NBH_DIST_FACTOR 7.5        /* radius in which to count neighbours */
#define GRAVITY 2500.0             /* gravity acting on all particles */
#define GRAVITY_X 0.0        /* horizontal gravity acting on all particles */
#define INCREASE_GRAVITY 0     /* set to 1 to increase gravity during the simulation */
#define GRAVITY_SCHEDULE 2     /* type of gravity schedule, see list in global_ljones.c */
#define GRAVITY_FACTOR 100.0    /* factor by which to increase gravity */
#define GRAVITY_INITIAL_TIME 200    /* time at start of simulation with constant gravity */
#define GRAVITY_RESTORE_TIME 700    /* time at end of simulation with gravity restored to initial value */

#define ROTATION 0           /* set to 1 to include rotation of particles */
#define COUPLE_ANGLE_TO_THERMOSTAT 0    /* set to 1 to couple angular degrees of freedom to thermostat */
#define DIMENSION_FACTOR 1.0  /* scaling factor taking into account number of degrees of freedom */  
#define KTORQUE 100.0          /* force constant in angular dynamics */
#define KTORQUE_B 10.0        /* force constant in angular dynamics */
#define KTORQUE_DIFF 150.0    /* force constant in angular dynamics for different particles */
#define DRAW_SPIN 0           /* set to 1 to draw spin vectors of particles */
#define DRAW_SPIN_B 0         /* set to 1 to draw spin vectors of particles */
#define DRAW_CROSS 1          /* set to 1 to draw cross on particles of second type */
#define SPIN_RANGE 7.0       /* range of spin-spin interaction */
#define SPIN_RANGE_B 5.0     /* range of spin-spin interaction for second type of particle */
#define QUADRUPOLE_RATIO 0.6  /* anisotropy in quadrupole potential */ 

#define INCREASE_BETA 0  /* set to 1 to increase BETA during simulation */
#define BETA_FACTOR 0.025   /* factor by which to change BETA during simulation */
#define N_TOSCILLATIONS 1.5   /* number of temperature oscillations in BETA schedule */
#define NO_OSCILLATION 1        /* set to 1 to have exponential BETA change only */
#define MIDDLE_CONSTANT_PHASE 2000   /* final phase in which temperature is constant */
#define FINAL_DECREASE_PHASE 1300    /* final phase in which temperature decreases */ 
#define FINAL_CONSTANT_PHASE -1     /* final phase in which temperature is constant */

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
#define PARTIAL_THERMO_REGION 1     /* region for partial thermostat coupling (see list in global_ljones.c) */
#define PARTIAL_THERMO_SHIFT 0.2    /* distance from obstacle at the right of which particles are coupled to thermostat */
#define PARTIAL_THERMO_WIDTH 0.5    /* vertical size of partial thermostat coupling */
#define PARTIAL_THERMO_HEIGHT 0.2   /* vertical size of partial thermostat coupling */

#define INCREASE_KREPEL 0   /* set to 1 to increase KREPEL during simulation */
#define KREPEL_FACTOR 1000.0   /* factor by which to change KREPEL during simulation */

#define PART_AT_BOTTOM 0     /* set to 1 to include "seed" particles at bottom */
#define MASS_PART_BOTTOM 10000.0 /* mass of particles at bottom */
#define NPART_BOTTOM 100     /* number of particles at the bottom */

#define ADD_PARTICLES 0    /* set to 1 to add particles */
#define ADD_TIME 0       /* time at which to add first particle */
#define ADD_PERIOD 10000       /* time interval between adding further particles */
#define FINAL_NOADD_PERIOD 200  /* final period where no particles are added */
#define SAFETY_FACTOR 2.0  /* no particles are added at distance less than MU*SAFETY_FACTOR of other particles */

#define TRACER_PARTICLE 0   /* set to 1 to have a tracer particle */
#define N_TRACER_PARTICLES 3    /* number of tracer particles */
#define TRAJECTORY_LENGTH 8000   /* length of recorded trajectory */
#define TRACER_PARTICLE_MASS 4.0    /* relative mass of tracer particle */
#define TRAJECTORY_WIDTH 3      /* width of tracer particle trajectory */

#define ROTATE_BOUNDARY 0           /* set to 1 to rotate the repelling segments */
#define SMOOTH_ROTATION 1           /* set to 1 to update segments at each time step (rather than at each movie frame) */
#define PERIOD_ROTATE_BOUNDARY 1000  /* period of rotating boundary */
#define ROTATE_INITIAL_TIME 0       /* initial time without rotation */
#define ROTATE_FINAL_TIME 100       /* final time without rotation */
#define ROTATE_CHANGE_TIME 0.33     /* relative duration of acceleration/deceleration phases */
#define OMEGAMAX 100.0              /* maximal rotation speed */
#define PRINT_OMEGA 0               /* set to 1 to print angular speed */
#define PRINT_PARTICLE_SPEEDS 0     /* set to 1 to print average speeds/momenta of particles */
#define PRINT_SEGMENTS_SPEEDS 1     /* set to 1 to print velocity of moving segments */

#define MOVE_BOUNDARY 0        /* set to 1 to move repelling segments, due to force from particles */
#define SEGMENTS_MASS 40.0     /* mass of collection of segments */
#define DEACTIVATE_SEGMENT 1    /* set to 1 to deactivate last segment after a certain time */
#define SEGMENT_DEACTIVATION_TIME 200   /* time at which to deactivate last segment */
#define RELEASE_ROCKET_AT_DEACTIVATION 1    /* set to 1 to limit segments velocity before segment release */
#define SEGMENTS_X0 1.5        /* initial position of segments */
#define SEGMENTS_Y0 0.0        /* initial position of segments */
#define SEGMENTS_VX0 0.0       /* initial velocity of segments */
#define SEGMENTS_VY0 0.0      /* initial velocity of segments */
#define DAMP_SEGS_AT_NEGATIVE_Y 0   /* set to 1 to dampen segments when y coordinate is negative */

#define MOVE_SEGMENT_GROUPS 1       /* set to 1 to group segments into moving units */
#define SEGMENT_GROUP_MASS 1000.0   /* mass of segment group */
#define SEGMENT_GROUP_I 1000.0      /* moment of inertia of segment group */
#define SEGMENT_GROUP_DAMPING 0.0   /* damping of segment groups */
#define GROUP_REPULSION 1           /* set to 1 for groups of segments to repel each other */
#define KSPRING_GROUPS 1.0e11       /* harmonic potential between segment groups */
#define GROUP_WIDTH 0.05            /* interaction width of groups */
#define GROUP_G_REPEL 1             /* set to 1 to add repulsion between centers of mass of groups */
#define GROUP_G_REPEL_RADIUS 1.2    /* radius within which centers of mass of groups repel each other */
#define TRACK_SEGMENT_GROUPS 1      /* set to 1 for view to track group of segments */
#define TRACK_X_PADDING 2.0         /* distance from x boundary where tracking starts */

#define POSITION_DEPENDENT_TYPE 0   /* set to 1 to make particle type depend on initial position */
#define POSITION_Y_DEPENDENCE 0     /* set to 1 for the separation between particles to be horizontal */
#define PRINT_ENTROPY 0     /* set to 1 to compute entropy */

#define REACTION_DIFFUSION 0    /* set to 1 to simulate a chemical reaction (particles may change type) */
#define RD_TYPES 3              /* number of types in reaction-diffusion equation */
#define REACTION_PROB 0.0045      /* probability controlling reaction term */ 

#define PRINT_PARTICLE_NUMBER 0     /* set to 1 to print total number of particles */
#define PRINT_LEFT 0        /* set to 1 to print certain parameters at the top left instead of right */
#define PLOT_SPEEDS 0       /* set to 1 to add a plot of obstacle speeds (e.g. for rockets) */
#define PLOT_TRAJECTORIES 0     /* set to 1 to add a plot of obstacle trajectories (e.g. for rockets) */
#define VMAX_PLOT_SPEEDS 0.6    /* vertical scale of plot of obstacle speeds */

#define EHRENFEST_COPY 0    /* set to 1 to add equal number of larger particles (for Ehrenfest model) */

#define LID_MASS 1000.0     /* mass of lid for BC_RECTANGLE_LID b.c. */
#define LID_WIDTH 0.1       /* width of lid for BC_RECTANGLE_LID b.c. */
#define WALL_MASS 2000.0    /* mass of wall for BC_RECTANGLE_WALL b.c. */
#define WALL_FRICTION 0.0   /* friction on wall for BC_RECTANGLE_WALL b.c. */
#define WALL_WIDTH 0.1      /* width of wall for BC_RECTANGLE_WALL b.c. */
#define WALL_VMAX 100.0     /* max speed of wall */
#define WALL_TIME 0         /* time during which to keep wall */

#define NXMAZE 10      /* width of maze */
#define NYMAZE 10      /* height of maze */
#define MAZE_MAX_NGBH 4     /* max number of neighbours of maze cell */
#define RAND_SHIFT 200      /* seed of random number generator */
#define MAZE_XSHIFT 0.0     /* horizontal shift of maze */

#define FLOOR_FORCE 1      /* set to 1 to limit force on particle to FMAX */
#define FMAX 1.0e12         /* maximal force */
#define FLOOR_OMEGA 0      /* set to 1 to limit particle momentum to PMAX */
#define PMAX 1000.0        /* maximal force */

#define HASHX 180   /* size of hashgrid in x direction */
#define HASHY 90    /* size of hashgrid in y direction */
#define HASHMAX 100  /* maximal number of particles per hashgrid cell */
#define HASHGRID_PADDING 0.1    /* padding of hashgrid outside simulation window */

#define DRAW_COLOR_SCHEME 0     /* set to 1 to plot the color scheme */
#define COLORBAR_RANGE 8.0    /* scale of color scheme bar */
#define COLORBAR_RANGE_B 12.0    /* scale of color scheme bar for 2nd part */
#define ROTATE_COLOR_SCHEME 0   /* set to 1 to draw color scheme horizontally */

```

### 9 December 22 - Shining a rotating laser beam into a 16 by 16 maze with slanted corners ###

**Program:** `particle_trajectory.c` 

**Initial condition in function `animation()`:** 
```
    x = -1.65;
    y = 0.03;
    alpha = angle_schedule(time);
    period = compute_trajectories_xy(x, y, alpha, alpha + DPI, configs, trajectory, traj_length, &xmax); 
```

```
#define MOVIE 1         /* set to 1 to generate movie */

#define WINWIDTH 	1280  /* window width */
#define WINHEIGHT 	720   /* window height */

#define XMIN -2.0
#define XMAX 2.0	/* x interval */
#define YMIN -1.125
#define YMAX 1.125	/* y interval for 9/16 aspect ratio */

#define SCALING_FACTOR 1.0       /* scaling factor of drawing, needed for flower billiards, otherwise set to 1.0 */

/* Choice of the billiard table, see global_particles.c */

#define B_DOMAIN 30     /* choice of domain shape */

#define CIRCLE_PATTERN 6   /* pattern of circles */
#define POLYLINE_PATTERN 11  /* pattern of polyline */

#define ABSORBING_CIRCLES 1 /* set to 1 for circular scatterers to be absorbing */

#define NMAXCIRCLES 50000        /* total number of circles (must be at least NCX*NCY for square grid) */
#define NMAXPOLY 50000        /* total number of sides of polygonal line */   
#define NCX 9             /* number of circles in x direction */
#define NCY 20            /* number of circles in y direction */
#define NPOISSON 500        /* number of points for Poisson C_RAND_POISSON arrangement */
#define NGOLDENSPIRAL 2000  /* max number of points for C_GOLDEN_SPIRAL arrandement */
#define SDEPTH 2            /* Sierpinski gastket depth */

#define LAMBDA -1.1	/* parameter controlling shape of domain */
#define MU 0.1          /* second parameter controlling shape of billiard */
#define FOCI 1          /* set to 1 to draw focal points of ellipse */
#define NPOLY 3             /* number of sides of polygon */
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
#define TRAJ_LENGTH 10000 /* length of trajectory */
#define LMAX 0.01       /* minimal segment length triggering resampling */ 
#define DMIN 0.02       /* minimal distance to boundary for triggering resampling */ 
#define CYCLE 1         /* set to 1 for closed curve (start in all directions) */
#define SHOWTRAILS 0    /* set to 1 to keep trails of the particles */
#define SHOWZOOM 0      /* set to 1 to show zoom on specific area */
#define TEST_ACTIVE 1   /* set to 1 to test whether particle is in billiard */
#define PRINT_TRAJECTORY_LENGTH 1   /* set to 1 to print length of trajectory 0 */
#define PRINT_TRAJECTORY_ANGLE 1    /* set to 1 to print angle of trajectory 0 */
#define PRINT_TRAJECTORY_PERIOD 0   /* set to 1 to print period of trajectory 0 */
#define DRAW_LENGTHS_PLOT 1         /* set to 1 to plot trajectory lengths */
#define LENGTHS_LOG_SCALE 1         /* set to 1 to use log scale for plot of lengths */
#define MIN_ANGLE 10.0         /* range of angles of trajectory */
#define MAX_ANGLE 50.0         /* range of angles of trajectory */
#define TEST_INITIAL_COND 1     /* set to 1 to allow only initial conditions that pass a test */

#define SLOW_AT_LONG_TRAJ 0     /* set to 1 to slow down movie for long trajectories */
#define ADD_SUCCESS_GALLERY 1   /* set to 1 to add gallery of successful trajectories at end of movie */
#define SUCCESS_GALLERY_FRAMES 25   /* number of frames per success */

#define NSTEPS 2500      /* number of frames of movie */
#define TIME 2500         /* time between movie frames, for fluidity of real-time simulation */ 
#define DPHI 0.00001     /* integration step */
#define NVID 150         /* number of iterations between images displayed on screen */

/* Decreasing TIME accelerates the animation and the movie                               */
/* For constant speed of movie, TIME*DPHI should be kept constant                        */
/* However, increasing DPHI too much deterioriates quality of simulation                 */
/* NVID tells how often a picture is drawn in the animation, increase it for faster anim */
/* For a good quality movie, take for instance TIME = 400, DPHI = 0.00005, NVID = 100    */

/* Colors and other graphical parameters */

#define COLOR_PALETTE 0     /* Color palette, see list in global_pdes.c  */

#define NCOLORS 500       /* number of colors */
#define COLORSHIFT 0     /* hue of initial color */ 
#define RAINBOW_COLOR 0  /* set to 1 to use different colors for all particles */
#define FLOWER_COLOR 0   /* set to 1 to adapt initial colors to flower billiard (tracks vs core) */
#define NSEG 100         /* number of segments of boundary */
#define LENGTH 0.03       /* length of velocity vectors */
#define BILLIARD_WIDTH 2    /* width of billiard */
#define PARTICLE_WIDTH 2    /* width of particles */
#define FRONT_WIDTH 3       /* width of wave front */
#define GRAPH_HUE 30.0      /* color hue of graph */
#define SUCCESS_HUE 90.0   /* color hue of success */

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

#define NXMAZE 16      /* width of maze */
#define NYMAZE 16     /* height of maze */
#define MAZE_MAX_NGBH 4     /* max number of neighbours of maze cell */
#define RAND_SHIFT 1       /* seed of random number generator */
#define MAZE_XSHIFT -0.5     /* horizontal shift of maze */

```

### 8 December 22 - Energy flux of a wave refracted on a honeycomb lattice ###

**Program:** `wave_billiard.c` 

**Initial condition in function `animation()`:** `init_circular_wave(0.0, 0.0, phi, psi, xy_in);`

```
#define MOVIE 1         /* set to 1 to generate movie */
#define DOUBLE_MOVIE 1  /* set to 1 to produce movies for wave height and energy simultaneously */
#define SAVE_MEMORY 0   /* set to 1 to save memory when writing tiff images */

/* General geometrical parameters */

#define WINWIDTH 	1920  /* window width */
#define WINHEIGHT 	1000  /* window height */
#define NX 3840          /* number of grid points on x axis */
#define NY 2000          /* number of grid points on y axis */

#define XMIN -2.0
#define XMAX 2.0	/* x interval  */
#define YMIN -1.041666667
#define YMAX 1.041666667	/* y interval for 9/16 aspect ratio */

#define HIGHRES 1       /* set to 1 if resolution of grid is double that of displayed image */

#define JULIA_SCALE 1.0 /* scaling for Julia sets */

/* Choice of the billiard table */

#define B_DOMAIN 57        /* choice of domain shape, see list in global_pdes.c */

#define CIRCLE_PATTERN 1   /* pattern of circles or polygons, see list in global_pdes.c */

#define COMPARISON 0        /* set to 1 to compare two different patterns (beta) */
#define B_DOMAIN_B 20       /* second domain shape, for comparisons */
#define CIRCLE_PATTERN_B 0  /* second pattern of circles or polygons */

#define P_PERCOL 0.25       /* probability of having a circle in C_RAND_PERCOL arrangement */
#define NPOISSON 300        /* number of points for Poisson C_RAND_POISSON arrangement */
#define RANDOM_POLY_ANGLE 1 /* set to 1 to randomize angle of polygons */

#define LAMBDA 0.25	    /* parameter controlling the dimensions of domain */
#define MU 0.0              /* parameter controlling the dimensions of domain */
#define NPOLY 6             /* number of sides of polygon */
#define APOLY 0.0           /* angle by which to turn polygon, in units of Pi/2 */ 
#define MDEPTH 6            /* depth of computation of Menger gasket */
#define MRATIO 3            /* ratio defining Menger gasket */
#define MANDELLEVEL 1000    /* iteration level for Mandelbrot set */
#define MANDELLIMIT 10.0    /* limit value for approximation of Mandelbrot set */
#define FOCI 1              /* set to 1 to draw focal points of ellipse */
#define NGRIDX 12           /* number of grid point for grid of disks */
#define NGRIDY 15           /* number of grid point for grid of disks */

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
#define OSCILLATION_SCHEDULE 1  /* oscillation schedule, see list in global_pdes.c */

#define OMEGA 0.0005       /* frequency of periodic excitation */
#define AMPLITUDE 0.8      /* amplitude of periodic excitation */ 
#define ACHIRP 0.25        /* acceleration coefficient in chirp */
#define DAMPING 0.0        /* damping of periodic excitation */
#define COURANT 0.05        /* Courant number */
#define COURANTB 0.1      /* Courant number in medium B */
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
#define OSCILLATING_SOURCE_PERIOD 25    /* period of oscillating source */
#define ALTERNATE_OSCILLATING_SOURCE 1  /* set to 1 to alternate sign of oscillating source */

/* Boundary conditions, see list in global_pdes.c  */

#define B_COND 2

/* Parameters for length and speed of simulation */

#define NSTEPS 1800       /* number of frames of movie */
#define NVID 12          /* number of iterations between images displayed on screen */
#define NSEG 1000         /* number of segments of boundary */
#define INITIAL_TIME 0      /* time after which to start saving frames */
#define BOUNDARY_WIDTH 2    /* width of billiard boundary */
#define PRINT_SPEED 0       /* print speed of moving source */

#define PAUSE 200       /* number of frames after which to pause */
#define PSLEEP 1         /* sleep time during pause */
#define SLEEP1  1        /* initial sleeping time */
#define SLEEP2  1        /* final sleeping time */
#define MID_FRAMES 20    /* number of still frames between parts of two-part movie */
#define END_FRAMES 100    /* number of still frames at end of movie */
#define FADE 1           /* set to 1 to fade at end of movie */

/* Parameters of initial condition */

#define INITIAL_AMP 0.35         /* amplitude of initial condition */
#define INITIAL_VARIANCE 0.0002  /* variance of initial condition */
#define INITIAL_WAVELENGTH  0.01  /* wavelength of initial condition */

/* Plot type, see list in global_pdes.c  */

#define PLOT 6

#define PLOT_B 7        /* plot type for second movie */

/* Color schemes */

#define COLOR_PALETTE 13     /* Color palette, see list in global_pdes.c  */
#define COLOR_PALETTE_B 18     /* Color palette, see list in global_pdes.c  */

#define BLACK 1          /* background */

#define COLOR_SCHEME 3   /* choice of color scheme, see list in global_pdes.c  */

#define SCALE 0          /* set to 1 to adjust color scheme to variance of field */
#define SLOPE 1.0        /* sensitivity of color on wave amplitude */
#define PHASE_FACTOR 1.0       /* factor in computation of phase in color scheme P_3D_PHASE */
#define PHASE_SHIFT 0.0      /* shift of phase in color scheme P_3D_PHASE */
#define ATTENUATION 0.0  /* exponential attenuation coefficient of contrast with time */
#define E_SCALE 100.0     /* scaling factor for energy representation */
#define LOG_SCALE 0.25     /* scaling factor for energy log representation */
#define LOG_SHIFT 0.0     /* shift of colors on log scale */
#define FLUX_SCALE 1.0e3    /* scaling factor for enegy flux represtnation */
#define RESCALE_COLOR_IN_CENTER 0   /* set to 1 to decrease color intentiy in the center (for wave escaping ring) */

#define COLORHUE 260     /* initial hue of water color for scheme C_LUM */
#define COLORDRIFT 0.0   /* how much the color hue drifts during the whole simulation */
#define LUMMEAN 0.5      /* amplitude of luminosity variation for scheme C_LUM */
#define LUMAMP 0.3       /* amplitude of luminosity variation for scheme C_LUM */
#define HUEMEAN 180.0    /* mean value of hue for color scheme C_HUE */
#define HUEAMP -180.0      /* amplitude of variation of hue for color scheme C_HUE */

#define DRAW_COLOR_SCHEME 1    /* set to 1 to plot the color scheme */
#define COLORBAR_RANGE 3.0     /* scale of color scheme bar */
#define COLORBAR_RANGE_B 1.0    /* scale of color scheme bar for 2nd part */
#define ROTATE_COLOR_SCHEME 0   /* set to 1 to draw color scheme horizontally */

#define SAVE_TIME_SERIES 0      /* set to 1 to save wave time series at a point */

#define NXMAZE 7      /* width of maze */
#define NYMAZE 7      /* height of maze */
#define MAZE_MAX_NGBH 4     /* max number of neighbours of maze cell */
#define RAND_SHIFT 24        /* seed of random number generator */
#define MAZE_XSHIFT 0.0     /* horizontal shift of maze */

/* for compatibility with sub_wave and sub_maze */
#define ADD_POTENTIAL 0
#define POT_MAZE 7
#define POTENTIAL 0
/* end of constants only used by sub_wave and sub_maze */


/* For debugging purposes only */
#define FLOOR 0         /* set to 1 to limit wave amplitude to VMAX */
#define VMAX 10.0       /* max value of wave amplitude */

```

### 7 December 22 - Van Gogh meets Munch: An airfoil in a channel ###

**Program:** `rde.c` 

**Initial condition in function `animation()`:** `init_shear_flow(-1.0, 0.0, 0.1, 1, 1, 0.0, phi, xy_in);`

```
#define MOVIE 1         /* set to 1 to generate movie */
#define DOUBLE_MOVIE 1  /* set to 1 to produce movies for wave height and energy simultaneously */
#define SAVE_MEMORY 0           /* set to 1 to save memory when writing tiff images */

/* General geometrical parameters */

#define WINWIDTH 	1920  /* window width */
#define WINHEIGHT 	1000  /* window height */
#define NX 960          /* number of grid points on x axis */
#define NY 500          /* number of grid points on y axis */
 
#define XMIN -2.0
#define XMAX 2.0	/* x interval */
#define YMIN -1.041666667
#define YMAX 1.041666667	/* y interval for 9/16 aspect ratio */

/* Choice of simulated equation */

#define RDE_EQUATION 6  /* choice of reaction term, see list in global_3d.c */
#define NFIELDS 2       /* number of fields in reaction-diffusion equation */
#define NLAPLACIANS 1   /* number of fields for which to compute Laplacian */

#define ADD_POTENTIAL 0 /* set to 1 to add a potential (for Schrodinger equation) */
#define ADD_MAGNETIC_FIELD 0    /* set to 1 to add a magnetic field (for Schrodinger equation) - then set POTENTIAL 1 */
#define POTENTIAL 7     /* type of potential or vector potential, see list in global_3d.c  */

#define ANTISYMMETRIZE_WAVE_FCT 0   /* set tot 1 to make wave function antisymmetric */

#define JULIA_SCALE 0.5 /* scaling for Julia sets */

/* Choice of the billiard table */

#define B_DOMAIN 197       /* choice of domain shape, see list in global_pdes.c  */

#define CIRCLE_PATTERN 99    /* pattern of circles, see list in global_pdes.c */

#define P_PERCOL 0.25       /* probability of having a circle in C_RAND_PERCOL arrangement */
#define NPOISSON 300        /* number of points for Poisson C_RAND_POISSON arrangement */
#define RANDOM_POLY_ANGLE 0 /* set to 1 to randomize angle of polygons */

#define LAMBDA 0.7	    /* parameter controlling the dimensions of domain */
#define MU 0.15	            /* parameter controlling the dimensions of domain */
#define NPOLY 5             /* number of sides of polygon */
#define APOLY 2.0          /* angle by which to turn polygon, in units of Pi/2 */
#define MDEPTH 7            /* depth of computation of Menger gasket */
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

#define DT 0.0000012

#define VISCOSITY 2.0

#define RPSA 0.75         /* parameter in Rock-Paper-Scissors-type interaction */
#define RPSLZB 0.75       /* second parameter in Rock-Paper-Scissors-Lizard-Spock type interaction */

#define EPSILON 0.8     /* time scale separation */
#define DELTA 0.1       /* time scale separation */
#define FHNA 1.0        /* parameter in FHN equation */
#define FHNC -0.01      /* parameter in FHN equation */
#define K_HARMONIC 1.0  /* spring constant of harmonic potential */
#define K_COULOMB 0.5   /* constant in Coulomb potential */
#define V_MAZE 0.4      /* potential in walls of maze */
#define BZQ 0.0008      /* parameter in BZ equation */
#define BZF 1.2         /* parameter in BZ equation */
#define B_FIELD 10.0    /* magnetic field */
#define AB_RADIUS 0.2   /* radius of region with magnetic field for Aharonov-Bohm effect */
#define K_EULER 50.0    /* constant in stream function integration of Euler equation */

#define SMOOTHEN_VORTICITY 1    /* set to 1 to smoothen vorticity field in Euler equation */
#define SMOOTHEN_PERIOD 10      /* period between smoothenings */
#define SMOOTH_FACTOR 0.03      /* factor by which to smoothen */

#define ADD_TRACERS 1    /* set to 1 to add tracer particles (for Euler equations) */
#define N_TRACERS 1000    /* number of tracer particles */

#define T_OUT 2.0       /* outside temperature */
#define T_IN 0.0        /* inside temperature */
#define SPEED 0.0       /* speed of drift to the right */

#define ADD_NOISE 0     /* set to 1 to add noise, set to 2 to add noise in right half */
#define NOISE_INTENSITY 0.005      /* noise intensity */
#define CHANGE_NOISE 1      /* set to 1 to increase noise intensity */
#define NOISE_FACTOR 40.0   /* factor by which to increase noise intensity */
#define NOISE_INITIAL_TIME 100  /* initial time during which noise remains constant */

#define CHANGE_VISCOSITY 0      /* set to 1 to change the viscosity in the course of the simulation */
#define ADJUST_INTSTEP 0       /* set to 1 to decrease integration step when viscosity increases */
#define VISCOSITY_INITIAL_TIME 10  /* initial time during which viscosity remains constant */
#define VISCOSITY_FACTOR 100.0   /* factor by which to change viscosity */
#define VISCOSITY_MAX 2.0        /* max value of viscosity beyond which NVID is increased and integration step is decrase, 
                                    for numerical stability */
                                        
#define CHANGE_RPSLZB 0         /* set to 1 to change second parameter in Rock-Paper-Scissors-Lizard-Spock equation */
#define RPSLZB_CHANGE 0.75      /* factor by which to rpslzb parameter */
#define RPSLZB_INITIAL_TIME 0   /* initial time during which rpslzb remains constant */
#define RPSLZB_FINAL_TIME 500   /* final time during which rpslzb remains constant */
                                      

/* Boundary conditions, see list in global_pdes.c  */

#define B_COND 1

#define EULER_GRADIENT_YSHIFT 0.0    /* y-shift in computation of gradient in Euler equation */

/* Parameters for length and speed of simulation */

#define NSTEPS 2250       /* number of frames of movie */
#define NVID 120          /* number of iterations between images displayed on screen */
#define ACCELERATION_FACTOR 1.0 /* factor by which to increase NVID in course of simulation */
#define DT_ACCELERATION_FACTOR 1.0 /* factor by which to increase time step in course of simulation  */
#define MAX_DT 0.024     /* maximal value of integration step */
#define NSEG 100         /* number of segments of boundary */
#define BOUNDARY_WIDTH 5    /* width of billiard boundary */

#define PAUSE 100       /* number of frames after which to pause */
#define PSLEEP 2         /* sleep time during pause */
#define SLEEP1  2        /* initial sleeping time */
#define SLEEP2  1        /* final sleeping time */
#define INITIAL_TIME 0  /* initial still time */
#define MID_FRAMES 50    /* number of still frames between parts of two-part movie */
#define END_FRAMES 50    /* number of still frames at end of movie */
#define FADE 1           /* set to 1 to fade at end of movie */

/* Visualisation */

#define PLOT_3D 0    /* controls whether plot is 2D or 3D */

#define ROTATE_VIEW 0       /* set to 1 to rotate position of observer */
#define ROTATE_ANGLE 360.0  /* total angle of rotation during simulation */

#define DRAW_PERIODICISED 1     /* set to 1 to repeat wave periodically in x and y directions */

/* Plot type - color scheme */

#define CPLOT 52
#define CPLOT_B 51

/* Plot type - height of 3D plot */

#define ZPLOT 52     /* z coordinate in 3D plot */
#define ZPLOT_B 51    /* z coordinate in second 3D plot */

#define AMPLITUDE_HIGH_RES 1    /* set to 1 to increase resolution of P_3D_AMPLITUDE plot */
#define SHADE_3D 1              /* set to 1 to change luminosity according to normal vector */
#define NON_DIRICHLET_BC 0      /* set to 1 to draw only facets in domain, if field is not zero on boundary */
#define WRAP_ANGLE 1            /* experimental: wrap angle to [0, 2Pi) for interpolation in angle schemes */
#define FADE_IN_OBSTACLE 0      /* set to 1 to fade color inside obstacles */
#define DRAW_OUTSIDE_GRAY 0     /* experimental - draw outside of billiard in gray */
#define ADD_POTENTIAL_TO_Z 1    /* set to 1 to add the external potential to z-coordinate of plot */
#define ADD_POT_CONSTANT 0.35   /* constant in front of added potential */

#define PLOT_SCALE_ENERGY 0.05      /* vertical scaling in energy plot */

#define PRINT_TIME 0        /* set to 1 to print running time */
#define PRINT_VISCOSITY 0   /* set to 1 to print viscosity */
#define PRINT_RPSLZB 0      /* set to 1 to print rpslzb parameter */
#define PRINT_PROBABILITIES 0   /* set to 1 to print probabilities (for Ehrenfest urn configuration) */
#define PRINT_NOISE 0       /* set to 1 to print noise intensity */

#define DRAW_FIELD_LINES 0  /* set to 1 to draw field lines */
#define FIELD_LINE_WIDTH 1  /* width of field lines */
#define N_FIELD_LINES 120   /* number of field lines */
#define FIELD_LINE_FACTOR 120 /* factor controlling precision when computing origin of field lines */
#define DRAW_BILLIARD 1     /* set to 1 to draw boundary */
#define DRAW_BILLIARD_FRONT 1     /* set to 1 to draw boundary */
#define FILL_BILLIARD_COMPLEMENT 1  /* set to 1 to fill complement of billiard (for certain shapes only) */

/* 3D representation */

#define REPRESENTATION_3D 1     /* choice of 3D representation */ 

#define REP_AXO_3D 0        /* linear projection (axonometry) */
#define REP_PROJ_3D 1       /* projection on plane orthogonal to observer line of sight */

/* Color schemes, see list in global_pdes.c  */

#define COLOR_PALETTE 14       /* Color palette, see list in global_pdes.c  */
#define COLOR_PALETTE_B 13     /* Color palette, see list in global_pdes.c  */

#define BLACK 1          /* black background */

#define COLOR_SCHEME 3   /* choice of color scheme */

#define COLOR_PHASE_SHIFT 0.0   /* phase shift of color scheme, in units of Pi */

#define SCALE 0          /* set to 1 to adjust color scheme to variance of field */
#define SLOPE 1.0        /* sensitivity of color on wave amplitude */
#define VSCALE_AMPLITUDE 1.5      /* additional scaling factor for color scheme P_3D_AMPLITUDE */
#define ATTENUATION 0.0  /* exponential attenuation coefficient of contrast with time */
#define CURL_SCALE 0.000015   /* scaling factor for curl representation */
#define RESCALE_COLOR_IN_CENTER 0   /* set to 1 to decrease color intentiy in the center (for wave escaping ring) */
#define SLOPE_SCHROD_LUM 50.0       /* sensitivity of luminosity on module, for color scheme Z_ARGUMENT */
#define MIN_SCHROD_LUM 0.2       /* minimal luminosity in color scheme Z_ARGUMENT*/

#define COLORHUE 260     /* initial hue of water color for scheme C_LUM */
#define COLORDRIFT 0.0   /* how much the color hue drifts during the whole simulation */
#define LUMMEAN 0.5      /* amplitude of luminosity variation for scheme C_LUM */
#define LUMAMP 0.3       /* amplitude of luminosity variation for scheme C_LUM */
#define HUEMEAN 359.0    /* mean value of hue for color scheme C_HUE */
#define HUEAMP -359.0    /* amplitude of variation of hue for color scheme C_HUE */
#define E_SCALE 100.0    /* scaling factor for energy representation */
#define FLUX_SCALE 100.0 /* scaling factor for energy representation */
#define LOG_SCALE 0.5    /* scaling factor for energy log representation */
#define LOG_SHIFT 1.0   
#define LOG_MIN 1.0e-3   /* floor value for log vorticity plot */

#define NXMAZE 7      /* width of maze */
#define NYMAZE 7      /* height of maze */
#define MAZE_MAX_NGBH 4     /* max number of neighbours of maze cell */
#define RAND_SHIFT 0        /* seed of random number generator */
#define MAZE_XSHIFT 0.0     /* horizontal shift of maze */

#define DRAW_COLOR_SCHEME 1     /* set to 1 to plot the color scheme */
#define COLORBAR_RANGE 3.0      /* scale of color scheme bar */
#define COLORBAR_RANGE_B 3.0   /* scale of color scheme bar for 2nd part */
#define ROTATE_COLOR_SCHEME 0   /* set to 1 to draw color scheme horizontally */

/* only for compatibility with wave_common.c */
#define TWOSPEEDS 0          /* set to 1 to replace hardcore boundary by medium with different speed */
#define OMEGA 0.005        /* frequency of periodic excitation */
#define COURANT 0.08       /* Courant number */
#define COURANTB 0.03      /* Courant number in medium B */
#define INITIAL_AMP 0.5         /* amplitude of initial condition */
#define INITIAL_VARIANCE 0.0002  /* variance of initial condition */
#define INITIAL_WAVELENGTH  0.1  /* wavelength of initial condition */
#define VSCALE_ENERGY 200.0       /* additional scaling factor for color scheme P_3D_ENERGY */
#define PHASE_FACTOR 20.0       /* factor in computation of phase in color scheme P_3D_PHASE */
#define PHASE_SHIFT 0.0      /* shift of phase in color scheme P_3D_PHASE */
#define OSCILLATION_SCHEDULE 0  /* oscillation schedule, see list in global_pdes.c */
#define AMPLITUDE 0.8      /* amplitude of periodic excitation */ 
#define ACHIRP 0.2        /* acceleration coefficient in chirp */
#define DAMPING 0.0        /* damping of periodic excitation */
#define COMPARISON 0        /* set to 1 to compare two different patterns (beta) */
#define B_DOMAIN_B 20       /* second domain shape, for comparisons */
#define CIRCLE_PATTERN_B 0  /* second pattern of circles or polygons */
/* end of constants added only for compatibility with wave_common.c */


double u_3d[2] = {0.75, -0.45};     /* projections of basis vectors for REP_AXO_3D representation */
double v_3d[2] = {-0.75, -0.45};
double w_3d[2] = {0.0, 0.015};
double light[3] = {0.816496581, -0.40824829, 0.40824829};      /* vector of "light" direction for P_3D_ANGLE color scheme */
double observer[3] = {8.0, 8.0, 8.0};    /* location of observer for REP_PROJ_3D representation */ 
int reset_view = 0;         /* switch to reset 3D view parameters (for option ROTATE_VIEW) */

#define Z_SCALING_FACTOR 0.08  /* overall scaling factor of z axis for REP_PROJ_3D representation */
#define XY_SCALING_FACTOR 1.7  /* overall scaling factor for on-screen (x,y) coordinates after projection */
#define ZMAX_FACTOR 1.0        /* max value of z coordinate for REP_PROJ_3D representation */
#define XSHIFT_3D -0.1         /* overall x shift for REP_PROJ_3D representation */
#define YSHIFT_3D 0.1          /* overall y shift for REP_PROJ_3D representation */
#define BORDER_PADDING 0       /* distance from boundary at which to plot points, to avoid boundary effects due to gradient */

/* For debugging purposes only */
#define FLOOR 1         /* set to 1 to limit wave amplitude to VMAX */
#define VMAX 10.0        /* max value of wave amplitude */
#define TEST_GRADIENT 0 /* print norm squared of gradient */

```

### 6 December 22 - Deep impact 6: Comparing a triangular and a rectangular obstacle ###

**Program:** `lennardjones.c` 

```
#define MOVIE 1         /* set to 1 to generate movie */
#define DOUBLE_MOVIE 1  /* set to 1 to produce movies for wave height and energy simultaneously */

#define TIME_LAPSE 0     /* set to 1 to add a time-lapse movie at the end */
                         /* so far incompatible with double movie */
#define TIME_LAPSE_FACTOR 3    /* factor of time-lapse movie */
#define TIME_LAPSE_FIRST 0  /* set to 1 to show time-lapse version first */

/* General geometrical parameters */

#define WINWIDTH 	1280  /* window width */
#define WINHEIGHT 	720   /* window height */

#define XMIN -2.0
#define XMAX 2.0	/* x interval */
#define YMIN -1.125
#define YMAX 1.125	/* y interval for 9/16 aspect ratio */

#define INITXMIN -1.95
#define INITXMAX 1.95	/* x interval for initial condition */
#define INITYMIN -1.1
#define INITYMAX 0.15	/* y interval for initial condition */

#define BCXMIN -2.0
#define BCXMAX 2.0	/* x interval for boundary condition */
#define BCYMIN -1.125
#define BCYMAX 1.125	/* y interval for boundary condition */

#define OBSXMIN -2.0
#define OBSXMAX 2.0     /* x interval for motion of obstacle */

#define CIRCLE_PATTERN 8   /* pattern of circles, see list in global_ljones.c */

#define ADD_FIXED_OBSTACLES 0   /* set to 1 do add fixed circular obstacles */
#define OBSTACLE_PATTERN 3  /* pattern of obstacles, see list in global_ljones.c */

#define ADD_FIXED_SEGMENTS 1    /* set to 1 to add fixed segments as obstacles */
#define SEGMENT_PATTERN 19    /* pattern of repelling segments, see list in global_ljones.c */
#define ROCKET_SHAPE 2        /* shape of rocket combustion chamber, see list in global_ljones.c */
#define ROCKET_SHAPE_B 2      /* shape of second rocket */
#define NOZZLE_SHAPE 1        /* shape of nozzle, see list in global_ljones.c */
#define NOZZLE_SHAPE_B 0      /* shape of nozzle for second rocket, see list in global_ljones.c */

#define TWO_TYPES 0         /* set to 1 to have two types of particles */
#define TPYE_PROPORTION 0.5 /* proportion of particles of first type */
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
#define PDISC_DISTANCE 3.0  /* minimal distance in Poisson disc process, controls density of particles */
#define PDISC_CANDIDATES 100 /* number of candidates in construction of Poisson disc process */
#define RANDOM_POLY_ANGLE 0 /* set to 1 to randomize angle of polygons */

#define LAMBDA 0.3	    /* parameter controlling the dimensions of domain */
#define MU 0.009 	    /* parameter controlling radius of particles */
#define MU_B 0.018          /* parameter controlling radius of particles of second type */
#define NPOLY 3             /* number of sides of polygon */
#define APOLY 1.0           /* angle by which to turn polygon, in units of Pi/2 */ 
#define MDEPTH 4            /* depth of computation of Menger gasket */
#define MRATIO 3            /* ratio defining Menger gasket */
#define MANDELLEVEL 1000    /* iteration level for Mandelbrot set */
#define MANDELLIMIT 10.0    /* limit value for approximation of Mandelbrot set */
#define FOCI 1              /* set to 1 to draw focal points of ellipse */
#define NGRIDX 90           /* number of grid point for grid of disks */
#define NGRIDY 100           /* number of grid point for grid of disks */
#define EHRENFEST_RADIUS 0.9    /* radius of container for Ehrenfest urn configuration */
#define EHRENFEST_WIDTH 0.035     /* width of tube for Ehrenfest urn configuration */
#define TWO_CIRCLES_RADIUS_RATIO 0.8    /* ratio of radii for S_TWO_CIRCLES_EXT segment configuration */
#define DAM_WIDTH 0.05       /* width of dam for S_DAM segment configuration */

#define X_SHOOTER -0.2
#define Y_SHOOTER -0.6
#define X_TARGET 0.4
#define Y_TARGET 0.7        /* shooter and target positions in laser fight */

/* Parameters for length and speed of simulation */

#define NSTEPS 1500      /* number of frames of movie */
#define NVID 100         /* number of iterations between images displayed on screen */
#define NSEG 250         /* number of segments of boundary */
#define INITIAL_TIME 250     /* time after which to start saving frames */
#define OBSTACLE_INITIAL_TIME 350     /* time after which to start moving obstacle */
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
#define PLOT_B 8        /* plot type for second movie */

#define DRAW_BONDS 0    /* set to 1 to draw bonds between neighbours */
#define COLOR_BONDS 1   /* set to 1 to color bonds according to length */
#define FILL_TRIANGLES 1    /* set to 1 to fill triangles between neighbours */
#define ALTITUDE_LINES 0    /* set to 1 to add horizontal lines to show altitude */
#define COLOR_SEG_GROUPS 1  /* set to 1 to collor segment groups differently */

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

#define ENERGY_HUE_MIN 330.0        /* color of original particle */
#define ENERGY_HUE_MAX 50.0         /* color of saturated particle */
#define PARTICLE_HUE_MIN 359.0      /* color of original particle */
#define PARTICLE_HUE_MAX 0.0        /* color of saturated particle */
#define PARTICLE_EMAX 3.0e2         /* energy of particle with hottest color */
#define HUE_TYPE0 70.0     /* hue of particles of type 0 */
#define HUE_TYPE1 280.0      /* hue of particles of type 1 */
#define HUE_TYPE2 .0      /* hue of particles of type 2 */
#define HUE_TYPE3 210.0     /* hue of particles of type 3 */

#define RANDOM_RADIUS 0     /* set to 1 for random circle radius */
#define DT_PARTICLE 3.0e-6    /* time step for particle displacement */
#define KREPEL 12.0          /* constant in repelling force between particles */
#define EQUILIBRIUM_DIST 4.0    /* Lennard-Jones equilibrium distance */
#define EQUILIBRIUM_DIST_B 3.5  /* Lennard-Jones equilibrium distance for second type of particle */
#define REPEL_RADIUS 20.0    /* radius in which repelling force acts (in units of particle radius) */
#define DAMPING 10.0          /* damping coefficient of particles */
#define INITIAL_DAMPING 35.0  /* damping coefficient of particles during initial phase */
#define PARTICLE_MASS 1.0    /* mass of particle of radius MU */
#define PARTICLE_MASS_B 1.0   /* mass of particle of radius MU */
#define PARTICLE_INERTIA_MOMENT 0.2     /* moment of inertia of particle */
#define PARTICLE_INERTIA_MOMENT_B 0.02     /* moment of inertia of second type of particle */
#define V_INITIAL 0.0        /* initial velocity range */
#define OMEGA_INITIAL 10.0        /* initial angular velocity range */

#define THERMOSTAT 0        /* set to 1 to switch on thermostat */
#define VARY_THERMOSTAT 0   /* set to 1 for time-dependent thermostat schedule */
#define SIGMA 5.0           /* noise intensity in thermostat */
#define BETA 0.02          /* initial inverse temperature */
#define MU_XI 0.01           /* friction constant in thermostat */
#define KSPRING_BOUNDARY 1.0e7    /* confining harmonic potential outside simulation region */
#define KSPRING_OBSTACLE 1.0e11    /* harmonic potential of obstacles */
#define NBH_DIST_FACTOR 7.5       /* radius in which to count neighbours */
#define GRAVITY 1500.0            /* gravity acting on all particles */
#define GRAVITY_X 0.0        /* horizontal gravity acting on all particles */
#define INCREASE_GRAVITY 0     /* set to 1 to increase gravity during the simulation */
#define GRAVITY_SCHEDULE 2     /* type of gravity schedule, see list in global_ljones.c */
#define GRAVITY_FACTOR 100.0    /* factor by which to increase gravity */
#define GRAVITY_INITIAL_TIME 200    /* time at start of simulation with constant gravity */
#define GRAVITY_RESTORE_TIME 700    /* time at end of simulation with gravity restored to initial value */

#define ROTATION 0           /* set to 1 to include rotation of particles */
#define COUPLE_ANGLE_TO_THERMOSTAT 0    /* set to 1 to couple angular degrees of freedom to thermostat */
#define DIMENSION_FACTOR 1.0  /* scaling factor taking into account number of degrees of freedom */  
#define KTORQUE 100.0          /* force constant in angular dynamics */
#define KTORQUE_B 10.0        /* force constant in angular dynamics */
#define KTORQUE_DIFF 150.0    /* force constant in angular dynamics for different particles */
#define DRAW_SPIN 0           /* set to 1 to draw spin vectors of particles */
#define DRAW_SPIN_B 0         /* set to 1 to draw spin vectors of particles */
#define DRAW_CROSS 1          /* set to 1 to draw cross on particles of second type */
#define SPIN_RANGE 7.0       /* range of spin-spin interaction */
#define SPIN_RANGE_B 5.0     /* range of spin-spin interaction for second type of particle */
#define QUADRUPOLE_RATIO 0.6  /* anisotropy in quadrupole potential */ 

#define INCREASE_BETA 0  /* set to 1 to increase BETA during simulation */
#define BETA_FACTOR 0.025   /* factor by which to change BETA during simulation */
#define N_TOSCILLATIONS 1.5   /* number of temperature oscillations in BETA schedule */
#define NO_OSCILLATION 1        /* set to 1 to have exponential BETA change only */
#define MIDDLE_CONSTANT_PHASE 1600   /* final phase in which temperature is constant */
#define FINAL_DECREASE_PHASE 1500    /* final phase in which temperature decreases */ 
#define FINAL_CONSTANT_PHASE -1     /* final phase in which temperature is constant */

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
#define PARTIAL_THERMO_REGION 1     /* region for partial thermostat coupling (see list in global_ljones.c) */
#define PARTIAL_THERMO_SHIFT 0.2    /* distance from obstacle at the right of which particles are coupled to thermostat */
#define PARTIAL_THERMO_WIDTH 0.5    /* vertical size of partial thermostat coupling */
#define PARTIAL_THERMO_HEIGHT 0.2   /* vertical size of partial thermostat coupling */

#define INCREASE_KREPEL 0   /* set to 1 to increase KREPEL during simulation */
#define KREPEL_FACTOR 1000.0   /* factor by which to change KREPEL during simulation */

#define PART_AT_BOTTOM 0     /* set to 1 to include "seed" particles at bottom */
#define MASS_PART_BOTTOM 10000.0 /* mass of particles at bottom */
#define NPART_BOTTOM 100     /* number of particles at the bottom */

#define ADD_PARTICLES 0    /* set to 1 to add particles */
#define ADD_TIME 0       /* time at which to add first particle */
#define ADD_PERIOD 10000       /* time interval between adding further particles */
#define FINAL_NOADD_PERIOD 200  /* final period where no particles are added */
#define SAFETY_FACTOR 2.0  /* no particles are added at distance less than MU*SAFETY_FACTOR of other particles */

#define TRACER_PARTICLE 0   /* set to 1 to have a tracer particle */
#define N_TRACER_PARTICLES 3    /* number of tracer particles */
#define TRAJECTORY_LENGTH 8000   /* length of recorded trajectory */
#define TRACER_PARTICLE_MASS 4.0    /* relative mass of tracer particle */
#define TRAJECTORY_WIDTH 3      /* width of tracer particle trajectory */

#define ROTATE_BOUNDARY 0           /* set to 1 to rotate the repelling segments */
#define SMOOTH_ROTATION 1           /* set to 1 to update segments at each time step (rather than at each movie frame) */
#define PERIOD_ROTATE_BOUNDARY 1000  /* period of rotating boundary */
#define ROTATE_INITIAL_TIME 0       /* initial time without rotation */
#define ROTATE_FINAL_TIME 100       /* final time without rotation */
#define ROTATE_CHANGE_TIME 0.33     /* relative duration of acceleration/deceleration phases */
#define OMEGAMAX 100.0              /* maximal rotation speed */
#define PRINT_OMEGA 0               /* set to 1 to print angular speed */
#define PRINT_PARTICLE_SPEEDS 0     /* set to 1 to print average speeds/momenta of particles */
#define PRINT_SEGMENTS_SPEEDS 0     /* set to 1 to print velocity of moving segments */

#define MOVE_BOUNDARY 0        /* set to 1 to move repelling segments, due to force from particles */
#define SEGMENTS_MASS 40.0     /* mass of collection of segments */
#define DEACTIVATE_SEGMENT 1    /* set to 1 to deactivate last segment after a certain time */
#define SEGMENT_DEACTIVATION_TIME 200   /* time at which to deactivate last segment */
#define RELEASE_ROCKET_AT_DEACTIVATION 1    /* set to 1 to limit segments velocity before segment release */
#define SEGMENTS_X0 1.0        /* initial position of segments */
#define SEGMENTS_Y0 0.8        /* initial position of segments */
#define SEGMENTS_VX0 0.0       /* initial velocity of segments */
#define SEGMENTS_VY0 0.0      /* initial velocity of segments */
#define DAMP_SEGS_AT_NEGATIVE_Y 0   /* set to 1 to dampen segments when y coordinate is negative */

#define MOVE_SEGMENT_GROUPS 1       /* set to 1 to group segments into moving units */
#define SEGMENT_GROUP_MASS 750.0    /* mass of segment group */
#define SEGMENT_GROUP_I 500.0       /* moment of inertia of segment group */
#define SEGMENT_GROUP_DAMPING 0.0   /* damping of segment groups */
#define GROUP_REPULSION 1           /* set to 1 for groups of segments to repel each other */
#define KSPRING_GROUPS 1.0e11       /* harmonic potential between segment groups */
#define GROUP_WIDTH 0.05            /* interaction width of groups */
#define GROUP_G_REPEL 0             /* set to 1 to add repulsion between centers of mass of groups */
#define GROUP_G_REPEL_RADIUS 1.2    /* radius within which centers of mass of groups repel each other */
#define TRACK_SEGMENT_GROUPS 0      /* set to 1 for view to track group of segments */
#define TRACK_X_PADDING 2.0         /* distance from x boundary where tracking starts */

#define POSITION_DEPENDENT_TYPE 0   /* set to 1 to make particle type depend on initial position */
#define POSITION_Y_DEPENDENCE 0     /* set to 1 for the separation between particles to be horizontal */
#define PRINT_ENTROPY 0     /* set to 1 to compute entropy */

#define REACTION_DIFFUSION 0    /* set to 1 to simulate a chemical reaction (particles may change type) */
#define RD_TYPES 3              /* number of types in reaction-diffusion equation */
#define REACTION_PROB 0.0045      /* probability controlling reaction term */ 

#define PRINT_PARTICLE_NUMBER 0     /* set to 1 to print total number of particles */
#define PRINT_LEFT 0        /* set to 1 to print certain parameters at the top left instead of right */
#define PLOT_SPEEDS 0       /* set to 1 to add a plot of obstacle speeds (e.g. for rockets) */
#define PLOT_TRAJECTORIES 0     /* set to 1 to add a plot of obstacle trajectories (e.g. for rockets) */
#define VMAX_PLOT_SPEEDS 0.6    /* vertical scale of plot of obstacle speeds */

#define EHRENFEST_COPY 0    /* set to 1 to add equal number of larger particles (for Ehrenfest model) */

#define LID_MASS 1000.0     /* mass of lid for BC_RECTANGLE_LID b.c. */
#define LID_WIDTH 0.1       /* width of lid for BC_RECTANGLE_LID b.c. */
#define WALL_MASS 2000.0    /* mass of wall for BC_RECTANGLE_WALL b.c. */
#define WALL_FRICTION 0.0   /* friction on wall for BC_RECTANGLE_WALL b.c. */
#define WALL_WIDTH 0.1      /* width of wall for BC_RECTANGLE_WALL b.c. */
#define WALL_VMAX 100.0     /* max speed of wall */
#define WALL_TIME 0         /* time during which to keep wall */

#define NXMAZE 10      /* width of maze */
#define NYMAZE 10      /* height of maze */
#define MAZE_MAX_NGBH 4     /* max number of neighbours of maze cell */
#define RAND_SHIFT 200      /* seed of random number generator */
#define MAZE_XSHIFT 0.0     /* horizontal shift of maze */

#define FLOOR_FORCE 1      /* set to 1 to limit force on particle to FMAX */
#define FMAX 1.0e12         /* maximal force */
#define FLOOR_OMEGA 0      /* set to 1 to limit particle momentum to PMAX */
#define PMAX 1000.0        /* maximal force */

#define HASHX 120   /* size of hashgrid in x direction */
#define HASHY 60    /* size of hashgrid in y direction */
#define HASHMAX 100  /* maximal number of particles per hashgrid cell */
#define HASHGRID_PADDING 0.1    /* padding of hashgrid outside simulation window */

#define DRAW_COLOR_SCHEME 0     /* set to 1 to plot the color scheme */
#define COLORBAR_RANGE 8.0    /* scale of color scheme bar */
#define COLORBAR_RANGE_B 12.0    /* scale of color scheme bar for 2nd part */
#define ROTATE_COLOR_SCHEME 0   /* set to 1 to draw color scheme horizontally */

```

### 5 December 22 - Shooting laser beams at a 16 by 16 maze with slanted corners ###

**Program:** `particle_billiard.c` 

**Initial condition in function `animation()`:** `init_drop_config(-1.15, 0.01, 0.00*PI, 0.3*PI, configs); `

```
#define MOVIE 1         /* set to 1 to generate movie */

#define WINWIDTH 	1280  /* window width */
#define WINHEIGHT 	720   /* window height */

#define XMIN -2.0
#define XMAX 2.0	/* x interval */
#define YMIN -1.125
#define YMAX 1.125	/* y interval for 9/16 aspect ratio */

#define SCALING_FACTOR 1.0       /* scaling factor of drawing, needed for flower billiards, otherwise set to 1.0 */

/* Choice of the billiard table, see global_particles.c */

#define B_DOMAIN 30     /* choice of domain shape */

#define CIRCLE_PATTERN 1    /* pattern of circles */
#define POLYLINE_PATTERN 11  /* pattern of polyline */

#define ABSORBING_CIRCLES 0 /* set to 1 for circular scatterers to be absorbing */

#define NMAXCIRCLES 100000     /* total number of circles (must be at least NCX*NCY for square grid) */
#define NMAXPOLY 100000        /* total number of sides of polygonal line */   
#define NCX 30            /* number of circles in x direction */
#define NCY 20            /* number of circles in y direction */
#define NPOISSON 500        /* number of points for Poisson C_RAND_POISSON arrangement */
#define NGOLDENSPIRAL 2000  /* max number of points for C_GOLDEN_SPIRAL arrandement */
#define SDEPTH 1            /* Sierpinski gastket depth */

#define LAMBDA 1.5	/* parameter controlling shape of domain */
#define MU 0.5          /* second parameter controlling shape of billiard */
#define FOCI 1          /* set to 1 to draw focal points of ellipse */
#define NPOLY 6             /* number of sides of polygon */
#define APOLY 0.0           /* angle by which to turn polygon, in units of Pi/2 */ 
#define PENROSE_RATIO 2.5    /* parameter controlling the shape of small ellipses in Penrose room */

#define DRAW_BILLIARD 1     /* set to 1 to draw billiard */
#define DRAW_CONSTRUCTION_LINES 0   /* set to 1 to draw additional construction lines for billiard */
#define PERIODIC_BC 0       /* set to 1 to enforce periodic boundary conditions when drawing particles */

#define RESAMPLE 0      /* set to 1 if particles should be added when dispersion too large */
#define DEBUG 0         /* draw trajectories, for debugging purposes */

/* Simulation parameters */

#define NPART 10000     /* number of particles */
#define NPARTMAX 100000	/* maximal number of particles after resampling */
#define LMAX 0.01       /* minimal segment length triggering resampling */ 
#define DMIN 0.02       /* minimal distance to boundary for triggering resampling */ 
#define CYCLE 1         /* set to 1 for closed curve (start in all directions) */
#define SHOWTRAILS 1    /* set to 1 to keep trails of the particles */
#define SHOWZOOM 0      /* set to 1 to show zoom on specific area */
#define PRINT_PARTICLE_NUMBER 0 /* set to 1 to print number of particles */
#define PRINT_LEFT_RIGHT_PARTICLE_NUMBER 1 /* set to 1 to print number of particles on left and right side */
#define PRINT_COLLISION_NUMBER 0 /* set to 1 to print number of collisions */
#define TEST_ACTIVE 1   /* set to 1 to test whether particle is in billiard */

#define TEST_INITIAL_COND 0     /* set to 1 to allow only initial conditions that pass a test */

#define NSTEPS 4500     /* number of frames of movie */
#define TIME 1500        /* time between movie frames, for fluidity of real-time simulation */ 
#define DPHI 0.00002     /* integration step */
#define NVID 25          /* number of iterations between images displayed on screen */
#define END_FRAMES 25    /* number of still frames at the end of the movie */

/* Decreasing TIME accelerates the animation and the movie                               */
/* For constant speed of movie, TIME*DPHI should be kept constant                        */
/* However, increasing DPHI too much deterioriates quality of simulation                 */
/* NVID tells how often a picture is drawn in the animation, increase it for faster anim */
/* For a good quality movie, take for instance TIME = 400, DPHI = 0.00005, NVID = 100    */

/* Colors and other graphical parameters */

#define COLOR_PALETTE 10     /* Color palette, see list in global_pdes.c  */

#define NCOLORS 128      /* number of colors */
#define COLORSHIFT 0     /* hue of initial color */ 
#define RAINBOW_COLOR 0  /* set to 1 to use different colors for all particles */
#define FLOWER_COLOR 0   /* set to 1 to adapt initial colors to flower billiard (tracks vs core) */
#define NSEG 100         /* number of segments of boundary */
#define LENGTH 0.04       /* length of velocity vectors */
#define BILLIARD_WIDTH 2    /* width of billiard */
#define PARTICLE_WIDTH 2    /* width of particles */
#define FRONT_WIDTH 3       /* width of wave front */

#define BLACK 1             /* set to 1 for black background */
#define COLOR_OUTSIDE 0     /* set to 1 for colored outside */ 
#define OUTER_COLOR 270.0   /* color outside billiard */
#define PAINT_INT 0         /* set to 1 to paint interior in other color (for polygon/Reuleaux) */
#define PAINT_EXT 1         /* set to 1 to paint exterior */

#define PAUSE 1000       /* number of frames after which to pause */
#define PSLEEP 2         /* sleep time during pause */
#define SLEEP1  1        /* initial sleeping time */
#define SLEEP2  1       /* final sleeping time */

#define NXMAZE 16      /* width of maze */
#define NYMAZE 16      /* height of maze */
#define MAZE_MAX_NGBH 4     /* max number of neighbours of maze cell */
#define RAND_SHIFT 1        /* seed of random number generator */
#define MAZE_XSHIFT 0.0     /* horizontal shift of maze */

```

### 4 December 22 - 3D view of a wave refracted on a honeycomb lattice ###

**Program:** `wave_3d.c` 

**Initial condition in function `animation()`:** 
```
    init_circular_wave_mod(0.0, 0.0, phi, psi, xy_in);
    add_circular_wave_mod(sign, 0.0, 0.0, phi, psi, xy_in);
```

```
#define MOVIE 1         /* set to 1 to generate movie */
#define DOUBLE_MOVIE 1  /* set to 1 to produce movies for wave height and energy simultaneously */
#define SAVE_MEMORY 0   /* set to 1 to save memory when writing tiff images */

/* General geometrical parameters */

#define WINWIDTH 	1920  /* window width */
#define WINHEIGHT 	1000  /* window height */
#define NX 2500          /* number of grid points on x axis */
#define NY 1250          /* number of grid points on y axis */

#define XMIN -2.0
#define XMAX 2.0	/* x interval  */
#define YMIN -1.041666667
#define YMAX 1.041666667	/* y interval for 9/16 aspect ratio */

#define HIGHRES 1        /* set to 1 if resolution of grid is double that of displayed image */

#define JULIA_SCALE 0.8 /* scaling for Julia sets */

/* Choice of the billiard table */

#define B_DOMAIN 57        /* choice of domain shape, see list in global_pdes.c */

#define CIRCLE_PATTERN 2   /* pattern of circles or polygons, see list in global_pdes.c */

#define COMPARISON 0        /* set to 1 to compare two different patterns */
#define B_DOMAIN_B 20       /* second domain shape, for comparisons */
#define CIRCLE_PATTERN_B 0  /* second pattern of circles or polygons */

#define VARIABLE_IOR 0      /* set to 1 for a variable index of refraction */
#define IOR 2               /* choice of index of refraction, see list in global_pdes.c */
#define MANDEL_IOR_SCALE -0.05   /* parameter controlling dependence of IoR on Mandelbrot escape speed */

#define P_PERCOL 0.25       /* probability of having a circle in C_RAND_PERCOL arrangement */
#define NPOISSON 300        /* number of points for Poisson C_RAND_POISSON arrangement */
#define RANDOM_POLY_ANGLE 1 /* set to 1 to randomize angle of polygons */

#define LAMBDA 0.25	    /* parameter controlling the dimensions of domain */
#define MU 0.0              /* parameter controlling the dimensions of domain */
#define NPOLY 3             /* number of sides of polygon */
#define APOLY 2.0           /* angle by which to turn polygon, in units of Pi/2 */ 
#define MDEPTH 7            /* depth of computation of Menger gasket */
#define MRATIO 3            /* ratio defining Menger gasket */
#define MANDELLEVEL 1000    /* iteration level for Mandelbrot set */
#define MANDELLIMIT 10.0    /* limit value for approximation of Mandelbrot set */
#define FOCI 1              /* set to 1 to draw focal points of ellipse */
#define NGRIDX 12           /* number of grid point for grid of disks */
#define NGRIDY 16           /* number of grid point for grid of disks */

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
#define COURANT 0.05      /* Courant number */
#define COURANTB 0.1      /* Courant number in medium B */
#define GAMMA 0.0          /* damping factor in wave equation */
#define GAMMAB 0.0            /* damping factor in wave equation */
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
#define OSCILLATING_SOURCE_PERIOD 100   /* period of oscillating source */
#define ALTERNATE_OSCILLATING_SOURCE 1  /* set to 1 to alternate sign of oscillating source */

/* Boundary conditions, see list in global_pdes.c  */

#define B_COND 2

#define PRECOMPUTE_BC 0     /* set to 1 to compute neighbours for Laplacian in advance */

/* Parameters for length and speed of simulation */

#define NSTEPS 1400        /* number of frames of movie */
#define NVID 5            /* number of iterations between images displayed on screen */
#define NSEG 1000         /* number of segments of boundary */
#define INITIAL_TIME 0      /* time after which to start saving frames */
#define BOUNDARY_WIDTH 2    /* width of billiard boundary */
#define PRINT_SPEED 0       /* set to 1 to print speed of moving source */

#define PAUSE 100       /* number of frames after which to pause */
#define PSLEEP 2         /* sleep time during pause */
#define SLEEP1  1        /* initial sleeping time */
#define SLEEP2  1        /* final sleeping time */
#define MID_FRAMES 100    /* number of still frames between parts of two-part movie */
#define END_FRAMES 100   /* number of still frames at end of movie */
#define FADE 1           /* set to 1 to fade at end of movie */
#define ROTATE_VIEW_WHILE_FADE 1    /* set to 1 to keep rotating viewpoint during fade */

/* Parameters of initial condition */

#define INITIAL_AMP 0.25         /* amplitude of initial condition */
#define INITIAL_VARIANCE 0.0005  /* variance of initial condition */
#define INITIAL_WAVELENGTH  0.02  /* wavelength of initial condition */


/* Plot type, see list in global_pdes.c  */

#define ZPLOT 103     /* wave height */
#define CPLOT 103     /* color scheme */

#define ZPLOT_B 107 
#define CPLOT_B 107        /* plot type for second movie */

#define AMPLITUDE_HIGH_RES 1    /* set to 1 to increase resolution of plot */
#define SHADE_3D 1              /* set to 1 to change luminosity according to normal vector */
#define NON_DIRICHLET_BC 0      /* set to 1 to draw only facets in domain, if field is not zero on boundary */
#define FLOOR_ZCOORD 1          /* set to 1 to draw only facets with z not too negative */
#define DRAW_BILLIARD 1         /* set to 1 to draw boundary */
#define DRAW_BILLIARD_FRONT 1   /* set to 1 to draw front of boundary after drawing wave */
#define DRAW_CONSTRUCTION_LINES 0   /* set to 1 to draw construction lines of certain domains */
#define FADE_IN_OBSTACLE 1      /* set to 1 to fade color inside obstacles */
#define DRAW_OUTSIDE_GRAY 0     /* experimental, draw outside of billiard in gray */

#define PLOT_SCALE_ENERGY 0.01         /* vertical scaling in energy plot */
#define PLOT_SCALE_LOG_ENERGY 0.25       /* vertical scaling in log energy plot */

/* 3D representation */

#define REPRESENTATION_3D 1     /* choice of 3D representation */ 

#define REP_AXO_3D 0        /* linear projection (axonometry) */
#define REP_PROJ_3D 1       /* projection on plane orthogonal to observer line of sight */

#define ROTATE_VIEW 1       /* set to 1 to rotate position of observer */
#define ROTATE_ANGLE 360.0  /* total angle of rotation during simulation */

/* Color schemes */

#define COLOR_PALETTE 10      /* Color palette, see list in global_pdes.c  */
#define COLOR_PALETTE_B 14    /* Color palette, see list in global_pdes.c  */

#define BLACK 1          /* background */

#define COLOR_SCHEME 3   /* choice of color scheme, see list in global_pdes.c  */

#define SCALE 0          /* set to 1 to adjust color scheme to variance of field */
#define SLOPE 1.0        /* sensitivity of color on wave amplitude */
#define VSCALE_AMPLITUDE 1.0   /* additional scaling factor for color scheme P_3D_AMPLITUDE */
#define VSCALE_ENERGY 100.0     /* additional scaling factor for color scheme P_3D_ENERGY */
#define PHASE_FACTOR 20.0       /* factor in computation of phase in color scheme P_3D_PHASE */
#define PHASE_SHIFT 0.0      /* shift of phase in color scheme P_3D_PHASE */
#define ATTENUATION 0.0    /* exponential attenuation coefficient of contrast with time */
#define E_SCALE 100.0      /* scaling factor for energy representation */
#define LOG_SCALE 0.75     /* scaling factor for energy log representation */
#define LOG_SHIFT 0.5      /* shift of colors on log scale */
#define LOG_ENERGY_FLOOR -10.0    /* floor value for log of (total) energy */
#define LOG_MEAN_ENERGY_SHIFT 1.0   /* additional shift for log of mean energy */
#define FLUX_SCALE 1.0e4    /* scaling factor for enegy flux represtnation */
#define RESCALE_COLOR_IN_CENTER 0   /* set to 1 to decrease color intentiy in the center (for wave escaping ring) */

#define COLORHUE 260     /* initial hue of water color for scheme C_LUM */
#define COLORDRIFT 0.0   /* how much the color hue drifts during the whole simulation */
#define LUMMEAN 0.5      /* amplitude of luminosity variation for scheme C_LUM */
#define LUMAMP 0.3       /* amplitude of luminosity variation for scheme C_LUM */
#define HUEMEAN 240.0    /* mean value of hue for color scheme C_HUE */
#define HUEAMP -200.0    /* amplitude of variation of hue for color scheme C_HUE */

#define NXMAZE 7      /* width of maze */
#define NYMAZE 7      /* height of maze */
#define MAZE_MAX_NGBH 4     /* max number of neighbours of maze cell */
#define RAND_SHIFT 24        /* seed of random number generator */
#define MAZE_XSHIFT 0.0     /* horizontal shift of maze */

#define DRAW_COLOR_SCHEME 1       /* set to 1 to plot the color scheme */
#define COLORBAR_RANGE 3.5      /* scale of color scheme bar */
#define COLORBAR_RANGE_B 3.0    /* scale of color scheme bar for 2nd part */
#define ROTATE_COLOR_SCHEME 0     /* set to 1 to draw color scheme horizontally */

#define SAVE_TIME_SERIES 0      /* set to 1 to save wave time series at a point */

/* for compatibility with sub_wave and sub_maze */
#define ADD_POTENTIAL 0
#define POTENTIAL 0
/* end of constants only used by sub_wave and sub_maze */


/* For debugging purposes only */
#define FLOOR 1         /* set to 1 to limit wave amplitude to VMAX */
#define VMAX 10.0       /* max value of wave amplitude */

/* Parameters controlling 3D projection */

double u_3d[2] = {0.75, -0.45};     /* projections of basis vectors for REP_AXO_3D representation */
double v_3d[2] = {-0.75, -0.45};
double w_3d[2] = {0.0, 0.015};
double light[3] = {0.816496581, -0.40824829, 0.40824829};      /* vector of "light" direction for P_3D_ANGLE color scheme */
double observer[3] = {8.0, 8.0, 8.0};    /* location of observer for REP_PROJ_3D representation */ 
int reset_view = 0;         /* switch to reset 3D view parameters (for option ROTATE_VIEW) */

#define Z_SCALING_FACTOR 0.3    /* overall scaling factor of z axis for REP_PROJ_3D representation */
#define XY_SCALING_FACTOR 2.4   /* overall scaling factor for on-screen (x,y) coordinates after projection */
#define ZMAX_FACTOR 1.0         /* max value of z coordinate for REP_PROJ_3D representation */
#define XSHIFT_3D -0.1          /* overall x shift for REP_PROJ_3D representation */
#define YSHIFT_3D 0.1           /* overall y shift for REP_PROJ_3D representation */

```

### 3 December 22 - A shear flow with an airfoil-shaped obstacle ###

**Program:** `rde.c` 

**Initial condition in function `animation()`:** `init_shear_flow(-1.0, 0.0, 0.1, 1, 1, 0.0, phi, xy_in);`

```
#define MOVIE 1         /* set to 1 to generate movie */
#define DOUBLE_MOVIE 1  /* set to 1 to produce movies for wave height and energy simultaneously */
#define SAVE_MEMORY 0           /* set to 1 to save memory when writing tiff images */

/* General geometrical parameters */

#define WINWIDTH 	1920  /* window width */
#define WINHEIGHT 	1000  /* window height */
#define NX 960          /* number of grid points on x axis */
#define NY 500          /* number of grid points on y axis */
 
#define XMIN -2.0
#define XMAX 2.0	/* x interval */
#define YMIN -1.041666667
#define YMAX 1.041666667	/* y interval for 9/16 aspect ratio */

/* Choice of simulated equation */

#define RDE_EQUATION 6  /* choice of reaction term, see list in global_3d.c */
#define NFIELDS 2       /* number of fields in reaction-diffusion equation */
#define NLAPLACIANS 1   /* number of fields for which to compute Laplacian */

#define ADD_POTENTIAL 0 /* set to 1 to add a potential (for Schrodinger equation) */
#define ADD_MAGNETIC_FIELD 0    /* set to 1 to add a magnetic field (for Schrodinger equation) - then set POTENTIAL 1 */
#define POTENTIAL 7     /* type of potential or vector potential, see list in global_3d.c  */

#define ANTISYMMETRIZE_WAVE_FCT 0   /* set tot 1 to make wave function antisymmetric */

#define JULIA_SCALE 0.5 /* scaling for Julia sets */

/* Choice of the billiard table */

#define B_DOMAIN 198       /* choice of domain shape, see list in global_pdes.c  */

#define CIRCLE_PATTERN 99    /* pattern of circles, see list in global_pdes.c */

#define P_PERCOL 0.25       /* probability of having a circle in C_RAND_PERCOL arrangement */
#define NPOISSON 300        /* number of points for Poisson C_RAND_POISSON arrangement */
#define RANDOM_POLY_ANGLE 0 /* set to 1 to randomize angle of polygons */

#define LAMBDA 0.7	    /* parameter controlling the dimensions of domain */
#define MU 0.15	            /* parameter controlling the dimensions of domain */
#define NPOLY 5             /* number of sides of polygon */
#define APOLY 2.0          /* angle by which to turn polygon, in units of Pi/2 */
#define MDEPTH 7            /* depth of computation of Menger gasket */
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

#define DT 0.0000012

#define VISCOSITY 2.0

#define RPSA 0.75         /* parameter in Rock-Paper-Scissors-type interaction */
#define RPSLZB 0.75       /* second parameter in Rock-Paper-Scissors-Lizard-Spock type interaction */

#define EPSILON 0.8     /* time scale separation */
#define DELTA 0.1       /* time scale separation */
#define FHNA 1.0        /* parameter in FHN equation */
#define FHNC -0.01      /* parameter in FHN equation */
#define K_HARMONIC 1.0  /* spring constant of harmonic potential */
#define K_COULOMB 0.5   /* constant in Coulomb potential */
#define V_MAZE 0.4      /* potential in walls of maze */
#define BZQ 0.0008      /* parameter in BZ equation */
#define BZF 1.2         /* parameter in BZ equation */
#define B_FIELD 10.0    /* magnetic field */
#define AB_RADIUS 0.2   /* radius of region with magnetic field for Aharonov-Bohm effect */
#define K_EULER 50.0    /* constant in stream function integration of Euler equation */

#define SMOOTHEN_VORTICITY 1    /* set to 1 to smoothen vorticity field in Euler equation */
#define SMOOTHEN_PERIOD 10      /* period between smoothenings */
#define SMOOTH_FACTOR 0.03      /* factor by which to smoothen */

#define ADD_TRACERS 1    /* set to 1 to add tracer particles (for Euler equations) */
#define N_TRACERS 1000    /* number of tracer particles */

#define T_OUT 2.0       /* outside temperature */
#define T_IN 0.0        /* inside temperature */
#define SPEED 0.0       /* speed of drift to the right */

#define ADD_NOISE 0     /* set to 1 to add noise, set to 2 to add noise in right half */
#define NOISE_INTENSITY 0.005      /* noise intensity */
#define CHANGE_NOISE 1      /* set to 1 to increase noise intensity */
#define NOISE_FACTOR 40.0   /* factor by which to increase noise intensity */
#define NOISE_INITIAL_TIME 100  /* initial time during which noise remains constant */

#define CHANGE_VISCOSITY 0      /* set to 1 to change the viscosity in the course of the simulation */
#define ADJUST_INTSTEP 0       /* set to 1 to decrease integration step when viscosity increases */
#define VISCOSITY_INITIAL_TIME 10  /* initial time during which viscosity remains constant */
#define VISCOSITY_FACTOR 100.0   /* factor by which to change viscosity */
#define VISCOSITY_MAX 2.0        /* max value of viscosity beyond which NVID is increased and integration step is decrase, 
                                    for numerical stability */
                                        
#define CHANGE_RPSLZB 0         /* set to 1 to change second parameter in Rock-Paper-Scissors-Lizard-Spock equation */
#define RPSLZB_CHANGE 0.75      /* factor by which to rpslzb parameter */
#define RPSLZB_INITIAL_TIME 0   /* initial time during which rpslzb remains constant */
#define RPSLZB_FINAL_TIME 500   /* final time during which rpslzb remains constant */
                                      
/* Boundary conditions, see list in global_pdes.c  */

#define B_COND 1

#define EULER_GRADIENT_YSHIFT 0.0    /* y-shift in computation of gradient in Euler equation */

/* Parameters for length and speed of simulation */

#define NSTEPS 2650       /* number of frames of movie */
#define NVID 120          /* number of iterations between images displayed on screen */
#define ACCELERATION_FACTOR 1.0 /* factor by which to increase NVID in course of simulation */
#define DT_ACCELERATION_FACTOR 1.0 /* factor by which to increase time step in course of simulation  */
#define MAX_DT 0.024     /* maximal value of integration step */
#define NSEG 100         /* number of segments of boundary */
#define BOUNDARY_WIDTH 5    /* width of billiard boundary */

#define PAUSE 100       /* number of frames after which to pause */
#define PSLEEP 2         /* sleep time during pause */
#define SLEEP1  2        /* initial sleeping time */
#define SLEEP2  1        /* final sleeping time */
#define INITIAL_TIME 0  /* initial still time */
#define MID_FRAMES 50    /* number of still frames between parts of two-part movie */
#define END_FRAMES 50    /* number of still frames at end of movie */
#define FADE 1           /* set to 1 to fade at end of movie */

/* Visualisation */

#define PLOT_3D 0    /* controls whether plot is 2D or 3D */

#define ROTATE_VIEW 0       /* set to 1 to rotate position of observer */
#define ROTATE_ANGLE 360.0  /* total angle of rotation during simulation */

#define DRAW_PERIODICISED 1     /* set to 1 to repeat wave periodically in x and y directions */

/* Plot type - color scheme */

#define CPLOT 52
#define CPLOT_B 51

/* Plot type - height of 3D plot */

#define ZPLOT 52     /* z coordinate in 3D plot */
#define ZPLOT_B 51    /* z coordinate in second 3D plot */

#define AMPLITUDE_HIGH_RES 1    /* set to 1 to increase resolution of P_3D_AMPLITUDE plot */
#define SHADE_3D 1              /* set to 1 to change luminosity according to normal vector */
#define NON_DIRICHLET_BC 0      /* set to 1 to draw only facets in domain, if field is not zero on boundary */
#define WRAP_ANGLE 1            /* experimental: wrap angle to [0, 2Pi) for interpolation in angle schemes */
#define FADE_IN_OBSTACLE 0      /* set to 1 to fade color inside obstacles */
#define DRAW_OUTSIDE_GRAY 0     /* experimental - draw outside of billiard in gray */
#define ADD_POTENTIAL_TO_Z 1    /* set to 1 to add the external potential to z-coordinate of plot */
#define ADD_POT_CONSTANT 0.35   /* constant in front of added potential */

#define PLOT_SCALE_ENERGY 0.05      /* vertical scaling in energy plot */

#define PRINT_TIME 0        /* set to 1 to print running time */
#define PRINT_VISCOSITY 0   /* set to 1 to print viscosity */
#define PRINT_RPSLZB 0      /* set to 1 to print rpslzb parameter */
#define PRINT_PROBABILITIES 0   /* set to 1 to print probabilities (for Ehrenfest urn configuration) */
#define PRINT_NOISE 0       /* set to 1 to print noise intensity */

#define DRAW_FIELD_LINES 0  /* set to 1 to draw field lines */
#define FIELD_LINE_WIDTH 1  /* width of field lines */
#define N_FIELD_LINES 120   /* number of field lines */
#define FIELD_LINE_FACTOR 120 /* factor controlling precision when computing origin of field lines */
#define DRAW_BILLIARD 1     /* set to 1 to draw boundary */
#define DRAW_BILLIARD_FRONT 1     /* set to 1 to draw boundary */
#define FILL_BILLIARD_COMPLEMENT 1  /* set to 1 to fill complement of billiard (for certain shapes only) */

/* 3D representation */

#define REPRESENTATION_3D 1     /* choice of 3D representation */ 

#define REP_AXO_3D 0        /* linear projection (axonometry) */
#define REP_PROJ_3D 1       /* projection on plane orthogonal to observer line of sight */

/* Color schemes, see list in global_pdes.c  */

#define COLOR_PALETTE 10       /* Color palette, see list in global_pdes.c  */
#define COLOR_PALETTE_B 17     /* Color palette, see list in global_pdes.c  */

#define BLACK 1          /* black background */

#define COLOR_SCHEME 3   /* choice of color scheme */

#define COLOR_PHASE_SHIFT 0.0   /* phase shift of color scheme, in units of Pi */

#define SCALE 0          /* set to 1 to adjust color scheme to variance of field */
#define SLOPE 1.0        /* sensitivity of color on wave amplitude */
#define VSCALE_AMPLITUDE 1.5      /* additional scaling factor for color scheme P_3D_AMPLITUDE */
#define ATTENUATION 0.0  /* exponential attenuation coefficient of contrast with time */
#define CURL_SCALE 0.000015   /* scaling factor for curl representation */
#define RESCALE_COLOR_IN_CENTER 0   /* set to 1 to decrease color intentiy in the center (for wave escaping ring) */
#define SLOPE_SCHROD_LUM 50.0       /* sensitivity of luminosity on module, for color scheme Z_ARGUMENT */
#define MIN_SCHROD_LUM 0.2       /* minimal luminosity in color scheme Z_ARGUMENT*/

#define COLORHUE 260     /* initial hue of water color for scheme C_LUM */
#define COLORDRIFT 0.0   /* how much the color hue drifts during the whole simulation */
#define LUMMEAN 0.5      /* amplitude of luminosity variation for scheme C_LUM */
#define LUMAMP 0.3       /* amplitude of luminosity variation for scheme C_LUM */
#define HUEMEAN 359.0    /* mean value of hue for color scheme C_HUE */
#define HUEAMP -359.0    /* amplitude of variation of hue for color scheme C_HUE */
#define E_SCALE 100.0    /* scaling factor for energy representation */
#define FLUX_SCALE 100.0 /* scaling factor for energy representation */
#define LOG_SCALE 0.75   /* scaling factor for energy log representation */
#define LOG_SHIFT 1.0 
#define LOG_MIN 1.0e-3   /* floor value for log vorticity plot */

#define NXMAZE 7      /* width of maze */
#define NYMAZE 7      /* height of maze */
#define MAZE_MAX_NGBH 4     /* max number of neighbours of maze cell */
#define RAND_SHIFT 0        /* seed of random number generator */
#define MAZE_XSHIFT 0.0     /* horizontal shift of maze */

#define DRAW_COLOR_SCHEME 1     /* set to 1 to plot the color scheme */
#define COLORBAR_RANGE 3.0      /* scale of color scheme bar */
#define COLORBAR_RANGE_B 3.0   /* scale of color scheme bar for 2nd part */
#define ROTATE_COLOR_SCHEME 0   /* set to 1 to draw color scheme horizontally */

/* only for compatibility with wave_common.c */
#define TWOSPEEDS 0          /* set to 1 to replace hardcore boundary by medium with different speed */
#define OMEGA 0.005        /* frequency of periodic excitation */
#define COURANT 0.08       /* Courant number */
#define COURANTB 0.03      /* Courant number in medium B */
#define INITIAL_AMP 0.5         /* amplitude of initial condition */
#define INITIAL_VARIANCE 0.0002  /* variance of initial condition */
#define INITIAL_WAVELENGTH  0.1  /* wavelength of initial condition */
#define VSCALE_ENERGY 200.0       /* additional scaling factor for color scheme P_3D_ENERGY */
#define PHASE_FACTOR 20.0       /* factor in computation of phase in color scheme P_3D_PHASE */
#define PHASE_SHIFT 0.0      /* shift of phase in color scheme P_3D_PHASE */
#define OSCILLATION_SCHEDULE 0  /* oscillation schedule, see list in global_pdes.c */
#define AMPLITUDE 0.8      /* amplitude of periodic excitation */ 
#define ACHIRP 0.2        /* acceleration coefficient in chirp */
#define DAMPING 0.0        /* damping of periodic excitation */
#define COMPARISON 0        /* set to 1 to compare two different patterns (beta) */
#define B_DOMAIN_B 20       /* second domain shape, for comparisons */
#define CIRCLE_PATTERN_B 0  /* second pattern of circles or polygons */
/* end of constants added only for compatibility with wave_common.c */


double u_3d[2] = {0.75, -0.45};     /* projections of basis vectors for REP_AXO_3D representation */
double v_3d[2] = {-0.75, -0.45};
double w_3d[2] = {0.0, 0.015};
double light[3] = {0.816496581, -0.40824829, 0.40824829};      /* vector of "light" direction for P_3D_ANGLE color scheme */
double observer[3] = {8.0, 8.0, 8.0};    /* location of observer for REP_PROJ_3D representation */ 
int reset_view = 0;         /* switch to reset 3D view parameters (for option ROTATE_VIEW) */

#define Z_SCALING_FACTOR 0.08  /* overall scaling factor of z axis for REP_PROJ_3D representation */
#define XY_SCALING_FACTOR 1.7  /* overall scaling factor for on-screen (x,y) coordinates after projection */
#define ZMAX_FACTOR 1.0        /* max value of z coordinate for REP_PROJ_3D representation */
#define XSHIFT_3D -0.1         /* overall x shift for REP_PROJ_3D representation */
#define YSHIFT_3D 0.1          /* overall y shift for REP_PROJ_3D representation */
#define BORDER_PADDING 0       /* distance from boundary at which to plot points, to avoid boundary effects due to gradient */

/* For debugging purposes only */
#define FLOOR 1         /* set to 1 to limit wave amplitude to VMAX */
#define VMAX 10.0        /* max value of wave amplitude */
#define TEST_GRADIENT 0 /* print norm squared of gradient */

```

### 2 December 22 - Deep inpact 5: Comparing a round and a rectangular obstacle ###

**Program:** `lennardjones.c` 

```
#define MOVIE 1         /* set to 1 to generate movie */
#define DOUBLE_MOVIE 1  /* set to 1 to produce movies for wave height and energy simultaneously */

#define TIME_LAPSE 0     /* set to 1 to add a time-lapse movie at the end */
                         /* so far incompatible with double movie */
#define TIME_LAPSE_FACTOR 3    /* factor of time-lapse movie */
#define TIME_LAPSE_FIRST 0  /* set to 1 to show time-lapse version first */

/* General geometrical parameters */

#define WINWIDTH 	1280  /* window width */
#define WINHEIGHT 	720   /* window height */

#define XMIN -2.0
#define XMAX 2.0	/* x interval */
#define YMIN -1.125
#define YMAX 1.125	/* y interval for 9/16 aspect ratio */

#define INITXMIN -1.95
#define INITXMAX 1.95	/* x interval for initial condition */
#define INITYMIN -1.1
#define INITYMAX 0.15	/* y interval for initial condition */

#define BCXMIN -2.0
#define BCXMAX 2.0	/* x interval for boundary condition */
#define BCYMIN -1.125
#define BCYMAX 1.125	/* y interval for boundary condition */

#define OBSXMIN -2.0
#define OBSXMAX 2.0     /* x interval for motion of obstacle */

#define CIRCLE_PATTERN 8   /* pattern of circles, see list in global_ljones.c */

#define ADD_FIXED_OBSTACLES 0   /* set to 1 do add fixed circular obstacles */
#define OBSTACLE_PATTERN 3  /* pattern of obstacles, see list in global_ljones.c */

#define ADD_FIXED_SEGMENTS 1    /* set to 1 to add fixed segments as obstacles */
#define SEGMENT_PATTERN 19    /* pattern of repelling segments, see list in global_ljones.c */
#define ROCKET_SHAPE 2        /* shape of rocket combustion chamber, see list in global_ljones.c */
#define ROCKET_SHAPE_B 2      /* shape of second rocket */
#define NOZZLE_SHAPE 1        /* shape of nozzle, see list in global_ljones.c */
#define NOZZLE_SHAPE_B 0      /* shape of nozzle for second rocket, see list in global_ljones.c */

#define TWO_TYPES 0         /* set to 1 to have two types of particles */
#define TPYE_PROPORTION 0.5 /* proportion of particles of first type */
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
#define PDISC_DISTANCE 3.0  /* minimal distance in Poisson disc process, controls density of particles */
#define PDISC_CANDIDATES 100 /* number of candidates in construction of Poisson disc process */
#define RANDOM_POLY_ANGLE 0 /* set to 1 to randomize angle of polygons */

#define LAMBDA 0.3	    /* parameter controlling the dimensions of domain */
#define MU 0.009 	    /* parameter controlling radius of particles */
#define MU_B 0.018          /* parameter controlling radius of particles of second type */
#define NPOLY 25            /* number of sides of polygon */
#define APOLY 0.666666666   /* angle by which to turn polygon, in units of Pi/2 */ 
#define MDEPTH 4            /* depth of computation of Menger gasket */
#define MRATIO 3            /* ratio defining Menger gasket */
#define MANDELLEVEL 1000    /* iteration level for Mandelbrot set */
#define MANDELLIMIT 10.0    /* limit value for approximation of Mandelbrot set */
#define FOCI 1              /* set to 1 to draw focal points of ellipse */
#define NGRIDX 90           /* number of grid point for grid of disks */
#define NGRIDY 100           /* number of grid point for grid of disks */
#define EHRENFEST_RADIUS 0.9    /* radius of container for Ehrenfest urn configuration */
#define EHRENFEST_WIDTH 0.035     /* width of tube for Ehrenfest urn configuration */
#define TWO_CIRCLES_RADIUS_RATIO 0.8    /* ratio of radii for S_TWO_CIRCLES_EXT segment configuration */
#define DAM_WIDTH 0.05       /* width of dam for S_DAM segment configuration */

#define X_SHOOTER -0.2
#define Y_SHOOTER -0.6
#define X_TARGET 0.4
#define Y_TARGET 0.7        /* shooter and target positions in laser fight */

/* Parameters for length and speed of simulation */

#define NSTEPS 1700      /* number of frames of movie */
#define NVID 100         /* number of iterations between images displayed on screen */
#define NSEG 250         /* number of segments of boundary */
#define INITIAL_TIME 250     /* time after which to start saving frames */
#define OBSTACLE_INITIAL_TIME 350     /* time after which to start moving obstacle */
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
#define PLOT_B 8        /* plot type for second movie */

#define DRAW_BONDS 0    /* set to 1 to draw bonds between neighbours */
#define COLOR_BONDS 1   /* set to 1 to color bonds according to length */
#define FILL_TRIANGLES 1    /* set to 1 to fill triangles between neighbours */
#define ALTITUDE_LINES 0    /* set to 1 to add horizontal lines to show altitude */
#define COLOR_SEG_GROUPS 1  /* set to 1 to collor segment groups differently */

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

#define ENERGY_HUE_MIN 330.0        /* color of original particle */
#define ENERGY_HUE_MAX 50.0         /* color of saturated particle */
#define PARTICLE_HUE_MIN 359.0      /* color of original particle */
#define PARTICLE_HUE_MAX 0.0        /* color of saturated particle */
#define PARTICLE_EMAX 3.0e2         /* energy of particle with hottest color */
#define HUE_TYPE0 70.0     /* hue of particles of type 0 */
#define HUE_TYPE1 280.0      /* hue of particles of type 1 */
#define HUE_TYPE2 .0      /* hue of particles of type 2 */
#define HUE_TYPE3 210.0     /* hue of particles of type 3 */

#define RANDOM_RADIUS 0     /* set to 1 for random circle radius */
#define DT_PARTICLE 3.0e-6    /* time step for particle displacement */
#define KREPEL 12.0          /* constant in repelling force between particles */
#define EQUILIBRIUM_DIST 4.0    /* Lennard-Jones equilibrium distance */
#define EQUILIBRIUM_DIST_B 3.5  /* Lennard-Jones equilibrium distance for second type of particle */
#define REPEL_RADIUS 20.0    /* radius in which repelling force acts (in units of particle radius) */
#define DAMPING 10.0          /* damping coefficient of particles */
#define INITIAL_DAMPING 35.0  /* damping coefficient of particles during initial phase */
#define PARTICLE_MASS 1.0    /* mass of particle of radius MU */
#define PARTICLE_MASS_B 1.0   /* mass of particle of radius MU */
#define PARTICLE_INERTIA_MOMENT 0.2     /* moment of inertia of particle */
#define PARTICLE_INERTIA_MOMENT_B 0.02     /* moment of inertia of second type of particle */
#define V_INITIAL 0.0        /* initial velocity range */
#define OMEGA_INITIAL 10.0        /* initial angular velocity range */

#define THERMOSTAT 0        /* set to 1 to switch on thermostat */
#define VARY_THERMOSTAT 0   /* set to 1 for time-dependent thermostat schedule */
#define SIGMA 5.0           /* noise intensity in thermostat */
#define BETA 0.02          /* initial inverse temperature */
#define MU_XI 0.01           /* friction constant in thermostat */
#define KSPRING_BOUNDARY 1.0e7    /* confining harmonic potential outside simulation region */
#define KSPRING_OBSTACLE 1.0e11    /* harmonic potential of obstacles */
#define NBH_DIST_FACTOR 7.5       /* radius in which to count neighbours */
#define GRAVITY 1500.0            /* gravity acting on all particles */
#define GRAVITY_X 0.0        /* horizontal gravity acting on all particles */
#define INCREASE_GRAVITY 0     /* set to 1 to increase gravity during the simulation */
#define GRAVITY_SCHEDULE 2     /* type of gravity schedule, see list in global_ljones.c */
#define GRAVITY_FACTOR 100.0    /* factor by which to increase gravity */
#define GRAVITY_INITIAL_TIME 200    /* time at start of simulation with constant gravity */
#define GRAVITY_RESTORE_TIME 700    /* time at end of simulation with gravity restored to initial value */

#define ROTATION 0           /* set to 1 to include rotation of particles */
#define COUPLE_ANGLE_TO_THERMOSTAT 0    /* set to 1 to couple angular degrees of freedom to thermostat */
#define DIMENSION_FACTOR 1.0  /* scaling factor taking into account number of degrees of freedom */  
#define KTORQUE 100.0          /* force constant in angular dynamics */
#define KTORQUE_B 10.0        /* force constant in angular dynamics */
#define KTORQUE_DIFF 150.0    /* force constant in angular dynamics for different particles */
#define DRAW_SPIN 0           /* set to 1 to draw spin vectors of particles */
#define DRAW_SPIN_B 0         /* set to 1 to draw spin vectors of particles */
#define DRAW_CROSS 1          /* set to 1 to draw cross on particles of second type */
#define SPIN_RANGE 7.0       /* range of spin-spin interaction */
#define SPIN_RANGE_B 5.0     /* range of spin-spin interaction for second type of particle */
#define QUADRUPOLE_RATIO 0.6  /* anisotropy in quadrupole potential */ 

#define INCREASE_BETA 0  /* set to 1 to increase BETA during simulation */
#define BETA_FACTOR 0.025   /* factor by which to change BETA during simulation */
#define N_TOSCILLATIONS 1.5   /* number of temperature oscillations in BETA schedule */
#define NO_OSCILLATION 1        /* set to 1 to have exponential BETA change only */
#define MIDDLE_CONSTANT_PHASE 1600   /* final phase in which temperature is constant */
#define FINAL_DECREASE_PHASE 1500    /* final phase in which temperature decreases */ 
#define FINAL_CONSTANT_PHASE -1     /* final phase in which temperature is constant */

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
#define PARTIAL_THERMO_REGION 1     /* region for partial thermostat coupling (see list in global_ljones.c) */
#define PARTIAL_THERMO_SHIFT 0.2    /* distance from obstacle at the right of which particles are coupled to thermostat */
#define PARTIAL_THERMO_WIDTH 0.5    /* vertical size of partial thermostat coupling */
#define PARTIAL_THERMO_HEIGHT 0.2   /* vertical size of partial thermostat coupling */

#define INCREASE_KREPEL 0   /* set to 1 to increase KREPEL during simulation */
#define KREPEL_FACTOR 1000.0   /* factor by which to change KREPEL during simulation */

#define PART_AT_BOTTOM 0     /* set to 1 to include "seed" particles at bottom */
#define MASS_PART_BOTTOM 10000.0 /* mass of particles at bottom */
#define NPART_BOTTOM 100     /* number of particles at the bottom */

#define ADD_PARTICLES 0    /* set to 1 to add particles */
#define ADD_TIME 0       /* time at which to add first particle */
#define ADD_PERIOD 10000       /* time interval between adding further particles */
#define FINAL_NOADD_PERIOD 200  /* final period where no particles are added */
#define SAFETY_FACTOR 2.0  /* no particles are added at distance less than MU*SAFETY_FACTOR of other particles */

#define TRACER_PARTICLE 0   /* set to 1 to have a tracer particle */
#define N_TRACER_PARTICLES 3    /* number of tracer particles */
#define TRAJECTORY_LENGTH 8000   /* length of recorded trajectory */
#define TRACER_PARTICLE_MASS 4.0    /* relative mass of tracer particle */
#define TRAJECTORY_WIDTH 3      /* width of tracer particle trajectory */

#define ROTATE_BOUNDARY 0           /* set to 1 to rotate the repelling segments */
#define SMOOTH_ROTATION 1           /* set to 1 to update segments at each time step (rather than at each movie frame) */
#define PERIOD_ROTATE_BOUNDARY 1000  /* period of rotating boundary */
#define ROTATE_INITIAL_TIME 0       /* initial time without rotation */
#define ROTATE_FINAL_TIME 100       /* final time without rotation */
#define ROTATE_CHANGE_TIME 0.33     /* relative duration of acceleration/deceleration phases */
#define OMEGAMAX 100.0              /* maximal rotation speed */
#define PRINT_OMEGA 0               /* set to 1 to print angular speed */
#define PRINT_PARTICLE_SPEEDS 0     /* set to 1 to print average speeds/momenta of particles */
#define PRINT_SEGMENTS_SPEEDS 0     /* set to 1 to print velocity of moving segments */

#define MOVE_BOUNDARY 0        /* set to 1 to move repelling segments, due to force from particles */
#define SEGMENTS_MASS 40.0     /* mass of collection of segments */
#define DEACTIVATE_SEGMENT 1    /* set to 1 to deactivate last segment after a certain time */
#define SEGMENT_DEACTIVATION_TIME 200   /* time at which to deactivate last segment */
#define RELEASE_ROCKET_AT_DEACTIVATION 1    /* set to 1 to limit segments velocity before segment release */
#define SEGMENTS_X0 1.0        /* initial position of segments */
#define SEGMENTS_Y0 0.8        /* initial position of segments */
#define SEGMENTS_VX0 0.0       /* initial velocity of segments */
#define SEGMENTS_VY0 0.0      /* initial velocity of segments */
#define DAMP_SEGS_AT_NEGATIVE_Y 0   /* set to 1 to dampen segments when y coordinate is negative */

#define MOVE_SEGMENT_GROUPS 1       /* set to 1 to group segments into moving units */
#define SEGMENT_GROUP_MASS 750.0    /* mass of segment group */
#define SEGMENT_GROUP_I 500.0       /* moment of inertia of segment group */
#define SEGMENT_GROUP_DAMPING 0.0   /* damping of segment groups */
#define GROUP_REPULSION 1           /* set to 1 for groups of segments to repel each other */
#define KSPRING_GROUPS 1.0e11       /* harmonic potential between segment groups */
#define GROUP_WIDTH 0.05            /* interaction width of groups */
#define GROUP_G_REPEL 0             /* set to 1 to add repulsion between centers of mass of groups */
#define GROUP_G_REPEL_RADIUS 1.2    /* radius within which centers of mass of groups repel each other */
#define TRACK_SEGMENT_GROUPS 0      /* set to 1 for view to track group of segments */
#define TRACK_X_PADDING 2.0         /* distance from x boundary where tracking starts */

#define POSITION_DEPENDENT_TYPE 0   /* set to 1 to make particle type depend on initial position */
#define POSITION_Y_DEPENDENCE 0     /* set to 1 for the separation between particles to be horizontal */
#define PRINT_ENTROPY 0     /* set to 1 to compute entropy */

#define REACTION_DIFFUSION 0    /* set to 1 to simulate a chemical reaction (particles may change type) */
#define RD_TYPES 3              /* number of types in reaction-diffusion equation */
#define REACTION_PROB 0.0045      /* probability controlling reaction term */ 

#define PRINT_PARTICLE_NUMBER 0     /* set to 1 to print total number of particles */
#define PRINT_LEFT 0        /* set to 1 to print certain parameters at the top left instead of right */
#define PLOT_SPEEDS 0       /* set to 1 to add a plot of obstacle speeds (e.g. for rockets) */
#define PLOT_TRAJECTORIES 0     /* set to 1 to add a plot of obstacle trajectories (e.g. for rockets) */
#define VMAX_PLOT_SPEEDS 0.6    /* vertical scale of plot of obstacle speeds */

#define EHRENFEST_COPY 0    /* set to 1 to add equal number of larger particles (for Ehrenfest model) */

#define LID_MASS 1000.0     /* mass of lid for BC_RECTANGLE_LID b.c. */
#define LID_WIDTH 0.1       /* width of lid for BC_RECTANGLE_LID b.c. */
#define WALL_MASS 2000.0    /* mass of wall for BC_RECTANGLE_WALL b.c. */
#define WALL_FRICTION 0.0   /* friction on wall for BC_RECTANGLE_WALL b.c. */
#define WALL_WIDTH 0.1      /* width of wall for BC_RECTANGLE_WALL b.c. */
#define WALL_VMAX 100.0     /* max speed of wall */
#define WALL_TIME 0         /* time during which to keep wall */

#define NXMAZE 10      /* width of maze */
#define NYMAZE 10      /* height of maze */
#define MAZE_MAX_NGBH 4     /* max number of neighbours of maze cell */
#define RAND_SHIFT 200      /* seed of random number generator */
#define MAZE_XSHIFT 0.0     /* horizontal shift of maze */

#define FLOOR_FORCE 1      /* set to 1 to limit force on particle to FMAX */
#define FMAX 1.0e12         /* maximal force */
#define FLOOR_OMEGA 0      /* set to 1 to limit particle momentum to PMAX */
#define PMAX 1000.0        /* maximal force */

#define HASHX 120   /* size of hashgrid in x direction */
#define HASHY 60    /* size of hashgrid in y direction */
#define HASHMAX 100  /* maximal number of particles per hashgrid cell */
#define HASHGRID_PADDING 0.1    /* padding of hashgrid outside simulation window */

#define DRAW_COLOR_SCHEME 0     /* set to 1 to plot the color scheme */
#define COLORBAR_RANGE 8.0    /* scale of color scheme bar */
#define COLORBAR_RANGE_B 12.0    /* scale of color scheme bar for 2nd part */
#define ROTATE_COLOR_SCHEME 0   /* set to 1 to draw color scheme horizontally */

```

### 1 December 22 - 45 degrees are not always optimal: Solving a 16 by 32 maze with laser beams ###

**Program:** `particle_billiard.c` 

**Initial condition in function `animation()`:** `init_drop_config(-1.15, 0.01, 0.00*PI, 0.3*PI, configs);`

```
#define MOVIE 1         /* set to 1 to generate movie */

#define WINWIDTH 	1280  /* window width */
#define WINHEIGHT 	720   /* window height */

#define XMIN -2.0
#define XMAX 2.0	/* x interval */
#define YMIN -1.125
#define YMAX 1.125	/* y interval for 9/16 aspect ratio */

#define SCALING_FACTOR 1.0       /* scaling factor of drawing, needed for flower billiards, otherwise set to 1.0 */

/* Choice of the billiard table, see global_particles.c */

#define B_DOMAIN 30     /* choice of domain shape */

#define CIRCLE_PATTERN 1    /* pattern of circles */
#define POLYLINE_PATTERN 10  /* pattern of polyline */

#define ABSORBING_CIRCLES 1 /* set to 1 for circular scatterers to be absorbing */

#define NMAXCIRCLES 100000     /* total number of circles (must be at least NCX*NCY for square grid) */
#define NMAXPOLY 100000        /* total number of sides of polygonal line */ 
#define NCX 30            /* number of circles in x direction */
#define NCY 20            /* number of circles in y direction */
#define NPOISSON 500        /* number of points for Poisson C_RAND_POISSON arrangement */
#define NGOLDENSPIRAL 2000  /* max number of points for C_GOLDEN_SPIRAL arrandement */
#define SDEPTH 1            /* Sierpinski gastket depth */

#define LAMBDA 1.5	/* parameter controlling shape of domain */
#define MU 0.5          /* second parameter controlling shape of billiard */
#define FOCI 1          /* set to 1 to draw focal points of ellipse */
#define NPOLY 6             /* number of sides of polygon */
#define APOLY 0.0           /* angle by which to turn polygon, in units of Pi/2 */ 
#define PENROSE_RATIO 2.5    /* parameter controlling the shape of small ellipses in Penrose room */

#define DRAW_BILLIARD 1     /* set to 1 to draw billiard */
#define DRAW_CONSTRUCTION_LINES 0   /* set to 1 to draw additional construction lines for billiard */
#define PERIODIC_BC 0       /* set to 1 to enforce periodic boundary conditions when drawing particles */

#define RESAMPLE 0      /* set to 1 if particles should be added when dispersion too large */
#define DEBUG 0         /* draw trajectories, for debugging purposes */

/* Simulation parameters */

#define NPART 10000     /* number of particles */
#define NPARTMAX 100000	/* maximal number of particles after resampling */
#define LMAX 0.01       /* minimal segment length triggering resampling */ 
#define DMIN 0.02       /* minimal distance to boundary for triggering resampling */ 
#define CYCLE 1         /* set to 1 for closed curve (start in all directions) */
#define SHOWTRAILS 1    /* set to 1 to keep trails of the particles */
#define SHOWZOOM 0      /* set to 1 to show zoom on specific area */
#define PRINT_PARTICLE_NUMBER 0 /* set to 1 to print number of particles */
#define PRINT_LEFT_RIGHT_PARTICLE_NUMBER 1 /* set to 1 to print number of particles on left and right side */
#define PRINT_COLLISION_NUMBER 0 /* set to 1 to print number of collisions */
#define TEST_ACTIVE 1   /* set to 1 to test whether particle is in billiard */

#define TEST_INITIAL_COND 0     /* set to 1 to allow only initial conditions that pass a test */

#define NSTEPS 6000     /* number of frames of movie */
#define TIME 1500        /* time between movie frames, for fluidity of real-time simulation */ 
#define DPHI 0.00002     /* integration step */
#define NVID 25          /* number of iterations between images displayed on screen */
#define END_FRAMES 25    /* number of still frames at the end of the movie */

/* Decreasing TIME accelerates the animation and the movie                               */
/* For constant speed of movie, TIME*DPHI should be kept constant                        */
/* However, increasing DPHI too much deterioriates quality of simulation                 */
/* NVID tells how often a picture is drawn in the animation, increase it for faster anim */
/* For a good quality movie, take for instance TIME = 400, DPHI = 0.00005, NVID = 100    */

/* Colors and other graphical parameters */

#define COLOR_PALETTE 10     /* Color palette, see list in global_pdes.c  */

#define NCOLORS 500      /* number of colors */
#define COLORSHIFT 0     /* hue of initial color */ 
#define RAINBOW_COLOR 1  /* set to 1 to use different colors for all particles */
#define FLOWER_COLOR 0   /* set to 1 to adapt initial colors to flower billiard (tracks vs core) */
#define NSEG 100         /* number of segments of boundary */
#define LENGTH 0.04       /* length of velocity vectors */
#define BILLIARD_WIDTH 2    /* width of billiard */
#define PARTICLE_WIDTH 2    /* width of particles */
#define FRONT_WIDTH 3       /* width of wave front */

#define BLACK 1             /* set to 1 for black background */
#define COLOR_OUTSIDE 0     /* set to 1 for colored outside */ 
#define OUTER_COLOR 270.0   /* color outside billiard */
#define PAINT_INT 0         /* set to 1 to paint interior in other color (for polygon/Reuleaux) */
#define PAINT_EXT 1         /* set to 1 to paint exterior */

#define PAUSE 1000       /* number of frames after which to pause */
#define PSLEEP 2         /* sleep time during pause */
#define SLEEP1  1        /* initial sleeping time */
#define SLEEP2  1       /* final sleeping time */

#define NXMAZE 16      /* width of maze */
#define NYMAZE 32      /* height of maze */
#define MAZE_MAX_NGBH 4     /* max number of neighbours of maze cell */
#define RAND_SHIFT 1        /* seed of random number generator */
#define MAZE_XSHIFT 0.0     /* horizontal shift of maze */

```

