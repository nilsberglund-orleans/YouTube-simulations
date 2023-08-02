 ### Parameter values for YouTube simulations ###

Created by **Nils Berglund** and optimized by **Marco Mancini**

C code for videos on YouTube Channel https://www.youtube.com/c/NilsBerglund

Below are parameter values used for different simulations, as well as initial conditions used in 
function animation. Some simulations use variants of the published code. The list is going to be 
updated gradually. 



### 2 August 23 - Faraday waves with 14-fold symmetry ###

**Program:** `rde.c` 

**Initial condition in function `animation()`:** `init_circular_vibration(0.0, 0.0, 0.015, LAMBDA, 0.001, 0.1, 14.0, SWATER_MIN_HEIGHT, phi, xy_in);`

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

#define RDE_EQUATION 8  /* choice of reaction term, see list in global_3d.c */
#define NFIELDS 3       /* number of fields in reaction-diffusion equation */
#define NLAPLACIANS 0   /* number of fields for which to compute Laplacian */

#define ADD_POTENTIAL 0 /* set to 1 to add a potential (for Schrodinger equation) */
#define ADD_MAGNETIC_FIELD 0    /* set to 1 to add a magnetic field (for Schrodinger equation) - then set POTENTIAL 1 */
#define ADD_FORCE_FIELD 1   /* set to 1 to add a foce field (for compressible Euler equation) */
#define POTENTIAL 7         /* type of potential or vector potential, see list in global_3d.c  */
#define FORCE_FIELD 5       /* type of force field, see list in global_3d.c  */
#define ADD_CORIOLIS_FORCE 0    /* set to 1 to add Coriolis force (quasigeostrophic Euler equations) */
#define VARIABLE_DEPTH 0    /* set to 1 for variable depth in shallow water equation */
#define SWATER_DEPTH 4      /* variable depth in shallow water equation */

#define ANTISYMMETRIZE_WAVE_FCT 0   /* set tot 1 to make wave function antisymmetric */
#define ADAPT_STATE_TO_BC 1     /* to smoothly adapt initial state to obstacles */
#define OBSTACLE_GEOMETRY 1     /* geometry of obstacles, as in B_DOMAIN */
#define BC_STIFFNESS 50.0        /* controls region of boundary condition control */
#define CHECK_INTEGRAL 1     /* set to 1 to check integral of first field */

#define JULIA_SCALE 0.5 /* scaling for Julia sets */

/* Choice of the billiard table */

#define B_DOMAIN 999          /* choice of domain shape, see list in global_pdes.c  */

#define CIRCLE_PATTERN 8     /* pattern of circles, see list in global_pdes.c */

#define P_PERCOL 0.25       /* probability of having a circle in C_RAND_PERCOL arrangement */
#define NPOISSON 300        /* number of points for Poisson C_RAND_POISSON arrangement */
#define RANDOM_POLY_ANGLE 0 /* set to 1 to randomize angle of polygons */

#define LAMBDA 1.0	    /* parameter controlling the dimensions of domain */
#define MU 0.06	            /* parameter controlling the dimensions of domain */
#define NPOLY 5             /* number of sides of polygon */
#define APOLY 2.0          /* angle by which to turn polygon, in units of Pi/2 */
#define MDEPTH 7            /* depth of computation of Menger gasket */
#define MRATIO 5            /* ratio defining Menger gasket */
#define MANDELLEVEL 1000      /* iteration level for Mandelbrot set */
#define MANDELLIMIT 10.0     /* limit value for approximation of Mandelbrot set */
#define FOCI 1              /* set to 1 to draw focal points of ellipse */
#define NGRIDX 6            /* number of grid point for grid of disks */
#define NGRIDY 8            /* number of grid point for grid of disks */
#define REVERSE_TESLA_VALVE 1   /* set to 1 to orient Tesla valve in blocking configuration */

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

/* Physical parameters of wave equation */

#define DT 0.00000025

#define VISCOSITY 1.5e-5
#define DISSIPATION 1.0e-8

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
#define G_FIELD 0.75     /* gravity */
#define BC_FIELD 1.0e-5     /* constant in repulsive field from obstacles */
#define AB_RADIUS 0.2   /* radius of region with magnetic field for Aharonov-Bohm effect */
#define K_EULER 50.0    /* constant in stream function integration of Euler equation */
#define K_EULER_INC 0.5    /* constant in incompressible Euler equation */

#define SMOOTHEN_VORTICITY 0    /* set to 1 to smoothen vorticity field in Euler equation */
#define SMOOTHEN_VELOCITY 1     /* set to 1 to smoothen velocity field in Euler equation */
#define SMOOTHEN_PERIOD 10      /* period between smoothenings */
#define SMOOTH_FACTOR 0.15       /* factor by which to smoothen */

#define ADD_TRACERS 1    /* set to 1 to add tracer particles (for Euler equations) */
#define N_TRACERS 1000    /* number of tracer particles */
#define TRACERS_STEP 0.005  /* step size in tracer evolution */

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
#define IN_OUT_FLOW_BC 0          /* type of in-flow/out-flow boundary conditions for Euler equation, 0 for no b.c. */
#define IN_OUT_BC_FACTOR 0.001    /* factor of convex combination between old and new flow */
#define BC_FLOW_TYPE 1            /* type of initial condition */
                                  /* see list in global_pdes.c */
#define IN_OUT_FLOW_MIN_AMP 0.25   /* amplitude of in-flow/out-flow boundary conditions (for Euler equation) - min value */
#define IN_OUT_FLOW_AMP 0.25       /* amplitude of in-flow/out-flow boundary conditions (for Euler equation) - max value */
#define LAMINAR_FLOW_MODULATION 0.01     /* asymmetry of laminar flow */
#define LAMINAR_FLOW_YPERIOD 1.0    /* period of laminar flow in y direction */
#define PRESSURE_GRADIENT 0.2       /* amplitude of pressure gradient for Euler equation */

#define SWATER_MIN_HEIGHT 0.5      /* min height of initial condition for shallow water equation */
#define DEPTH_FACTOR 0.015        /* proportion of min height in variable depth */
#define TANH_FACTOR 1.0         /* steepness of variable depth */

#define EULER_GRADIENT_YSHIFT 0.0    /* y-shift in computation of gradient in Euler equation */

/* Boundary conditions, see list in global_pdes.c  */

#define B_COND 0

#define B_COND_LEFT 0
#define B_COND_RIGHT 0
#define B_COND_TOP 0
#define B_COND_BOTTOM 0


/* Parameters for length and speed of simulation */

#define NSTEPS 1400         /* number of frames of movie */
#define NVID 50            /* number of iterations between images displayed on screen */
#define ACCELERATION_FACTOR 1.0 /* factor by which to increase NVID in course of simulation */
#define DT_ACCELERATION_FACTOR 1.0 /* factor by which to increase time step in course of simulation  */
#define MAX_DT 0.024     /* maximal value of integration step */
#define NSEG 999         /* number of segments of boundary */
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

#define PLOT_3D 1    /* controls whether plot is 2D or 3D */

#define ROTATE_VIEW 0       /* set to 1 to rotate position of observer */
#define ROTATE_ANGLE 360.0  /* total angle of rotation during simulation */

#define DRAW_PERIODICISED 0     /* set to 1 to repeat wave periodically in x and y directions */

/* Plot type - color scheme */

#define CPLOT 70
#define CPLOT_B 74

/* Plot type - height of 3D plot */

#define ZPLOT 70     /* z coordinate in 3D plot */
#define ZPLOT_B 71    /* z coordinate in second 3D plot */

#define AMPLITUDE_HIGH_RES 1    /* set to 1 to increase resolution of P_3D_AMPLITUDE plot */
#define SHADE_3D 1              /* set to 1 to change luminosity according to normal vector */
#define NON_DIRICHLET_BC 0      /* set to 1 to draw only facets in domain, if field is not zero on boundary */
#define WRAP_ANGLE 1            /* experimental: wrap angle to [0, 2Pi) for interpolation in angle schemes */
#define FADE_IN_OBSTACLE 0      /* set to 1 to fade color inside obstacles */
#define FADE_WATER_DEPTH 0      /* set to 1 to make wave color depth-dependent */
#define DRAW_OUTSIDE_GRAY 0     /* experimental - draw outside of billiard in gray */
#define ADD_POTENTIAL_TO_Z 0    /* set to 1 to add the external potential to z-coordinate of plot */
#define ADD_POT_CONSTANT 0.35   /* constant added to wave height */
#define ADD_HEIGHT_CONSTANT -0.5   /* constant added to wave height */
#define DRAW_DEPTH 0            /* set to 1 to draw water depth */
#define DEPTH_SCALE 0.75         /* vertical scaling of depth plot */
#define DEPTH_SHIFT -0.015      /* vertical shift of depth plot */

#define PLOT_SCALE_ENERGY 0.05      /* vertical scaling in energy plot */

#define PRINT_TIME 0        /* set to 1 to print running time */
#define PRINT_VISCOSITY 0   /* set to 1 to print viscosity */
#define PRINT_RPSLZB 0      /* set to 1 to print rpslzb parameter */
#define PRINT_PROBABILITIES 0   /* set to 1 to print probabilities (for Ehrenfest urn configuration) */
#define PRINT_NOISE 0       /* set to 1 to print noise intensity */
#define PRINT_FLOW_SPEED 0      /* set to 1 to print speed of flow */
#define PRINT_AVERAGE_SPEED 0   /* set to 1 to print average speed of flow */
#define PRINT_LEFT 1            /* set to 1 to print parameters at left side */

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

#define COLOR_PALETTE 13       /* Color palette, see list in global_pdes.c  */
#define COLOR_PALETTE_B 0      /* Color palette, see list in global_pdes.c  */

#define BLACK 1          /* black background */

#define COLOR_SCHEME 3   /* choice of color scheme */

#define COLOR_PHASE_SHIFT 0.5   /* phase shift of color scheme, in units of Pi */

#define SCALE 0          /* set to 1 to adjust color scheme to variance of field */
#define SLOPE 1.0        /* sensitivity of color on wave amplitude */
#define VSCALE_AMPLITUDE 10.0      /* additional scaling factor for color scheme P_3D_AMPLITUDE */
#define ATTENUATION 0.0  /* exponential attenuation coefficient of contrast with time */
#define CURL_SCALE 0.000015   /* scaling factor for curl representation */
#define RESCALE_COLOR_IN_CENTER 0   /* set to 1 to decrease color intentiy in the center (for wave escaping ring) */
#define SLOPE_SCHROD_LUM 200.0       /* sensitivity of luminosity on module, for color scheme Z_ARGUMENT */
#define MIN_SCHROD_LUM 0.2       /* minimal luminosity in color scheme Z_ARGUMENT*/
#define VSCALE_PRESSURE 2.0      /* additional scaling factor for color scheme Z_EULER_PRESSURE */
#define PRESSURE_SHIFT 10.0        /* shift for color scheme Z_EULER_PRESSURE */
#define PRESSURE_LOG_SHIFT -2.5     /* shift for color scheme Z_EULER_PRESSURE */
#define VSCALE_WATER_HEIGHT 0.4     /* vertical scaling of water height */

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
#define VSCALE_SPEED 200.0      /* additional scaling factor for color scheme Z_EULER_SPEED */
#define VMEAN_SPEED 0.0       /* mean value around which to scale for color scheme Z_EULER_SPEED */
#define SHIFT_DENSITY 8.5        /* shift for color scheme Z_EULER_DENSITY */
#define VSCALE_DENSITY 3.0        /* additional scaling factor for color scheme Z_EULER_DENSITY */
#define VSCALE_VORTICITY 20.0     /* additional scaling factor for color scheme Z_EULERC_VORTICITY */
#define VORTICITY_SHIFT 0.0     /* vertical shift of vorticity */
#define ZSCALE_SPEED 0.3        /* additional scaling factor for z-coord Z_EULER_SPEED and Z_SWATER_SPEED */
#define VSCALE_SWATER 300.0        /* additional scaling factor for color scheme Z_EULER_DENSITY */

#define NXMAZE 7      /* width of maze */
#define NYMAZE 7      /* height of maze */
#define MAZE_MAX_NGBH 4     /* max number of neighbours of maze cell */
#define RAND_SHIFT 3        /* seed of random number generator */
#define MAZE_XSHIFT 0.0     /* horizontal shift of maze */
#define MAZE_WIDTH 0.04     /* half width of maze walls */

#define DRAW_COLOR_SCHEME 1     /* set to 1 to plot the color scheme */
#define COLORBAR_RANGE 2.5      /* scale of color scheme bar */
#define COLORBAR_RANGE_B 2.5    /* scale of color scheme bar for 2nd part */
#define ROTATE_COLOR_SCHEME 0   /* set to 1 to draw color scheme horizontally */
#define CIRC_COLORBAR 0         /* set to 1 to draw circular color scheme */
#define CIRC_COLORBAR_B 1       /* set to 1 to draw circular color scheme */

/* only for compatibility with wave_common.c */
#define TWOSPEEDS 0          /* set to 1 to replace hardcore boundary by medium with different speed */
#define VARIABLE_IOR 0      /* set to 1 for a variable index of refraction */
#define IOR 4               /* choice of index of refraction, see list in global_pdes.c */
#define IOR_TOTAL_TURNS 1.5 /* total angle of rotation for IOR_PERIODIC_WELLS_ROTATING */
#define MANDEL_IOR_SCALE -0.05   /* parameter controlling dependence of IoR on Mandelbrot escape speed */
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
#define FLUX_WINDOW 20      /* averaging window for energy flux */
#define ADD_WAVE_PACKET_SOURCES 1       /* set to 1 to add several sources emitting wave packets */
#define WAVE_PACKET_SOURCE_TYPE 1       /* type of wave packet sources */
#define N_WAVE_PACKETS 15               /* number of wave packets */
#define WAVE_PACKET_RADIUS 20            /* radius of wave packets */
#define OSCIL_LEFT_YSHIFT 25.0   /* y-dependence of left oscillation (for non-horizontal waves) */
/* end of constants added only for compatibility with wave_common.c */


double u_3d[2] = {0.75, -0.45};     /* projections of basis vectors for REP_AXO_3D representation */
double v_3d[2] = {-0.75, -0.45};
double w_3d[2] = {0.0, 0.015};
double light[3] = {0.816496581, 0.40824829, 0.40824829};      /* vector of "light" direction for P_3D_ANGLE color scheme */
double observer[3] = {8.0, 8.0, 7.0};    /* location of observer for REP_PROJ_3D representation */ 
int reset_view = 0;         /* switch to reset 3D view parameters (for option ROTATE_VIEW) */

#define Z_SCALING_FACTOR 75.0   /* overall scaling factor of z axis for REP_PROJ_3D representation */
#define XY_SCALING_FACTOR 2.4  /* overall scaling factor for on-screen (x,y) coordinates after projection */
#define ZMAX_FACTOR 1.0        /* max value of z coordinate for REP_PROJ_3D representation */
#define XSHIFT_3D 0.0          /* overall x shift for REP_PROJ_3D representation */
#define YSHIFT_3D 0.1         /* overall y shift for REP_PROJ_3D representation */
#define BORDER_PADDING 0       /* distance from boundary at which to plot points, to avoid boundary effects due to gradient */

/* For debugging purposes only */
#define FLOOR 1         /* set to 1 to limit wave amplitude to VMAX */
#define VMAX 100.0        /* max value of wave amplitude */
#define TEST_GRADIENT 0 /* print norm squared of gradient */

```

### 1 August 23 - Bragg diffraction ###

**Program:** `wave_billiard.c` 

**Initial condition in function `animation()`:** `init_wave_flat(phi, psi, xy_in);`

```
#define MOVIE 1         /* set to 1 to generate movie */
#define DOUBLE_MOVIE 1  /* set to 1 to produce movies for wave height and energy simultaneously */
#define SAVE_MEMORY 1   /* set to 1 to save memory when writing tiff images */
#define NO_EXTRA_BUFFER_SWAP 1    /* some OS require one less buffer swap when recording images */

#define VARIABLE_IOR 0      /* set to 1 for a variable index of refraction */
#define IOR 9               /* choice of index of refraction, see list in global_pdes.c */
#define IOR_TOTAL_TURNS 1.5 /* total angle of rotation for IOR_PERIODIC_WELLS_ROTATING */
#define MANDEL_IOR_SCALE -0.05   /* parameter controlling dependence of IoR on Mandelbrot escape speed */


/* General geometrical parameters */

#define WINWIDTH 	1920  /* window width */
#define WINHEIGHT 	1150  /* window height */
#define NX 3840          /* number of grid points on x axis */
#define NY 2300          /* number of grid points on y axis */

#define XMIN -2.0
#define XMAX 2.0	/* x interval */
#define YMIN -0.797916667
#define YMAX 1.597916667	/* y interval for 9/16 aspect ratio */

#define HIGHRES 1       /* set to 1 if resolution of grid is double that of displayed image */

#define JULIA_SCALE 1.0 /* scaling for Julia sets */

/* Choice of the billiard table */

#define B_DOMAIN 20        /* choice of domain shape, see list in global_pdes.c */

#define CIRCLE_PATTERN 101   /* pattern of circles or polygons, see list in global_pdes.c */

#define COMPARISON 0        /* set to 1 to compare two different patterns (beta) */
#define B_DOMAIN_B 20       /* second domain shape, for comparisons */
#define CIRCLE_PATTERN_B 0  /* second pattern of circles or polygons */

#define P_PERCOL 0.25       /* probability of having a circle in C_RAND_PERCOL arrangement */
#define NPOISSON 1000        /* number of points for Poisson C_RAND_POISSON arrangement */
#define RANDOM_POLY_ANGLE 1 /* set to 1 to randomize angle of polygons */

#define LAMBDA 0.1197916667	    /* parameter controlling the dimensions of domain */
#define MU 0.01             /* parameter controlling the dimensions of domain */
#define NPOLY 6             /* number of sides of polygon */
#define APOLY 0.0           /* angle by which to turn polygon, in units of Pi/2 */ 
#define MDEPTH 6            /* depth of computation of Menger gasket */
#define MRATIO 3            /* ratio defining Menger gasket */
#define MANDELLEVEL 1000    /* iteration level for Mandelbrot set */
#define MANDELLIMIT 10.0    /* limit value for approximation of Mandelbrot set */
#define FOCI 1              /* set to 1 to draw focal points of ellipse */
#define NGRIDX 40            /* number of grid point for grid of disks */
#define NGRIDY 8            /* number of grid point for grid of disks */

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
#define OSCILLATE_LEFT 1     /* set to 1 to add oscilating boundary condition on the left */
#define OSCILLATE_TOPBOT 0   /* set to 1 to enforce a planar wave on top and bottom boundary */
#define OSCILLATION_SCHEDULE 0  /* oscillation schedule, see list in global_pdes.c */

#define OMEGA 0.012        /* frequency of periodic excitation */
#define AMPLITUDE 1.0      /* amplitude of periodic excitation */ 
#define ACHIRP 0.25        /* acceleration coefficient in chirp */
#define DAMPING 0.0        /* damping of periodic excitation */
#define COURANT 0.1        /* Courant number */
#define COURANTB 0.01      /* Courant number in medium B */
#define GAMMA 0.0          /* damping factor in wave equation */
#define GAMMAB 0.0         /* damping factor in wave equation */
#define GAMMA_SIDES 1.0e-4      /* damping factor on boundary */
#define GAMMA_TOPBOT 1.0e-7     /* damping factor on boundary */
#define KAPPA 0.0           /* "elasticity" term enforcing oscillations */
#define KAPPA_SIDES 5.0e-4  /* "elasticity" term on absorbing boundary */
#define KAPPA_TOPBOT 0.0    /* "elasticity" term on absorbing boundary */
#define OSCIL_LEFT_YSHIFT -200.0   /* y-dependence of left oscillation (for non-horizontal waves) */
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

/* Parameters for length and speed of simulation */

#define NSTEPS 2800       /* number of frames of movie */
#define NVID 7            /* number of iterations between images displayed on screen */
#define NSEG 1000         /* number of segments of boundary */
#define INITIAL_TIME 0      /* time after which to start saving frames */
#define BOUNDARY_WIDTH 2    /* width of billiard boundary */
#define PRINT_SPEED 0       /* print speed of moving source */
#define PRINT_FREQUENCY 0       /* print frequency (for phased array) */

#define PAUSE 200       /* number of frames after which to pause */
#define PSLEEP 1         /* sleep time during pause */
#define SLEEP1  1        /* initial sleeping time */
#define SLEEP2  1        /* final sleeping time */
#define MID_FRAMES 20    /* number of still frames between parts of two-part movie */
#define END_FRAMES 100    /* number of still frames at end of movie */
#define FADE 1           /* set to 1 to fade at end of movie */

/* Parameters of initial condition */

#define INITIAL_AMP 2.0            /* amplitude of initial condition */
#define INITIAL_VARIANCE 0.000025    /* variance of initial condition */
#define INITIAL_WAVELENGTH  0.05  /* wavelength of initial condition */

/* Plot type, see list in global_pdes.c  */

#define PLOT 0

#define PLOT_B 5        /* plot type for second movie */

/* Color schemes */

#define COLOR_PALETTE 17      /* Color palette, see list in global_pdes.c  */
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
#define LOG_SHIFT 1.5     /* shift of colors on log scale */
#define FLUX_SCALE 5.0e3    /* scaling factor for enegy flux represtnation */
#define RESCALE_COLOR_IN_CENTER 0   /* set to 1 to decrease color intentiy in the center (for wave escaping ring) */

#define COLORHUE 260     /* initial hue of water color for scheme C_LUM */
#define COLORDRIFT 0.0   /* how much the color hue drifts during the whole simulation */
#define LUMMEAN 0.5      /* amplitude of luminosity variation for scheme C_LUM */
#define LUMAMP 0.3       /* amplitude of luminosity variation for scheme C_LUM */
#define HUEMEAN 180.0    /* mean value of hue for color scheme C_HUE */
#define HUEAMP -180.0      /* amplitude of variation of hue for color scheme C_HUE */

#define DRAW_COLOR_SCHEME 1    /* set to 1 to plot the color scheme */
#define COLORBAR_RANGE 2.0     /* scale of color scheme bar */
#define COLORBAR_RANGE_B 0.8   /* scale of color scheme bar for 2nd part */
#define ROTATE_COLOR_SCHEME 0   /* set to 1 to draw color scheme horizontally */
#define CIRC_COLORBAR 0         /* set to 1 to draw circular color scheme */
#define CIRC_COLORBAR_B 0       /* set to 1 to draw circular color scheme */

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

