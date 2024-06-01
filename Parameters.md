 ### Parameter values for YouTube simulations ###

Created by **Nils Berglund** and optimized by **Marco Mancini**

C code for videos on YouTube Channel https://www.youtube.com/c/NilsBerglund

Below are parameter values used for different simulations, as well as initial conditions used in function animation. Some simulations use variants of the published code. The list is going to be updated gradually. 


### 01 June 2024 - More stable weather with 16 pressure systems, with westerlies and trade winds ###

**Program:** `rde.c` 

**Initial condition in function `animation()`:** 

```
    init_laminar_flow_earth(0.04, phi, xy_in);

    /* high pressure systems, northern hemisphere */
    init_vortex_state_sphere_mod(1, 0.5, 2.5*PID, PID + 0.6, 0.03, 0.04, phi, xy_in, wsphere);
    init_vortex_state_sphere_mod(1, 0.6, 3.7*PID, PID + 0.7, 0.035, 0.04, phi, xy_in, wsphere);
    init_vortex_state_sphere_mod(1, 0.2, 2.5*PID, PID + 1.3, 0.01, 0.01, phi, xy_in, wsphere);
    init_vortex_state_sphere_mod(1, 0.5, 1.0*PID, PID + 1.0, 0.03, 0.02, phi, xy_in, wsphere);
    /* low pressure systems, northern hemisphere */
    init_vortex_state_sphere_mod(1, -0.5, 3.0*PID, PID + 0.6, 0.03, -0.04, phi, xy_in, wsphere);
    init_vortex_state_sphere_mod(1, -0.2, 2.3*PID, PID + 1.1, 0.01, -0.02, phi, xy_in, wsphere);
    init_vortex_state_sphere_mod(1, -0.5, 1.2*PID, PID + 0.4, 0.01, -0.04, phi, xy_in, wsphere);
    init_vortex_state_sphere_mod(1, -0.8, 0.3*PID, PID + 0.6, 0.01, -0.04, phi, xy_in, wsphere);
    /* high pressure systems, southern hemisphere */
    init_vortex_state_sphere_mod(1, -0.6, 2.4*PID, PID - 0.7, 0.03, 0.04, phi, xy_in, wsphere);
    init_vortex_state_sphere_mod(1, -0.6, 1.6*PID, PID - 0.7, 0.02, 0.04, phi, xy_in, wsphere);
    init_vortex_state_sphere_mod(1, -0.4, 0.5*PID, PID - 0.7, 0.02, 0.04, phi, xy_in, wsphere);
   /* low pressure systems, southern hemisphere */
    init_vortex_state_sphere_mod(1, 0.5, 3.4*PID, PID - 0.7, 0.03, -0.04, phi, xy_in, wsphere);
    init_vortex_state_sphere_mod(1, 0.2, 2.0*PID, PID - 0.1, 0.04, -0.04, phi, xy_in, wsphere);
    init_vortex_state_sphere_mod(1, 0.3, 1.0*PID, PID - 0.6, 0.03, -0.04, phi, xy_in, wsphere);
    init_vortex_state_sphere_mod(1, 0.4, 0.1*PID, PID - 0.5, 0.03, -0.04, phi, xy_in, wsphere);
    init_vortex_state_sphere_mod(1, 0.2, 0.1*PID, 0.1, 0.1, -0.01, phi, xy_in, wsphere);

```

**3D parts:** 

```
#define MOVIE 1         /* set to 1 to generate movie */
#define DOUBLE_MOVIE 1  /* set to 1 to produce movies for wave height and energy simultaneously */
#define SAVE_MEMORY 1           /* set to 1 to save memory when writing tiff images */
#define NO_EXTRA_BUFFER_SWAP 1    /* some OS require one less buffer swap when recording images */

/* General geometrical parameters */

#define WINWIDTH 	1920  /* window width */
#define WINHEIGHT 	1150  /* window height */
#define NX 780          /* number of grid points on x axis */
#define NY 400          /* number of grid points on y axis */
#define HRES 2          /* factor for high resolution plots */

#define XMIN -2.0
#define XMAX 2.0	/* x interval */
#define YMIN -1.041666667
#define YMAX 1.041666667	/* y interval for 9/16 aspect ratio */

/* Choice of simulated equation */

#define RDE_EQUATION 7  /* choice of reaction term, see list in global_3d.c */
#define NFIELDS 3       /* number of fields in reaction-diffusion equation */
#define NLAPLACIANS 0   /* number of fields for which to compute Laplacian */

#define SPHERE 1        /* set to 1 to simulate equation on sphere */
#define DPOLE 1         /* safety distance to poles */
#define DSMOOTH 1       /* size of neighbourhood of poles that are smoothed */
#define SMOOTHPOLE 0.01  /* smoothing coefficient at poles */
#define SMOOTHCOTPOLE 0.01  /* smoothing coefficient of cotangent at poles */
#define PHISHIFT 0.0    /* shift of phi in 2D plot (in degrees) */
#define SMOOTHBLOCKS 1  /* set to 1 to use blocks of points near the poles */
#define BLOCKDIST 64    /* distance to poles where points are blocked */ 
#define ZERO_MERIDIAN 190.0     /* choice of zero meridian (will be at left/right boundary of 2d plot) */

#define ADD_POTENTIAL 0 /* set to 1 to add a potential (for Schrodinger equation) */
#define ADD_MAGNETIC_FIELD 0    /* set to 1 to add a magnetic field (for Schrodinger equation) - then set POTENTIAL 1 */
#define ADD_FORCE_FIELD 1   /* set to 1 to add a foce field (for compressible Euler equation) */
#define POTENTIAL 7         /* type of potential or vector potential, see list in global_3d.c  */
#define FORCE_FIELD 6       /* type of force field, see list in global_3d.c  */
#define ADD_CORIOLIS_FORCE 1    /* set to 1 to add Coriolis force (quasigeostrophic Euler equations) */
#define VARIABLE_DEPTH 0    /* set to 1 for variable depth in shallow water equation */
#define SWATER_DEPTH 4      /* variable depth in shallow water equation */

#define ANTISYMMETRIZE_WAVE_FCT 0   /* set tot 1 to make wave function antisymmetric */
#define ADAPT_STATE_TO_BC 1     /* to smoothly adapt initial state to obstacles */
#define OBSTACLE_GEOMETRY 84     /* geometry of obstacles, as in B_DOMAIN */
#define BC_STIFFNESS 50.0        /* controls region of boundary condition control */
#define CHECK_INTEGRAL 1     /* set to 1 to check integral of first field */

#define JULIA_SCALE 0.5 /* scaling for Julia sets */
#define JULIA_ROT -20.0       /* rotation of Julia set, in degrees */
#define JULIA_RE 0.5    
#define JULIA_IM 0.462    /* parameters for Julia sets */

/* Choice of the billiard table */

#define B_DOMAIN 999          /* choice of domain shape, see list in global_pdes.c  */

#define CIRCLE_PATTERN 8     /* pattern of circles, see list in global_pdes.c */

#define P_PERCOL 0.25       /* probability of having a circle in C_RAND_PERCOL arrangement */
#define NPOISSON 300        /* number of points for Poisson C_RAND_POISSON arrangement */
#define RANDOM_POLY_ANGLE 0 /* set to 1 to randomize angle of polygons */

#define LAMBDA 1.0	    /* parameter controlling the dimensions of domain */
#define MU 1.0	            /* parameter controlling the dimensions of domain */
#define NPOLY 5             /* number of sides of polygon */
#define APOLY 2.0           /* angle by which to turn polygon, in units of Pi/2 */
#define MDEPTH 7            /* depth of computation of Menger gasket */
#define MRATIO 5            /* ratio defining Menger gasket */
#define MANDELLEVEL 1000      /* iteration level for Mandelbrot set */
#define MANDELLIMIT 10.0     /* limit value for approximation of Mandelbrot set */
#define FOCI 1              /* set to 1 to draw focal points of ellipse */
#define NGRIDX 6            /* number of grid point for grid of disks */
#define NGRIDY 8            /* number of grid point for grid of disks */
#define REVERSE_TESLA_VALVE 1   /* set to 1 to orient Tesla valve in blocking configuration */
#define WALL_WIDTH 0.05      /* width of wall separating lenses */

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

#define DT 0.00000015

#define VISCOSITY 0.02
#define POISSON_STIFFNESS 1.0   /* stiffness of Poisson equation solver for incompressible Euler */
#define DISSIPATION 0.0

#define RPSA 0.75         /* parameter in Rock-Paper-Scissors-type interaction */
#define RPSLZB 0.0       /* second parameter in Rock-Paper-Scissors-Lizard-Spock type interaction */

#define K_AC 0.1        /* force constant in Allen-Cahn equation */
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
#define G_FIELD 0.03    /* gravity/Coriolis force */
#define BC_FIELD 1.0e-5     /* constant in repulsive field from obstacles */
#define AB_RADIUS 0.2   /* radius of region with magnetic field for Aharonov-Bohm effect */
#define K_EULER 50.0    /* constant in stream function integration of Euler equation */
#define K_EULER_INC 0.5    /* constant in incompressible Euler equation */
#define C_EULER_COMP 0.1   /* constant in compressible Euler equation */

#define SMOOTHEN_VORTICITY 0    /* set to 1 to smoothen vorticity field in Euler equation */
#define SMOOTHEN_VELOCITY 1     /* set to 1 to smoothen velocity field in Euler equation */
#define SMOOTHEN_PERIOD 10      /* period between smoothenings */
#define SMOOTH_FACTOR 0.05      /* factor by which to smoothen */

#define ADD_OSCILLATING_SOURCE 0        /* set to 1 to add an oscillating wave source */
#define OSCILLATING_SOURCE_PERIOD 1     /* period of oscillating source */
#define OSCILLATING_SOURCE_OMEGA 0.2    /* frequency of oscillating source */

#define ADD_TRACERS 1    /* set to 1 to add tracer particles (for Euler equations) */
#define N_TRACERS 2000    /* number of tracer particles */
#define TRACERS_STEP 0.005  /* step size in tracer evolution */

#define T_OUT 2.0       /* outside temperature */
#define T_IN 0.0        /* inside temperature */
#define SPEED 0.0       /* speed of drift to the right */

#define ADD_NOISE 0     /* set to 1 to add noise, set to 2 to add noise in right half */
#define NOISE_INTENSITY 0.01      /* noise intensity */
#define CHANGE_NOISE 0      /* set to 1 to increase noise intensity */
#define NOISE_FACTOR 40.0   /* factor by which to increase noise intensity */
#define NOISE_INITIAL_TIME 100  /* initial time during which noise remains constant */

#define CHANGE_VISCOSITY 0      /* set to 1 to change the viscosity in the course of the simulation */
#define ADJUST_INTSTEP 0       /* set to 1 to decrease integration step when viscosity increases */
#define VISCOSITY_INITIAL_TIME 10  /* initial time during which viscosity remains constant */
#define VISCOSITY_FACTOR 100.0   /* factor by which to change viscosity */
#define VISCOSITY_MAX 2.0        /* max value of viscosity beyond which NVID is increased and integration step is decrase, for numerical stability */
                                        
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

#define B_COND 1

#define B_COND_LEFT 0
#define B_COND_RIGHT 0
#define B_COND_TOP 0
#define B_COND_BOTTOM 0

/* Parameters for length and speed of simulation */

#define NSTEPS 1900           /* number of frames of movie */
#define NVID 100          /* number of iterations between images displayed on screen */
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
#define MID_FRAMES 100    /* number of still frames between parts of two-part movie */
#define END_FRAMES 250    /* number of still frames at end of movie */
#define FADE 1           /* set to 1 to fade at end of movie */

/* Visualisation */

#define PLOT_3D 1    /* controls whether plot is 2D or 3D */
#define PLOT_SPHERE 1   /* draws fields on a sphere */

#define ROTATE_VIEW 1       /* set to 1 to rotate position of observer */
#define ROTATE_ANGLE 360.0  /* total angle of rotation during simulation */
#define SHADE_3D 1              /* set to 1 to change luminosity according to normal vector */
#define SHADE_2D 0              /* set to 1 to change luminosity according to normal vector */
#define VIEWPOINT_TRAJ 0    /* type of viewpoint trajectory */
#define MAX_LATITUDE 45.0   /* maximal latitude for viewpoint trajectory VP_ORBIT2 */

#define DRAW_PERIODICISED 0     /* set to 1 to repeat wave periodically in x and y directions */

/* Plot type - color scheme */

#define CPLOT 60
#define CPLOT_B 61

/* Plot type - height of 3D plot */

#define ZPLOT 60     /* z coordinate in 3D plot */
#define ZPLOT_B 61    /* z coordinate in second 3D plot */

#define AMPLITUDE_HIGH_RES 1    /* set to 1 to increase resolution of P_3D_AMPLITUDE plot */
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

#define COLOR_PALETTE 11       /* Color palette, see list in global_pdes.c  */
#define COLOR_PALETTE_B 16     /* Color palette, see list in global_pdes.c  */

#define BLACK 1          /* black background */
#define COLOR_OUT_R 1.0    /* color outside domain */
#define COLOR_OUT_G 1.0    
#define COLOR_OUT_B 1.0    

#define COLOR_SCHEME 3   /* choice of color scheme */

#define COLOR_PHASE_SHIFT 0.25   /* phase shift of color scheme, in units of Pi */

#define SCALE 0          /* set to 1 to adjust color scheme to variance of field */
#define SLOPE 1.0        /* sensitivity of color on wave amplitude */
#define VSCALE_AMPLITUDE 100.0      /* additional scaling factor for color scheme P_3D_AMPLITUDE */
#define ATTENUATION 0.0  /* exponential attenuation coefficient of contrast with time */
#define CURL_SCALE 1.0   /* scaling factor for curl representation */
#define RESCALE_COLOR_IN_CENTER 0   /* set to 1 to decrease color intentiy in the center (for wave escaping ring) */
#define SLOPE_SCHROD_LUM 100.0       /* sensitivity of luminosity on module, for color scheme Z_ARGUMENT */
#define MIN_SCHROD_LUM 0.1       /* minimal luminosity in color scheme Z_ARGUMENT*/
#define VSCALE_PRESSURE 2.0      /* additional scaling factor for color scheme Z_EULER_PRESSURE */
#define PRESSURE_SHIFT 10.0        /* shift for color scheme Z_EULER_PRESSURE */
#define PRESSURE_LOG_SHIFT -2.5     /* shift for color scheme Z_EULER_PRESSURE */
#define VSCALE_WATER_HEIGHT 0.4     /* vertical scaling of water height */
#define SHADE_SCALE_2D 0.25     /* controls "depth" of 2D shading */

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
#define VSCALE_SPEED 50.0      /* additional scaling factor for color scheme Z_EULER_SPEED */
#define VMEAN_SPEED 0.0       /* mean value around which to scale for color scheme Z_EULER_SPEED */
#define SHIFT_DENSITY 1.0        /* shift for color scheme Z_EULER_DENSITY */
#define VSCALE_DENSITY 30.0      /* additional scaling factor for color scheme Z_EULER_DENSITY */
#define VSCALE_VORTICITY 15.0     /* additional scaling factor for color scheme Z_EULERC_VORTICITY */
#define VORTICITY_SHIFT 0.0     /* vertical shift of vorticity */
#define ZSCALE_SPEED 0.5        /* additional scaling factor for z-coord Z_EULER_SPEED and Z_SWATER_SPEED */
#define ZSHIFT_SPEED 0.0        /* additional shift of z-coord Z_EULER_SPEED and Z_SWATER_SPEED */
#define ZSCALE_NORMGRADIENT -0.0001  /* vertical scaling for Z_NORM_GRADIENT */
#define VSCALE_SWATER 250.0        /* additional scaling factor for color scheme Z_EULER_DENSITY */

#define NXMAZE 7      /* width of maze */
#define NYMAZE 7      /* height of maze */
#define MAZE_MAX_NGBH 4     /* max number of neighbours of maze cell */
#define RAND_SHIFT 3        /* seed of random number generator */
#define MAZE_XSHIFT 0.0     /* horizontal shift of maze */
#define MAZE_WIDTH 0.04     /* half width of maze walls */

#define DRAW_COLOR_SCHEME 1     /* set to 1 to plot the color scheme */
#define COLORBAR_RANGE 3.0      /* scale of color scheme bar */
#define COLORBAR_RANGE_B 2.5    /* scale of color scheme bar for 2nd part */
#define ROTATE_COLOR_SCHEME 0   /* set to 1 to draw color scheme horizontally */
#define CIRC_COLORBAR 0         /* set to 1 to draw circular color scheme */
#define CIRC_COLORBAR_B 0       /* set to 1 to draw circular color scheme */

/* only for compatibility with wave_common.c */
#define TWOSPEEDS 0          /* set to 1 to replace hardcore boundary by medium with different speed */
#define VARIABLE_IOR 0      /* set to 1 for a variable index of refraction */
#define IOR 4               /* choice of index of refraction, see list in global_pdes.c */
#define IOR_TOTAL_TURNS 1.5 /* total angle of rotation for IOR_PERIODIC_WELLS_ROTATING */
#define MANDEL_IOR_SCALE -0.05   /* parameter controlling dependence of IoR on Mandelbrot escape speed */
#define OMEGA 0.005        /* frequency of periodic excitation */
#define OSCIL_YMAX 0.2      /* defines oscillation range */
#define COURANT 0.08       /* Courant number */
#define COURANTB 0.03      /* Courant number in medium B */
#define INITIAL_AMP 0.5         /* amplitude of initial condition */
#define INITIAL_VARIANCE 0.0002  /* variance of initial condition */
#define INITIAL_WAVELENGTH  0.1  /* wavelength of initial condition */
#define VSCALE_ENERGY 200.0       /* additional scaling factor for color scheme P_3D_ENERGY */
#define PHASE_FACTOR 20.0       /* factor in computation of phase in color scheme P_3D_PHASE */
#define PHASE_SHIFT 0.0      /* shift of phase in color scheme P_3D_PHASE */
#define AVRG_E_FACTOR 0.99   /* controls time window size in P_AVERAGE_ENERGY scheme */
#define OSCILLATION_SCHEDULE 0  /* oscillation schedule, see list in global_pdes.c */
#define AMPLITUDE 0.8      /* amplitude of periodic excitation */ 
#define ACHIRP 0.2        /* acceleration coefficient in chirp */
#define DAMPING 0.0        /* damping of periodic excitation */
#define COMPARISON 0        /* set to 1 to compare two different patterns (beta) */
#define B_DOMAIN_B 20       /* second domain shape, for comparisons */
#define CIRCLE_PATTERN_B 0  /* second pattern of circles or polygons */
#define FLUX_WINDOW 20      /* averaging window for energy flux */
#define ADD_WAVE_PACKET_SOURCES 0       /* set to 1 to add several sources emitting wave packets */
#define WAVE_PACKET_SOURCE_TYPE 1       /* type of wave packet sources */
#define N_SOURCES 1                     /* number of wave sources */
#define N_WAVE_PACKETS 15               /* number of wave packets */
#define WAVE_PACKET_RADIUS 20            /* radius of wave packets */
#define OSCIL_LEFT_YSHIFT 25.0   /* y-dependence of left oscillation (for non-horizontal waves) */
#define INITIAL_SHIFT 20.0          /* time shift of initial wave packet (in oscillation periods) */
#define WAVE_PACKET_SHIFT 200.0     /* time shift between wave packets (in oscillation periods) */
#define DRAW_WAVE_PROFILE 0     /* set to 1 to draw a profile of the wave */
#define HORIZONTAL_WAVE_PROFILE 0 /* set to 1 to draw wave profile vertically */
#define VERTICAL_WAVE_PROFILE 0 /* set to 1 to draw wave profile vertically */
#define WAVE_PROFILE_X 1.9      /* value of x to sample wave profile */
#define WAVE_PROFILE_Y -1.0      /* value of y to sample wave profile */
#define PROFILE_AT_BOTTOM 1     /* draw wave profile at bottom instead of top */
#define AVERAGE_WAVE_PROFILE 1  /* set to 1 to draw time-average of wave profile squared*/
#define MU_B 1.0           /* parameter controlling the dimensions of domain */
#define GAMMA 0.0          /* damping factor in wave equation */
#define GAMMAB 0.0         /* damping factor in wave equation */
#define VERTICAL_WAVE_PROFILE 0 /* set to 1 to draw wave profile vertically */
#define DRAW_WAVE_TIMESERIES 0  /* set to 1 to draw a time series of the wave */
#define TIMESERIES_NVALUES 400  /* number of values plotted in time series */
#define DRAW_WAVE_SOURCE 0      /* set to 1 to draw source of wave at (wave_source_x, wave_source_y) */
#define MESSAGE_LDASH 1         /* length of dash for Morse code message */
#define MESSAGE_LDOT 1          /* length of dot for Morse code message */
#define MESSAGE_LINTERVAL 1     /* length of interval between dashes/dots for Morse code message */
#define MESSAGE_LINTERLETTER 1  /* length of interval between letters for Morse code message */
#define MESSAGE_LSPACE 1        /* length of space for Morse code message */
#define MESSAGE_INITIAL_TIME 1  /* initial time before starting message for Morse code message */    
/* end of constants added only for compatibility with wave_common.c */

double u_3d[2] = {0.75, -0.45};     /* projections of basis vectors for REP_AXO_3D representation */
double v_3d[2] = {-0.75, -0.45};
double w_3d[2] = {0.0, 0.015};
double light[3] = {-0.816496581, -0.40824829, 0.40824829};      /* vector of "light" direction for P_3D_ANGLE color scheme */
double observer[3] = {-8.0, -4.0, 4.0};    /* location of observer for REP_PROJ_3D representation */ 
int reset_view = 0;         /* switch to reset 3D view parameters (for option ROTATE_VIEW) */

/* constants for simulations on planets */
#define ADD_DEM 1               /* add DEM (digital elevation model) */
#define ADD_NEGATIVE_DEM 0      /* add DEM with bathymetric data */
#define RSCALE_DEM 0.1          /* scaling factor of radial component for DEM */
#define SMOOTH_DEM 0            /* set to 1 to smoothen DEM (to make altitude less constant) */
#define DEM_SMOOTH_STEPS 1      /* number of smoothening steps */
#define DEM_SMOOTH_HEIGHT 2.0   /* relative height below which to smoothen */
#define DEM_MAXHEIGHT 9000.0     /* max height of DEM (estimated from Everest/Olympus Mons) */
#define DEM_MAXDEPTH -10000     /* max depth of DEM */
#define PLANET_SEALEVEL 0.0      /* sea level for flooded planet */
#define VENUS_NODATA_FACTOR 0.5     /* altitude to assign to DEM points without data (fraction of mean altitude) */

#define Z_SCALING_FACTOR 0.8   /* overall scaling factor of z axis for REP_PROJ_3D representation */
#define XY_SCALING_FACTOR 2.0  /* overall scaling factor for on-screen (x,y) coordinates after projection */
#define ZMAX_FACTOR 1.0        /* max value of z coordinate for REP_PROJ_3D representation */
#define XSHIFT_3D 0.0          /* overall x shift for REP_PROJ_3D representation */
#define YSHIFT_3D 0.0          /* overall y shift for REP_PROJ_3D representation */
#define BORDER_PADDING 0       /* distance from boundary at which to plot points, to avoid boundary effects due to gradient */
#define DRAW_ARROW 0           /* set to 1 to draw arrow above sphere */

#define RSCALE 0.01            /* scaling factor of radial component */
#define RSHIFT -0.01           /* shift in radial component */
#define RMAX 2.0               /* max value of radial component */
#define RMIN 0.5               /* min value of radial component */
#define COS_VISIBLE -1.1       /* limit on cosine of normal to shown facets */

/* For debugging purposes only */
#define FLOOR 1         /* set to 1 to limit wave amplitude to VMAX */
#define VMAX 100.0        /* max value of wave amplitude */
#define TEST_GRADIENT 0 /* print norm squared of gradient */

```

**2D parts** 

```
#define MOVIE 1         /* set to 1 to generate movie */
#define DOUBLE_MOVIE 1  /* set to 1 to produce movies for wave height and energy simultaneously */
#define SAVE_MEMORY 1           /* set to 1 to save memory when writing tiff images */
#define NO_EXTRA_BUFFER_SWAP 1    /* some OS require one less buffer swap when recording images */

/* General geometrical parameters */

#define WINWIDTH 	1920  /* window width */
#define WINHEIGHT 	1150  /* window height */
#define NX 780          /* number of grid points on x axis */
#define NY 400          /* number of grid points on y axis */
#define HRES 2          /* factor for high resolution plots */

#define XMIN -2.0
#define XMAX 2.0	/* x interval */
#define YMIN -1.041666667
#define YMAX 1.041666667	/* y interval for 9/16 aspect ratio */

/* Choice of simulated equation */

#define RDE_EQUATION 7  /* choice of reaction term, see list in global_3d.c */
#define NFIELDS 3       /* number of fields in reaction-diffusion equation */
#define NLAPLACIANS 0   /* number of fields for which to compute Laplacian */

#define SPHERE 1        /* set to 1 to simulate equation on sphere */
#define DPOLE 1         /* safety distance to poles */
#define DSMOOTH 1       /* size of neighbourhood of poles that are smoothed */
#define SMOOTHPOLE 0.01  /* smoothing coefficient at poles */
#define SMOOTHCOTPOLE 0.01  /* smoothing coefficient of cotangent at poles */
#define PHISHIFT 0.0    /* shift of phi in 2D plot (in degrees) */
#define SMOOTHBLOCKS 1  /* set to 1 to use blocks of points near the poles */
#define BLOCKDIST 64    /* distance to poles where points are blocked */ 
#define ZERO_MERIDIAN 190.0     /* choice of zero meridian (will be at left/right boundary of 2d plot) */

#define ADD_POTENTIAL 0 /* set to 1 to add a potential (for Schrodinger equation) */
#define ADD_MAGNETIC_FIELD 0    /* set to 1 to add a magnetic field (for Schrodinger equation) - then set POTENTIAL 1 */
#define ADD_FORCE_FIELD 1   /* set to 1 to add a foce field (for compressible Euler equation) */
#define POTENTIAL 7         /* type of potential or vector potential, see list in global_3d.c  */
#define FORCE_FIELD 6       /* type of force field, see list in global_3d.c  */
#define ADD_CORIOLIS_FORCE 1    /* set to 1 to add Coriolis force (quasigeostrophic Euler equations) */
#define VARIABLE_DEPTH 0    /* set to 1 for variable depth in shallow water equation */
#define SWATER_DEPTH 4      /* variable depth in shallow water equation */

#define ANTISYMMETRIZE_WAVE_FCT 0   /* set tot 1 to make wave function antisymmetric */
#define ADAPT_STATE_TO_BC 1     /* to smoothly adapt initial state to obstacles */
#define OBSTACLE_GEOMETRY 84     /* geometry of obstacles, as in B_DOMAIN */
#define BC_STIFFNESS 50.0        /* controls region of boundary condition control */
#define CHECK_INTEGRAL 1     /* set to 1 to check integral of first field */

#define JULIA_SCALE 0.5 /* scaling for Julia sets */
#define JULIA_ROT -20.0       /* rotation of Julia set, in degrees */
#define JULIA_RE 0.5    
#define JULIA_IM 0.462    /* parameters for Julia sets */

/* Choice of the billiard table */

#define B_DOMAIN 999          /* choice of domain shape, see list in global_pdes.c  */

#define CIRCLE_PATTERN 8     /* pattern of circles, see list in global_pdes.c */

#define P_PERCOL 0.25       /* probability of having a circle in C_RAND_PERCOL arrangement */
#define NPOISSON 300        /* number of points for Poisson C_RAND_POISSON arrangement */
#define RANDOM_POLY_ANGLE 0 /* set to 1 to randomize angle of polygons */

#define LAMBDA 1.0	    /* parameter controlling the dimensions of domain */
#define MU 1.0	            /* parameter controlling the dimensions of domain */
#define NPOLY 5             /* number of sides of polygon */
#define APOLY 2.0           /* angle by which to turn polygon, in units of Pi/2 */
#define MDEPTH 7            /* depth of computation of Menger gasket */
#define MRATIO 5            /* ratio defining Menger gasket */
#define MANDELLEVEL 1000      /* iteration level for Mandelbrot set */
#define MANDELLIMIT 10.0     /* limit value for approximation of Mandelbrot set */
#define FOCI 1              /* set to 1 to draw focal points of ellipse */
#define NGRIDX 6            /* number of grid point for grid of disks */
#define NGRIDY 8            /* number of grid point for grid of disks */
#define REVERSE_TESLA_VALVE 1   /* set to 1 to orient Tesla valve in blocking configuration */
#define WALL_WIDTH 0.05      /* width of wall separating lenses */

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

#define DT 0.00000015

#define VISCOSITY 0.02
#define POISSON_STIFFNESS 1.0   /* stiffness of Poisson equation solver for incompressible Euler */
#define DISSIPATION 0.0

#define RPSA 0.75         /* parameter in Rock-Paper-Scissors-type interaction */
#define RPSLZB 0.0       /* second parameter in Rock-Paper-Scissors-Lizard-Spock type interaction */

#define K_AC 0.1        /* force constant in Allen-Cahn equation */
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
#define G_FIELD 0.02    /* gravity/Coriolis force */
#define BC_FIELD 1.0e-5     /* constant in repulsive field from obstacles */
#define AB_RADIUS 0.2   /* radius of region with magnetic field for Aharonov-Bohm effect */
#define K_EULER 50.0    /* constant in stream function integration of Euler equation */
#define K_EULER_INC 0.5    /* constant in incompressible Euler equation */
#define C_EULER_COMP 0.05   /* constant in compressible Euler equation */

#define SMOOTHEN_VORTICITY 0    /* set to 1 to smoothen vorticity field in Euler equation */
#define SMOOTHEN_VELOCITY 1     /* set to 1 to smoothen velocity field in Euler equation */
#define SMOOTHEN_PERIOD 10      /* period between smoothenings */
#define SMOOTH_FACTOR 0.05      /* factor by which to smoothen */

#define ADD_OSCILLATING_SOURCE 0        /* set to 1 to add an oscillating wave source */
#define OSCILLATING_SOURCE_PERIOD 1     /* period of oscillating source */
#define OSCILLATING_SOURCE_OMEGA 0.2    /* frequency of oscillating source */

#define ADD_TRACERS 1    /* set to 1 to add tracer particles (for Euler equations) */
#define N_TRACERS 2000    /* number of tracer particles */
#define TRACERS_STEP 0.005  /* step size in tracer evolution */

#define T_OUT 2.0       /* outside temperature */
#define T_IN 0.0        /* inside temperature */
#define SPEED 0.0       /* speed of drift to the right */

#define ADD_NOISE 0     /* set to 1 to add noise, set to 2 to add noise in right half */
#define NOISE_INTENSITY 0.01      /* noise intensity */
#define CHANGE_NOISE 0      /* set to 1 to increase noise intensity */
#define NOISE_FACTOR 40.0   /* factor by which to increase noise intensity */
#define NOISE_INITIAL_TIME 100  /* initial time during which noise remains constant */

#define CHANGE_VISCOSITY 0      /* set to 1 to change the viscosity in the course of the simulation */
#define ADJUST_INTSTEP 0       /* set to 1 to decrease integration step when viscosity increases */
#define VISCOSITY_INITIAL_TIME 10  /* initial time during which viscosity remains constant */
#define VISCOSITY_FACTOR 100.0   /* factor by which to change viscosity */
#define VISCOSITY_MAX 2.0        /* max value of viscosity beyond which NVID is increased and integration step is decrase, for numerical stability */
                                        
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

#define B_COND 1

#define B_COND_LEFT 0
#define B_COND_RIGHT 0
#define B_COND_TOP 0
#define B_COND_BOTTOM 0

/* Parameters for length and speed of simulation */

#define NSTEPS 1900           /* number of frames of movie */
#define NVID 100          /* number of iterations between images displayed on screen */
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
#define MID_FRAMES 100    /* number of still frames between parts of two-part movie */
#define END_FRAMES 250    /* number of still frames at end of movie */
#define FADE 1           /* set to 1 to fade at end of movie */

/* Visualisation */

#define PLOT_3D 0    /* controls whether plot is 2D or 3D */
#define PLOT_SPHERE 1   /* draws fields on a sphere */

#define ROTATE_VIEW 1       /* set to 1 to rotate position of observer */
#define ROTATE_ANGLE 360.0  /* total angle of rotation during simulation */
#define SHADE_3D 0              /* set to 1 to change luminosity according to normal vector */
#define SHADE_2D 1              /* set to 1 to change luminosity according to normal vector */
#define VIEWPOINT_TRAJ 0    /* type of viewpoint trajectory */
#define MAX_LATITUDE 45.0   /* maximal latitude for viewpoint trajectory VP_ORBIT2 */

#define DRAW_PERIODICISED 0     /* set to 1 to repeat wave periodically in x and y directions */

/* Plot type - color scheme */

#define CPLOT 60
#define CPLOT_B 61

/* Plot type - height of 3D plot */

#define ZPLOT 60     /* z coordinate in 3D plot */
#define ZPLOT_B 61    /* z coordinate in second 3D plot */

#define AMPLITUDE_HIGH_RES 1    /* set to 1 to increase resolution of P_3D_AMPLITUDE plot */
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

#define COLOR_PALETTE 11       /* Color palette, see list in global_pdes.c  */
#define COLOR_PALETTE_B 16     /* Color palette, see list in global_pdes.c  */

#define BLACK 1          /* black background */
#define COLOR_OUT_R 1.0    /* color outside domain */
#define COLOR_OUT_G 1.0    
#define COLOR_OUT_B 1.0    

#define COLOR_SCHEME 3   /* choice of color scheme */

#define COLOR_PHASE_SHIFT 0.25   /* phase shift of color scheme, in units of Pi */

#define SCALE 0          /* set to 1 to adjust color scheme to variance of field */
#define SLOPE 1.0        /* sensitivity of color on wave amplitude */
#define VSCALE_AMPLITUDE 100.0      /* additional scaling factor for color scheme P_3D_AMPLITUDE */
#define ATTENUATION 0.0  /* exponential attenuation coefficient of contrast with time */
#define CURL_SCALE 1.0   /* scaling factor for curl representation */
#define RESCALE_COLOR_IN_CENTER 0   /* set to 1 to decrease color intentiy in the center (for wave escaping ring) */
#define SLOPE_SCHROD_LUM 100.0       /* sensitivity of luminosity on module, for color scheme Z_ARGUMENT */
#define MIN_SCHROD_LUM 0.1       /* minimal luminosity in color scheme Z_ARGUMENT*/
#define VSCALE_PRESSURE 2.0      /* additional scaling factor for color scheme Z_EULER_PRESSURE */
#define PRESSURE_SHIFT 10.0        /* shift for color scheme Z_EULER_PRESSURE */
#define PRESSURE_LOG_SHIFT -2.5     /* shift for color scheme Z_EULER_PRESSURE */
#define VSCALE_WATER_HEIGHT 0.4     /* vertical scaling of water height */
#define SHADE_SCALE_2D 0.25     /* controls "depth" of 2D shading */

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
#define VSCALE_SPEED 100.0      /* additional scaling factor for color scheme Z_EULER_SPEED */
#define VMEAN_SPEED 1.04       /* mean value around which to scale for color scheme Z_EULER_SPEED */
#define SHIFT_DENSITY 1.0        /* shift for color scheme Z_EULER_DENSITY */
#define VSCALE_DENSITY 30.0      /* additional scaling factor for color scheme Z_EULER_DENSITY */
#define VSCALE_VORTICITY 15.0     /* additional scaling factor for color scheme Z_EULERC_VORTICITY */
#define VORTICITY_SHIFT 0.0     /* vertical shift of vorticity */
#define ZSCALE_SPEED 0.5        /* additional scaling factor for z-coord Z_EULER_SPEED and Z_SWATER_SPEED */
#define ZSHIFT_SPEED 1.0        /* additional shift of z-coord Z_EULER_SPEED and Z_SWATER_SPEED */
#define ZSCALE_NORMGRADIENT -0.0001  /* vertical scaling for Z_NORM_GRADIENT */
#define VSCALE_SWATER 250.0        /* additional scaling factor for color scheme Z_EULER_DENSITY */

#define NXMAZE 7      /* width of maze */
#define NYMAZE 7      /* height of maze */
#define MAZE_MAX_NGBH 4     /* max number of neighbours of maze cell */
#define RAND_SHIFT 3        /* seed of random number generator */
#define MAZE_XSHIFT 0.0     /* horizontal shift of maze */
#define MAZE_WIDTH 0.04     /* half width of maze walls */

#define DRAW_COLOR_SCHEME 0     /* set to 1 to plot the color scheme */
#define COLORBAR_RANGE 3.0      /* scale of color scheme bar */
#define COLORBAR_RANGE_B 2.5    /* scale of color scheme bar for 2nd part */
#define ROTATE_COLOR_SCHEME 0   /* set to 1 to draw color scheme horizontally */
#define CIRC_COLORBAR 0         /* set to 1 to draw circular color scheme */
#define CIRC_COLORBAR_B 0       /* set to 1 to draw circular color scheme */

/* only for compatibility with wave_common.c */
#define TWOSPEEDS 0          /* set to 1 to replace hardcore boundary by medium with different speed */
#define VARIABLE_IOR 0      /* set to 1 for a variable index of refraction */
#define IOR 4               /* choice of index of refraction, see list in global_pdes.c */
#define IOR_TOTAL_TURNS 1.5 /* total angle of rotation for IOR_PERIODIC_WELLS_ROTATING */
#define MANDEL_IOR_SCALE -0.05   /* parameter controlling dependence of IoR on Mandelbrot escape speed */
#define OMEGA 0.005        /* frequency of periodic excitation */
#define OSCIL_YMAX 0.2      /* defines oscillation range */
#define COURANT 0.08       /* Courant number */
#define COURANTB 0.03      /* Courant number in medium B */
#define INITIAL_AMP 0.5         /* amplitude of initial condition */
#define INITIAL_VARIANCE 0.0002  /* variance of initial condition */
#define INITIAL_WAVELENGTH  0.1  /* wavelength of initial condition */
#define VSCALE_ENERGY 200.0       /* additional scaling factor for color scheme P_3D_ENERGY */
#define PHASE_FACTOR 20.0       /* factor in computation of phase in color scheme P_3D_PHASE */
#define PHASE_SHIFT 0.0      /* shift of phase in color scheme P_3D_PHASE */
#define AVRG_E_FACTOR 0.99   /* controls time window size in P_AVERAGE_ENERGY scheme */
#define OSCILLATION_SCHEDULE 0  /* oscillation schedule, see list in global_pdes.c */
#define AMPLITUDE 0.8      /* amplitude of periodic excitation */ 
#define ACHIRP 0.2        /* acceleration coefficient in chirp */
#define DAMPING 0.0        /* damping of periodic excitation */
#define COMPARISON 0        /* set to 1 to compare two different patterns (beta) */
#define B_DOMAIN_B 20       /* second domain shape, for comparisons */
#define CIRCLE_PATTERN_B 0  /* second pattern of circles or polygons */
#define FLUX_WINDOW 20      /* averaging window for energy flux */
#define ADD_WAVE_PACKET_SOURCES 0       /* set to 1 to add several sources emitting wave packets */
#define WAVE_PACKET_SOURCE_TYPE 1       /* type of wave packet sources */
#define N_SOURCES 1                     /* number of wave sources */
#define N_WAVE_PACKETS 15               /* number of wave packets */
#define WAVE_PACKET_RADIUS 20            /* radius of wave packets */
#define OSCIL_LEFT_YSHIFT 25.0   /* y-dependence of left oscillation (for non-horizontal waves) */
#define INITIAL_SHIFT 20.0          /* time shift of initial wave packet (in oscillation periods) */
#define WAVE_PACKET_SHIFT 200.0     /* time shift between wave packets (in oscillation periods) */
#define DRAW_WAVE_PROFILE 0     /* set to 1 to draw a profile of the wave */
#define HORIZONTAL_WAVE_PROFILE 0 /* set to 1 to draw wave profile vertically */
#define VERTICAL_WAVE_PROFILE 0 /* set to 1 to draw wave profile vertically */
#define WAVE_PROFILE_X 1.9      /* value of x to sample wave profile */
#define WAVE_PROFILE_Y -1.0      /* value of y to sample wave profile */
#define PROFILE_AT_BOTTOM 1     /* draw wave profile at bottom instead of top */
#define AVERAGE_WAVE_PROFILE 1  /* set to 1 to draw time-average of wave profile squared*/
#define MU_B 1.0           /* parameter controlling the dimensions of domain */
#define GAMMA 0.0          /* damping factor in wave equation */
#define GAMMAB 0.0         /* damping factor in wave equation */
#define VERTICAL_WAVE_PROFILE 0 /* set to 1 to draw wave profile vertically */
#define DRAW_WAVE_TIMESERIES 0  /* set to 1 to draw a time series of the wave */
#define TIMESERIES_NVALUES 400  /* number of values plotted in time series */
#define DRAW_WAVE_SOURCE 0      /* set to 1 to draw source of wave at (wave_source_x, wave_source_y) */
#define MESSAGE_LDASH 1         /* length of dash for Morse code message */
#define MESSAGE_LDOT 1          /* length of dot for Morse code message */
#define MESSAGE_LINTERVAL 1     /* length of interval between dashes/dots for Morse code message */
#define MESSAGE_LINTERLETTER 1  /* length of interval between letters for Morse code message */
#define MESSAGE_LSPACE 1        /* length of space for Morse code message */
#define MESSAGE_INITIAL_TIME 1  /* initial time before starting message for Morse code message */    
/* end of constants added only for compatibility with wave_common.c */

double u_3d[2] = {0.75, -0.45};     /* projections of basis vectors for REP_AXO_3D representation */
double v_3d[2] = {-0.75, -0.45};
double w_3d[2] = {0.0, 0.015};
double light[3] = {0.816496581, -0.40824829, 0.40824829};      /* vector of "light" direction for P_3D_ANGLE color scheme */
double observer[3] = {-8.0, -4.0, 4.0};    /* location of observer for REP_PROJ_3D representation */ 
int reset_view = 0;         /* switch to reset 3D view parameters (for option ROTATE_VIEW) */

/* constants for simulations on planets */
#define ADD_DEM 1               /* add DEM (digital elevation model) */
#define ADD_NEGATIVE_DEM 0      /* add DEM with bathymetric data */
#define RSCALE_DEM 0.1          /* scaling factor of radial component for DEM */
#define SMOOTH_DEM 0            /* set to 1 to smoothen DEM (to make altitude less constant) */
#define DEM_SMOOTH_STEPS 1      /* number of smoothening steps */
#define DEM_SMOOTH_HEIGHT 2.0   /* relative height below which to smoothen */
#define DEM_MAXHEIGHT 9000.0     /* max height of DEM (estimated from Everest/Olympus Mons) */
#define DEM_MAXDEPTH -10000     /* max depth of DEM */
#define PLANET_SEALEVEL 0.0      /* sea level for flooded planet */
#define VENUS_NODATA_FACTOR 0.5     /* altitude to assign to DEM points without data (fraction of mean altitude) */

#define Z_SCALING_FACTOR 0.8   /* overall scaling factor of z axis for REP_PROJ_3D representation */
#define XY_SCALING_FACTOR 2.0  /* overall scaling factor for on-screen (x,y) coordinates after projection */
#define ZMAX_FACTOR 1.0        /* max value of z coordinate for REP_PROJ_3D representation */
#define XSHIFT_3D 0.0          /* overall x shift for REP_PROJ_3D representation */
#define YSHIFT_3D 0.0          /* overall y shift for REP_PROJ_3D representation */
#define BORDER_PADDING 0       /* distance from boundary at which to plot points, to avoid boundary effects due to gradient */
#define DRAW_ARROW 0           /* set to 1 to draw arrow above sphere */

#define RSCALE 0.25            /* scaling factor of radial component */
#define RSHIFT -0.25           /* shift in radial component */
#define RMAX 2.0               /* max value of radial component */
#define RMIN 0.5               /* min value of radial component */
#define COS_VISIBLE -1.1       /* limit on cosine of normal to shown facets */

/* For debugging purposes only */
#define FLOOR 1         /* set to 1 to limit wave amplitude to VMAX */
#define VMAX 100.0        /* max value of wave amplitude */
#define TEST_GRADIENT 0 /* print norm squared of gradient */

```

