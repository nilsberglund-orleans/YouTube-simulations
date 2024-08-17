/*********************************************************************************/
/*                                                                               */
/*  Animation of reaction-diffusion equation in a planar domain                  */
/*                                                                               */
/*  N. Berglund, January 2022                                                    */
/*                                                                               */
/*  Feel free to reuse, but if doing so it would be nice to drop a               */
/*  line to nils.berglund@univ-orleans.fr - Thanks!                              */
/*                                                                               */
/*  compile with                                                                 */
/*  gcc -o rde rde.c                                                             */
/* -L/usr/X11R6/lib -ltiff -lm -lGL -lGLU -lX11 -lXmu -lglut -O3 -fopenmp        */
/*                                                                               */
/*  OMP acceleration may be more effective after executing                       */
/*  export OMP_NUM_THREADS=2 in the shell before running the program             */
/*                                                                               */
/*  To make a video, set MOVIE to 1 and create subfolder tif_bz                  */
/*  It may be possible to increase parameter PAUSE                               */
/*                                                                               */
/*  create movie using                                                           */
/*  ffmpeg -i wave.%05d.tif -vcodec libx264 wave.mp4                             */
/*                                                                               */
/*********************************************************************************/

/*********************************************************************************/
/*                                                                               */
/* NB: The algorithm used to simulate the wave equation is highly paralellizable */
/* One could make it much faster by using a GPU                                  */
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

#define MOVIE 0         /* set to 1 to generate movie */
#define DOUBLE_MOVIE 1  /* set to 1 to produce movies for wave height and energy simultaneously */
#define SAVE_MEMORY 1           /* set to 1 to save memory when writing tiff images */
#define NO_EXTRA_BUFFER_SWAP 1    /* some OS require one less buffer swap when recording images */

/* General geometrical parameters */

#define WINWIDTH 	1920  /* window width */
#define WINHEIGHT 	1150  /* window height */
#define NX 960          /* number of grid points on x axis */
#define NY 575          /* number of grid points on y axis */
#define HRES 1          /* factor for high resolution plots */

#define XMIN -2.0
#define XMAX 2.0	/* x interval */
#define YMIN -1.041666667
#define YMAX 1.041666667	/* y interval for 9/16 aspect ratio */

/* Choice of simulated equation */

#define RDE_EQUATION 8  /* choice of reaction term, see list in global_3d.c */
#define NFIELDS 3       /* number of fields in reaction-diffusion equation */
#define NLAPLACIANS 0   /* number of fields for which to compute Laplacian */

#define SPHERE 1        /* set to 1 to simulate equation on sphere */
#define DPOLE 1         /* safety distance to poles */
#define DSMOOTH 1       /* size of neighbourhood of poles that are smoothed */
#define SMOOTHPOLE 0.05  /* smoothing coefficient at poles */
#define SMOOTHCOTPOLE 0.05  /* smoothing coefficient of cotangent at poles */
#define PHISHIFT 0.0    /* shift of phi in 2D plot (in degrees) */
#define SMOOTHBLOCKS 1  /* set to 1 to use blocks of points near the poles */
#define BLOCKDIST 64    /* distance to poles where points are blocked */ 
#define ZERO_MERIDIAN 190.0     /* choice of zero meridian (will be at left/right boundary of 2d plot) */
#define POLE_NODRAW 10  /* distance around poles where wave is not drawn */

#define ADD_POTENTIAL 0 /* set to 1 to add a potential (for Schrodinger equation) */
#define ADD_MAGNETIC_FIELD 0    /* set to 1 to add a magnetic field (for Schrodinger equation) - then set POTENTIAL 1 */
#define ADD_FORCE_FIELD 0   /* set to 1 to add a foce field (for compressible Euler equation) */
#define POTENTIAL 7         /* type of potential or vector potential, see list in global_3d.c  */
#define FORCE_FIELD 6       /* type of force field, see list in global_3d.c  */
#define ADD_CORIOLIS_FORCE 1    /* set to 1 to add Coriolis force (quasigeostrophic Euler equations) */
#define VARIABLE_DEPTH 1    /* set to 1 for variable depth in shallow water equation */
#define SWATER_DEPTH 5      /* variable depth in shallow water equation */

#define ANTISYMMETRIZE_WAVE_FCT 0   /* set tot 1 to make wave function antisymmetric */
#define ADAPT_STATE_TO_BC 0     /* to smoothly adapt initial state to obstacles */
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
#define PDISC_FACTOR 3.25    /* controls density of Poisson disc process (default: 3.25) */
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
#define RADIUS_FACTOR 0.3   /* controls inner radius for C_RING arrangements */

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

#define VISCOSITY 0.02
#define POISSON_STIFFNESS 1.0   /* stiffness of Poisson equation solver for incompressible Euler */
#define DISSIPATION 1.0e-8

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
#define G_FIELD 0.002   /* gravity/Coriolis force */
#define BC_FIELD 1.0e-5     /* constant in repulsive field from obstacles */
#define AB_RADIUS 0.2   /* radius of region with magnetic field for Aharonov-Bohm effect */
#define K_EULER 50.0    /* constant in stream function integration of Euler equation */
#define K_EULER_INC 0.5    /* constant in incompressible Euler equation */
#define C_EULER_COMP 0.1   /* constant in compressible Euler equation */

#define SMOOTHEN_VORTICITY 0    /* set to 1 to smoothen vorticity field in Euler equation */
#define SMOOTHEN_VELOCITY 1     /* set to 1 to smoothen velocity field in Euler equation */
#define SMOOTHEN_PERIOD 7       /* period between smoothenings */
#define SMOOTH_FACTOR 0.15      /* factor by which to smoothen */

#define ADD_OSCILLATING_SOURCE 0        /* set to 1 to add an oscillating wave source */
#define OSCILLATING_SOURCE_PERIOD 1     /* period of oscillating source */
#define OSCILLATING_SOURCE_OMEGA 0.2    /* frequency of oscillating source */

#define ADD_TRACERS 1    /* set to 1 to add tracer particles (for Euler equations) */
#define N_TRACERS 2000    /* number of tracer particles */
#define TRACERS_STEP 0.1  /* step size in tracer evolution */

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

#define SWATER_MIN_HEIGHT 0.025      /* min height of initial condition for shallow water equation */
#define DEPTH_FACTOR 0.075       /* proportion of min height in variable depth */
#define TANH_FACTOR 1.0         /* steepness of variable depth */

#define EULER_GRADIENT_YSHIFT 0.0    /* y-shift in computation of gradient in Euler equation */

/* Boundary conditions, see list in global_pdes.c  */

#define B_COND 1

#define B_COND_LEFT 0
#define B_COND_RIGHT 0
#define B_COND_TOP 0
#define B_COND_BOTTOM 0

/* Parameters for length and speed of simulation */

#define NSTEPS 950           /* number of frames of movie */
#define NVID 85           /* number of iterations between images displayed on screen */
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
#define ROTATE_ANGLE 360.0   /* total angle of rotation during simulation */
#define SHADE_3D 1              /* set to 1 to change luminosity according to normal vector */
#define SHADE_2D 0              /* set to 1 to change luminosity according to normal vector */
#define VIEWPOINT_TRAJ 1    /* type of viewpoint trajectory */
#define MAX_LATITUDE 45.0   /* maximal latitude for viewpoint trajectory VP_ORBIT2 */

#define DRAW_PERIODICISED 0     /* set to 1 to repeat wave periodically in x and y directions */

/* Plot type - color scheme */

#define CPLOT 70
#define CPLOT_B 74

/* Plot type - height of 3D plot */

#define ZPLOT 70     /* z coordinate in 3D plot */
#define ZPLOT_B 71    /* z coordinate in second 3D plot */

#define AMPLITUDE_HIGH_RES 1    /* set to 1 to increase resolution of P_3D_AMPLITUDE plot */
#define NON_DIRICHLET_BC 0      /* set to 1 to draw only facets in domain, if field is not zero on boundary */
#define WRAP_ANGLE 1            /* experimental: wrap angle to [0, 2Pi) for interpolation in angle schemes */
#define FADE_IN_OBSTACLE 0      /* set to 1 to fade color inside obstacles */
#define FADE_WATER_DEPTH 0      /* set to 1 to make wave color depth-dependent */
#define DRAW_OUTSIDE_GRAY 0     /* experimental - draw outside of billiard in gray */
#define ADD_POTENTIAL_TO_Z 0    /* set to 1 to add the external potential to z-coordinate of plot */
#define ADD_POT_CONSTANT 0.35   /* constant added to wave height */
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

#define COLOR_PALETTE 10       /* Color palette, see list in global_pdes.c  */
#define COLOR_PALETTE_B 12      /* Color palette, see list in global_pdes.c  */

#define BLACK 1          /* black background */
#define COLOR_OUT_R 1.0    /* color outside domain */
#define COLOR_OUT_G 1.0    
#define COLOR_OUT_B 1.0    

#define COLOR_SCHEME 3   /* choice of color scheme */

#define PHASE_SHIFT -0.25   /* phase shift of color scheme, in units of Pi (formerly COLOR_PHASE_SHIFT) */

#define SCALE 0          /* set to 1 to adjust color scheme to variance of field */
#define SLOPE 1.0        /* sensitivity of color on wave amplitude */
#define VSHIFT_AMPLITUDE 0.0   /* additional shift for wave amplitude */
#define VSCALE_AMPLITUDE 15.0      /* additional scaling factor for color scheme P_3D_AMPLITUDE */
#define ATTENUATION 0.0  /* exponential attenuation coefficient of contrast with time */
#define CURL_SCALE 1.0   /* scaling factor for curl representation */
#define RESCALE_COLOR_IN_CENTER 0   /* set to 1 to decrease color intentiy in the center (for wave escaping ring) */
#define SLOPE_SCHROD_LUM 10.0      /* sensitivity of luminosity on module, for color scheme Z_ARGUMENT */
#define MIN_SCHROD_LUM 0.1       /* minimal luminosity in color scheme Z_ARGUMENT*/
#define VSCALE_PRESSURE 2.0      /* additional scaling factor for color scheme Z_EULER_PRESSURE */
#define PRESSURE_SHIFT 10.0        /* shift for color scheme Z_EULER_PRESSURE */
#define PRESSURE_LOG_SHIFT -2.5     /* shift for color scheme Z_EULER_PRESSURE */
#define VSCALE_WATER_HEIGHT 30.0     /* vertical scaling of water height */
#define ADD_HEIGHT_CONSTANT -0.025   /* constant added to wave height */
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
#define ZSCALE_SPEED 300.0        /* additional scaling factor for z-coord Z_EULER_SPEED and Z_SWATER_SPEED */
#define ZSHIFT_SPEED 0.0        /* additional shift of z-coord Z_EULER_SPEED and Z_SWATER_SPEED */
#define ZSCALE_NORMGRADIENT -0.0001  /* vertical scaling for Z_NORM_GRADIENT */
#define VSCALE_SWATER 50.0        /* additional scaling factor for color scheme Z_EULER_DENSITY */

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
#define OSCIL_YMAX 0.2      /* defines oscillation range */
#define COURANT 0.08       /* Courant number */
#define COURANTB 0.03      /* Courant number in medium B */
#define INITIAL_AMP 0.5         /* amplitude of initial condition */
#define INITIAL_VARIANCE 0.0002  /* variance of initial condition */
#define INITIAL_WAVELENGTH  0.1  /* wavelength of initial condition */
#define VSCALE_ENERGY 200.0       /* additional scaling factor for color scheme P_3D_ENERGY */
#define PHASE_FACTOR 20.0       /* factor in computation of phase in color scheme P_3D_PHASE */
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
double light[3] = {0.40824829, -0.816496581, 0.40824829};      /* vector of "light" direction for P_3D_ANGLE color scheme */
double observer[3] = {8.0, 4.0, 2.5};    /* location of observer for REP_PROJ_3D representation */ 
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
#define PLANET_SEALEVEL 2500.0      /* sea level for flooded planet */
#define VENUS_NODATA_FACTOR 0.5     /* altitude to assign to DEM points without data (fraction of mean altitude) */

#define Z_SCALING_FACTOR 0.8   /* overall scaling factor of z axis for REP_PROJ_3D representation */
#define XY_SCALING_FACTOR 2.0  /* overall scaling factor for on-screen (x,y) coordinates after projection */
#define ZMAX_FACTOR 1.0        /* max value of z coordinate for REP_PROJ_3D representation */
#define XSHIFT_3D 0.0          /* overall x shift for REP_PROJ_3D representation */
#define YSHIFT_3D 0.0          /* overall y shift for REP_PROJ_3D representation */
#define BORDER_PADDING 0       /* distance from boundary at which to plot points, to avoid boundary effects due to gradient */
#define DRAW_ARROW 1           /* set to 1 to draw arrow above sphere */

#define RSCALE 0.2               /* scaling factor of radial component */
#define RSHIFT -0.01             /* shift in radial component */
#define RMAX 2.0                 /* max value of radial component */
#define RMIN 0.5                 /* min value of radial component */
#define COS_VISIBLE -0.75         /* limit on cosine of normal to shown facets */

/* For debugging purposes only */
#define FLOOR 1         /* set to 1 to limit wave amplitude to VMAX */
#define VMAX 100.0        /* max value of wave amplitude */
#define TEST_GRADIENT 0 /* print norm squared of gradient */


#define REFRESH_B (ZPLOT_B != ZPLOT)||(CPLOT_B != CPLOT)    /* to save computing time, to be improved */
#define COMPUTE_WRAP_ANGLE ((WRAP_ANGLE)&&((cplot == Z_ANGLE_GRADIENT)||(cplot == Z_ANGLE_GRADIENTX)||(cplot == Z_ARGUMENT)||(cplot == Z_ANGLE_GRADIENTX)||(cplot == Z_EULER_DIRECTION_SPEED)||(cplot == Z_SWATER_DIRECTION_SPEED)))
#define PRINT_PARAMETERS ((PRINT_TIME)||(PRINT_VISCOSITY)||(PRINT_RPSLZB)||(PRINT_PROBABILITIES)||(PRINT_NOISE)||(PRINT_FLOW_SPEED)||(PRINT_AVERAGE_SPEED))
#define COMPUTE_PRESSURE ((ZPLOT == Z_EULER_PRESSURE)||(CPLOT == Z_EULER_PRESSURE)||(ZPLOT_B == Z_EULER_PRESSURE)||(CPLOT_B == Z_EULER_PRESSURE))

#define ASYM_SPEED_COLOR (VMEAN_SPEED == 0.0)
#define XYIN_INITIALISED (B_DOMAIN == D_IMAGE)

int block_sizes[NY];      /* table of block sizes for blocking around poles */
int block_numbers[NY];    /* table of block numbers for blocking around poles */

#include "global_pdes.c"
#include "global_3d.c"          /* constants and global variables */

#include "sub_maze.c"
#include "sub_wave.c"
#include "wave_common.c"        /* common functions for wave_billiard, wave_comparison, etc */

#include "sub_wave_3d_rde.c"    /* should be later replaced by sub_wave_rde.c */

#include "sub_rde.c"    


double f_aharonov_bohm(double r2)
/* radial part of Aharonov-Bohm vector potential */
{
    double r02 = AB_RADIUS*AB_RADIUS;
    
    if (r2 > r02) return(-0.25*r02/r2);
    else return(0.25*(r2 - 2.0*r02)/r02);

//     if (r2 > r02) return(1.0/r2);
//     else return((2.0*r02 - r2)/(r02*r02));    
}

double potential(int i, int j)
/* compute potential (e.g. for Schrödinger equation), or potential part if there is a magnetic field */
{
    double x, y, xy[2], r, small = 1.0e-1, kx, ky, lx = XMAX - XMIN, r1, r2, r3, f;
    int rect;
    
    ij_to_xy(i, j, xy);
    x = xy[0];
    y = xy[1];
    
    switch (POTENTIAL) {
        case (POT_HARMONIC):
        {
            return (K_HARMONIC*(x*x + y*y));
        }
        case (POT_COULOMB):
        {
//             r = module2(x, y);
            r = sqrt(x*x + y*y + small*small);
//             if (r < small) r = small;
            return (-K_COULOMB/r);
        }
        case (POT_PERIODIC):
        {
            kx = 4.0*DPI/(XMAX - XMIN);
            ky = 2.0*DPI/(YMAX - YMIN);
            return(-K_HARMONIC*cos(kx*x)*cos(ky*y));
        }
        case (POT_FERMIONS):
        {
            r = sqrt((x-y)*(x-y) + small*small);
            return (-K_COULOMB/r);
        }
        case (POT_FERMIONS_PERIODIC):
        {
            r1 = sqrt((x-y)*(x-y) + small*small);
            r2 = sqrt((x-lx-y)*(x-lx-y) + small*small);
            r3 = sqrt((x+lx-y)*(x+lx-y) + small*small);
//             r = r/3.0;
            return (-0.5*K_COULOMB*(1.0/r1 + 1.0/r2 + 1.0/r3));
        }
        case (VPOT_CONSTANT_FIELD):
        {
            return (K_HARMONIC*(x*x + y*y));        /* magnetic field strength b is chosen such that b^2 = K_HARMONIC */
        }
        case (VPOT_AHARONOV_BOHM):
        {
            r2 = x*x + y*y;
            f = f_aharonov_bohm(r2);
            return (B_FIELD*B_FIELD*f*f*r2);    /* magnetic field strength b is chosen such that b^2 = K_HARMONIC */
//             return (K_HARMONIC*f);    /* magnetic field strength b is chosen such that b^2 = K_HARMONIC */
        }
        case (POT_MAZE):
        {
            for (rect=0; rect<npolyrect; rect++)
                if (ij_in_polyrect(i, j, polyrect[rect])) return(V_MAZE);
            return(0.0);    
        }
        default:
        {
            return(0.0);
        }
    }
}   


void compute_vector_potential(int i, int j, double *ax, double *ay)
/* initialize the vector potential, for Schrodinger equation in a magnetic field */
{
    double x, y, xy[2], r2, f;
    
    ij_to_xy(i, j, xy);
    x = xy[0];
    y = xy[1];
    
    switch (POTENTIAL) {
        case (VPOT_CONSTANT_FIELD):
        {
            *ax = B_FIELD*y;
            *ay = -B_FIELD*x;
            break;
        }
        case (VPOT_AHARONOV_BOHM):
        {
            r2 = x*x + y*y;
            f = f_aharonov_bohm(r2);
            *ax = B_FIELD*y*f;
            *ay = -B_FIELD*x*f;
            break;
        }
        default:
        {
            *ax = 0.0;
            *ay = 0.0;
        }
    }
}

void compute_gfield(int i, int j, double *gx, double *gy)
/* initialize the exterior field, for the compressible Euler equation */
{
    double x, y, xy[2], r, f, a = 0.4, x1, y1, hx, hy, h;
    
    ij_to_xy(i, j, xy);
    x = xy[0];
    y = xy[1];
    
    switch (FORCE_FIELD) {
        case (GF_VERTICAL):
        {
            *gx = 0.0;
            *gy = -G_FIELD;
            break;
        }
        case (GF_CIRCLE):
        {
            r = module2(x,y) + 1.0e-2;
            f = 0.5*(1.0 - tanh(BC_STIFFNESS*(r - LAMBDA))); 
            *gx = BC_FIELD*f*x/r;
            *gy = BC_FIELD*f*y/r;
            break;
        }
        case (GF_ELLIPSE):
        {
            r = module2(x/LAMBDA,y/MU) + 1.0e-2;
            f = 0.5*(1.0 - tanh(BC_STIFFNESS*(r - 1.0))); 
            *gx = BC_FIELD*f*x/(LAMBDA*LAMBDA);
            *gy = BC_FIELD*f*y/(MU*MU);
            break;
        }
        case (GF_AIRFOIL):
        {
            y1 = y + a*x*x;
            r = module2(x/LAMBDA,y1/MU) + 1.0e-2;
            f = 0.5*(1.0 - tanh(BC_STIFFNESS*(r - 1.0))); 
            *gx = BC_FIELD*f*(x/(LAMBDA*LAMBDA) + a*y1/(MU*MU));
            *gy = BC_FIELD*f*y1/(MU*MU);
            break;
        }
        case (GF_WING):
        {
            if (x >= LAMBDA)
            {
                *gx = 0.0;
                *gy = 0.0;
            }
            else
            {
                x1 = 1.0 - x/LAMBDA;
                if (x1 < 0.1) x1 = 0.1;
                y1 = y + a*x*x;
                r = module2(x/LAMBDA,y1/(MU*x1)) + 1.0e-2;
                f = 0.5*(1.0 - tanh(BC_STIFFNESS*(r - 1.0))); 
                *gx = BC_FIELD*f*(x/(LAMBDA*LAMBDA) + 2.0*a*x*y1/(MU*MU*x1*x1) - y1*y1/(MU*MU*x1*x1*x1));
                *gy = BC_FIELD*f*y1/(MU*MU*x1*x1);
//                 *gx = 0.1*G_FIELD*f*(x/(LAMBDA*LAMBDA) + 2.0*a*x*y1/(MU*MU*x1*x1) - y1*y1/(MU*MU*x1*x1*x1));
//                 *gy = 0.1*G_FIELD*f*y1/(MU*MU*x1*x1);
//                 hx = x/(LAMBDA*LAMBDA) + 2.0*a*x*y1/(MU*MU*x1*x1) - y1*y1/(MU*MU*x1*x1*x1);
//                 hy = y1/(MU*MU*x1*x1);
//                 h = module2(hx, hy) + 1.0e-2;
//                 *gx = G_FIELD*f*hx/h;
//                 *gy = G_FIELD*f*hy/h;
            }
            break;
        }
        default:
        {
            *gx = 0.0;
            *gy = 0.0;
        }
    }
}
void initialize_potential(double potential_field[NX*NY])
/* initialize the potential field, e.g. for the Schrödinger equation */
{
    int i, j;
    
    #pragma omp parallel for private(i,j)
    for (i=0; i<NX; i++){
        for (j=0; j<NY; j++){
            potential_field[i*NY+j] = potential(i,j);
        }
    }
}

void initialize_vector_potential(double vpotential_field[2*NX*NY])
/* initialize the potential field, e.g. for the Schrödinger equation */
{
    int i, j;
    
    #pragma omp parallel for private(i,j)
    for (i=0; i<NX; i++){
        for (j=0; j<NY; j++){
            compute_vector_potential(i, j, &vpotential_field[i*NY+j], &vpotential_field[NX*NY+i*NY+j]);
        }
    }
}

void initialize_gfield(double gfield[2*NX*NY], double bc_field[NX*NY], double bc_field2[NX*NY], t_wave_sphere *wsphere, t_wave_sphere *wsphere_hr)
/* initialize the exterior field, e.g. for the compressible Euler equation */
{
    int i, j;
    double dx, dy;
        
    switch (FORCE_FIELD) 
    {
        case (GF_COMPUTE_FROM_BC):
        {
            dx = (XMAX - XMIN)/(double)NX;
            dy = (YMAX - YMIN)/(double)NY;

            #pragma omp parallel for private(i,j)
            for (i=1; i<NX-1; i++){
                for (j=1; j<NY-1; j++){
                    gfield[i*NY+j] = BC_FIELD*(bc_field2[(i+1)*NY+j] - bc_field2[(i-1)*NY+j])/dx;
                    gfield[NX*NY+i*NY+j] = BC_FIELD*(bc_field2[i*NY+j+1] - bc_field2[i*NY+j-1])/dy;
//                 printf("gfield at (%i,%i): (%.3lg, %.3lg)\n", i, j, gfield[i*NY+j], gfield[NX*NY+i*NY+j]);
                }
            }
        
            /* boundaries */
            for (i=0; i<NX; i++)
            {
                gfield[i*NY] = 0.0;
                gfield[NX*NY+i*NY] = 0.0;
                gfield[i*NY+NY-1] = 0.0;
                gfield[NX*NY+i*NY+NY-1] = 0.0;
            }
            for (j=0; j<NY; j++)
            {
                gfield[j] = 0.0;
                gfield[NX*NY+j] = 0.0;
                gfield[(NX-1)*NY+j] = 0.0;
                gfield[NX*NY+(NX-1)*NY+j] = 0.0;            
            }
            break;
        }
        case (GF_EARTH):
        {
            dx = (XMAX - XMIN)/(double)NX;
            dy = (YMAX - YMIN)/(double)NY;
            init_earth_map_rde(wsphere, 1);
            init_earth_map_rde(wsphere_hr, HRES);
        
            #pragma omp parallel for private(i,j)
            for (i=1; i<NX-1; i++){
                for (j=1; j<NY-1; j++){
                    gfield[i*NY+j] = BC_FIELD*(wsphere[(i+1)*NY+j].altitude - wsphere[(i-1)*NY+j].altitude)/dx;
                    gfield[NX*NY+i*NY+j] = BC_FIELD*(wsphere[i*NY+j+1].altitude - wsphere[i*NY+j-1].altitude)/dy;
                }
            }
        
            /* boundaries TODO */
            for (i=0; i<NX; i++)
            {
                gfield[i*NY] = 0.0;
                gfield[NX*NY+i*NY] = 0.0;
                gfield[i*NY+NY-1] = 0.0;
                gfield[NX*NY+i*NY+NY-1] = 0.0;
            }
            for (j=0; j<NY; j++)
            {
                gfield[j] = 0.0;
                gfield[NX*NY+j] = 0.0;
                gfield[(NX-1)*NY+j] = 0.0;
                gfield[NX*NY+(NX-1)*NY+j] = 0.0;            
            }
            break;
        }
        case (GF_MARS):
        {
            dx = (XMAX - XMIN)/(double)NX;
            dy = (YMAX - YMIN)/(double)NY;
            init_planet_map_rde(wsphere, D_SPHERE_MARS, 1);
            init_planet_map_rde(wsphere_hr, D_SPHERE_MARS, HRES);
        
            #pragma omp parallel for private(i,j)
            for (i=1; i<NX-1; i++){
                for (j=1; j<NY-1; j++){
                    gfield[i*NY+j] = BC_FIELD*(wsphere[(i+1)*NY+j].altitude - wsphere[(i-1)*NY+j].altitude)/dx;
                    gfield[NX*NY+i*NY+j] = BC_FIELD*(wsphere[i*NY+j+1].altitude - wsphere[i*NY+j-1].altitude)/dy;
                }
            }
        
            /* boundaries TODO */
            for (i=0; i<NX; i++)
            {
                gfield[i*NY] = 0.0;
                gfield[NX*NY+i*NY] = 0.0;
                gfield[i*NY+NY-1] = 0.0;
                gfield[NX*NY+i*NY+NY-1] = 0.0;
            }
            for (j=0; j<NY; j++)
            {
                gfield[j] = 0.0;
                gfield[NX*NY+j] = 0.0;
                gfield[(NX-1)*NY+j] = 0.0;
                gfield[NX*NY+(NX-1)*NY+j] = 0.0;            
            }
            break;
        }
        case (GF_VENUS):
        {
            dx = (XMAX - XMIN)/(double)NX;
            dy = (YMAX - YMIN)/(double)NY;
            init_planet_map_rde(wsphere, D_SPHERE_VENUS, 1);
            init_planet_map_rde(wsphere_hr, D_SPHERE_VENUS, HRES);
        
            #pragma omp parallel for private(i,j)
            for (i=1; i<NX-1; i++){
                for (j=1; j<NY-1; j++){
                    gfield[i*NY+j] = BC_FIELD*(wsphere[(i+1)*NY+j].altitude - wsphere[(i-1)*NY+j].altitude)/dx;
                    gfield[NX*NY+i*NY+j] = BC_FIELD*(wsphere[i*NY+j+1].altitude - wsphere[i*NY+j-1].altitude)/dy;
                }
            }
        
            /* boundaries TODO */
            for (i=0; i<NX; i++)
            {
                gfield[i*NY] = 0.0;
                gfield[NX*NY+i*NY] = 0.0;
                gfield[i*NY+NY-1] = 0.0;
                gfield[NX*NY+i*NY+NY-1] = 0.0;
            }
            for (j=0; j<NY; j++)
            {
                gfield[j] = 0.0;
                gfield[NX*NY+j] = 0.0;
                gfield[(NX-1)*NY+j] = 0.0;
                gfield[NX*NY+(NX-1)*NY+j] = 0.0;            
            }
            break;
        }
        default:
        {
            #pragma omp parallel for private(i,j)
            for (i=0; i<NX; i++){
                for (j=0; j<NY; j++){
                    compute_gfield(i, j, &gfield[i*NY+j], &gfield[NX*NY+i*NY+j]);
                }
            }
        }
    }
}

void compute_water_depth(int i, int j, t_rde *rde, t_wave_sphere wsphere[NX*NY], int ncircles)
/* initialize the vector potential, for Schrodinger equation in a magnetic field */
{
    double x, y, xy[2], r0, r, z, h, h0, hh, x1, y1, z1, d, rmax, rmin, phi;
    int n, nx, ny;
    static double phi1, phi2, phi3, sq3, rico;
    static int first = 1;
    
    if (first)
    {
        phi = 0.5*(sqrt(5.0)+1.0);
        phi1 = phi - 1.0;
        phi2 = phi - 2.0;
        phi3 = 2.0*phi - 3.0;
        sq3 = 0.5/sqrt(3.0);
        rico = 0.5/(2.0 + phi);
        first = 0;
    }
    
    ij_to_xy(i, j, xy);
    x = xy[0];
    y = xy[1];
    
    switch (SWATER_DEPTH) {
        case (SH_CIRCLE):
        {
            r0 = module2(x,y);
            r = r0/LAMBDA;
            h = tanh(TANH_FACTOR*(r-1.0));
            hh = 1.0 - h*h;
            z = SWATER_MIN_HEIGHT*DEPTH_FACTOR*TANH_FACTOR*hh/(r0*LAMBDA);
            rde->depth = SWATER_MIN_HEIGHT*DEPTH_FACTOR*h;
            rde->gradx = x*z;
            rde->grady = y*z;
            break;
        }
        case (SH_CIRCLES):
        {
            h = 0.0;
            z = 0.0;
            r = 10.0;
            /* compute minimal distance to circles */
            for (n = 0; n < ncircles; n++)
                for (nx = -1; nx < 2; nx++)
                    for (ny = -1; ny < 2; ny++)
                    {
                        x1 = circles[n].xc + (double)nx*(XMAX-XMIN);
                        y1 = circles[n].yc + (double)ny*(YMAX-YMIN);
                        r0 = module2(x - x1,y - y1)/circles[n].radius;
                        if (r0 < r) r = r0;
                    }
            h = tanh(TANH_FACTOR*(r-1.0));
            h = 0.5*(h + 1.0);
            rde->depth = SWATER_MIN_HEIGHT*DEPTH_FACTOR*h;
            break;
        }
        case (SH_COAST):
        {
            x1 = PI*(x + 0.1 - 0.1*cos(PI*y/YMAX))/XMAX;
            d = 3.0*sin(x1);
            h = tanh(-TANH_FACTOR*d);
            h = 0.5*(h + 1.0);
            rde->depth = SWATER_MIN_HEIGHT*DEPTH_FACTOR*h;
            break;
        }
        case (SH_COAST_MONOTONE):
        {
            x1 = PI*(x + 0.1 - 0.2*cos(PI*y/YMAX))/XMAX;
            h = tanh(-TANH_FACTOR*x1);
            h = 0.5*(h + 1.0);
            rde->depth = SWATER_MIN_HEIGHT*DEPTH_FACTOR*h;
            break;
        }
        case (SH_SPHERE_CUBE):
        {
            rmax = 0.0;
            /* compute distance from cube to origin, given by 0.5/rmax */
            if (vabs(wsphere[i*NY+j].x) > rmax) rmax = vabs(wsphere[i*NY+j].x);
            if (vabs(wsphere[i*NY+j].y) > rmax) rmax = vabs(wsphere[i*NY+j].y);
            if (vabs(wsphere[i*NY+j].z) > rmax) rmax = vabs(wsphere[i*NY+j].z);
            h = 1.0 - 0.5/rmax;
//             printf("h = %.3lg\n", h);
            if (h < 0.0) h = 0.0;
            rde->depth = SWATER_MIN_HEIGHT*DEPTH_FACTOR*h;
            break;
        }
        case (SH_SPHERE_OCTAHEDRON):
        {
            rmax = 0.0;
            /* compute distance from octahedron to origin, given by 0.5/rmax */
            rmax = 1.0/(vabs(wsphere[i*NY+j].x) + vabs(wsphere[i*NY+j].y) + vabs(wsphere[i*NY+j].z));
            h = 1.0 - 0.5/rmax;
            printf("h = %.3lg\n", h);
            if (h < 0.0) h = 0.0;
            rde->depth = SWATER_MIN_HEIGHT*DEPTH_FACTOR*h;
            break;
        }
        case (SH_SPHERE_DODECAHEDRON):
        {
            rmax = 0.0;
            x1 = wsphere[i*NY+j].x;
            y1 = wsphere[i*NY+j].y;
            z1 = wsphere[i*NY+j].z;
            /* compute distance from dodecahedron */
            d = vabs(phi1*y1 - phi2*z1) - sq3;
            if (d > rmax) rmax = d; 
            d = vabs(-phi2*x1 + phi1*z1) - sq3;
            if (d > rmax) rmax = d; 
            d = vabs(phi1*x1 - phi2*y1) - sq3;
            if (d > rmax) rmax = d; 
            d = vabs(phi1*y1 + phi2*z1) - sq3;
            if (d > rmax) rmax = d; 
            d = vabs(phi2*x1 + phi1*z1) - sq3;
            if (d > rmax) rmax = d; 
            d = vabs(phi1*x1 + phi2*y1) - sq3;
            if (d > rmax) rmax = d; 
            h = rmax;
            printf("h = %.3lg\n", h);
            if (h < 0.0) h = 0.0;
            rde->depth = SWATER_MIN_HEIGHT*DEPTH_FACTOR*h;
            break;
        }
        case (SH_SPHERE_ICOSAHEDRON):
        {
            rmax = 0.0;
            x1 = wsphere[i*NY+j].x;
            y1 = wsphere[i*NY+j].y;
            z1 = wsphere[i*NY+j].z;
            /* compute distance from dodecahedron */
            d = vabs(phi2*(x1+y1+z1)) - rico;
            if (d > rmax) rmax = d; 
            d = vabs(phi2*(-x1+y1+z1)) - rico;
            if (d > rmax) rmax = d; 
            d = vabs(phi2*(x1-y1+z1)) - rico;
            if (d > rmax) rmax = d; 
            d = vabs(phi2*(x1+y1-z1)) - rico;
            if (d > rmax) rmax = d; 
            d = vabs(phi3*x1 + phi1*z1) - rico;
            if (d > rmax) rmax = d; 
            d = vabs(phi3*z1 + phi1*y1) - rico;
            if (d > rmax) rmax = d; 
            d = vabs(phi3*y1 + phi1*x1) - rico;
            if (d > rmax) rmax = d; 
            d = vabs(phi3*x1 - phi1*z1) - rico;
            if (d > rmax) rmax = d; 
            d = vabs(phi3*z1 - phi1*y1) - rico;
            if (d > rmax) rmax = d; 
            d = vabs(phi3*y1 - phi1*x1) - rico;
            if (d > rmax) rmax = d; 
            h = rmax;
            printf("h = %.3lg\n", h);
            if (h < 0.0) h = 0.0;
            rde->depth = SWATER_MIN_HEIGHT*DEPTH_FACTOR*h;
            break;
        }
    }
}

double initialize_water_depth(t_rde rde[NX*NY], t_wave_sphere wsphere[NX*NY])
{
    int i, j, ncircles;
    double dx, dy, min, max, pscal, norm, vz = 0.01;
    
    if (SWATER_DEPTH == SH_CIRCLES) ncircles = init_circle_config_pattern(circles, CIRCLE_PATTERN);
    
    #pragma omp parallel for private(i,j)
    for (i=0; i<NX; i++){
        for (j=0; j<NY; j++){
            compute_water_depth(i, j, &rde[i*NY + j], wsphere, ncircles);
            if (rde[i*NY + j].depth < 0.0) rde[i*NY + j].depth = 0.0;
        }
    }
    
    dx = (XMAX - XMIN)/(double)NX;
    dy = (YMAX - YMIN)/(double)NY;
    
    /* compute x gradient */
    #pragma omp parallel for private(i,j)
    for (i=1; i<NX-1; i++){
        for (j=0; j<NY; j++){
            rde[i*NY + j].gradx = (rde[(i+1)*NY + j].depth - rde[(i-1)*NY + j].depth)/dx;
        }
    }
    
    /* left boundary */
    for (j=0; j<NY; j++){
        switch (B_COND_LEFT) {
            case (BC_PERIODIC): 
            {
                rde[j].gradx = (rde[NY + j].depth - rde[(NX-1)*NY + j].depth)/dx;
                break;
            }
            case (BC_DIRICHLET):
            {
                rde[j].gradx = (rde[NY + j].depth - rde[j].depth)*2.0/dx;
                break;
            }
        }
    }
    
    /* right boundary */
    for (j=0; j<NY; j++){
        switch (B_COND_RIGHT) {
            case (BC_PERIODIC): 
            {
                rde[(NX-1)*NY + j].gradx = (rde[j].depth - rde[(NX-2)*NY + j].depth)/dx;
                break;
            }
            case (BC_DIRICHLET):
            {
                rde[(NX-1)*NY + j].gradx = (rde[(NX-1)*NY + j].depth - rde[(NX-2)*NY + j].depth)*2.0/dx;
                break;
            }
        }
    }

    /* compute y gradient */
    #pragma omp parallel for private(i,j)
    for (i=0; i<NX; i++){
        for (j=1; j<NY-1; j++){
            rde[i*NY + j].grady = (rde[i*NY + j+1].depth - rde[i*NY + j-1].depth)/dy;
        }
    }
    
    /* bottom boundary */
    for (i=0; i<NX; i++){
        switch (B_COND_BOTTOM) {
            case (BC_PERIODIC): 
            {
                rde[i*NY].grady = (rde[i*NY + 1].depth - rde[i*NY + NY-1].depth)/dy;
                break;
            }
            case (BC_DIRICHLET):
            {
                rde[i*NY].grady = (rde[i*NY + 1].depth - rde[i*NY].depth)*2.0/dy;
                break;
            }
        }
    }
    
    /* top boundary */
    for (i=0; i<NX; i++){
        switch (B_COND_TOP) {
            case (BC_PERIODIC): 
            {
                rde[i*NY + NY-1].grady = (rde[i*NY].depth - rde[i*NY + NY-2].depth)/dy;
                break;
            }
            case (BC_DIRICHLET):
            {
                rde[i*NY + NY-1].grady = (rde[i*NY + NY-1].depth - rde[i*NY + NY-2].depth)*2.0/dy;
                break;
            }
        }
    }

    /* compute light angle */
    #pragma omp parallel for private(i,j)
    for (i=0; i<NX; i++){
        for (j=1; j<NY-1; j++){
            norm = sqrt(vz*vz + rde[i*NY + j].gradx*rde[i*NY + j].gradx + rde[i*NY + j].grady*rde[i*NY + j].grady);
            pscal = -rde[i*NY + j].gradx*light[0] - rde[i*NY + j].grady*light[1] + vz;
            rde[i*NY+j].cos_depth_angle = pscal/norm;
        }
    }
    
    /* for monitoring */
    min = 0.0;
    max = 0.0;
    #pragma omp parallel for private(i,j)
    for (i=0; i<NX; i++){
        for (j=0; j<NY; j++){
            if (rde[i*NY + j].gradx > max) max = rde[i*NY + j].gradx;
            if (rde[i*NY + j].gradx < min) min = rde[i*NY + j].gradx;
        }
    }
    printf("gradx min = %.3lg, max = %.3lg\n", min, max);

    min = 0.0;
    max = 0.0;
    #pragma omp parallel for private(i,j)
    for (i=0; i<NX; i++){
        for (j=0; j<NY; j++){
            if (rde[i*NY + j].grady > max) max = rde[i*NY + j].grady;
            if (rde[i*NY + j].grady < min) min = rde[i*NY + j].grady;
        }
    }
    printf("grady min = %.3lg, max = %.3lg\n", min, max);
    
    min = 0.0;
    max = 0.0;
    #pragma omp parallel for private(i,j)
    for (i=0; i<NX; i++){
        for (j=0; j<NY; j++){
            if (rde[i*NY + j].depth > max) max = rde[i*NY + j].depth;
            if (rde[i*NY + j].depth < min) min = rde[i*NY + j].depth;
        }
    }
    printf("Depth min = %.3lg, max = %.3lg\n", min, max);
    
    
//     for (j=0; j<NY; j++)
//     {
//         for (i=0; i<2; i++){
//             printf("point (%03i,%03i): depth %.03lg, gradient (%.03lg, %.03lg)\n", i, j, rde[i*NY + j].depth, rde[i*NY + j].gradx, rde[i*NY + j].grady);
//         }
//         for (i=NX-2; i<NX; i++){
//             printf("point (%03i,%03i): depth %.03lg, gradient (%.03lg, %.03lg)\n", i, j, rde[i*NY + j].depth, rde[i*NY + j].gradx, rde[i*NY + j].grady);
//         }
//     }
    
    return(max);
}

void evolve_wave_half(double *phi_in[NFIELDS], double *phi_out[NFIELDS], short int xy_in[NX*NY], double potential_field[NX*NY], double vector_potential_field[2*NX*NY], 
double gfield[2*NX*NY], t_rde rde[NX*NY], t_wave_sphere wsphere[NX*NY])
/* time step of field evolution */
{
    int i, j, k, iplus, iminus, jplus, jminus, ropening, w;
    double x, y, z, deltax, deltay, deltaz, rho, rhox, rhoy, pot, u, v, ux, uy, vx, vy, test = 0.0, dx, dy, xy[2], padding, a, eta, etax, etay, sum;
    double *delta_phi[NLAPLACIANS], *nabla_phi, *nabla_psi, *nabla_omega, *delta_vorticity, *delta_pressure, *delta_p, *delta_u, *delta_v, *nabla_rho, *nabla_u, *nabla_v, *nabla_eta;
//     double u_bc[NY], v_bc[NY]; 
    static double invsqr3 = 0.577350269;    /* 1/sqrt(3) */
    static int smooth = 0, y_channels, y_channels1, imin, imax, first = 1;
    
    if (first)  /* for D_MAZE_CHANNELS boundary conditions in Euler equation */
    {
        ropening = (NYMAZE+1)/2;
        padding = 0.02;
        dy = (YMAX - YMIN - 2.0*padding)/(double)(NYMAZE);
        y = YMIN + 0.02 + dy*((double)ropening);
        x = YMAX - padding + MAZE_XSHIFT;
        xy_to_pos(x, y, xy);
        if ((B_DOMAIN == D_MAZE_CHANNELS)||(OBSTACLE_GEOMETRY == D_MAZE_CHANNELS)||(OBSTACLE_GEOMETRY == D_MAZE_CHANNELS_INT))
        {
            imax = xy[0] + 2;
            x = YMIN + padding + MAZE_XSHIFT;
            xy_to_pos(x, y, xy);
            imin = xy[0] - 2;
            if (imin < 5) imin = 5;
        }
        else
        {
            imin = 0;
            imax = NX;
        }
        
        /* average values at poles on sphere */
        if (SPHERE) for (k=0; k<NFIELDS; k++)
        {
            sum = 0.0;
            for (i=0; i<NX; i++) for (j=0; j<DPOLE; j++) 
                sum += phi_in[k][i*NY+j];
            sum = sum/(double)(NX*DPOLE);
            if (sum > 1.0) sum = 1.0;
            if (sum < -1.0) sum = -1.0;
            for (i=0; i<NX; i++) for (j=0; j<DPOLE; j++)
                phi_in[k][i*NY+j] = sum;
        
            sum = 0.0;
            for (i=0; i<NX; i++) for (j=NY-DPOLE-1; j<NY; j++)
                sum += phi_in[k][i*NY+j];
            sum = sum/(double)(NX*DPOLE);
            if (sum > 1.0) sum = 1.0;
            if (sum < -1.0) sum = -1.0;
            for (i=0; i<NX; i++) for (j=NY-DPOLE-1; j<NY; j++)
                phi_in[k][i*NY+j] = sum;
        }
        first = 0;
    }
    
    /* smooth vector fields at poles */
    if ((SPHERE)&&(!SMOOTHBLOCKS)) for (k=0; k<NFIELDS; k++) smooth_poles(phi_in[k]);
    else for (k=0; k<NFIELDS; k++) block_poles(phi_in[k]);
        
    for (i=0; i<NLAPLACIANS; i++) delta_phi[i] = (double *)malloc(NX*NY*sizeof(double));
    
    if (COMPUTE_PRESSURE) 
    {
        delta_pressure = (double *)malloc(NX*NY*sizeof(double));
        delta_p = (double *)malloc(NX*NY*sizeof(double));
    }
    
    /* compute the Laplacian of phi */
    for (i=0; i<NLAPLACIANS; i++) compute_laplacian_rde(phi_in[i], delta_phi[i], xy_in, wsphere);
    
    if (COMPUTE_PRESSURE) compute_laplacian_rde(phi_in[2], delta_pressure, xy_in, wsphere);
    
    /* compute the gradient of phi if there is a magnetic field */
    if (ADD_MAGNETIC_FIELD) 
    {
        nabla_phi = (double *)malloc(2*NX*NY*sizeof(double));
        nabla_psi = (double *)malloc(2*NX*NY*sizeof(double));
        compute_gradient_xy(phi_in[0], nabla_phi);
        compute_gradient_xy(phi_in[1], nabla_psi);
    }
    
    switch (RDE_EQUATION){
        /* compute gradients of stream function and vorticity for Euler equation */
        case (E_EULER_INCOMP):
        {
            nabla_psi = (double *)malloc(2*NX*NY*sizeof(double));
            nabla_omega = (double *)malloc(2*NX*NY*sizeof(double));
            compute_gradient_euler(phi_in[0], nabla_psi, wsphere, EULER_GRADIENT_YSHIFT);
            compute_gradient_euler(phi_in[1], nabla_omega, wsphere, 0.0);
        
            if (COMPUTE_PRESSURE) compute_pressure_laplacian(phi_in, delta_p);
        
            dx = (XMAX-XMIN)/((double)NX);
            dy = (YMAX-YMIN)/((double)NY);
        
            if (SMOOTHEN_VORTICITY)     /* beta: try to reduce formation of ripples */
            {
                if (smooth == 0)
                {
                    delta_vorticity = (double *)malloc(NX*NY*sizeof(double));
                    compute_laplacian_rde(phi_in[1], delta_vorticity, xy_in, wsphere); 
//                 #pragma omp parallel for private(i,delta_vorticity)
                    for (i=0; i<NX*NY; i++) phi_in[1][i] += intstep*SMOOTH_FACTOR*delta_vorticity[i];
                    free(delta_vorticity);
                }
                smooth++;
                if (smooth >= SMOOTHEN_PERIOD) smooth = 0;
            }
            break;
        }
    
        /* compute gradients of fields for compressible Euler equation */
        case (E_EULER_COMP):
        {
            nabla_rho = (double *)malloc(2*NX*NY*sizeof(double));
            compute_gradient_euler_test(phi_in[0], nabla_rho, xy_in, wsphere);
            compute_velocity_gradients(phi_in, rde, xy_in, wsphere);
        
            if (SMOOTHEN_VELOCITY)     /* beta: try to reduce formation of ripples */
            {
                if (smooth == 0)
                {
                    delta_u = (double *)malloc(NX*NY*sizeof(double));
                    delta_v = (double *)malloc(NX*NY*sizeof(double));
                    compute_laplacian_rde(phi_in[1], delta_u, xy_in, wsphere); 
                    compute_laplacian_rde(phi_in[2], delta_v, xy_in, wsphere); 
                    #pragma omp parallel for private(i)
                    for (i=0; i<NX*NY; i++) phi_in[1][i] += intstep*SMOOTH_FACTOR*delta_u[i];
                    #pragma omp parallel for private(i)
                    for (i=0; i<NX*NY; i++) phi_in[2][i] += intstep*SMOOTH_FACTOR*delta_v[i];
                    free(delta_u);
                    free(delta_v);
                }
                smooth++;
                if (smooth >= SMOOTHEN_PERIOD) smooth = 0;
            }
            break;
        }
        
        case (E_SHALLOW_WATER):
        {
            nabla_eta = (double *)malloc(2*NX*NY*sizeof(double));
            compute_gradient_euler_test(phi_in[0], nabla_eta, xy_in, wsphere);
            compute_velocity_gradients(phi_in, rde, xy_in, wsphere);
            
            if (VISCOSITY > 0.0)     
            {
                delta_u = (double *)malloc(NX*NY*sizeof(double));
                delta_v = (double *)malloc(NX*NY*sizeof(double));
                compute_laplacian_rde(phi_in[1], delta_u, xy_in, wsphere); 
                compute_laplacian_rde(phi_in[2], delta_v, xy_in, wsphere); 
            }
            break;
        }
        
        default: 
        {
            /* do nothing */
        }
    }
    
    if (TEST_GRADIENT) {
        test = 0.0;
        for (i=0; i<2*NX*NY; i++){
            test += nabla_v[i]*nabla_v[i];
//             test += nabla_omega[i]*nabla_omega[i];
//             test += nabla_psi[i]*nabla_psi[i];
        }
        printf("nabla square = %.5lg\n", test/((double)NX*NY));
    }
    
    
    
    #pragma omp parallel for private(i,j,k,x,y,z,deltax,deltay,deltaz,rho)
    for (i=imin; i<imax; i++){
        for (j=0; j<NY; j++){
            if (xy_in[i*NY+j]) switch (RDE_EQUATION){
                case (E_HEAT):
                {
                    deltax = viscosity*delta_phi[0][i*NY+j];
                    phi_out[0][i*NY+j] = phi_in[0][i*NY+j] + intstep*deltax;
                    break;
                }
                case (E_ALLEN_CAHN):
                {
                    x = phi_in[0][i*NY+j];
                    deltax = viscosity*delta_phi[0][i*NY+j];
                    phi_out[0][i*NY+j] = phi_in[0][i*NY+j] + intstep*(deltax + K_AC*x*(1.0-x*x));
                    break;
                }
                case (E_CAHN_HILLIARD):
                {
                    /* TO DO */
                    break;
                }
                case (E_FHN):
                {
                    x = phi_in[0][i*NY+j];
                    y = phi_in[1][i*NY+j];
                    deltax = viscosity*delta_phi[0][i*NY+j];
                    phi_out[0][i*NY+j] = phi_in[0][i*NY+j] + intstep*(deltax + x*(1.0-x*x) + y);
                    phi_out[1][i*NY+j] = phi_in[0][i*NY+j] + intstep*EPSILON*(- invsqr3 - FHNC - FHNA*x);
                    break;
                }
                case (E_RPS):
                {
                    x = phi_in[0][i*NY+j];
                    y = phi_in[1][i*NY+j];
                    z = phi_in[2][i*NY+j];
                    rho = x + y + z;
                    deltax = viscosity*delta_phi[0][i*NY+j];
                    deltay = viscosity*delta_phi[1][i*NY+j];
                    deltaz = viscosity*delta_phi[2][i*NY+j];
                
                    phi_out[0][i*NY+j] = x + intstep*(deltax + x*(1.0 - rho - RPSA*y));
                    phi_out[1][i*NY+j] = y + intstep*(deltay + y*(1.0 - rho - RPSA*z));
                    phi_out[2][i*NY+j] = z + intstep*(deltaz + z*(1.0 - rho - RPSA*x));
                    break;
                }
                case (E_RPSLZ):
                {
                    rho = 0.0;
                    for (k=0; k<5; k++) rho += phi_in[k][i*NY+j];
                    
                    for (k=0; k<5; k++) 
                    {
                        x = phi_in[k][i*NY+j];
                        y = phi_in[(k+1)%5][i*NY+j];
                        z = phi_in[(k+3)%5][i*NY+j];
                        phi_out[k][i*NY+j] = x + intstep*(delta_phi[k][i*NY+j] + x*(1.0 - rho - RPSA*y - rpslzb*z));
                    }
                    break;
                }
                case (E_SCHRODINGER):
                {
                    phi_out[0][i*NY+j] = phi_in[0][i*NY+j] - intstep*delta_phi[1][i*NY+j];
                    phi_out[1][i*NY+j] = phi_in[1][i*NY+j] + intstep*delta_phi[0][i*NY+j];
                    if ((ADD_POTENTIAL)||(ADD_MAGNETIC_FIELD))
                    {
                        pot = potential_field[i*NY+j];
                        phi_out[0][i*NY+j] += intstep*pot*phi_in[1][i*NY+j];
                        phi_out[1][i*NY+j] -= intstep*pot*phi_in[0][i*NY+j];
                    }
                    if (ADD_MAGNETIC_FIELD)
                    {
                        vx = vector_potential_field[i*NY+j];
                        vy = vector_potential_field[NX*NY+i*NY+j];
                        phi_out[0][i*NY+j] -= 2.0*intstep*(vx*nabla_phi[i*NY+j] + vy*nabla_phi[NX*NY+i*NY+j]);
                        phi_out[1][i*NY+j] -= 2.0*intstep*(vx*nabla_psi[i*NY+j] + vy*nabla_psi[NX*NY+i*NY+j]);
                    }
                    break;
                }
                case (E_EULER_INCOMP):
                {
                    phi_out[0][i*NY+j] = phi_in[0][i*NY+j] + intstep*POISSON_STIFFNESS*(delta_phi[0][i*NY+j] + phi_in[1][i*NY+j]*dx*dx);
//                     phi_out[0][i*NY+j] += intstep*EULER_GRADIENT_YSHIFT;
                    phi_out[1][i*NY+j] = phi_in[1][i*NY+j] - intstep*K_EULER*(nabla_omega[i*NY+j]*nabla_psi[NX*NY+i*NY+j]);
                    phi_out[1][i*NY+j] += intstep*K_EULER*(nabla_omega[NX*NY+i*NY+j]*nabla_psi[i*NY+j]);
                        
                    if (COMPUTE_PRESSURE)
                    {
                        phi_out[2][i*NY+j] = phi_in[2][i*NY+j] + intstep*POISSON_STIFFNESS*(delta_pressure[i*NY+j] - delta_p[i*NY+j]);
                        phi_out[2][i*NY+j] *= exp(-2.0e-3);
                    }
                    break;
                }
                case (E_EULER_COMP):
                {
                    rho = phi_in[0][i*NY+j];
                    if (rho == 0.0) rho = 1.0e-1;
                    u = phi_in[1][i*NY+j];
                    v = phi_in[2][i*NY+j];
                    rhox = nabla_rho[i*NY+j];
                    rhoy = nabla_rho[NX*NY+i*NY+j];
                    
                    ux = rde[i*NY+j].dxu;
                    uy = rde[i*NY+j].dyu;
                    vx = rde[i*NY+j].dxv;
                    vy = rde[i*NY+j].dyv;
                    
                    phi_out[0][i*NY+j] = rho - intstep*C_EULER_COMP*(u*rhox + v*rhoy + rho*(ux + vy));
                    phi_out[1][i*NY+j] = u - intstep*(u*ux + v*uy + K_EULER_INC*rhox/rho);
                    phi_out[2][i*NY+j] = v - intstep*(u*vx + v*vy + K_EULER_INC*rhoy/rho);
                    
                    if (ADD_FORCE_FIELD)
                    {
                        phi_out[1][i*NY+j] += intstep*gfield[i*NY+j];
                        phi_out[2][i*NY+j] += intstep*gfield[NX*NY+i*NY+j];
                    }
                    if (ADD_CORIOLIS_FORCE)
                    {
                        if (SPHERE)
                        {
                            phi_out[1][i*NY+j] += intstep*G_FIELD*v*wsphere[i*NY+j].ctheta;
                            phi_out[2][i*NY+j] -= intstep*G_FIELD*u*wsphere[i*NY+j].reg_cottheta;                            
//                             phi_out[1][i*NY+j] += intstep*G_FIELD*v;
//                             phi_out[2][i*NY+j] -= intstep*G_FIELD*u;                            
//                             phi_out[1][i*NY+j] -= intstep*G_FIELD*v;
//                             phi_out[2][i*NY+j] += intstep*G_FIELD*u;                            
//                             phi_out[1][i*NY+j] -= intstep*G_FIELD*v*wsphere[i*NY+j].ctheta;
//                             phi_out[2][i*NY+j] += intstep*G_FIELD*u*wsphere[i*NY+j].ctheta;                            
                        }
                        else
                        {
                            phi_out[1][i*NY+j] += intstep*G_FIELD*v;
                            phi_out[2][i*NY+j] -= intstep*G_FIELD*u;
                        }
                    }
                    break;
                }
                case (E_SHALLOW_WATER):
                {
                    eta = phi_in[0][i*NY+j];
                    if (eta == 0.0) eta = 1.0e-6;
                    u = phi_in[1][i*NY+j];
                    v = phi_in[2][i*NY+j];
                    etax = nabla_eta[i*NY+j];
                    etay = nabla_eta[NX*NY+i*NY+j];
                    
                    ux = rde[i*NY+j].dxu;
                    uy = rde[i*NY+j].dyu;
                    vx = rde[i*NY+j].dxv;
                    vy = rde[i*NY+j].dyv;
                    
//                     printf("etax = %.3lg, etay = %.3lg, ux = %.3lg, uy = %.3lg, vx = %.3lg, vy = %.3lg\n", etax, etay, ux, uy, vx, vy);
                    
                    if (SPHERE)
                    {
                        /* TODO */
                        if ((j > DSMOOTH)&&(j < NY-DSMOOTH))
                        {
                            phi_out[0][i*NY+j] = eta - intstep*(u*etax + v*etay + eta*(ux + vy));
                            phi_out[1][i*NY+j] = u - intstep*(u*ux + v*uy + G_FIELD*etax);
                            phi_out[2][i*NY+j] = v - intstep*(u*vx + v*vy + G_FIELD*etax);
                        }   
                        else
                        {
                            phi_out[0][i*NY+j] = SWATER_MIN_HEIGHT;
                            phi_out[1][i*NY+j] = 0.0;
                            phi_out[2][i*NY+j] = 0.0;
                        }
                    }
                    else
                    {
                        phi_out[0][i*NY+j] = eta - intstep*(u*etax + v*etay + eta*(ux + vy));
                        phi_out[1][i*NY+j] = u - intstep*(u*ux + v*uy + G_FIELD*etax);
                        phi_out[2][i*NY+j] = v - intstep*(u*vx + v*vy + G_FIELD*etay);
                    }
                    
                    if (VISCOSITY > 0.0)
                    {
                        phi_out[1][i*NY+j] += intstep*VISCOSITY*delta_u[i*NY+j];
                        phi_out[2][i*NY+j] += intstep*VISCOSITY*delta_v[i*NY+j];
                    }
                    if (DISSIPATION > 0.0)
                    {
                        phi_out[1][i*NY+j] -= intstep*DISSIPATION*u;
                        phi_out[2][i*NY+j] -= intstep*DISSIPATION*v;
                    }
                    if (ADD_FORCE_FIELD)
                    {
                        phi_out[1][i*NY+j] += intstep*gfield[i*NY+j];
                        phi_out[2][i*NY+j] += intstep*gfield[NX*NY+i*NY+j];
                    }
                    if (VARIABLE_DEPTH)
                    {
                        phi_out[0][i*NY+j] -= intstep*rde[i*NY+j].depth*(ux + vy);
                        phi_out[0][i*NY+j] -= intstep*rde[i*NY+j].gradx*u;
                        phi_out[0][i*NY+j] -= intstep*rde[i*NY+j].grady*v;
                    }   
                    
                    break;
                }
            }
        }
    }
    
    /* in-flow/out-flow b.c. for incompressible Euler equation */
    if (((RDE_EQUATION == E_EULER_INCOMP)||(RDE_EQUATION == E_EULER_COMP))&&(IN_OUT_FLOW_BC > 0))
        set_in_out_flow_bc(phi_out, xy_in, flow_speed);
        
//     if (TEST_GRADIENT) {
//         test = 0.0;
//         for (i=0; i<NX*NY; i++){
//             test += delta_phi[0][i] + phi_out[1][i]*dx*dx;
//         }
//         printf("Delta psi + omega = %.5lg\n", test/((double)NX*NY));
//     }
                
    if (FLOOR) for (i=0; i<NX; i++){
        for (j=0; j<NY; j++){
            if (xy_in[i*NY+j] != 0) for (k=0; k<NFIELDS; k++)
            {
                if (phi_out[k][i*NY+j] > VMAX) phi_out[k][i*NY+j] = VMAX;
                if (phi_out[k][i*NY+j] < -VMAX) phi_out[k][i*NY+j] = -VMAX;
            }
        }
    }
    
    for (i=0; i<NLAPLACIANS; i++) free(delta_phi[i]);
    
    if (ADD_MAGNETIC_FIELD) 
    {
        free(nabla_phi);
        free(nabla_psi);
    }
    
    if (RDE_EQUATION == E_EULER_INCOMP) 
    {
        free(nabla_psi);
        free(nabla_omega);
    }
    else if (RDE_EQUATION == E_EULER_COMP)
    {
        free(nabla_rho);
    }
    else if (RDE_EQUATION == E_SHALLOW_WATER)
    {
        free(nabla_eta);
        if (VISCOSITY > 0.0)
        {
            free(delta_u);
            free(delta_v);
        }
    }
    
    if (COMPUTE_PRESSURE) 
    {
        free(delta_pressure);
        free(delta_p);
    }
}

void evolve_wave(double *phi[NFIELDS], double *phi_tmp[NFIELDS], short int xy_in[NX*NY], 
                 double potential_field[NX*NY], double vector_potential_field[2*NX*NY], 
                 double gfield[2*NX*NY], t_rde rde[NX*NY], t_wave_sphere wsphere[NX*NY])
/* time step of field evolution */
{
    evolve_wave_half(phi, phi_tmp, xy_in, potential_field, vector_potential_field, gfield, rde, wsphere);
    evolve_wave_half(phi_tmp, phi, xy_in, potential_field, vector_potential_field, gfield, rde, wsphere);
}

void update_tracer_table(double tracers[2*N_TRACERS*NSTEPS], t_rde rde[NX*NY], int time)
/* update tracer information in rde */
{
    int tracer, t, t1, maxtime, i, j, n, ij[2], length = 50, cell, oldcell; 
    double x, y; 
    
    #pragma omp parallel for private(cell)
    for (cell=0; cell<NX*NY; cell++) 
    {
        rde[cell].tracer = 0;
        rde[cell].prev_cell = cell;
        rde[cell].n_tracer_pts = 0;
    }
    
    maxtime = length;
    if (maxtime > time) maxtime = time;
    
    #pragma omp parallel for private(tracer)
    for (tracer = 0; tracer < N_TRACERS; tracer++)
    {
        for (t = 0; t < maxtime; t++)
        {
            t1 = time - t;
            
            x = tracers[t1*2*N_TRACERS + 2*tracer];
            y = tracers[t1*2*N_TRACERS + 2*tracer + 1];
            
            xy_to_ij(x, y, ij);
            cell = ij[0]*NY + ij[1];
            
            n = rde[cell].n_tracer_pts;
            if (n < NMAX_TRACER_PTS)
            {
                rde[cell].tracerx[n] = x;
                rde[cell].tracery[n] = y;
                rde[cell].n_tracer_pts++;
                rde[cell].tracer = length - t;
            }
//             else printf("More than %i tracer points per cell\n", NMAX_TRACER_PTS);
            
            if ((cell != oldcell)&&(t > 0))
                rde[cell].prev_cell = oldcell;
            
            oldcell = cell;
        }
    }
}

void evolve_tracers(double *phi[NFIELDS], double tracers[2*N_TRACERS*NSTEPS], t_rde rde[NX*NY], int time, int nsteps, double step)
/* time steps of tracer particle evolution (for Euler equation) */
{
    int tracer, i, j, n, t, ij[2], iplus, jplus, prev_cell, new_cell;
    double x, y, xy[2], vx, vy; 
    
    step = TRACERS_STEP;
    
    for (tracer = 0; tracer < N_TRACERS; tracer++)
    {
        x = tracers[time*2*N_TRACERS + 2*tracer];
        y = tracers[time*2*N_TRACERS + 2*tracer + 1];
        
//         printf("Tracer %i position (%.2f, %.2f)\n", tracer, x, y);
        
        for (t=0; t<nsteps; t++) 
        {
            xy_to_ij_safe(x, y, ij);
            i = ij[0];
            j = ij[1];
            
            switch (RDE_EQUATION) {
                case (E_EULER_INCOMP): 
                {
                    iplus = i + 1;  if (iplus == NX) iplus = 0;
                    jplus = j + 1;  if (jplus == NY) jplus = 0;
        
                    vx = phi[0][i*NY+jplus] - phi[0][i*NY+j];
                    vy = -(phi[0][iplus*NY+j] - phi[0][i*NY+j]);
            
                    if (j == 0) vx += EULER_GRADIENT_YSHIFT;
                    else if (j == NY-1) vx -= EULER_GRADIENT_YSHIFT;
                    break;
                }
                case (E_EULER_COMP):
                {
                    vx = phi[1][i*NY+j];
                    vy = phi[2][i*NY+j];
                    break;
                }
                case (E_SHALLOW_WATER):
                {
                    vx = phi[1][i*NY+j];
                    vy = phi[2][i*NY+j];
                    break;
                }
            }
            
//             v = module2(vx, vy);
//             if ((v > 0.0)&&(v < 0.1)) 
//             {
//                 vx = vx*0.1/v;
//                 vy = vy*0.1/v;
//             }
            
//             printf("(i, j) = (%i, %i), Tracer %i velocity (%.6f, %.6f)\n", i, j, tracer, vx, vy);
            
            x += vx*step;
            y += vy*step;
        }
//         printf("Tracer %i velocity (%.2f, %.2f)\n", tracer, vx, vy);
        
        if (x > XMAX) x += (XMIN - XMAX);
        if (x < XMIN) x += (XMAX - XMIN);
        if (y > YMAX) y += (YMIN - YMAX);
        if (y < YMIN) y += (YMAX - YMIN);
        
        if (time+1 < NSTEPS)
        {
            tracers[(time+1)*2*N_TRACERS + 2*tracer] = x;
            tracers[(time+1)*2*N_TRACERS + 2*tracer + 1] = y;
        }
    }
    
    if ((PLOT_3D)&&(time+1 < NSTEPS)) update_tracer_table(tracers, rde, time);
}


void print_level(int level)
{
    double pos[2];
    char message[50];
    
    glColor3f(1.0, 1.0, 1.0);
    sprintf(message, "Level %i", level);
    xy_to_pos(XMIN + 0.1, YMAX - 0.2, pos);
    write_text(pos[0], pos[1], message);
}



void print_parameters(double *phi[NFIELDS], t_rde rde[NX*NY], short int xy_in[NX*NY], double time, short int left, double viscosity, double noise)
{
    char message[100], message2[100];
    double density, hue, rgb[3], logratio, x, y, pos[2], probas[2], speed1, speed2;
    static double xbox, xtext, boxwidth, boxheight;
    static int first = 1;
    
    if (first)
    {
        if (WINWIDTH > 1280)
        {
            boxheight = 0.035;
            boxwidth = 0.21;
            if (left)
            {
                xbox = XMIN + 0.4;
                xtext = XMIN + 0.2;
            }
            else
            {
                xbox = XMAX - 0.49;
                xtext = XMAX - 0.65;
//                 xbox = XMAX - 0.39;
//                 xtext = XMAX - 0.55;
            }
        }
        else
        {
            boxwidth = 0.3;
            boxheight = 0.05;
            if (left)
            {
                xbox = XMIN + 0.4;
                xtext = XMIN + 0.1;
            }
            else
            {
                xbox = XMAX - 0.49;
                xtext = XMAX - 0.71;
//                 xbox = XMAX - 0.39;
//                 xtext = XMAX - 0.61;
            }
        }
         
        first = 0;
    }
    
    if (PRINT_PROBABILITIES)
    {
        compute_probabilities(rde, xy_in, probas);
        printf("pleft = %.3lg, pright = %.3lg\n", probas[0], probas[1]);
        
        x = XMIN + 0.15*(XMAX - XMIN);
        y = YMIN + 0.3*(YMAX - YMIN);
        erase_area_hsl(x, y, boxwidth, boxheight, 0.0, 0.9, 0.0);
        glColor3f(1.0, 1.0, 1.0);
        sprintf(message, "Proba %.3f", probas[0]);
        write_text(x, y, message);
        
        x = XMIN + 0.72*(XMAX - XMIN);
        y = YMIN + 0.68*(YMAX - YMIN);
        erase_area_hsl(x, y, boxwidth, boxheight, 0.0, 0.9, 0.0);
        glColor3f(1.0, 1.0, 1.0);
        sprintf(message, "Proba %.3f", probas[1]);
        write_text(x, y, message);
    }
    else
    {
        y = YMAX - 0.1;
        erase_area_hsl(xbox, y + 0.02, boxwidth, boxheight, 0.0, 0.9, 0.0);
        glColor3f(1.0, 1.0, 1.0);
        if (PRINT_TIME) sprintf(message, "Time %.3f", time);
        else if (PRINT_VISCOSITY) sprintf(message, "Viscosity %.3f", viscosity);
        else if (PRINT_RPSLZB) sprintf(message, "b = %.3f", rpslzb);
        else if (PRINT_NOISE) sprintf(message, "noise %.3f", noise);
        else if (PRINT_FLOW_SPEED) sprintf(message, "Speed %.3f", flow_speed);
        else if (PRINT_AVERAGE_SPEED) 
        {
            compute_average_speeds(phi, rde, &speed1, &speed2);
            sprintf(message, "Average vx %.3f", speed2);
            sprintf(message2, "Average vx %.3f", speed1);
        }
        if (PLOT_3D) write_text(xtext, y, message);
        else
        {
            xy_to_pos(xtext, y, pos);
            write_text(pos[0], pos[1], message);
            if (PRINT_AVERAGE_SPEED)
            {
                y = YMIN + 0.1;
                erase_area_hsl(xbox, y + 0.02, boxwidth, boxheight, 0.0, 0.9, 0.0);
                glColor3f(1.0, 1.0, 1.0);
                xy_to_pos(xtext, y, pos);
                write_text(pos[0], pos[1], message2);
            }
        }
    }
}

void draw_color_bar_palette(int plot, double range, int palette, int circular, int fade, double fade_value)
{
    double width = 0.14;
//     double width = 0.2;
    
    if (PLOT_3D)
    {
        if (ROTATE_COLOR_SCHEME) 
            draw_color_scheme_palette_3d(XMIN + 0.3, YMIN + 0.1, XMAX - 0.3, YMIN + 0.1 + width, plot, -range, range, palette, fade, fade_value);
        else if (circular)
            draw_circular_color_scheme_palette_3d(XMAX - 2.0*width, YMAX - 2.0*width, 1.0*width, plot, -range, range, palette, fade, fade_value);
        else 
            draw_color_scheme_palette_3d(XMAX - 1.5*width, YMIN + 0.1, XMAX - 0.5*width, YMAX - 0.1, plot, -range, range, palette, fade, fade_value);
    }
    else
    {
        if (circular)
            draw_circular_color_scheme_palette_fade(XMAX - 2.0*width, YMAX - 2.0*width, 1.0*width, plot, -range, range, palette, fade, fade_value);
        else if (ROTATE_COLOR_SCHEME) 
            draw_color_scheme_palette_fade(XMIN + 0.8, YMIN + 0.05, XMAX - 0.8, YMIN + 0.05 + width, plot, -range, range, palette, fade, fade_value);
        else 
            draw_color_scheme_palette_fade(XMAX - 1.5*width, YMIN + 0.1, XMAX - 0.5*width, YMAX - 0.1, plot, -range, range, palette, fade, fade_value);
    }
}

double noise_schedule(int i)
{
    double ratio;
    
    if (i < NOISE_INITIAL_TIME) return (NOISE_INTENSITY);
    else 
    {
        ratio = (double)(i - NOISE_INITIAL_TIME)/(double)(NSTEPS - NOISE_INITIAL_TIME);
        return (NOISE_INTENSITY*(1.0 + ratio*(NOISE_FACTOR - 1.0)));
    }
}


double viscosity_schedule(int i)
{
    double ratio;
    
    if (i < VISCOSITY_INITIAL_TIME) return (VISCOSITY);
    else 
    {
        ratio = (double)(i - VISCOSITY_INITIAL_TIME)/(double)(NSTEPS - VISCOSITY_INITIAL_TIME);
        return (VISCOSITY*(1.0 + ratio*(VISCOSITY_FACTOR - 1.0)));
    }
}

double rpslzb_schedule(int i)
{
    double ratio;
    
    if (i < RPSLZB_INITIAL_TIME) return (RPSLZB);
    else if (i > NSTEPS - RPSLZB_FINAL_TIME) return(RPSLZB - RPSLZB_CHANGE);
    else 
    {
        ratio = (double)(i - RPSLZB_INITIAL_TIME)/(double)(NSTEPS - RPSLZB_INITIAL_TIME - RPSLZB_FINAL_TIME);
        return (RPSLZB - ratio*RPSLZB_CHANGE);
    }
}

double flow_speed_schedule(int i)
{
    double ratio;
    
    ratio = (double)i/(double)NSTEPS;
    return (IN_OUT_FLOW_MIN_AMP + (IN_OUT_FLOW_AMP - IN_OUT_FLOW_MIN_AMP)*ratio);
}


void viewpoint_schedule(int i)
/* change position of observer */
{
    int j;
    double angle, ca, sa, r1, interpolate, rho;
    static double observer_initial[3], r, ratio, rho0, zmax;
    static int first = 1;
    
    if (first)
    {
        for (j=0; j<3; j++) observer_initial[j] = observer[j];
        r1 = observer[0]*observer[0] + observer[1]*observer[1];
        r = sqrt(r1 + observer[2]*observer[2]);
        ratio = r/sqrt(r1);
        rho0 = module2(observer[0], observer[1]);
        if (vabs(rho0) < 0.001) rho0 = 0.001; 
        zmax = r*sin(MAX_LATITUDE*PI/180.0);
        first = 0;
    }
    
    interpolate = (double)i/(double)NSTEPS;
    angle = (ROTATE_ANGLE*DPI/360.0)*interpolate;
//     printf("i = %i, interpolate = %.3lg, angle = %.3lg\n", i, interpolate, angle);
    ca = cos(angle);
    sa = sin(angle);
    switch (VIEWPOINT_TRAJ)
    {
        case (VP_HORIZONTAL):
        {
            observer[0] = ca*observer_initial[0] - sa*observer_initial[1];
            observer[1] = sa*observer_initial[0] + ca*observer_initial[1];
            break;
        }
        case (VP_ORBIT):
        {
            observer[0] = ca*observer_initial[0] - sa*observer_initial[1]*ratio;
            observer[1] = ca*observer_initial[1] + sa*observer_initial[0]*ratio;
            observer[2] = ca*observer_initial[2];
            break;
        }
        case (VP_ORBIT2):
        {
            observer[0] = ca*observer_initial[0] - sa*observer_initial[1]*ratio;
            observer[1] = ca*observer_initial[1] + sa*observer_initial[0]*ratio;
            observer[2] = sa*zmax;
            break;
        }
        case (VP_POLAR):
        {
            rho = -sa*observer_initial[2] + ca*rho0;
            observer[0] = observer_initial[0]*rho/rho0;
            observer[1] = observer_initial[1]*rho/rho0;
            observer[2] = ca*observer_initial[2] + sa*rho0;
            break; 
        }
    }
    
    printf("Angle %.3lg, Observer position (%.3lg, %.3lg, %.3lg)\n", angle, observer[0], observer[1], observer[2]);
}


void animation()
{
    double time = 0.0, scale, dx, var, jangle, cosj, sinj, sqrintstep, phishift, thetashift, amp, 
        intstep0, viscosity_printed, fade_value, noise = NOISE_INTENSITY, x, y, sign, phase;
    double *phi[NFIELDS], *phi_tmp[NFIELDS], *potential_field, *vector_potential_field, *tracers, *gfield, *bc_field, *bc_field2;
    short int *xy_in;
    int i, j, k, s, nvid, field;
    static int counter = 0;
    t_rde *rde;
    t_wave_sphere *wsphere, *wsphere_hr;

    /* Since NX and NY are big, it seemed wiser to use some memory allocation here */
    for (i=0; i<NFIELDS; i++)
    {
        phi[i] = (double *)malloc(NX*NY*sizeof(double));
        phi_tmp[i] = (double *)malloc(NX*NY*sizeof(double));
    }

    xy_in = (short int *)malloc(NX*NY*sizeof(short int));
    rde = (t_rde *)malloc(NX*NY*sizeof(t_rde));
    
    if (SPHERE) 
    {
        wsphere = (t_wave_sphere *)malloc(NX*NY*sizeof(t_wave_sphere));
        init_wave_sphere_rde(wsphere,1);
        /* high resolution version for planet simulations */
        wsphere_hr = (t_wave_sphere *)malloc(HRES*HRES*NX*NY*sizeof(t_wave_sphere));
        init_wave_sphere_rde(wsphere_hr,HRES);
    }
    
    npolyline = init_polyline(MDEPTH, polyline);
    for (i=0; i<npolyline; i++) printf("vertex %i: (%.3f, %.3f)\n", i, polyline[i].x, polyline[i].y);
    
    npolyrect = init_polyrect(polyrect);
    for (i=0; i<npolyrect; i++) printf("polyrect vertex %i: (%.3f, %.3f) - (%.3f, %.3f)\n", i, polyrect[i].x1, polyrect[i].y1, polyrect[i].x2, polyrect[i].y2);

    if (ADD_POTENTIAL) 
    {
        potential_field = (double *)malloc(NX*NY*sizeof(double));
        initialize_potential(potential_field);
    }
    else if (ADD_MAGNETIC_FIELD)
    {
        potential_field = (double *)malloc(NX*NY*sizeof(double));
        vector_potential_field = (double *)malloc(2*NX*NY*sizeof(double));
        initialize_potential(potential_field);
        initialize_vector_potential(vector_potential_field);
    }
    
    if (VARIABLE_DEPTH) 
    {
//         water_depth = (t_swater_depth *)malloc(NX*NY*sizeof(t_swater_depth));
        max_depth = initialize_water_depth(rde, wsphere);
        printf("Max depth = %.3lg\n", max_depth);
    }
    
    if (ADAPT_STATE_TO_BC)
    {
        bc_field = (double *)malloc(NX*NY*sizeof(double));
        bc_field2 = (double *)malloc(NX*NY*sizeof(double));
        
        initialize_bcfield(bc_field, bc_field2, polyrect, wsphere);
    }
    if (ADD_FORCE_FIELD)
    {
        if (!ADAPT_STATE_TO_BC)
        {
            printf("Error: if ADD_FORCE_FIELD = 1, then ADAPT_STATE_TO_BC should be 1\n");
            exit(1);
        }
        gfield = (double *)malloc(2*NX*NY*sizeof(double));
        initialize_gfield(gfield, bc_field, bc_field2, wsphere, wsphere_hr);
    }
        
    
//     if (ADD_TRACERS) tracers = (double *)malloc(2*NSTEPS*N_TRACERS*sizeof(double));
    if (ADD_TRACERS) tracers = (double *)malloc(4*NSTEPS*N_TRACERS*sizeof(double));

    dx = (XMAX-XMIN)/((double)NX);
    intstep = DT/(dx*dx);
    
    intstep0 = intstep;
    intstep1 = DT/dx;
    
    viscosity = VISCOSITY;
    
    sqrintstep = sqrt(intstep*(double)NVID);
        
    printf("Integration step %.3lg\n", intstep);

    /* initialize field */
//     init_random(0.5, 0.4, phi, xy_in);
//     init_random(0.5, 0.25, phi, xy_in, wsphere);
//     init_random_smoothed(1.0, 0.1, phi, xy_in, wsphere);
    
//     init_gaussian(x, y, mean, amplitude, scalex, phi, xy_in)
//     init_coherent_state(0.0, 0.0, 10.0, 0.0, 0.1, phi, xy_in);
    
//     init_coherent_state_sphere(0, 0.0, PID, 10.0, 5.0, 0.1, phi, xy_in, wsphere);
//     init_coherent_state_sphere(1, PI, PID, -10.0, 5.0, 0.1, phi, xy_in, wsphere);
//     init_coherent_state_sphere(1, PID, PID, 10.0, -5.0, 0.1, phi, xy_in, wsphere);
//     init_coherent_state_sphere(1, 3.0*PID, PID, -10.0, -5.0, 0.1, phi, xy_in, wsphere);
    
//     add_coherent_state(-0.75, -0.75, 0.0, 5.0, 0.1, phi, xy_in);
//     init_fermion_state(-0.5, 0.5, 2.0, 0.0, 0.1, phi, xy_in);
//     init_boson_state(-0.5, 0.5, 2.0, 0.0, 0.1, phi, xy_in);

//     init_laminar_flow(0.05, LAMINAR_FLOW_MODULATION, 0.015, 0.1, 1.0, 0.0, 0.1, phi, xy_in);
//     for (i=0; i<5; i++)
//         for (j=0; j<3; j++)
//         {
//             x = XMIN + ((double)i - 0.5)*(XMAX-XMIN)/4.0 + 0.15*gaussian();
//             y = YMIN + ((double)j - 0.5)*(YMAX-YMIN)/2.0 + 0.15*gaussian();
//             sign = (double)((i+j)%2)*2.0 - 1.0;
//             add_vortex_state(0.1*sign, x, y, 0.1, -0.35*sign + 0.1*gaussian(), phi, xy_in);
//         }
//     add_vortex_state(0.2, 0.75, 0.1, 0.3, -0.5, phi, xy_in);
//     add_vortex_state(-0.35, -0.75, -0.1, 0.4, 0.5, phi, xy_in);
//     add_vortex_state(0.1, -0.3, 0.7, 0.1, -0.5, phi, xy_in);
    
//     init_vortex_state(0.1, 0.4, 0.0, 0.3, -0.1, phi, xy_in);
//     add_vortex_state(0.1, -0.4, 0.0, 0.3, 0.1, phi, xy_in);

//     init_laminar_flow(double amp, double xmodulation, double ymodulation, double xperiod, double yperiod, double yshift, double density_mod, double *phi[NFIELDS], short int xy_in[NX*NY])
    
//     init_laminar_flow_earth(0.01, phi, xy_in);
//     
//     add_random_onefield_smoothed(0, 0.0, 0.03, phi, xy_in, wsphere);
//     add_random_onefield_smoothed(1, 0.0, 0.07, phi, xy_in, wsphere);
//     add_random_onefield_smoothed(2, 0.0, 0.07, phi, xy_in, wsphere);
    

//     init_pressure_gradient_flow(flow_speed_schedule(0), 1.0 + PRESSURE_GRADIENT, 1.0 - PRESSURE_GRADIENT, phi, xy_in, bc_field);

//     init_shear_flow_sphere(1.0, 0.25, 0.25, 6, 6, 0.75, phi, xy_in, wsphere);
    
//     phishift = 0.0;
//     thetashift = 0.2*PID;
//     amp = 1.0;
//     init_vortex_state_sphere(0, amp, phishift, PID + thetashift, 0.15, 0.3, phi, xy_in, wsphere);
    
    
//     init_gaussian_wave(0.7, 0.0, 0.0225, 0.1, SWATER_MIN_HEIGHT, phi, xy_in);
//     init_gaussian_wave_sphere(4.0, PID, 0.05, 0.05, SWATER_MIN_HEIGHT, phi, xy_in, wsphere);

    //     init_linear_wave(-1.0, 0.01, 5.0e-8, 0.25, SWATER_MIN_HEIGHT, phi, xy_in);
    
    init_swater_laminar_flow_sphere(0, -0.75, 0.0075, 6, 1, 0.0025, SWATER_MIN_HEIGHT, phi, xy_in, wsphere);
        
//     add_gaussian_wave(-1.6, -0.5, 0.015, 0.25, SWATER_MIN_HEIGHT, phi, xy_in);
    
    if (SMOOTHBLOCKS) for (k=0; k<NFIELDS; k++) block_poles(phi[k]);
    if (ADAPT_STATE_TO_BC) adapt_state_to_bc(phi, bc_field, xy_in);
    
    init_cfield_rde(phi, xy_in, CPLOT, rde, 0);
    if ((PLOT_3D)||(SHADE_2D)) init_zfield_rde(phi, xy_in, ZPLOT, rde, 0);

    if (DOUBLE_MOVIE)
    {
        init_cfield_rde(phi, xy_in, CPLOT_B, rde, 1);
        if  ((PLOT_3D)||(SHADE_2D)) init_zfield_rde(phi, xy_in, ZPLOT_B, rde, 1);
    }
    
    if (ADD_TRACERS) for (i=0; i<N_TRACERS; i++)
    {
        tracers[2*i] = XMIN + 0.05 + (XMAX - XMIN - 0.1)*rand()/RAND_MAX;
        tracers[2*i+1] = YMIN + 0.05 + (YMAX - YMIN - 0.1)*rand()/RAND_MAX;
    }
    
    blank();
    glColor3f(0.0, 0.0, 0.0);
    

    glutSwapBuffers();
    
    printf("Drawing wave\n");
    
    draw_wave_rde(0, phi, xy_in, rde, wsphere, wsphere_hr, potential_field, ZPLOT, CPLOT, COLOR_PALETTE, 0, 1.0, 1);
//     draw_billiard();
    if (PRINT_PARAMETERS) print_parameters(phi, rde, xy_in, time, PRINT_LEFT, VISCOSITY, noise);
    if (DRAW_COLOR_SCHEME) draw_color_bar_palette(CPLOT, COLORBAR_RANGE, COLOR_PALETTE, CIRC_COLORBAR, 0, 1.0);

    glutSwapBuffers();

    sleep(SLEEP1);
//     printf("Saving frame %i\n", i);
    if (MOVIE) for (i=0; i<INITIAL_TIME; i++) save_frame();

    for (i=0; i<=NSTEPS; i++)
    {
        nvid = NVID;
        if (CHANGE_VISCOSITY) 
        {
            viscosity = viscosity_schedule(i);
            viscosity_printed = viscosity;
            printf("Viscosity = %.3lg\n", viscosity); 
            if ((ADJUST_INTSTEP)&&(viscosity > VISCOSITY_MAX))
            {
                nvid = (int)((double)NVID*viscosity/VISCOSITY_MAX);
//                 viscosity = VISCOSITY_MAX;
                intstep = intstep0*VISCOSITY_MAX/viscosity;
                printf("Nvid = %i, intstep = %.3lg\n", nvid, intstep);
            }
        }
        if (CHANGE_RPSLZB) rpslzb = rpslzb_schedule(i);
        if (CHANGE_FLOW_SPEED) flow_speed = flow_speed_schedule(i); 
        else flow_speed = IN_OUT_FLOW_AMP;
        
        if (ROTATE_VIEW) 
        {
            viewpoint_schedule(i - INITIAL_TIME);
            reset_view = 1;
        }
        
        printf("Drawing wave %i\n", i);
        draw_wave_rde(0, phi, xy_in, rde, wsphere, wsphere_hr, potential_field, ZPLOT, CPLOT, COLOR_PALETTE, 0, 1.0, 1);
        
//         nvid = (int)((double)NVID*(1.0 + (ACCELERATION_FACTOR - 1.0)*(double)i/(double)NSTEPS));
        /* increase integration step */
//         intstep = intstep0*exp(log(DT_ACCELERATION_FACTOR)*(double)i/(double)NSTEPS);
//         if (intstep > MAX_DT)
//         {
//             nvid *= intstep/MAX_DT;
//             intstep = MAX_DT;
//         }
//         printf("Steps per frame: %i\n", nvid);
//         printf("Integration step %.5lg\n", intstep);
        
        printf("Evolving wave\n");
        for (j=0; j<nvid; j++) evolve_wave(phi, phi_tmp, xy_in, potential_field, vector_potential_field, gfield, rde, wsphere);

        if (ADAPT_STATE_TO_BC) adapt_state_to_bc(phi, bc_field, xy_in);
        
        if ((ADD_OSCILLATING_SOURCE)&&(i%OSCILLATING_SOURCE_PERIOD == 0))
        {
            phase = (double)i*OSCILLATING_SOURCE_OMEGA;
            printf("Adding vibration with phase %.3lg\n", phase); 
            add_elliptical_vibration(0.0, 0.0, 0.00007*cos(phase), LAMBDA, 1.0, 0.00004*cos(phase), 0.1, 1.0, SWATER_MIN_HEIGHT, phi, xy_in);
        }
                
        if (ADD_TRACERS)
        {
            printf("Evolving tracer particles\n");
            evolve_tracers(phi, tracers, rde, i, 10, 0.1);
//             for (j=0; j<N_TRACERS; j++) 
//                 printf("Tracer %i position (%.2f, %.2f)\n", j, tracers[2*N_TRACERS*i + 2*j], tracers[2*N_TRACERS*i + 2*j + 1]);
            if (!PLOT_3D) 
            {
                printf("Drawing tracer particles\n");
                draw_tracers(phi, tracers, i, 0, 1.0);
            }
        }
        
        if (ANTISYMMETRIZE_WAVE_FCT) antisymmetrize_wave_function(phi, xy_in);
        
        for (j=0; j<NFIELDS; j++) printf("field[%i] = %.3lg\t", j, phi[j][NY/2]);
//         for (j=0; j<NFIELDS; j++) printf("field[%i] = %.3lg\t", j, phi[j][NX*NY/2]);
//         for (j=0; j<NFIELDS; j++) printf("field[%i] = %.3lg\t", j, phi[j][0]);
        printf("\n");
        
        if (ADD_NOISE == 1) 
        {
//             #pragma omp parallel for private(field,j,k)
            for (field=0; field<NFIELDS; field++)
                for (j=0; j<NX; j++) 
                    for (k=0; k<NY; k++)
                        phi[field][j*NY+k] += sqrintstep*NOISE_INTENSITY*gaussian();
        }
        else if (ADD_NOISE == 2) 
        {
            if (CHANGE_NOISE)
            {
                noise = noise_schedule(i);
//                 #pragma omp parallel for private(field,j,k)
                for (field=0; field<NFIELDS; field++)
                    for (j=NX/2; j<NX; j++) 
                        for (k=0; k<NY; k++)
                            phi[field][j*NY+k] += sqrintstep*noise*gaussian();
            }
            else
            {
//                 #pragma omp parallel for private(field,j,k)
                for (field=0; field<NFIELDS; field++)
                    for (j=NX/2; j<NX; j++) 
                        for (k=0; k<NY; k++)
                            phi[field][j*NY+k] += sqrintstep*NOISE_INTENSITY*gaussian();
            }
        }
        time += nvid*intstep;
        
//         draw_billiard();
        if (PRINT_PARAMETERS) print_parameters(phi, rde, xy_in, time, PRINT_LEFT, viscosity_printed, noise);
        if (DRAW_COLOR_SCHEME) draw_color_bar_palette(CPLOT, COLORBAR_RANGE, COLOR_PALETTE, CIRC_COLORBAR, 0, 1.0); 
        
//         print_level(MDEPTH);
//         print_Julia_parameters(i);

        if (!((NO_EXTRA_BUFFER_SWAP)&&(MOVIE))) glutSwapBuffers();
        
//         glutSwapBuffers();
//         save_frame();
        
        /* modify Julia set */
//         set_Julia_parameters(i, phi, xy_in);

	if (MOVIE)
//         if (0)
        {
            printf("Saving frame %i\n", i);
//             if (NO_EXTRA_BUFFER_SWAP) glutSwapBuffers();
            save_frame();
            
            if ((i >= INITIAL_TIME)&&(DOUBLE_MOVIE))
            {
                draw_wave_rde(1, phi, xy_in, rde, wsphere, wsphere_hr, potential_field, ZPLOT_B, CPLOT_B, COLOR_PALETTE_B, 0, 1.0, REFRESH_B);
                if ((ADD_TRACERS)&&(!PLOT_3D)) draw_tracers(phi, tracers, i, 0, 1.0);
//                 draw_billiard();
                if (PRINT_PARAMETERS) print_parameters(phi, rde, xy_in, time, PRINT_LEFT, viscosity_printed, noise);
                if (DRAW_COLOR_SCHEME) draw_color_bar_palette(CPLOT_B, COLORBAR_RANGE_B, COLOR_PALETTE_B, CIRC_COLORBAR_B, 0, 1.0);  
                glutSwapBuffers();
//                 if (NO_EXTRA_BUFFER_SWAP) glutSwapBuffers();
                save_frame_counter(NSTEPS + MID_FRAMES + 1 + counter);
                counter++;
            }
            else if (NO_EXTRA_BUFFER_SWAP) glutSwapBuffers();
            
            /* TEST */
//              if (ADAPT_STATE_TO_BC) adapt_state_to_bc(phi, bc_field, xy_in);

            /* it seems that saving too many files too fast can cause trouble with the file system */
            /* so this is to make a pause from time to time - parameter PAUSE may need adjusting   */
            if (i % PAUSE == PAUSE - 1)
            {
                printf("Making a short pause\n");
                sleep(PSLEEP);
                s = system("mv wave*.tif tif_bz/");
            }
        }
        else printf("Computing frame %i\n", i);

    }
    
    if (MOVIE) 
    {
        if (DOUBLE_MOVIE) 
        {
            draw_wave_rde(0, phi, xy_in, rde, wsphere, wsphere_hr, potential_field, ZPLOT, CPLOT, COLOR_PALETTE, 0, 1.0, 1);
            if ((ADD_TRACERS)&&(!PLOT_3D)) draw_tracers(phi, tracers, NSTEPS, 0, 1.0);
//             draw_billiard();
            if (PRINT_PARAMETERS) print_parameters(phi, rde, xy_in, time, PRINT_LEFT, viscosity_printed, noise);
            if (DRAW_COLOR_SCHEME) draw_color_bar_palette(CPLOT, COLORBAR_RANGE, COLOR_PALETTE, CIRC_COLORBAR, 0, 1.0); 
//             if (!NO_EXTRA_BUFFER_SWAP) glutSwapBuffers();
            glutSwapBuffers();
            
            if (!FADE) for (i=0; i<MID_FRAMES; i++) save_frame();
            else for (i=0; i<MID_FRAMES; i++) 
            {
                fade_value = 1.0 - (double)i/(double)MID_FRAMES;
                draw_wave_rde(0, phi, xy_in, rde, wsphere, wsphere_hr, potential_field, ZPLOT, CPLOT, COLOR_PALETTE, 1, fade_value, 0);
                if ((ADD_TRACERS)&&(!PLOT_3D)) draw_tracers(phi, tracers, NSTEPS, 1, fade_value);
//                 draw_billiard();
                if (PRINT_PARAMETERS) print_parameters(phi, rde, xy_in, time, PRINT_LEFT, viscosity_printed, noise);
                if (DRAW_COLOR_SCHEME) draw_color_bar_palette(CPLOT, COLORBAR_RANGE, COLOR_PALETTE, CIRC_COLORBAR, 1, fade_value);   
                if (!NO_EXTRA_BUFFER_SWAP) glutSwapBuffers();
                save_frame_counter(NSTEPS + i + 1);
            }
            draw_wave_rde(1, phi, xy_in, rde, wsphere, wsphere_hr, potential_field, ZPLOT_B, CPLOT_B, COLOR_PALETTE_B, 0, 1.0, REFRESH_B);
            if ((ADD_TRACERS)&&(!PLOT_3D)) draw_tracers(phi, tracers, NSTEPS, 0, 1.0);
            if (PRINT_PARAMETERS) print_parameters(phi, rde, xy_in, time, PRINT_LEFT, viscosity_printed, noise);
            if (DRAW_COLOR_SCHEME) draw_color_bar_palette(CPLOT_B, COLORBAR_RANGE_B, COLOR_PALETTE_B, CIRC_COLORBAR_B, 0, 1.0); 
            glutSwapBuffers();
            
            if (!FADE) for (i=0; i<END_FRAMES; i++) save_frame_counter(NSTEPS + MID_FRAMES + 1 + counter + i);
            else for (i=0; i<END_FRAMES; i++) 
            {
                fade_value = 1.0 - (double)i/(double)END_FRAMES;
                draw_wave_rde(1, phi, xy_in, rde, wsphere, wsphere_hr, potential_field, ZPLOT_B, CPLOT_B, COLOR_PALETTE_B, 1, fade_value, 0);
                if ((ADD_TRACERS)&&(!PLOT_3D)) draw_tracers(phi, tracers, NSTEPS, 1, fade_value);
                if (PRINT_PARAMETERS) print_parameters(phi, rde, xy_in, time, PRINT_LEFT, viscosity_printed, noise);
                if (DRAW_COLOR_SCHEME) draw_color_bar_palette(CPLOT_B, COLORBAR_RANGE_B, COLOR_PALETTE_B, CIRC_COLORBAR_B, 1, fade_value);   
                glutSwapBuffers();
                save_frame_counter(NSTEPS + MID_FRAMES + 1 + counter + i);
            }
        }
        else
        {
            if (!FADE) for (i=0; i<END_FRAMES; i++) save_frame_counter(NSTEPS + MID_FRAMES + 1 + counter + i);
            else for (i=0; i<END_FRAMES; i++) 
            {
                fade_value = 1.0 - (double)i/(double)END_FRAMES;
                draw_wave_rde(0, phi, xy_in, rde, wsphere, wsphere_hr, potential_field, ZPLOT, CPLOT, COLOR_PALETTE, 1, fade_value, 0);
                if ((ADD_TRACERS)&&(!PLOT_3D)) draw_tracers(phi, tracers, NSTEPS, 1, fade_value);
                if (DRAW_COLOR_SCHEME) draw_color_bar_palette(CPLOT, COLORBAR_RANGE, COLOR_PALETTE, CIRC_COLORBAR, 1, fade_value); 
                glutSwapBuffers();
                save_frame_counter(NSTEPS + 1 + counter + i);
            }
        }
        
        s = system("mv wave*.tif tif_bz/");
    }

    for (i=0; i<NFIELDS; i++)
    {
        free(phi[i]);
        free(phi_tmp[i]);
    }
    free(xy_in);
    if (SPHERE) 
    {
        free(wsphere);
        free(wsphere_hr);
    }
    if (ADD_POTENTIAL) free(potential_field);
    else if (ADD_MAGNETIC_FIELD) 
    {
        free(potential_field);
        free(vector_potential_field);
    }
//     if (VARIABLE_DEPTH) free(water_depth);
    if (ADD_TRACERS) free(tracers);
    if (ADD_FORCE_FIELD) free(gfield);
    if (ADAPT_STATE_TO_BC) 
    {
        free(bc_field);
        free(bc_field2);
    }
    
    printf("Time %.5lg\n", time);

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
    glutCreateWindow("FitzHugh-Nagumo equation in a planar domain");

    if (PLOT_3D) init_3d(); 
    else init_hres(HRES);
//     else init();

    glutDisplayFunc(display);

    glutMainLoop();

    return 0;
}

