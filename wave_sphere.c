/*********************************************************************************/
/*                                                                               */
/*  Animation of wave equation on a sphere                                       */
/*                                                                               */
/*  N. Berglund, july 2023                                                       */
/*                                                                               */
/*  Feel free to reuse, but if doing so it would be nice to drop a               */
/*  line to nils.berglund@univ-orleans.fr - Thanks!                              */
/*                                                                               */
/*  compile with                                                                 */
/*  gcc -o wave_sphere wave_sphere.c                                             */
/* -L/usr/X11R6/lib -ltiff -lm -lGL -lGLU -lX11 -lXmu -lglut -O3 -fopenmp        */
/*                                                                               */
/*  OMP acceleration may be more effective after executing                       */
/*  export OMP_NUM_THREADS=2 in the shell before running the program             */
/*                                                                               */
/*  To make a video, set MOVIE to 1 and create subfolder tif_wave                */
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
#define DOUBLE_MOVIE 0  /* set to 1 to produce movies for wave height and energy simultaneously */
#define SAVE_MEMORY 1   /* set to 1 to save memory when writing tiff images */
#define NO_EXTRA_BUFFER_SWAP 1    /* some OS require one less buffer swap when recording images */

/* General geometrical parameters */

#define WINWIDTH 	1920  /* window width */
#define WINHEIGHT 	1150  /* window height */
#define NX 2560          /* number of grid points on x axis */
#define NY 1280           /* number of grid points on y axis */

#define DPOLE 15         /* safety distance to poles */
#define SMOOTHPOLE 0.1     /* smoothing coefficient at poles */
#define ZERO_MERIDIAN 0.0     /* choice of zero meridian (will be at left/right boundary of 2d plot) */

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

#define B_DOMAIN 81         /* choice of domain shape, see list in global_pdes.c */

#define CIRCLE_PATTERN 31   /* pattern of circles or polygons, see list in global_pdes.c */

#define COMPARISON 0        /* set to 1 to compare two different patterns */
#define B_DOMAIN_B 20       /* second domain shape, for comparisons */
#define CIRCLE_PATTERN_B 0  /* second pattern of circles or polygons */

#define VARIABLE_IOR 0      /* set to 1 for a variable index of refraction */
#define IOR 20              /* choice of index of refraction, see list in global_pdes.c */
#define IOR_TOTAL_TURNS 1.0 /* total angle of rotation for IOR_PERIODIC_WELLS_ROTATING */
#define MANDEL_IOR_SCALE -0.05   /* parameter controlling dependence of IoR on Mandelbrot escape speed */

#define P_PERCOL 0.25       /* probability of having a circle in C_RAND_PERCOL arrangement */
#define NPOISSON 1000        /* number of points for Poisson C_RAND_POISSON arrangement */
#define PDISC_FACTOR 3.25   /* controls density of Poisson disc process (default: 3.25) */
#define RANDOM_POLY_ANGLE 1 /* set to 1 to randomize angle of polygons */

#define LAMBDA 0.75	    /* parameter controlling the dimensions of domain */
#define MU 0.2              /* parameter controlling the dimensions of domain */
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
#define WALL_WIDTH 0.6      /* width of wall separating lenses */
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

#define OMEGA 0.001       /* frequency of periodic excitation */
#define AMPLITUDE 0.8     /* amplitude of periodic excitation */ 
#define ACHIRP 0.2        /* acceleration coefficient in chirp */
#define DAMPING 0.0       /* damping of periodic excitation */
#define COURANT 0.1       /* Courant number */
#define COURANTB 0.005    /* Courant number in medium B */
#define GAMMA 0.0         /* damping factor in wave equation */
#define GAMMAB 0.0        /* damping factor in wave equation */
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
#define OSCILLATING_SOURCE_PERIOD 25    /* period of oscillating source */
#define ALTERNATE_OSCILLATING_SOURCE 1  /* set to 1 to alternate sign of oscillating source */

#define ADD_WAVE_PACKET_SOURCES 0       /* set to 1 to add several sources emitting wave packets */
#define WAVE_PACKET_SOURCE_TYPE 1       /* type of wave packet sources */
#define N_WAVE_PACKETS 15               /* number of wave packets */
#define WAVE_PACKET_RADIUS 20            /* radius of wave packets */

#define ADD_FORCING 0       /* set to 1 to add periodic forcing */
#define FORCING_AMP 0.0     /* amplitude of periodic forcing */
#define FORCING_CONST_AMP 1.0e-10    /* amplitude of periodic forcing */
#define FORCING_PERIOD 2400  /* period of forcing */

#define DRIFT_WAVE 0        /* add drift of wave to the right (experimental) */
#define DRIFT_FREQ 5        /* frequency of drift adding (experimental) */

#define MOVING_FRAME 0      /* set to 1 to use wave equation in moving frame */
#define VOVERC 0.025          /* moving frame speed over wave speed */

/* Boundary conditions, see list in global_pdes.c  */

#define B_COND 2

#define PRECOMPUTE_BC 0     /* set to 1 to compute neighbours for Laplacian in advance */

/* Parameters for length and speed of simulation */

// #define NSTEPS 2800       /* number of frames of movie */
#define NSTEPS 1800         /* number of frames of movie */
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

#define INITIAL_AMP 2.5            /* amplitude of initial condition */
#define INITIAL_VARIANCE 0.003    /* variance of initial condition */
#define INITIAL_WAVELENGTH  0.1   /* wavelength of initial condition */

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
#define NON_DIRICHLET_BC 1      /* set to 1 to draw only facets in domain, if field is not zero on boundary */
#define FLOOR_ZCOORD 0          /* set to 1 to draw only facets with z not too negative */
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
#define FLOODING 0              /* set to 1 to draw waves above altitude (for Earth representations) */

#define REP_AXO_3D 0        /* linear projection (axonometry) */
#define REP_PROJ_3D 1       /* projection on plane orthogonal to observer line of sight */

#define ROTATE_VIEW 1       /* set to 1 to rotate position of observer */
#define ROTATE_ANGLE 360.0   /* total angle of rotation during simulation */
#define ROTATE_VIEW_WHILE_FADE 1    /* set to 1 to keep rotating viewpoint during fade */

#define VIEWPOINT_TRAJ 1    /* type of viewpoint trajectory */
#define MAX_LATITUDE 45.0   /* maximal latitude for viewpoint trajectory VP_ORBIT2 */

/* Color schemes */

#define COLOR_PALETTE 14      /* Color palette, see list in global_pdes.c  */
#define COLOR_PALETTE_B 16     /* Color palette, see list in global_pdes.c  */

#define BLACK 1          /* background */
#define COLOR_OUT_R 1.0    /* color outside domain */
#define COLOR_OUT_G 1.0    
#define COLOR_OUT_B 1.0    

#define COLOR_SCHEME 3   /* choice of color scheme, see list in global_pdes.c  */

#define SCALE 0          /* set to 1 to adjust color scheme to variance of field */
#define SLOPE 1.0        /* sensitivity of color on wave amplitude */
#define VSCALE_AMPLITUDE 0.5  /* additional scaling factor for color scheme P_3D_AMPLITUDE */
#define VSHIFT_AMPLITUDE 0.0   /* additional shift for wave amplitude */
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

#define MESSAGE_LDASH 14         /* length of dash for Morse code message */
#define MESSAGE_LDOT 8           /* length of dot for Morse code message */
#define MESSAGE_LINTERVAL 54     /* length of interval between dashes/dots for Morse code message */
#define MESSAGE_LINTERLETTER 60  /* length of interval between letters for Morse code message */
#define MESSAGE_LSPACE 48        /* length of space for Morse code message */
#define MESSAGE_INITIAL_TIME 100 /* initial time before starting message for Morse code message */    

#define NXMAZE 8      /* width of maze */
#define NYMAZE 32      /* height of maze */
#define MAZE_MAX_NGBH 5     /* max number of neighbours of maze cell */
#define RAND_SHIFT 1        /* seed of random number generator */
#define MAZE_XSHIFT 0.0     /* horizontal shift of maze */
#define MAZE_WIDTH 0.02     /* half width of maze walls */

#define DRAW_COLOR_SCHEME 0       /* set to 1 to plot the color scheme */
#define COLORBAR_RANGE 1.5      /* scale of color scheme bar */
#define COLORBAR_RANGE_B 6.0    /* scale of color scheme bar for 2nd part */
#define ROTATE_COLOR_SCHEME 0     /* set to 1 to draw color scheme horizontally */
#define CIRC_COLORBAR 0         /* set to 1 to draw circular color scheme */
#define CIRC_COLORBAR_B 0       /* set to 1 to draw circular color scheme */
#define DRAW_MOON_POSITION 0    /* set to 1 to draw position of Moon (for tide simulation) */

#define DRAW_WAVE_PROFILE 0     /* set to 1 to draw a profile of the wave */
#define VERTICAL_WAVE_PROFILE 0 /* set to 1 to draw wave profile vertically */
#define DRAW_WAVE_TIMESERIES 0  /* set to 1 to draw a time series of the wave */
#define TIMESERIES_NVALUES 400  /* number of values plotted in time series */
#define SAVE_TIME_SERIES 0      /* set to 1 to save wave time series at a point */
#define OSCIL_YMAX 0.35      /* defines oscillation range */
#define AVRG_E_FACTOR 0.95   /* controls time window size in P_AVERAGE_ENERGY scheme */
#define HORIZONTAL_WAVE_PROFILE 0 /* set to 1 to draw wave profile vertically */
#define WAVE_PROFILE_X 2.1      /* value of x to sample wave profile */
#define WAVE_PROFILE_Y -1.0      /* value of y to sample wave profile */
#define PROFILE_AT_BOTTOM 0     /* draw wave profile at bottom instead of top */
#define AVERAGE_WAVE_PROFILE 0  /* set to 1 to draw time-average of wave profile squared*/
#define TIMESERIES_NVALUES 400  /* number of values plotted in time series */
#define DRAW_WAVE_SOURCE 0      /* set to 1 to draw source of wave at (wave_source_x, wave_source_y) */

#define ADD_POTENTIAL 0         /* set to 1 to add potential to z coordinate */
#define POTENTIAL 10
#define POT_FACT 20.0
#define HRES 1          /* dummy, only used by rde.c */
#define INITIAL_SHIFT 20.0          /* time shift of initial wave packet (in oscillation periods) */
#define WAVE_PACKET_SHIFT 200.0     /* time shift between wave packets (in oscillation periods) */
#define FADE_IN_OBSTACLE 0      /* set to 1 to fade color inside obstacles */
#define N_SOURCES 20                     /* number of sources, for option draw_sources */
#define XYIN_INITIALISED (B_DOMAIN == D_IMAGE)
/* end of constants only used by sub_wave and sub_maze */

/* For debugging purposes only */
#define FLOOR 1         /* set to 1 to limit wave amplitude to VMAX */
#define VMAX 5.0       /* max value of wave amplitude */

/* Parameters controlling 3D projection */

double u_3d[2] = {0.75, -0.45};     /* projections of basis vectors for REP_AXO_3D representation */
double v_3d[2] = {-0.75, -0.45};
double w_3d[2] = {0.0, 0.015};
double light[3] = {-0.816496581, -0.40824829, 0.40824829};      /* vector of "light" direction for P_3D_ANGLE color scheme */
double observer[3] = {-5.0, -5.0, 3.5};    /* location of observer for REP_PROJ_3D representation */ 
int reset_view = 0;         /* switch to reset 3D view parameters (for option ROTATE_VIEW) */

#define ADD_DEM 0               /* add DEM (digital elevation model) */
#define ADD_NEGATIVE_DEM 0      /* add DEM with bathymetric data */
#define RSCALE_DEM 0.1         /* scaling factor of radial component for DEM */
#define SMOOTH_DEM 5            /* set to 1 to smoothen DEM (to make altitude less constant) */
#define DEM_SMOOTH_STEPS 1      /* number of smoothening steps */
#define DEM_SMOOTH_HEIGHT 2.0   /* relative height below which to smoothen */
#define DEM_MAXHEIGHT 9000.0     /* max height of DEM (estimated from Everest/Olympus Mons) */
#define DEM_MAXDEPTH -10000     /* max depth of DEM */
#define PLANET_SEALEVEL -10000.0      /* sea level for flooded planet */
#define VENUS_NODATA_FACTOR 0.5     /* altitude to assign to DEM points without data (fraction of mean altitude) */
#define TRANSPARENT_WAVE 0      /* set to 1 for waves to be "transparent" */

#define RSCALE 1.0              /* scaling factor of radial component */
#define RMAX 1.1               /* max value of radial component */
#define RMIN 0.9               /* min value of radial component */
#define Z_SCALING_FACTOR 0.8   /* overall scaling factor of z axis for REP_PROJ_3D representation */
#define XY_SCALING_FACTOR 2.2   /* overall scaling factor for on-screen (x,y) coordinates after projection */
#define ZMAX_FACTOR 1.0         /* max value of z coordinate for REP_PROJ_3D representation */
#define XSHIFT_3D 0.0           /* overall x shift for REP_PROJ_3D representation */
#define YSHIFT_3D 0.0           /* overall y shift for REP_PROJ_3D representation */
#define COS_VISIBLE -0.02        /* limit on cosine of normal to shown facets */

#include "global_pdes.c"        /* constants and global variables */
#include "global_3d.c"          /* constants and global variables */
#include "sub_maze.c"           /* support for generating mazes */
#include "sub_wave.c"           /* common functions for wave_billiard, heat and schrodinger */
#include "wave_common.c"        /* common functions for wave_billiard, wave_comparison, etc */

#include "sub_wave_3d.c"        /* graphical functions specific to wave_3d */
#include "sub_sphere.c"         /* graphical functions specific to wave_sphere */

FILE *time_series_left, *time_series_right, *image_file;

double courant2, courantb2;  /* Courant parameters squared */

void compute_forcing_schedule(int t, t_wave_sphere wsphere[NX*NY])
/* compute periodic forcing */
{
    int i, j;
    double phase;
    
    phase = DPI*(double)t/(double)FORCING_PERIOD + ZERO_MERIDIAN*PI/180.0 + 0.5*PID;
    
    for (i=0; i<NX; i++) for (j=0.0; j<NY; j++)
    {
        wsphere[i*NY+j].force = FORCING_AMP*wsphere[i*NY+j].stheta*sin(2.0*(wsphere[i*NY+j].phi - phase));
        wsphere[i*NY+j].force += FORCING_CONST_AMP*wsphere[i*NY+j].stheta*cos(2.0*(wsphere[i*NY+j].phi));
    }
    
    moon_position = (int)((double)NX*(phase/DPI + 0.625));
    while (moon_position < 0) moon_position += NX;
    while (moon_position >= NX) moon_position -= NX;
    
//     i = NX/2;
    i = moon_position;
    printf("Phase = %.5lg, Forcing at i = %i: %.5lg, Moon position = %i\n", phase, i, wsphere[i*NY+NY/2].force, moon_position);
}


void shift_fields(double phi[NX*NY], short int xy_in[NX*NY])
/* EXPERIMENTAL */
{
    int i, j;
    double temp;
    
    #pragma omp parallel for private(j,i)
//     for (j=DPOLE; j<NY-DPOLE; j++)
    for (j=0; j<NY; j++)
    {
        temp = phi[j];
        for (i=0; i<NX-1; i++) if (xy_in[i*NY+j])
        {
            phi[i*NY+j] = phi[(i+1)*NY+j];
        }
        if (xy_in[j]) 
        {
            phi[(NX-1)*NY+j] = temp;
        }
    }
    
//     for (j=0; j<NY; j++)
//     {
//         temp = phi[(NX-1)*NY+j];
//         for (i=NX-1; i>1; i--) if (xy_in[i*NY+j])
//         {
//             phi[i*NY+j] = phi[(i-1)*NY+j];
//         }
//         if (xy_in[j]) 
//         {
//             phi[j] = temp;
//         }
//     }
    
    /* smoothing at poles */
    for (j=0; j<DPOLE; j++)
        for (i=0; i<NX; i++)
            phi[i*NY+j] = 0.0;
    
    for (j=NY-DPOLE; j<NY; j++)
        for (i=0; i<NX; i++)
            phi[i*NY+j] = 0.0;
}


void evolve_wave_half(double phi_in[NX*NY], double psi_in[NX*NY], double phi_out[NX*NY], 
                      short int xy_in[NX*NY], double tc[NX*NY], double tcc[NX*NY], double tgamma[NX*NY], t_wave_sphere wsphere[NX*NY], int t)
/* time step of field evolution */
/* phi is value of field at time t, psi at time t-1 */
/* this version of the function has been rewritten in order to minimize the number of if-branches */
{
    int i, j, iplus, iminus, jplus, jminus, jtop, jbot;
    double delta, x, y, c, cc, gamma, sintheta, cottheta, invstheta, sum, avrg, force;
    double tmp_phi[NY];
    static long time = 0;
    static short int first = 1;
    static double dphi, dtheta, cphiphi, ctheta, dt, vdrift;
    static int counter = 0;
    
    if (first)
    {
        dphi = DPI/(double)NX;
        dtheta = PI/(double)NY;
        cphiphi = dphi*dphi/(dtheta*dtheta);
        ctheta = dphi*dphi/(2.0*dtheta);
        dt = (double)NVID*COURANT*(XMAX-XMIN)/((double)NX);
        vdrift = COURANT*VOVERC;
        
        printf("dphi = %.5lg, dtheta = %.5lg, cphiphi = %.5lg, ctheta = %.5lg, dt = %.5lg\n", dphi, dtheta, cphiphi, ctheta, dt); 
        
        first = 0;
    }
    
    time++;
    
    #pragma omp parallel for private(i,j,iplus,iminus,jplus,jminus,delta,x,y)
    /* evolution in the bulk */
    for (j=DPOLE; j<NY-DPOLE; j++){
//         sintheta = sin(j*dtheta);
        /* NEW: used wsphere to accelerate computation */
        sintheta = wsphere[i*NY+j].stheta;
//         invstheta = 1.0/(sintheta*sintheta);
        invstheta = 1.0/(sintheta*sintheta + SMOOTHPOLE*SMOOTHPOLE);
//         cottheta = ctheta*cos(j*dtheta)/sintheta;
//         cottheta = ctheta*cos(j*dtheta)/(sintheta + SMOOTHPOLE);
        cottheta = ctheta*wsphere[i*NY+j].ctheta/(sintheta + SMOOTHPOLE);
        
        
        for (i=1; i<NX-1; i++){
            if ((TWOSPEEDS)||(xy_in[i*NY+j] != 0)){
                x = phi_in[i*NY+j];
		y = psi_in[i*NY+j];
                
                /* discretized Laplacian */
                /* 2nd phi derivative */
                delta = invstheta*(phi_in[(i+1)*NY+j] + phi_in[(i-1)*NY+j] - 2.0*x);
                
                /* 2nd theta derivative */
                delta += cphiphi*(phi_in[i*NY+j+1] + phi_in[i*NY+j-1] - 2.0*x);
                
                /* first theta derivative */
                delta += cottheta*(phi_in[i*NY+j+1] - phi_in[i*NY+j-1]);

                /* evolve phi */
                phi_out[i*NY+j] = -y + 2*x + tcc[i*NY+j]*delta - KAPPA*x - tgamma[i*NY+j]*(x-y);
                
                if (ADD_FORCING) phi_out[i*NY+j] += wsphere[i*NY+j].force*dt;
                
                if (MOVING_FRAME)
                {
//                     phi_out[i*NY+j] += vdrift*wsphere[i*NY+j].stheta*(phi_in[(i+1)*NY+j] - phi_in[(i-1)*NY+j] - psi_in[(i+1)*NY+j] + psi_in[(i-1)*NY+j]);
                    phi_out[i*NY+j] += VOVERC*tc[i*NY+j]*wsphere[i*NY+j].stheta*(phi_in[(i+1)*NY+j] - phi_in[(i-1)*NY+j] - psi_in[(i+1)*NY+j] + psi_in[(i-1)*NY+j]);
                    /* the .stheta factor is not really correct, but prevents instability */
                }
                
                /* optional: add const*delta to smoothen ? */
            }
        }
    }
    
    /* evolution at longitude zero */
    for (j=DPOLE; j<NY-DPOLE; j++){
//         sintheta = sin(j*dtheta);
        sintheta = wsphere[j].stheta;
        invstheta = 1.0/(sintheta*sintheta + SMOOTHPOLE*SMOOTHPOLE);
//         cottheta = ctheta*cos(j*dtheta)/(sintheta + SMOOTHPOLE);
        cottheta = ctheta*wsphere[j].ctheta/(sintheta + SMOOTHPOLE);
        
        
        /* i = 0 */
        if ((TWOSPEEDS)||(xy_in[j] != 0)){
            x = phi_in[j];
            y = psi_in[j];
                
            /* discretized Laplacian */
            /* 2nd phi derivative */
            delta = invstheta*(phi_in[NY+j] + phi_in[(NX-1)*NY+j] - 2.0*x);
                
            /* 2nd theta derivative */
            delta += cphiphi*(phi_in[j+1] + phi_in[j-1] - 2.0*x);
                
            /* first theta derivative */
            delta += cottheta*(phi_in[j+1] - phi_in[j-1]);

            /* evolve phi */
            phi_out[j] = -y + 2*x + tcc[j]*delta - KAPPA*x - tgamma[j]*(x-y);
            
            if (ADD_FORCING) phi_out[j] += wsphere[j].force*dt;
            
            if (MOVING_FRAME)
            {
//                 phi_out[j] += vdrift*wsphere[j].stheta*(phi_in[NY+j] - phi_in[(NX-1)*NY+j] - psi_in[NY+j] + psi_in[(NX-1)*NY+j]);
                phi_out[j] += VOVERC*tc[j]*wsphere[j].stheta*(phi_in[NY+j] - phi_in[(NX-1)*NY+j] - psi_in[NY+j] + psi_in[(NX-1)*NY+j]);
            }
//             if (j%100==0) 
//                 printf("tcc*delta = %.5lg, forcing = %.5lg\n", tcc[j]*delta, forcing_schedule(t, i, 0, wsphere)*dt);
        }
        
        /* i = NX-1 */
        if ((TWOSPEEDS)||(xy_in[(NX-1)*NY+j] != 0)){
            x = phi_in[(NX-1)*NY+j];
            y = psi_in[(NX-1)*NY+j];
                
            /* discretized Laplacian */
            /* 2nd phi derivative */
            delta = invstheta*(phi_in[j] + phi_in[(NX-2)*NY+j] - 2.0*x);
                
            /* 2nd theta derivative */
            delta += cphiphi*(phi_in[(NX-1)*NY+j+1] + phi_in[(NX-1)*NY+j-1] - 2.0*x);
                
            /* first theta derivative */
            delta += cottheta*(phi_in[(NX-1)*NY+j+1] - phi_in[(NX-1)*NY+j-1]);

            /* evolve phi */
            phi_out[(NX-1)*NY+j] = -y + 2*x + tcc[(NX-1)*NY+j]*delta - KAPPA*x - tgamma[(NX-1)*NY+j]*(x-y);
            
            if (ADD_FORCING) phi_out[(NX-1)*NY + j] += wsphere[(NX-1)*NY + j].force*dt;
            
            if (MOVING_FRAME)
            {
//                 phi_out[(NX-1)*NY+j] += vdrift*wsphere[(NX-1)*NY+j].stheta*(phi_in[j] - phi_in[(NX-2)*NY+j] - psi_in[j] + psi_in[(NX-2)*NY+j]);
                phi_out[(NX-1)*NY+j] += VOVERC*tc[i*NY+j]*wsphere[(NX-1)*NY+j].stheta*(phi_in[j] - phi_in[(NX-2)*NY+j] - psi_in[j] + psi_in[(NX-2)*NY+j]);
            }
        }
    }
    
    /* compute average at north pole */
    sum = 0.0;
//     for (i=0; i<NX; i++) sum += phi_out[i*NY + DPOLE];
    for (i=0; i<NX; i++) sum += phi_out[i*NY + DPOLE - 1];
    avrg = sum/(double)NX;
    for (i=0; i<NX; i++) for (j=0; j<DPOLE; j++)
        phi_out[i*NY + j] = avrg;
    
    /* compute average at south pole */
    sum = 0.0;
//     for (i=0; i<NX; i++) sum += phi_out[i*NY + NY-1-DPOLE];
    for (i=0; i<NX; i++) sum += phi_out[i*NY + NY-DPOLE];
    avrg = sum/(double)NX;
    for (i=0; i<NX; i++) for (j=NY-DPOLE; j<NY; j++) 
        phi_out[i*NY + j] = avrg;
    
    /* for debugging purposes/if there is a risk of blow-up */
    if (FLOOR) for (i=0; i<NX; i++){
        for (j=0; j<NY; j++){
            if (xy_in[i*NY+j] != 0) 
            {
                if (phi_out[i*NY+j] > VMAX) phi_out[i*NY+j] = VMAX;
                if (phi_out[i*NY+j] < -VMAX) phi_out[i*NY+j] = -VMAX;
            }
        }
    }
}


void evolve_wave(double phi[NX*NY], double psi[NX*NY], double tmp[NX*NY], short int xy_in[NX*NY],
    double tc[NX*NY], double tcc[NX*NY], double tgamma[NX*NY], t_laplacian laplace[NX*NY], t_laplacian laplace1[NX*NY], t_laplacian laplace2[NX*NY], t_wave_sphere wsphere[NX*NY], int t)
/* time step of field evolution */
/* phi is value of field at time t, psi at time t-1 */
{
    evolve_wave_half(phi, psi, tmp, xy_in, tc, tcc, tgamma, wsphere, t);
    evolve_wave_half(tmp, phi, psi, xy_in, tc, tcc, tgamma, wsphere, t);
    evolve_wave_half(psi, tmp, phi, xy_in, tc, tcc, tgamma, wsphere, t);
}


void draw_color_bar_palette(int plot, double range, int palette, int circular, int fade, double fade_value)
{
//     double width = 0.2;
    double width = 0.12;
//     double width = 0.2;
    
    width *= (double)NX/(double)WINWIDTH;
    
    if (PLOT_2D)
    {
        if (ROTATE_COLOR_SCHEME) 
            draw_color_scheme_palette_fade(XMIN + 0.6, YMIN + 0.1, XMAX - 0.6, YMIN + 0.1 + width, plot, -range, range, palette, fade, fade_value);
        else 
            draw_color_scheme_palette_fade(XMAX - 1.5*width, YMIN + 0.1, XMAX - 0.5*width, YMAX - 0.1, plot, -range, range, palette, fade, fade_value);
    }
    else
    {
        if (ROTATE_COLOR_SCHEME) 
            draw_color_scheme_palette_3d(-1.0, -0.8, XMAX - 0.1, -1.0, plot, -range, range, palette, fade, fade_value);
        else if (circular)
            draw_circular_color_scheme_palette_3d(XMAX - 2.0*width, YMIN + 2.0*width, 1.5*width, 1.3*width, plot, -range, range, palette, fade, fade_value);
        else 
            draw_color_scheme_palette_3d(XMAX - 1.5*width, YMIN + 0.1, XMAX - 0.5*width, YMAX - 0.1, plot, -range, range, palette, fade, fade_value);
    }

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
    double time, scale, ratio, startleft[2], startright[2], sign = -1.0, r2, xy[2], fade_value, yshift, speed = 0.0, a, b, c, angle = 0.0, lambda1, y, x1, sign1, omega, phase_shift, theta, amp; 
    double *phi, *psi, *tmp, *color_scale, *tc, *tcc, *tgamma;
//     double *total_energy;
    short int *xy_in;
    int i, j, s, sample_left[2], sample_right[2], period = 0, fade, source_counter = 0, k, p, q, drift_counter = 0;
    static int counter = 0, first_source = 1;
    long int wave_value;
    t_wave *wave;
    t_laplacian *laplace, *laplace1, *laplace2;
    t_wave_source wave_source[25];
    t_wave_sphere *wsphere;
    
    if (SAVE_TIME_SERIES)
    {
        time_series_left = fopen("wave_left.dat", "w");
        time_series_right = fopen("wave_right.dat", "w");
    }

    /* Since NX and NY are big, it seemed wiser to use some memory allocation here */
    xy_in = (short int *)malloc(NX*NY*sizeof(short int));
    phi = (double *)malloc(NX*NY*sizeof(double));
    psi = (double *)malloc(NX*NY*sizeof(double));
    tmp = (double *)malloc(NX*NY*sizeof(double));
//     total_energy = (double *)malloc(NX*NY*sizeof(double));
    color_scale = (double *)malloc(NX*NY*sizeof(double));
    tc = (double *)malloc(NX*NY*sizeof(double));
    tcc = (double *)malloc(NX*NY*sizeof(double));
    tgamma = (double *)malloc(NX*NY*sizeof(double));
    
    wave = (t_wave *)malloc(NX*NY*sizeof(t_wave));
    wsphere = (t_wave_sphere *)malloc(NX*NY*sizeof(t_wave_sphere));
    
    laplace = (t_laplacian *)malloc(NX*NY*sizeof(t_laplacian));
    laplace1 = (t_laplacian *)malloc(NX*NY*sizeof(t_laplacian));
    laplace2 = (t_laplacian *)malloc(NX*NY*sizeof(t_laplacian));
    
    /* initialise positions and radii of circles */
    if (COMPARISON)
    {
        if ((B_DOMAIN == D_CIRCLES)||(B_DOMAIN == D_CIRCLES_IN_RECT))
            ncircles = init_circle_config_pattern(circles, CIRCLE_PATTERN);
        else if (B_DOMAIN == D_POLYGONS) ncircles = init_polygon_config_pattern(polygons, CIRCLE_PATTERN);
        
        if ((B_DOMAIN_B == D_CIRCLES)||(B_DOMAIN_B == D_CIRCLES_IN_RECT))
            ncircles_b = init_circle_config_pattern(circles_b, CIRCLE_PATTERN_B);
        else if (B_DOMAIN_B == D_POLYGONS) ncircles_b = init_polygon_config_pattern(polygons_b, CIRCLE_PATTERN_B);
        /* TO DO: adapt to different polygon patterns */
    }
    else
    {
        if ((B_DOMAIN == D_CIRCLES)||(B_DOMAIN == D_CIRCLES_IN_RECT)) ncircles = init_circle_config(circles);
        else if (B_DOMAIN == D_POLYGONS) ncircles = init_polygon_config(polygons);
    }
    
    printf("Polygons initialized\n");
    
    /* initialise polyline for von Koch and similar domains */
    npolyline = init_polyline(MDEPTH, polyline);
    for (i=0; i<npolyline; i++) printf("vertex %i: (%.3f, %.3f)\n", i, polyline[i].x, polyline[i].y);
    if (COMPARISON) npolyline_b = init_polyline(MDEPTH, polyline);
    
    npolyrect = init_polyrect(polyrect);
    for (i=0; i<npolyrect; i++) printf("polyrect vertex %i: (%.3f, %.3f) - (%.3f, %.3f)\n", i, polyrect[i].x1, polyrect[i].y1, polyrect[i].x2, polyrect[i].y2);

    init_polyrect_arc(polyrectrot, polyarc, &npolyrect_rot, &npolyarc);
    printf("Rotated rectangles and arcs initialized\n");
    printf("%i rotated rectangles, %i arcs\n", npolyrect_rot, npolyarc);
    
    if ((B_DOMAIN == D_SPHERE_CIRCLES)||(B_DOMAIN_B == D_SPHERE_CIRCLES))
    {
        ncircles = init_circle_sphere(circ_sphere, CIRCLE_PATTERN); 
    }
    
    courant2 = COURANT*COURANT;
    courantb2 = COURANTB*COURANTB;
    c = COURANT*(XMAX - XMIN)/(double)NX;

    a = 0.05;
    b = 0.0016;
    
    /* initialize color scale, for option RESCALE_COLOR_IN_CENTER */
    if (RESCALE_COLOR_IN_CENTER)
    {
        for (i=0; i<NX; i++)
            for (j=0; j<NY; j++)
            {
                ij_to_xy(i, j, xy);
                r2 = xy[0]*xy[0] + xy[1]*xy[1];
                color_scale[i*NY+j] = 1.0 - exp(-4.0*r2/LAMBDA*LAMBDA);
            }
    }

    /* initialize wave with a drop at one point, zero elsewhere */
//     init_circular_wave(0.0, -LAMBDA, phi, psi, xy_in);

    
    init_wave_fields(wave);
    
    init_wave_sphere(wsphere); 
    
    /* initialize total energy table - no longer needed */
//     if ((ZPLOT == P_MEAN_ENERGY)||(ZPLOT_B == P_MEAN_ENERGY)||(ZPLOT == P_LOG_MEAN_ENERGY)||(ZPLOT_B == P_LOG_MEAN_ENERGY))
//         for (i=0; i<NX; i++)
//             for (j=0; j<NY; j++) 
//                 total_energy[i*NY+j] = 0.0;
    
    ratio = (XMAX - XMIN)/8.4;  /* for Tokarsky billiard */
    

//     init_circular_wave_mod(polyline[85].x, polyline[85].y, phi, psi, xy_in);
//     init_circular_wave_mod(LAMBDA*cos(APOLY*PID), LAMBDA*sin(APOLY*PID), phi, psi, xy_in);
    lambda1 = LAMBDA;
    angle = DPI/(double)NPOLY;
//     init_circular_wave_mod(lambda1*cos(0.5*angle), lambda1*sin(0.5*angle), phi, psi, xy_in);
//     for (j=1; j<NPOLY; j++)
//         add_circular_wave_mod(1.0, lambda1*cos(((double)j+0.5)*angle), lambda1*sin(((double)j+0.5)*angle), phi, psi, xy_in);
        
    /* Kilauea, Hawaii, USA */
//     init_circular_wave_sphere(155.2867*PI/180.0, 19.42109*PI/180.0, phi, psi, xy_in, wsphere);
    
//     init_moving_tidal_wave_sphere(0.0, 0.05, 0.025, phi, psi, xy_in, wsphere);

//     init_moving_tidal_wave_sphere(0.0, 2.5, 1.25, phi, psi, xy_in, wsphere);
    
//     init_circular_wave_sphere(0.0, -0.5, phi, psi, xy_in, wsphere);
    
     init_wave_flat_sphere(phi, psi, xy_in, wsphere);
    
//     init_wave_flat_sphere(phi, psi, xy_in, wsphere);
//     for (j=0; j<6; j++)
//         add_circular_wave_sphere(1.0, (double)j*PI/3.0, 0.925*PID, phi, psi, xy_in, wsphere);
    
//     add_circular_wave_sphere(1.0, 1.0 - PI/9.0 + DPI/3.0, 0.15, phi, psi, xy_in, wsphere);
//     add_circular_wave_sphere(1.0, 1.0 - PI/9.0 + 2.0*DPI/3.0, 0.15, phi, psi, xy_in, wsphere);
    
    for (j=0; j<5; j++)
    {
        a = circ_sphere[1].x + circ_sphere[2].x;
        b = circ_sphere[1].y + circ_sphere[2].y;
        c = 1.0 + circ_sphere[1].z + circ_sphere[2].z;
        theta = acos(c/sqrt(a*a + b*b + c*c));
        
        wave_source_x[4*j] = (double)j*DPI/5.0;
        wave_source_y[4*j] = PID - theta;
        wave_source_x[4*j+1] = (double)j*DPI/5.0 + PI/5.0;
        wave_source_y[4*j+1] = theta - PID;
        
        a = circ_sphere[1].x + circ_sphere[2].x + circ_sphere[6].x;
        b = circ_sphere[1].y + circ_sphere[2].y + circ_sphere[6].y;
        c = circ_sphere[1].z + circ_sphere[2].z + circ_sphere[6].z;
        theta = acos(c/sqrt(a*a + b*b + c*c));
        
        wave_source_x[4*j+2] = (double)j*DPI/5.0 + PI/5.0;
        wave_source_y[4*j+2] = theta - PID;
        wave_source_x[4*j+3] = (double)j*DPI/5.0;
        wave_source_y[4*j+3] = PID - theta;
        
    }
    
    
//     printf("Wave initialized\n");

    /* initialize table of wave speeds/dissipation */
    init_speed_dissipation_sphere(xy_in, tc, tcc, tgamma, wsphere);
    
    
    
    /* initialze potential to add to z coordinate */ 
    if (ADD_POTENTIAL)
    {
        if (POTENTIAL == POT_IOR) 
            for (i=0; i<NX*NY; i++)
                wave[i].potential = &tcc[i];
    }
    
    init_zfield(phi, psi, xy_in, ZPLOT, wave, 0);
    init_cfield(phi, psi, xy_in, CPLOT, wave, 0);
    
    if (DOUBLE_MOVIE)
    {
        init_zfield(phi, psi, xy_in, ZPLOT_B, wave, 1);
        init_cfield(phi, psi, xy_in, CPLOT_B, wave, 1);
    }

    blank();
    glColor3f(0.0, 0.0, 0.0);
    draw_wave_sphere(0, phi, psi, xy_in, wave, wsphere, ZPLOT, CPLOT, COLOR_PALETTE, 0, 1.0, 1);
//     draw_billiard();
    
    
    if (DRAW_COLOR_SCHEME) draw_color_bar_palette(CPLOT, COLORBAR_RANGE, COLOR_PALETTE, CIRC_COLORBAR, 0, 1.0);

    glutSwapBuffers();



    sleep(SLEEP1);

    for (i=0; i<=INITIAL_TIME + NSTEPS; i++)
    {
        global_time++;
        
	//printf("%d\n",i);
        /* compute the variance of the field to adjust color scheme */
        /* the color depends on the field divided by sqrt(1 + variance) */
        if (SCALE)
        {
            scale = sqrt(1.0 + compute_variance_mod(phi,psi, xy_in));
//             printf("Scaling factor: %5lg\n", scale);
        }
        else scale = 1.0;
        
        if (ROTATE_VIEW) 
        {
            viewpoint_schedule(i - INITIAL_TIME);
            reset_view = 1;
        }
        
        if (ADD_FORCING) compute_forcing_schedule(i, wsphere);
        
        draw_wave_sphere(0, phi, psi, xy_in, wave, wsphere, ZPLOT, CPLOT, COLOR_PALETTE, 0, 1.0, 1);
        
        for (j=0; j<NVID; j++) 
        {
            evolve_wave(phi, psi, tmp, xy_in, tc, tcc, tgamma, laplace, laplace1, laplace2, wsphere, i);
            if (SAVE_TIME_SERIES)
            {
                wave_value = (long int)(phi[sample_left[0]*NY+sample_left[1]]*1.0e16);
                fprintf(time_series_left, "%019ld\n", wave_value);
                wave_value = (long int)(phi[sample_right[0]*NY+sample_right[1]]*1.0e16);
                fprintf(time_series_right, "%019ld\n", wave_value);
                if ((j == 0)&&(i%10 == 0)) printf("Frame %i of %i\n", i, NSTEPS);
//                 fprintf(time_series_right, "%.15f\n", phi[sample_right[0]][sample_right[1]]);
                
                /* Experimental */
                if (DRIFT_WAVE)
                {
                    if (drift_counter == DRIFT_FREQ)
                    {
                        shift_fields(phi, xy_in);
                        shift_fields(psi, xy_in);
                        shift_fields(tmp, xy_in);    /* needed? */    
                        drift_counter = 0;
                    }
                    else drift_counter++;
                }
            }
//             if (i % 10 == 9) oscillate_linear_wave(0.2*scale, 0.15*(double)(i*NVID + j), -1.5, YMIN, -1.5, YMAX, phi, psi);
        }
        
        /* Experimental */
//         if (DRIFT_WAVE != 0)
//             if (i%DRIFT_FREQ == 0)
//             {
//                 shift_fields(phi, xy_in);
//                 shift_fields(psi, xy_in);
//                 shift_fields(tmp, xy_in);    /* needed? */        
//             }
        
//         draw_billiard();
        
        if (DRAW_COLOR_SCHEME) draw_color_bar_palette(CPLOT, COLORBAR_RANGE, COLOR_PALETTE, CIRC_COLORBAR, fade, fade_value); 
        
        /* add oscillating waves */
        if ((ADD_OSCILLATING_SOURCE)&&(i%OSCILLATING_SOURCE_PERIOD == 1))
        {
            if (ALTERNATE_OSCILLATING_SOURCE) sign = -sign;
            for (j=0; j<N_SOURCES; j++)
                add_circular_wave_sphere(sign,wave_source_x[j], wave_source_y[j], phi, psi, xy_in, wsphere);
        }
        if (PRINT_SPEED) print_speed_3d(speed, 0, 1.0);

	if (!((NO_EXTRA_BUFFER_SWAP)&&(MOVIE))) glutSwapBuffers();

	if (MOVIE)
        {
            if (i >= INITIAL_TIME) save_frame();
//             if (i >= INITIAL_TIME) save_frame_counter(i);
//             if (i >= INITIAL_TIME) save_frame_counter(NSTEPS + MID_FRAMES + 1 + counter);
            else printf("Initial phase time %i of %i\n", i, INITIAL_TIME);
            
            if ((i >= INITIAL_TIME)&&(DOUBLE_MOVIE))
            {
                draw_wave_sphere(1, phi, psi, xy_in, wave, wsphere, ZPLOT_B, CPLOT_B, COLOR_PALETTE_B, 0, 1.0, 1);
                if (DRAW_COLOR_SCHEME) draw_color_bar_palette(CPLOT_B, COLORBAR_RANGE_B, COLOR_PALETTE_B, CIRC_COLORBAR_B, 0, 1.0);  
                if (PRINT_SPEED) print_speed_3d(speed, 0, 1.0);
                glutSwapBuffers();
                save_frame_counter(NSTEPS + MID_FRAMES + 1 + counter);
//                 save_frame_counter(i);
                counter++;
            }
            else if (NO_EXTRA_BUFFER_SWAP) glutSwapBuffers();

            /* it seems that saving too many files too fast can cause trouble with the file system */
            /* so this is to make a pause from time to time - parameter PAUSE may need adjusting   */
            if (i % PAUSE == PAUSE - 1)
            {
                printf("Making a short pause\n");
                sleep(PSLEEP);
                s = system("mv wave*.tif tif_wave/");
            }
        }
        else printf("Computing frame %i\n", i);

    }

    if (MOVIE) 
    {
        if (DOUBLE_MOVIE) 
        {
            draw_wave_sphere(0, phi, psi, xy_in, wave, wsphere, ZPLOT, CPLOT, COLOR_PALETTE, 0, 1.0, 1);
            if (DRAW_COLOR_SCHEME) draw_color_bar_palette(CPLOT, COLORBAR_RANGE, COLOR_PALETTE, CIRC_COLORBAR, 0, 1.0);   
            if (PRINT_SPEED) print_speed_3d(speed, 0, 1.0);
            glutSwapBuffers();
            
            if (!FADE) for (i=0; i<MID_FRAMES; i++) save_frame();
            else for (i=0; i<MID_FRAMES; i++) 
            {
                if ((ROTATE_VIEW)&&(ROTATE_VIEW_WHILE_FADE))
                {
                    viewpoint_schedule(NSTEPS - INITIAL_TIME + i);
                    reset_view = 1;
                }
                fade_value = 1.0 - (double)i/(double)MID_FRAMES;
                draw_wave_sphere(0, phi, psi, xy_in, wave, wsphere, ZPLOT, CPLOT, COLOR_PALETTE, 1, fade_value, 1);
                if (DRAW_COLOR_SCHEME) draw_color_bar_palette(CPLOT, COLORBAR_RANGE, COLOR_PALETTE, CIRC_COLORBAR, 1, fade_value);   
                if (PRINT_SPEED) print_speed_3d(speed, 1, fade_value);
                if (!NO_EXTRA_BUFFER_SWAP) glutSwapBuffers();
                save_frame_counter(NSTEPS + i + 1);
            }
            
            if ((ROTATE_VIEW)&&(ROTATE_VIEW_WHILE_FADE)) 
            {
                viewpoint_schedule(NSTEPS - INITIAL_TIME);
                reset_view = 1;
            }
            draw_wave_sphere(1, phi, psi, xy_in, wave, wsphere, ZPLOT_B, CPLOT_B, COLOR_PALETTE_B, 0, 1.0, 1);
            if (DRAW_COLOR_SCHEME) draw_color_bar_palette(CPLOT_B, COLORBAR_RANGE_B, COLOR_PALETTE_B, CIRC_COLORBAR_B, 0, 1.0); 
            if (PRINT_SPEED) print_speed_3d(speed, 0, 1.0);
            glutSwapBuffers();
            
            if (!FADE) for (i=0; i<END_FRAMES; i++) save_frame_counter(NSTEPS + MID_FRAMES + 1 + counter + i);
            else for (i=0; i<END_FRAMES; i++) 
            {
                if ((ROTATE_VIEW)&&(ROTATE_VIEW_WHILE_FADE)) 
                {
                    viewpoint_schedule(NSTEPS - INITIAL_TIME + i);
                    reset_view = 1;
                }
                fade_value = 1.0 - (double)i/(double)END_FRAMES;
                draw_wave_sphere(1, phi, psi, xy_in, wave, wsphere, ZPLOT_B, CPLOT_B, COLOR_PALETTE_B, 1, fade_value, 1);
                if (DRAW_COLOR_SCHEME) draw_color_bar_palette(CPLOT_B, COLORBAR_RANGE_B, COLOR_PALETTE_B, CIRC_COLORBAR_B, 1, fade_value);   
                if (PRINT_SPEED) print_speed_3d(speed, 1, fade_value);
                glutSwapBuffers();
                save_frame_counter(NSTEPS + MID_FRAMES + 1 + counter + i);
            }
        }
        else
        {
            if (!FADE) for (i=0; i<END_FRAMES; i++) save_frame_counter(NSTEPS + MID_FRAMES + 1 + counter + i);
            else for (i=0; i<END_FRAMES; i++) 
            {
                if ((ROTATE_VIEW)&&(ROTATE_VIEW_WHILE_FADE))
                {
                    viewpoint_schedule(NSTEPS - INITIAL_TIME + i);
                    reset_view = 1;
                }
                fade_value = 1.0 - (double)i/(double)END_FRAMES;
                draw_wave_sphere(0, phi, psi, xy_in, wave, wsphere, ZPLOT, CPLOT, COLOR_PALETTE, 1, fade_value, 1);
                if (DRAW_COLOR_SCHEME) draw_color_bar_palette(CPLOT, COLORBAR_RANGE, COLOR_PALETTE, CIRC_COLORBAR, 1, fade_value); 
                if (PRINT_SPEED) print_speed_3d(speed, 1, fade_value);
                glutSwapBuffers();
                save_frame_counter(NSTEPS + 1 + counter + i);
            }
        }
        
        s = system("mv wave*.tif tif_wave/");
    }
    
    free(xy_in);
    free(phi);
    free(psi);
    free(tmp);
//     free(total_energy);
    free(color_scale);
    free(tc);
    free(tcc);
    free(tgamma);
    
    free(wave);
    free(laplace); 
    free(laplace1); 
    free(laplace2); 

    
    if (SAVE_TIME_SERIES)
    {
        fclose(time_series_left);
        fclose(time_series_right);
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
    glutCreateWindow("Wave equation in a planar domain");

    if (PLOT_2D) init_sphere_2D(); 
    else init_3d();

    glutDisplayFunc(display);

    glutMainLoop();

    return 0;
}

