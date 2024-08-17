/*********************************************************************************/
/*                                                                               */
/*  Animation of wave equation in a planar domain                                */
/*                                                                               */
/*  N. Berglund, december 2012, may  2021                                        */
/*                                                                               */
/*  UPDATE 24/04: distinction between damping and "elasticity" parameters        */
/*  UPDATE 27/04: new billiard shapes, bug in color scheme fixed                 */
/*  UPDATE 28/04: code made more efficient, with help of Marco Mancini           */
/*                                                                               */
/*  Feel free to reuse, but if doing so it would be nice to drop a               */
/*  line to nils.berglund@univ-orleans.fr - Thanks!                              */
/*                                                                               */
/*  compile with                                                                 */
/*  gcc -o wave_billiard wave_billiard.c                                         */
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
// #define NX 1920          /* number of grid points on x axis */
// #define NY 1000          /* number of grid points on y axis */
// #define YMID 500        /* mid point of display */
#define NX 3840          /* number of grid points on x axis */
#define NY 2300          /* number of grid points on y axis */
#define YMID 1150        /* mid point of display */
#define XMIN -2.0
#define XMAX 2.0	/* x interval */
#define YMIN -1.197916667
#define YMAX 1.197916667	/* y interval for 9/16 aspect ratio */

// #define WINWIDTH 	1280  /* window width */
// #define WINHEIGHT 	720   /* window height */
// #define NX 1280          /* number of grid points on x axis */
// #define NY 720          /* number of grid points on y axis */
// #define YMID 360        /* mid point of display */
// #define XMIN -2.0
// #define XMAX 2.0	/* x interval */
// #define YMIN -1.125
// #define YMAX 1.125	/* y interval for 9/16 aspect ratio */

#define HIGHRES 1       /* set to 1 if resolution of grid is double that of displayed image */

#define JULIA_SCALE 1.0 /* scaling for Julia sets */

/* Choice of the billiard table */

#define B_DOMAIN 63      /* choice of domain shape, see list in global_pdes.c */
#define B_DOMAIN_B 63    /* choice of domain shape, see list in global_pdes.c */

#define CIRCLE_PATTERN 13      /* pattern of circles, see list in global_pdes.c */
#define CIRCLE_PATTERN_B 13     /* pattern of circles, see list in global_pdes.c */

#define P_PERCOL 0.25       /* probability of having a circle in C_RAND_PERCOL arrangement */
#define NPOISSON 300        /* number of points for Poisson C_RAND_POISSON arrangement */
#define PDISC_FACTOR 3.25   /* controls density of Poisson disc process (default: 3.25) */
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

// #define NSTEPS 500      /* number of frames of movie */
#define NSTEPS 2600      /* number of frames of movie */
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
#define VSHIFT_AMPLITUDE 0.0   /* additional shift for wave amplitude */
#define VSCALE_AMPLITUDE 0.5   /* additional scaling factor for wave amplitude */
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
#define MU_B 1.0           /* parameter controlling the dimensions of domain */
#define VERTICAL_WAVE_PROFILE 0 /* set to 1 to draw wave profile vertically */
#define DRAW_WAVE_TIMESERIES 0  /* set to 1 to draw a time series of the wave */
#define WALL_WIDTH 0.1      /* width of wall separating lenses */
#define RADIUS_FACTOR 0.3   /* controls inner radius for C_RING arrangements */
#define OSCIL_YMAX 0.35      /* defines oscillation range */
#define MESSAGE_LDASH 14         /* length of dash for Morse code message */
#define MESSAGE_LDOT 8           /* length of dot for Morse code message */
#define MESSAGE_LINTERVAL 54     /* length of interval between dashes/dots for Morse code message */
#define MESSAGE_LINTERLETTER 60  /* length of interval between letters for Morse code message */
#define MESSAGE_LSPACE 48        /* length of space for Morse code message */
#define MESSAGE_INITIAL_TIME 100 /* initial time before starting message for Morse code message */    
#define AVRG_E_FACTOR 0.95   /* controls time window size in P_AVERAGE_ENERGY scheme */
#define HORIZONTAL_WAVE_PROFILE 0 /* set to 1 to draw wave profile vertically */
#define WAVE_PROFILE_X 2.1      /* value of x to sample wave profile */
#define WAVE_PROFILE_Y -1.0      /* value of y to sample wave profile */
#define PROFILE_AT_BOTTOM 1     /* draw wave profile at bottom instead of top */
#define AVERAGE_WAVE_PROFILE 1  /* set to 1 to draw time-average of wave profile squared*/
#define TIMESERIES_NVALUES 400  /* number of values plotted in time series */
#define DRAW_WAVE_SOURCE 0      /* set to 1 to draw source of wave at (wave_source_x, wave_source_y) */
#define HRES 1          /* dummy, only used by rde.c */
#define INITIAL_SHIFT 20.0          /* time shift of initial wave packet (in oscillation periods) */
#define WAVE_PACKET_SHIFT 200.0     /* time shift between wave packets (in oscillation periods) */
#define FADE_IN_OBSTACLE 0      /* set to 1 to fade color inside obstacles */
#define SHADE_2D 0       /* set to 1 to add pseudo-3d shading effect */ 
#define SHADE_SCALE_2D 0.05  /* lower value increases sensitivity of shading */
#define N_SOURCES 1                     /* number of sources, for option draw_sources */
#define XYIN_INITIALISED (B_DOMAIN == D_IMAGE)
double light[2] = {0.40824829, 0.816496581};   /* location of light source for SHADE_2D option*/
/* end of constants only used by sub_wave and sub_maze */


#include "global_pdes.c"        /* constants and global variables */
#include "sub_maze.c"           /* support for generating mazes */
#include "sub_wave.c"           /* common functions for wave_billiard, heat and schrodinger */
#include "wave_common.c"        /* common functions for wave_billiard, wave_comparison, etc */
#include "sub_wave_comp.c"      /* some functions specific to wave_comparison */

double courant2, courantb2;  /* Courant parameters squared */

FILE *monitor_sources;

/*********************/
/* animation part    */
/*********************/

void evolve_wave_half(double *phi_in[NX], double *psi_in[NX], double *phi_out[NX],
                      short int *xy_in[NX], double *tcc[NX])
/* time step of field evolution */
/* phi is value of field at time t, psi at time t-1 */
{
    int i, j, iplus, iminus, jplus, jminus, jmid = NY/2;
    double delta, x, y, c, cc, gamma;
    static long time = 0;
    static double tc[NX][NY], tgamma[NX][NY];
    static short int first = 1;
    
    time++;
    
    /* initialize tables with wave speeds and dissipation */
    if (first)
    {
        for (i=0; i<NX; i++){
            for (j=0; j<NY; j++){
                if (xy_in[i][j] != 0)
                {
                    tc[i][j] = COURANT;
                    if (!VARIABLE_IOR) tcc[i][j] = courant2;
                    if (xy_in[i][j] == 1) tgamma[i][j] = GAMMA;
                    else tgamma[i][j] = GAMMAB;
                }
                else if (TWOSPEEDS)
                {
                    tc[i][j] = COURANTB;
                    if (!VARIABLE_IOR) tcc[i][j] = courantb2;
                    tgamma[i][j] = GAMMAB;
                }
            }
        }
        first = 0;
    }

    #pragma omp parallel for private(i,j,iplus,iminus,jplus,jminus,delta,x,y,c,cc,gamma)
    /* evolution in the bulk */
    for (i=1; i<NX-1; i++){
        for (j=1; j<jmid-1; j++){
            if ((TWOSPEEDS)||(xy_in[i][j] != 0)){
                x = phi_in[i][j];
		y = psi_in[i][j];
                
                /* discretized Laplacian */
                delta = phi_in[i+1][j] + phi_in[i-1][j] + phi_in[i][j+1] + phi_in[i][j-1] - 4.0*x;

                /* evolve phi */
                phi_out[i][j] = -y + 2*x + tcc[i][j]*delta - KAPPA*x - tgamma[i][j]*(x-y);
            }
        }
        for (j=jmid+1; j<NY-1; j++){
            if ((TWOSPEEDS)||(xy_in[i][j] != 0)){
                x = phi_in[i][j];
		y = psi_in[i][j];
                
                /* discretized Laplacian */
                delta = phi_in[i+1][j] + phi_in[i-1][j] + phi_in[i][j+1] + phi_in[i][j-1] - 4.0*x;

                /* evolve phi */
                phi_out[i][j] = -y + 2*x + tcc[i][j]*delta - KAPPA*x - tgamma[i][j]*(x-y);
            }
        }
    }
    
    /* left boundary */
    if (OSCILLATE_LEFT) {
        for (j=1; j<jmid-1; j++) phi_out[0][j] = AMPLITUDE*cos((double)time*OMEGA);
        for (j=jmid+1; j<NY-1; j++) phi_out[0][j] = AMPLITUDE*cos((double)time*OMEGA);
    }
    else for (j=1; j<NY-1; j++) if ((j!=jmid-1)&&(j!=jmid)) {
        if ((TWOSPEEDS)||(xy_in[0][j] != 0)){
            x = phi_in[0][j];
            y = psi_in[0][j];
                    
            switch (B_COND) {
                case (BC_DIRICHLET):
                {
                    delta = phi_in[1][j] + phi_in[0][j+1] + phi_in[0][j-1] - 3.0*x;
                    phi_out[0][j] = -y + 2*x + tcc[0][j]*delta - KAPPA*x - tgamma[0][j]*(x-y);
                    break;
                }
                case (BC_PERIODIC):
                {
                    delta = phi_in[1][j] + phi_in[NX-1][j] + phi_in[0][j+1] + phi_in[0][j-1] - 4.0*x;
                    phi_out[0][j] = -y + 2*x + tcc[0][j]*delta - KAPPA*x - tgamma[0][j]*(x-y);
                    break;
                }
                case (BC_ABSORBING):
                {
                    delta = phi_in[1][j] + phi_in[0][j+1] + phi_in[0][j-1] - 3.0*x;
                    phi_out[0][j] = x - tc[0][j]*(x - phi_in[1][j]) - KAPPA_SIDES*x - GAMMA_SIDES*(x-y);
                    break;
                }
                case (BC_VPER_HABS):
                {
                    delta = phi_in[1][j] + phi_in[0][j+1] + phi_in[0][j-1] - 3.0*x;
                    phi_out[0][j] = x - tc[0][j]*(x - phi_in[1][j]) - KAPPA_SIDES*x - GAMMA_SIDES*(x-y);
                    break;
                }
                case (BC_ABS_REFLECT):
                {
                    delta = phi_in[1][j] + phi_in[0][j+1] + phi_in[0][j-1] - 3.0*x;
                    phi_out[0][j] = x - tc[0][j]*(x - phi_in[1][j]) - KAPPA_SIDES*x - GAMMA_SIDES*(x-y);
                    break;
                }
            }
        }
    }
    
    /* right boundary */
    for (j=1; j<NY-1; j++) if ((j!=jmid-1)&&(j!=jmid)) {
        if ((TWOSPEEDS)||(xy_in[NX-1][j] != 0)){
            x = phi_in[NX-1][j];
            y = psi_in[NX-1][j];
                    
            switch (B_COND) {
                case (BC_DIRICHLET):
                {
                    delta = phi_in[NX-2][j] + phi_in[NX-1][j+1] + phi_in[NX-1][j-1] - 3.0*x;
                    phi_out[NX-1][j] = -y + 2*x + tcc[NX-1][j]*delta - KAPPA*x - tgamma[NX-1][j]*(x-y);
                    break;
                }
                case (BC_PERIODIC):
                {
                    delta = phi_in[NX-2][j] + phi_in[0][j] + phi_in[NX-1][j+1] + phi_in[NX-1][j-1] - 4.0*x;
                    phi_out[NX-1][j] = -y + 2*x + tcc[NX-1][j]*delta - KAPPA*x - tgamma[NX-1][j]*(x-y);
                    break;
                }
                case (BC_ABSORBING):
                {
                    delta = phi_in[NX-2][j] + phi_in[NX-1][j+1] + phi_in[NX-1][j-1] - 3.0*x;
                    phi_out[NX-1][j] = x - tc[NX-1][j]*(x - phi_in[NX-2][j]) - KAPPA_SIDES*x - GAMMA_SIDES*(x-y);
                    break;
                }
                case (BC_VPER_HABS):
                {
                    delta = phi_in[NX-2][j] + phi_in[NX-1][j+1] + phi_in[NX-1][j-1] - 3.0*x;
                    phi_out[NX-1][j] = x - tc[NX-1][j]*(x - phi_in[NX-2][j]) - KAPPA_SIDES*x - GAMMA_SIDES*(x-y);
                    break;
                }
                case (BC_ABS_REFLECT):
                {
                    delta = phi_in[NX-2][j] + phi_in[NX-1][j+1] + phi_in[NX-1][j-1] - 3.0*x;
                    phi_out[NX-1][j] = x - tc[NX-1][j]*(x - phi_in[NX-2][j]) - KAPPA_SIDES*x - GAMMA_SIDES*(x-y);
                    break;
                }
            }
        }
    }
    
    /* top mid boundary */
    for (i=0; i<NX; i++){
        if ((TWOSPEEDS)||(xy_in[i][jmid-1] != 0)){
            x = phi_in[i][jmid-1];
            y = psi_in[i][jmid-1];
                    
            switch (B_COND) {
                case (BC_DIRICHLET):
                {
                    iplus = i+1;   if (iplus == NX) iplus = NX-1;
                    iminus = i-1;  if (iminus == -1) iminus = 0;
                    
                    delta = phi_in[iplus][jmid-1] + phi_in[iminus][jmid-1] + phi_in[i][jmid-2] - 3.0*x;
                    phi_out[i][jmid-1] = -y + 2*x + tcc[i][jmid-1]*delta - KAPPA*x - tgamma[i][jmid-1]*(x-y);
                    break;
                }
                case (BC_PERIODIC):
                {
                    iplus = (i+1) % NX;
                    iminus = (i-1) % NX;    if (iminus < 0) iminus += NX;
                    
                    delta = phi_in[iplus][jmid-1] + phi_in[iminus][jmid-1] + phi_in[i][jmid-2] + phi_in[i][0] - 4.0*x;
                    phi_out[i][jmid-1] = -y + 2*x + tcc[i][jmid-1]*delta - KAPPA*x - tgamma[i][jmid-1]*(x-y);
                    break;
                }
                case (BC_ABSORBING):
                {
                    iplus = (i+1);   if (iplus == NX) iplus = NX-1;
                    iminus = (i-1);  if (iminus == -1) iminus = 0;
                    
                    delta = phi_in[iplus][jmid-1] + phi_in[iminus][jmid-1] + phi_in[i][jmid-2] - 3.0*x;
                    phi_out[i][jmid-1] = x - tc[i][jmid-1]*(x - phi_in[i][jmid-2]) - KAPPA_TOPBOT*x - GAMMA_TOPBOT*(x-y);
                    break;
                }
                case (BC_VPER_HABS):
                {
                    iplus = (i+1);   if (iplus == NX) iplus = NX-1;
                    iminus = (i-1);  if (iminus == -1) iminus = 0;

                    delta = phi_in[iplus][jmid-1] + phi_in[iminus][jmid-1] + phi_in[i][jmid-2] + phi_in[i][0] - 4.0*x;
                    if (i==0) phi_out[0][jmid-1] = x - tc[0][jmid-1]*(x - phi_in[1][jmid-1]) - KAPPA_SIDES*x - GAMMA_SIDES*(x-y);
                    else phi_out[i][jmid-1] = -y + 2*x + tcc[i][jmid-1]*delta - KAPPA*x - tgamma[i][jmid-1]*(x-y);
                   break;
                }
                case (BC_ABS_REFLECT):
                {
                    iplus = i+1;   if (iplus == NX) iplus = NX-1;
                    iminus = i-1;  if (iminus == -1) iminus = 0;
                    
                    delta = phi_in[iplus][jmid-1] + phi_in[iminus][jmid-1] + phi_in[i][jmid-2] - 3.0*x;
                    phi_out[i][jmid-1] = -y + 2*x + tcc[i][jmid-1]*delta - KAPPA*x - tgamma[i][jmid-1]*(x-y);
                    break;
                }
            }
        }
    }
    
    /* bottom boundary */
    for (i=0; i<NX; i++){
        if ((TWOSPEEDS)||(xy_in[i][0] != 0)){
            x = phi_in[i][0];
            y = psi_in[i][0];
                    
            switch (B_COND) {
                case (BC_DIRICHLET):
                {
                    iplus = i+1;   if (iplus == NX) iplus = NX-1;
                    iminus = i-1;  if (iminus == -1) iminus = 0;
                    
                    delta = phi_in[iplus][0] + phi_in[iminus][0] + phi_in[i][1] - 3.0*x;
                    phi_out[i][0] = -y + 2*x + tcc[i][0]*delta - KAPPA*x - tgamma[i][0]*(x-y);
                    break;
                }
                case (BC_PERIODIC):
                {
                    iplus = (i+1) % NX;
                    iminus = (i-1) % NX;    if (iminus < 0) iminus += NX;
                    
                    delta = phi_in[iplus][0] + phi_in[iminus][0] + phi_in[i][1] + phi_in[i][jmid-1] - 4.0*x;
                    phi_out[i][0] = -y + 2*x + tcc[i][0]*delta - KAPPA*x - tgamma[i][0]*(x-y);
                    break;
                }
                case (BC_ABSORBING):
                {
                    iplus = (i+1);   if (iplus == NX) iplus = NX-1;
                    iminus = (i-1);  if (iminus == -1) iminus = 0;
                    
                    delta = phi_in[iplus][0] + phi_in[iminus][0] + phi_in[i][1] - 3.0*x;
                    phi_out[i][0] = x - tc[i][0]*(x - phi_in[i][1]) - KAPPA_TOPBOT*x - GAMMA_TOPBOT*(x-y);
                    break;
                }
                case (BC_VPER_HABS):
                {
                    iplus = (i+1);   if (iplus == NX) iplus = NX-1;
                    iminus = (i-1);  if (iminus == -1) iminus = 0;

                    delta = phi_in[iplus][0] + phi_in[iminus][0] + phi_in[i][1] + phi_in[i][jmid-1] - 4.0*x;
                    if (i==0) phi_out[0][0] = x - tc[0][0]*(x - phi_in[1][0]) - KAPPA_SIDES*x - GAMMA_SIDES*(x-y);
                    else phi_out[i][0] = -y + 2*x + tcc[i][0]*delta - KAPPA*x - tgamma[i][0]*(x-y);
                    break;
                }
                case (BC_ABS_REFLECT):
                {
                    iplus = (i+1);   if (iplus == NX) iplus = NX-1;
                    iminus = (i-1);  if (iminus == -1) iminus = 0;
                    
                    delta = phi_in[iplus][0] + phi_in[iminus][0] + phi_in[i][1] - 3.0*x;
                    phi_out[i][0] = x - tc[i][0]*(x - phi_in[i][1]) - KAPPA_TOPBOT*x - GAMMA_TOPBOT*(x-y);
                    break;
                }
            }
        }
    }
    
    /* top boundary */
    for (i=0; i<NX; i++){
        if ((TWOSPEEDS)||(xy_in[i][NY-1] != 0)){
            x = phi_in[i][NY-1];
            y = psi_in[i][NY-1];
                    
            switch (B_COND) {
                case (BC_DIRICHLET):
                {
                    iplus = i+1;   if (iplus == NX) iplus = NX-1;
                    iminus = i-1;  if (iminus == -1) iminus = 0;
                    
                    delta = phi_in[iplus][NY-1] + phi_in[iminus][NY-1] + phi_in[i][NY-2] - 3.0*x;
                    phi_out[i][NY-1] = -y + 2*x + tcc[i][NY-1]*delta - KAPPA*x - tgamma[i][NY-1]*(x-y);
                    break;
                }
                case (BC_PERIODIC):
                {
                    iplus = (i+1) % NX;
                    iminus = (i-1) % NX;    if (iminus < 0) iminus += NX;
                    
                    delta = phi_in[iplus][NY-1] + phi_in[iminus][NY-1] + phi_in[i][NY-2] + phi_in[i][jmid] - 4.0*x;
                    phi_out[i][NY-1] = -y + 2*x + tcc[i][NY-1]*delta - KAPPA*x - tgamma[i][NY-1]*(x-y);
                    break;
                }
                case (BC_ABSORBING):
                {
                    iplus = (i+1);   if (iplus == NX) iplus = NX-1;
                    iminus = (i-1);  if (iminus == -1) iminus = 0;
                    
                    delta = phi_in[iplus][NY-1] + phi_in[iminus][NY-1] + phi_in[i][NY-2] - 3.0*x;
                    phi_out[i][NY-1] = x - tc[i][NY-1]*(x - phi_in[i][NY-2]) - KAPPA_TOPBOT*x - GAMMA_TOPBOT*(x-y);
                    break;
                }
                case (BC_VPER_HABS):
                {
                    iplus = (i+1);   if (iplus == NX) iplus = NX-1;
                    iminus = (i-1);  if (iminus == -1) iminus = 0;

                    delta = phi_in[iplus][NY-1] + phi_in[iminus][NY-1] + phi_in[i][NY-2] + phi_in[i][jmid] - 4.0*x;
                    if (i==0) phi_out[0][NY-1] = x - tc[0][NY-1]*(x - phi_in[1][NY-1]) - KAPPA_SIDES*x - GAMMA_SIDES*(x-y);
                    else phi_out[i][NY-1] = -y + 2*x + tcc[i][NY-1]*delta - KAPPA*x - tgamma[i][NY-1]*(x-y);
                   break;
                }
                case (BC_ABS_REFLECT):
                {
                    iplus = (i+1);   if (iplus == NX) iplus = NX-1;
                    iminus = (i-1);  if (iminus == -1) iminus = 0;
                    
                    delta = phi_in[iplus][NY-1] + phi_in[iminus][NY-1] + phi_in[i][NY-2] - 3.0*x;
                    phi_out[i][NY-1] = x - tc[i][NY-1]*(x - phi_in[i][NY-2]) - KAPPA_TOPBOT*x - GAMMA_TOPBOT*(x-y);
                    break;
                }
            }
        }
    }
    
    /* bottom mid boundary */
    for (i=0; i<NX; i++){
        if ((TWOSPEEDS)||(xy_in[i][jmid] != 0)){
            x = phi_in[i][jmid];
            y = psi_in[i][jmid];
                    
            switch (B_COND) {
                case (BC_DIRICHLET):
                {
                    iplus = i+1;   if (iplus == NX) iplus = NX-1;
                    iminus = i-1;  if (iminus == -1) iminus = 0;
                    
                    delta = phi_in[iplus][jmid] + phi_in[iminus][jmid] + phi_in[i][jmid+1] - 3.0*x;
                    phi_out[i][jmid] = -y + 2*x + tcc[i][jmid]*delta - KAPPA*x - tgamma[i][jmid]*(x-y);
                    break;
                }
                case (BC_PERIODIC):
                {
                    iplus = (i+1) % NX;
                    iminus = (i-1) % NX;    if (iminus < 0) iminus += NX;
                    
                    delta = phi_in[iplus][jmid] + phi_in[iminus][jmid] + phi_in[i][jmid+1] + phi_in[i][NY-1] - 4.0*x;
                    phi_out[i][jmid] = -y + 2*x + tcc[i][jmid]*delta - KAPPA*x - tgamma[i][jmid]*(x-y);
                    break;
                }
                case (BC_ABSORBING):
                {
                    iplus = (i+1);   if (iplus == NX) iplus = NX-1;
                    iminus = (i-1);  if (iminus == -1) iminus = 0;
                    
                    delta = phi_in[iplus][jmid] + phi_in[iminus][jmid] + phi_in[i][jmid+1] - 3.0*x;
                    phi_out[i][jmid] = x - tc[i][jmid]*(x - phi_in[i][1]) - KAPPA_TOPBOT*x - GAMMA_TOPBOT*(x-y);
                    break;
                }
                case (BC_VPER_HABS):
                {
                    iplus = (i+1);   if (iplus == NX) iplus = NX-1;
                    iminus = (i-1);  if (iminus == -1) iminus = 0;

                    delta = phi_in[iplus][jmid] + phi_in[iminus][jmid] + phi_in[i][jmid+1] + phi_in[i][NY-1] - 4.0*x;
                    if (i==0) phi_out[0][jmid] = x - tc[0][jmid]*(x - phi_in[1][jmid]) - KAPPA_SIDES*x - GAMMA_SIDES*(x-y);
                    else phi_out[i][jmid] = -y + 2*x + tcc[i][jmid]*delta - KAPPA*x - tgamma[i][jmid]*(x-y);
                    break;
                }
                case (BC_ABS_REFLECT):
                {
                    iplus = i+1;   if (iplus == NX) iplus = NX-1;
                    iminus = i-1;  if (iminus == -1) iminus = 0;
                    
                    delta = phi_in[iplus][jmid] + phi_in[iminus][jmid] + phi_in[i][jmid+1] - 3.0*x;
                    phi_out[i][jmid] = -y + 2*x + tcc[i][jmid]*delta - KAPPA*x - tgamma[i][jmid]*(x-y);
                    break;
                }
            }
        }
    }
    
    /* add oscillating boundary condition on the left corners */
    if ((i == 0)&&(OSCILLATE_LEFT))
    {
        phi_out[i][0] = AMPLITUDE*cos((double)time*OMEGA);
        phi_out[i][jmid-1] = AMPLITUDE*cos((double)time*OMEGA);
        phi_out[i][jmid] = AMPLITUDE*cos((double)time*OMEGA);
        phi_out[i][NY-1] = AMPLITUDE*cos((double)time*OMEGA);
    }
    
    /* for debugging purposes/if there is a risk of blow-up */
    if (FLOOR) for (i=0; i<NX; i++){
        for (j=0; j<NY; j++){
            if (xy_in[i][j] != 0) 
            {
                if (phi_out[i][j] > VMAX) phi_out[i][j] = VMAX;
                if (phi_out[i][j] < -VMAX) phi_out[i][j] = -VMAX;
            }
        }
    }
//     printf("phi(0,0) = %.3lg, psi(0,0) = %.3lg\n", phi[NX/2][NY/2], psi[NX/2][NY/2]);
}


void evolve_wave(double *phi[NX], double *psi[NX], double *tmp[NX], short int *xy_in[NX], double *tcc_table[NX])
/* time step of field evolution */
/* phi is value of field at time t, psi at time t-1 */
{
    evolve_wave_half(phi, psi, tmp, xy_in, tcc_table);
    evolve_wave_half(tmp, phi, psi, xy_in, tcc_table);
    evolve_wave_half(psi, tmp, phi, xy_in, tcc_table);
}


// void draw_color_bar(int plot, double range)
// {
//     if (ROTATE_COLOR_SCHEME) draw_color_scheme(-1.0, -0.8, XMAX - 0.1, -1.0, plot, -range, range);
//     else draw_color_scheme(XMAX - 0.3, YMIN + 0.1, XMAX - 0.1, YMAX - 0.1, plot, -range, range);
// }

void draw_color_bar_palette(int plot, double range, int palette, int fade, double fade_value)
{
    double width = 0.14;
//     double width = 0.2;
    
    if (ROTATE_COLOR_SCHEME) 
        draw_color_scheme_palette_fade(-1.0, -0.8, XMAX - 0.1, -1.0, plot, -range, range, palette, fade, fade_value);
    else 
        draw_color_scheme_palette_fade(XMAX - 1.5*width, YMIN + 0.1, XMAX - 0.5*width, YMAX - 0.1, plot, -range, range, palette, fade, fade_value);
}

void animation()
{
    double time, scale, energies[6], top_energy, bottom_energy, omega, angle, fade_value, sign = 1.0, ior_angle = 0.0;
    double *phi[NX], *psi[NX], *tmp[NX], *total_energy[NX], *tcc_table[NX];
    short int *xy_in[NX];
    int i, j, s, counter = 0, k, first_source = 1, fade, resol = HIGHRES + 1;
    t_wave_source wave_source[NSOURCES];
    
    monitor_sources = fopen("monitor_sources.dat", "w");

    /* Since NX and NY are big, it seemed wiser to use some memory allocation here */
    for (i=0; i<NX; i++)
    {
        phi[i] = (double *)malloc(NY*sizeof(double));
        psi[i] = (double *)malloc(NY*sizeof(double));
        tmp[i] = (double *)malloc(NY*sizeof(double));
        total_energy[i] = (double *)malloc(NY*sizeof(double));
        xy_in[i] = (short int *)malloc(NY*sizeof(short int));
        tcc_table[i] = (double *)malloc(NX*sizeof(double));
    }
    
    /* initialise positions and radii of circles */
    printf("initializing circle configuration\n");
    if ((B_DOMAIN == D_CIRCLES)||(B_DOMAIN_B == D_CIRCLES)) init_circle_config_comp(circles);
    if ((B_DOMAIN == D_POLYGONS)|(B_DOMAIN_B == D_POLYGONS)) init_polygon_config_comp(polygons);
//     for (i=0; i<ncircles; i++) printf("polygon %i at (%.3f, %.3f) radius %.3f\n", i, polygons[i].xc, polygons[i].yc, polygons[i].radius);
    
    /* initialise polyline for von Koch and similar domains */
    npolyline = init_polyline(MDEPTH, polyline);
    for (i=0; i<npolyline; i++) printf("vertex %i: (%.3f, %.3f)\n", i, polyline[i].x, polyline[i].y);

    /* initialize total energy table */
    if ((PLOT == P_MEAN_ENERGY)||(PLOT_B == P_MEAN_ENERGY)||(PLOT == P_LOG_MEAN_ENERGY)||(PLOT_B == P_LOG_MEAN_ENERGY))
        for (i=0; i<NX; i++)
            for (j=0; j<NY; j++) 
                total_energy[i][j] = 0.0;

    courant2 = COURANT*COURANT;
    courantb2 = COURANTB*COURANTB;

    /* initialize wave with a drop at one point, zero elsewhere */
    init_wave_flat_comp(phi, psi, xy_in);
//     int_planar_wave_comp(XMIN + 0.015, 0.0, phi, psi, xy_in);
//     int_planar_wave_comp(XMIN + 0.5, 0.0, phi, psi, xy_in);
    printf("initializing wave\n");
//     int_planar_wave_comp(XMIN + 0.1, 0.0, phi, psi, xy_in);
//     int_planar_wave_comp(XMIN + 1.0, 0.0, phi, psi, xy_in);
//     init_wave(-1.5, 0.0, phi, psi, xy_in);
//     init_wave(0.0, 0.0, phi, psi, xy_in);

    /* add a drop at another point */
//     add_drop_to_wave(1.0, 0.7, 0.0, phi, psi);
//     add_drop_to_wave(1.0, -0.7, 0.0, phi, psi);
//     add_drop_to_wave(1.0, 0.0, -0.7, phi, psi);

    if (VARIABLE_IOR) init_ior_2d_comp(xy_in, tcc_table, ior_angle);
    
    /* initialize energies */
    if (COMPUTE_ENERGIES) 
    {
        printf("computing energies\n");
        compute_energy_tblr(phi, psi, xy_in, energies);
        top_energy = energies[0] + energies[1] + energies[2]; 
        bottom_energy = energies[3] + energies[4] + energies[5];
        printf("computed energies\n");
    }

    blank();
    glColor3f(0.0, 0.0, 0.0);
    printf("drawing wave\n");
    draw_wave_comp(phi, psi, xy_in, 1.0, 0, PLOT);

    printf("drawing billiard\n");
    draw_billiard_comp(0, 1.0);

    if (DRAW_COLOR_SCHEME) draw_color_bar_palette(PLOT, COLORBAR_RANGE, COLOR_PALETTE, fade, fade_value);

    glutSwapBuffers();



    sleep(SLEEP1);

    for (i=0; i<=INITIAL_TIME + NSTEPS; i++)
    {
	//printf("%d\n",i);
        /* compute the variance of the field to adjust color scheme */
        /* the color depends on the field divided by sqrt(1 + variance) */
        if (SCALE)
        {
            scale = sqrt(1.0 + compute_variance(phi,psi, xy_in));
//             printf("Scaling factor: %5lg\n", scale);
        }
        else scale = 1.0;

//         draw_wave_comp(phi, psi, xy_in, scale, i, PLOT);
        draw_wave_comp_highres_palette(resol, phi, psi, total_energy, xy_in, scale, i, PLOT, COLOR_PALETTE, 0, 1.0);
        
       for (j=0; j<NVID; j++) 
        {
            evolve_wave(phi, psi, tmp, xy_in, tcc_table);
//             if (i % 10 == 9) oscillate_linear_wave(0.2*scale, 0.15*(double)(i*NVID + j), -1.5, YMIN, -1.5, YMAX, phi, psi);
        }
        
        /* add oscillating waves */
        if ((ADD_OSCILLATING_SOURCE)&&(i%OSCILLATING_SOURCE_PERIOD == 1))
        {
            if (ALTERNATE_OSCILLATING_SOURCE) sign = -sign;
            add_circular_wave(sign, -1.75, 0.16, phi, psi, xy_in);
            add_circular_wave(sign, -1.75, -0.16, phi, psi, xy_in);
            
//             if (first_source) for (k=0; k<NSOURCES; k++)
//             {
//                 omega = DPI/(double)NSOURCES;
//                 wave_source[k].xc = 0.05*cos(((double)k + 0.5)*omega);;
//                 wave_source[k].yc = 0.05*sin(((double)k + 0.5)*omega);
//                 if (wave_source[k].yc < 0.0)
//                 {
//                     wave_source[k].xc *= 0.5;
//                     wave_source[k].yc *= 0.5;
//                 }
// //                 wave_source[k].phase = 0.99 - 1.4*sin(0.35*(1.0 + wave_source[k].xc/0.1));
// //                 wave_source[k].phase = 0.99 - 1.4*sin(0.7*(1.0 + wave_source[k].xc/0.1));
// //                 wave_source[k].phase = 0.99 - 1.4*sin(0.7*(1.0 + wave_source[k].xc/0.05));
//                 wave_source[k].phase = 0.99 - 1.4*sin(0.35*(1.0 + wave_source[k].xc/0.05));
//                 wave_source[k].amp = 1.0;
// //                 if (wave_source[k].phase) wave_source[k].sign = 1;
// //                 else wave_source[k].sign = -1;
//                 wave_source[k].sign = 1;
//                 
//                 first_source = 0;
//             }
//             
//             fprintf(monitor_sources, "Frame %i\n\n", i);
//             
//             for (k=0; k<NSOURCES; k++) /*if (wave_source[k].xc > 0.0) */
//             {
//                 wave_source[k].phase += 0.06;
// //                 if (k==1) printf("x = %.3lg, phase = %.3lg\n", wave_source[k].xc, wave_source[k].phase);
// //                 fprintf(monitor_sources, "x = %.3lg, y = %.3lg, phase = %.3lg\n", wave_source[k].xc, wave_source[k].yc, wave_source[k].phase);
//                 if (wave_source[k].phase > 1.0)
//                 {
//                     add_circular_wave_comp((double)wave_source[k].sign*wave_source[k].amp, wave_source[k].xc, wave_source[k].yc, phi, psi, xy_in, (wave_source[k].yc > 0));
//                     fprintf(monitor_sources, "x = %.3lg, y = %.3lg, phase = %.3lg\n", wave_source[k].xc, wave_source[k].yc, wave_source[k].phase);
//                     printf("Adding wave at (%.2lg, %.2lg)\n", wave_source[k].xc, wave_source[k].yc);
//                     fprintf(monitor_sources, "Adding wave at (%.2lg, %.2lg)\n", wave_source[k].xc, wave_source[k].yc);
//                     wave_source[k].phase -= 1.0;
//                     wave_source[k].sign *= -1;
//                 }
//             }
        }

        draw_billiard_comp(0, 1.0);
        
	if (COMPUTE_ENERGIES) 
        {
            compute_energy_tblr(phi, psi, xy_in, energies);
            if (i < INITIAL_TIME)
            {
                top_energy = energies[0] + energies[1] + energies[2]; 
                bottom_energy = energies[3] + energies[4] + energies[5];
            }
            print_energies(energies, top_energy, bottom_energy);
        }
        
        
        
        if (DRAW_COLOR_SCHEME) draw_color_bar_palette(PLOT, COLORBAR_RANGE, COLOR_PALETTE, fade, fade_value);

        if (!((NO_EXTRA_BUFFER_SWAP)&&(MOVIE))) glutSwapBuffers();

	if (MOVIE)
        {
            if (i >= INITIAL_TIME) save_frame();
            
            if ((i >= INITIAL_TIME)&&(DOUBLE_MOVIE))
            {
//                 save_frame();
                if ((TIME_LAPSE)&&((i - INITIAL_TIME)%TIME_LAPSE_FACTOR == 0))
                {
                    save_frame_counter(NSTEPS + END_FRAMES + (i - INITIAL_TIME)/TIME_LAPSE_FACTOR);
                    if (DRAW_COLOR_SCHEME) draw_color_bar_palette(PLOT_B, COLORBAR_RANGE_B, COLOR_PALETTE_B, 0, 1.0); 
                    counter++;
                }
                else /*if (DOUBLE_MOVIE)*/
                {
//                     draw_wave_comp(phi, psi, xy_in, scale, i, PLOT_B);
                    draw_wave_comp_highres_palette(resol, phi, psi, total_energy, xy_in, scale, i, PLOT_B, COLOR_PALETTE_B, 0, 1.0);
                    draw_billiard_comp(0, 1.0);
                    if (DRAW_COLOR_SCHEME) draw_color_bar_palette(PLOT_B, COLORBAR_RANGE_B, COLOR_PALETTE_B, 0, 1.0);
                    glutSwapBuffers();
                    save_frame_counter(NSTEPS + MID_FRAMES + 1 + counter);
                    counter++;
                }
            }
            else if (NO_EXTRA_BUFFER_SWAP) glutSwapBuffers();
            
//             else printf("Initial phase time %i of %i\n", i, INITIAL_TIME);
            
            

            /* it seems that saving too many files too fast can cause trouble with the file system */
            /* so this is to make a pause from time to time - parameter PAUSE may need adjusting   */
            if (i % PAUSE == PAUSE - 1)
            {
                printf("Making a short pause\n");
                sleep(PSLEEP);
                s = system("mv wave*.tif tif_wave/");
            }
        }

    }

    if (MOVIE) 
    {
        if (DOUBLE_MOVIE) 
        {
//             draw_wave_comp(phi, psi, xy_in, scale, i, PLOT);
            draw_wave_comp_highres_palette(resol, phi, psi, total_energy, xy_in, scale, NSTEPS, PLOT, COLOR_PALETTE, 0, 1.0);
            draw_billiard_comp(0, 1.0);
            if (DRAW_COLOR_SCHEME) draw_color_bar_palette(PLOT, COLORBAR_RANGE, COLOR_PALETTE, 0, 1.0);   
            glutSwapBuffers();
        }
        if (!FADE) for (i=0; i<MID_FRAMES; i++) save_frame();
        else for (i=0; i<MID_FRAMES; i++) 
        {
            fade_value = 1.0 - (double)i/(double)MID_FRAMES;
            draw_wave_comp_highres_palette(resol, phi, psi, total_energy, xy_in, scale, NSTEPS, PLOT, COLOR_PALETTE, 1, fade_value);
            draw_billiard_comp(1, fade_value);
            if (DRAW_COLOR_SCHEME) draw_color_bar_palette(PLOT, COLORBAR_RANGE, COLOR_PALETTE, 1, fade_value); 
            if (!NO_EXTRA_BUFFER_SWAP) glutSwapBuffers();
            save_frame_counter(NSTEPS + i + 1);
        }
        if (DOUBLE_MOVIE) 
        {
//             draw_wave_comp(phi, psi, xy_in, scale, i, PLOT_B);
            draw_wave_comp_highres_palette(resol, phi, psi, total_energy, xy_in, scale, NSTEPS, PLOT_B, COLOR_PALETTE_B, 0, 1.0);
            draw_billiard_comp(0, 1.0);
            if (DRAW_COLOR_SCHEME) draw_color_bar_palette(PLOT_B, COLORBAR_RANGE_B, COLOR_PALETTE_B, 0, 1.0); 
            glutSwapBuffers();
            
            if (!FADE) for (i=0; i<END_FRAMES; i++) save_frame_counter(NSTEPS + MID_FRAMES + 1 + counter + i);
            else for (i=0; i<END_FRAMES; i++) 
            {
                fade_value = 1.0 - (double)i/(double)END_FRAMES;
                draw_wave_comp_highres_palette(resol, phi, psi, total_energy, xy_in, scale, NSTEPS, PLOT_B, COLOR_PALETTE_B, 1, fade_value);
                draw_billiard_comp(1, fade_value);
                if (DRAW_COLOR_SCHEME) draw_color_bar_palette(PLOT_B, COLORBAR_RANGE_B, COLOR_PALETTE_B, 1, fade_value); 
                glutSwapBuffers();
                save_frame_counter(NSTEPS + MID_FRAMES + 1 + counter + i);
            }
        }
        
        
        if (TIME_LAPSE) for (i=0; i<END_FRAMES; i++) save_frame_counter(NSTEPS + END_FRAMES + NSTEPS/TIME_LAPSE_FACTOR + i);
        s = system("mv wave*.tif tif_wave/");
    }
    for (i=0; i<NX; i++)
    {
        free(phi[i]);
        free(psi[i]);
        free(tmp[i]);
        free(total_energy[i]);
        free(xy_in[i]);
        free(tcc_table[i]);
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

    init();

    glutDisplayFunc(display);

    glutMainLoop();

    return 0;
}
