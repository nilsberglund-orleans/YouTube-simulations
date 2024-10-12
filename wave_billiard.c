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
/**********************************************************************************/

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

#define VARIABLE_IOR 0      /* set to 1 for a variable index of refraction */
#define IOR 191             /* choice of index of refraction, see list in global_pdes.c */
#define IOR_TOTAL_TURNS 1.5 /* total angle of rotation for IOR_PERIODIC_WELLS_ROTATING */
#define MANDEL_IOR_SCALE -0.05   /* parameter controlling dependence of IoR on Mandelbrot escape speed */

/* General geometrical parameters */

#define WINWIDTH 	1920  /* window width */
#define WINHEIGHT 	1150  /* window height */
#define NX 3840          /* number of grid points on x axis */
#define NY 2300          /* number of grid points on y axis */

#define XMIN -2.0
#define XMAX 2.0	/* x interval */
#define YMIN -0.997916667
#define YMAX 1.397916667	/* y interval for 9/16 aspect ratio */

#define HIGHRES 1       /* set to 1 if resolution of grid is double that of displayed image */
#define HRES 1          /* dummy, only used by rde.c */

#define JULIA_SCALE 1.0 /* scaling for Julia sets */

/* Choice of the billiard table */

#define B_DOMAIN 792        /* choice of domain shape, see list in global_pdes.c */

#define CIRCLE_PATTERN 202   /* pattern of circles or polygons, see list in global_pdes.c */

#define COMPARISON 0        /* set to 1 to compare two different patterns (beta) */
#define B_DOMAIN_B 20      /* second domain shape, for comparisons */
#define CIRCLE_PATTERN_B 0  /* second pattern of circles or polygons */

#define P_PERCOL 0.15       /* probability of having a circle in C_RAND_PERCOL arrangement */
#define NPOISSON 1000       /* number of points for Poisson C_RAND_POISSON arrangement */
#define PDISC_FACTOR 3.5    /* controls density of Poisson disc process (default: 3.25) */
#define RANDOM_POLY_ANGLE 1 /* set to 1 to randomize angle of polygons */

#define LAMBDA 0.85	    /* parameter controlling the dimensions of domain */
#define MU 0.18             /* parameter controlling the dimensions of domain */
#define MU_B 1.0            /* parameter controlling the dimensions of domain */
#define NPOLY 3            /* number of sides of polygon */
#define APOLY -0.17           /* angle by which to turn polygon, in units of Pi/2 */ 
#define MDEPTH 6            /* depth of computation of Menger gasket */
#define MRATIO 3            /* ratio defining Menger gasket */
#define MANDELLEVEL 1000    /* iteration level for Mandelbrot set */
#define MANDELLIMIT 10.0    /* limit value for approximation of Mandelbrot set */
#define FOCI 1              /* set to 1 to draw focal points of ellipse */
#define NGRIDX 60           /* number of grid point for grid of disks */
#define NGRIDY 25           /* number of grid point for grid of disks */
#define WALL_WIDTH 0.05     /* width of wall separating lenses */
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

#define TWOSPEEDS 0         /* set to 1 to replace hardcore boundary by medium with different speed */
#define OSCILLATE_LEFT 0     /* set to 1 to add oscilating boundary condition on the left */
#define OSCILLATE_TOPBOT 0   /* set to 1 to enforce a planar wave on top and bottom boundary */
#define OSCILLATION_SCHEDULE 61  /* oscillation schedule, see list in global_pdes.c */
#define OSCIL_YMAX 0.35      /* defines oscillation range */
#define INITIAL_SHIFT 20.0          /* time shift of initial wave packet (in oscillation periods) */
#define WAVE_PACKET_SHIFT 200.0     /* time shift between wave packets (in oscillation periods) */

#define OMEGA 0.025         /* frequency of periodic excitation */
#define AMPLITUDE 0.5      /* amplitude of periodic excitation */ 
#define ACHIRP 0.25        /* acceleration coefficient in chirp */
#define DAMPING 0.0        /* damping of periodic excitation */
#define COURANT 0.12       /* Courant number */
// #define COURANT 0.07       /* Courant number */
#define COURANTB 0.08       /* Courant number in medium B */
// #define GAMMA 1.0e-7          /* damping factor in wave equation */
#define GAMMA 0.0          /* damping factor in wave equation */
#define GAMMAB 0.0         /* damping factor in wave equation */
#define GAMMA_SIDES 1.0e-4      /* damping factor on boundary */
#define GAMMA_TOPBOT 1.0e-7     /* damping factor on boundary */
#define KAPPA 0.0           /* "elasticity" term enforcing oscillations */
#define KAPPA_SIDES 5.0e-4  /* "elasticity" term on absorbing boundary */
#define KAPPA_TOPBOT 0.0    /* "elasticity" term on absorbing boundary */
#define OSCIL_LEFT_YSHIFT 40.0   /* y-dependence of left oscillation (for non-horizontal waves) */
/* The Courant number is given by c*DT/DX, where DT is the time step and DX the lattice spacing */
/* The physical damping coefficient is given by GAMMA/(DT)^2 */
/* Increasing COURANT speeds up the simulation, but decreases accuracy */
/* For similar wave forms, COURANT^2*GAMMA should be kept constant */

#define ADD_OSCILLATING_SOURCE 1        /* set to 1 to add an oscillating wave source */
#define OSCILLATING_SOURCE_PERIOD 16    /* period of oscillating source */
#define ALTERNATE_OSCILLATING_SOURCE 1  /* set to 1 to alternate sign of oscillating source */
#define N_SOURCES 3                     /* number of sources, for option draw_sources */
#define ALTERNATE_SOURCE_PHASES 0       /* set to 1 to alternate initial phases of sources */

#define ADD_WAVE_PACKET_SOURCES 0       /* set to 1 to add several sources emitting wave packets */
#define WAVE_PACKET_SOURCE_TYPE 3       /* type of wave packet sources */
#define N_WAVE_PACKETS 5                /* number of wave packets */
#define WAVE_PACKET_RADIUS 50           /* radius of wave packets */

#define USE_INPUT_TIMESERIES 0          /* set to 1 to use a time series (Morse code) as input * /

/* Boundary conditions, see list in global_pdes.c  */

#define B_COND 2

/* Parameters for length and speed of simulation */

#define NSTEPS 2300         /* number of frames of movie */
#define NVID 8              /* number of iterations between images displayed on screen */
#define NSEG 1000           /* number of segments of boundary */
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

#define INITIAL_AMP 1.5            /* amplitude of initial condition */
#define INITIAL_VARIANCE 0.00001    /* variance of initial condition */
#define INITIAL_WAVELENGTH  0.025   /* wavelength of initial condition */

/* Plot type, see list in global_pdes.c  */

#define PLOT 0

#define PLOT_B 8        /* plot type for second movie */

/* Color schemes */

#define COLOR_PALETTE 13     /* Color palette, see list in global_pdes.c  */
#define COLOR_PALETTE_B 17   /* Color palette, see list in global_pdes.c  */

#define BLACK 1          /* background */

#define COLOR_SCHEME 3   /* choice of color scheme, see list in global_pdes.c  */

#define SCALE 0          /* set to 1 to adjust color scheme to variance of field */
#define SLOPE 1.0        /* sensitivity of color on wave amplitude */
#define PHASE_FACTOR 1.0       /* factor in computation of phase in color scheme P_3D_PHASE */
#define PHASE_SHIFT 0.0      /* shift of phase in color scheme P_3D_PHASE */
#define ATTENUATION 0.0   /* exponential attenuation coefficient of contrast with time */
#define VSHIFT_AMPLITUDE 0.0   /* additional shift for wave amplitude */
#define VSCALE_AMPLITUDE 0.5   /* additional scaling factor for wave amplitude */
#define E_SCALE 10.0       /* scaling factor for energy representation */
#define LOG_SCALE 0.75     /* scaling factor for energy log representation */
#define LOG_SHIFT 0.75     /* shift of colors on log scale */
#define FLUX_SCALE 250.0    /* scaling factor for energy flux represtnation */
#define AVRG_E_FACTOR 0.8   /* controls time window size in P_AVERAGE_ENERGY scheme */
#define RESCALE_COLOR_IN_CENTER 0   /* set to 1 to decrease color intentiy in the center (for wave escaping ring) */
#define FADE_IN_OBSTACLE 1      /* set to 1 to fade color inside obstacles */
#define SHADE_2D 0       /* set to 1 to add pseudo-3d shading effect */ 
#define SHADE_SCALE_2D 0.25  /* lower value increases sensitivity of shading */

#define COLORHUE 260     /* initial hue of water color for scheme C_LUM */
#define COLORDRIFT 0.0   /* how much the color hue drifts during the whole simulation */
#define LUMMEAN 0.5      /* amplitude of luminosity variation for scheme C_LUM */
#define LUMAMP 0.3       /* amplitude of luminosity variation for scheme C_LUM */
#define HUEMEAN 180.0    /* mean value of hue for color scheme C_HUE */
#define HUEAMP -180.0      /* amplitude of variation of hue for color scheme C_HUE */

#define DRAW_COLOR_SCHEME 1     /* set to 1 to plot the color scheme */
#define COLORBAR_RANGE 2.0      /* scale of color scheme bar */
#define COLORBAR_RANGE_B 0.2    /* scale of color scheme bar for 2nd part */
#define ROTATE_COLOR_SCHEME 0   /* set to 1 to draw color scheme horizontally */
#define CIRC_COLORBAR 0         /* set to 1 to draw circular color scheme */
#define CIRC_COLORBAR_B 0       /* set to 1 to draw circular color scheme */

#define DRAW_WAVE_PROFILE 0     /* set to 1 to draw a profile of the wave */
#define HORIZONTAL_WAVE_PROFILE 0 /* set to 1 to draw wave profile vertically */
#define VERTICAL_WAVE_PROFILE 1 /* set to 1 to draw wave profile vertically */
#define WAVE_PROFILE_X 1.9      /* value of x to sample wave profile */
#define WAVE_PROFILE_Y 0.12     /* value of y to sample wave profile */
#define PROFILE_AT_BOTTOM 1     /* draw wave profile at bottom instead of top */
#define AVERAGE_WAVE_PROFILE 0  /* set to 1 to draw time-average of wave profile squared*/
#define DRAW_WAVE_TIMESERIES 0  /* set to 1 to draw a time series of the wave, 2 to also draw it at the top */
#define TIMESERIES_NVALUES 400  /* number of values plotted in time series */
#define SAVE_TIME_SERIES 0      /* set to 1 to save wave time series at a point */
#define DRAW_WAVE_SOURCE 0      /* set to 1 to draw source of wave at (wave_source_x, wave_source_y), set to 2 to also draw focus */

#define MESSAGE_LDASH 14         /* length of dash for Morse code message */
#define MESSAGE_LDOT 8           /* length of dot for Morse code message */
#define MESSAGE_LINTERVAL 54     /* length of interval between dashes/dots for Morse code message */
#define MESSAGE_LINTERLETTER 60  /* length of interval between letters for Morse code message */
#define MESSAGE_LSPACE 48        /* length of space for Morse code message */
#define MESSAGE_INITIAL_TIME 100 /* initial time before starting message for Morse code message */    

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

#define MEAN_FLUX (PLOT == P_TOTAL_ENERGY_FLUX)||(PLOT_B == P_TOTAL_ENERGY_FLUX)
#define REFRESH_IOR ((IOR == IOR_PERIODIC_WELLS_ROTATING)||(IOR == IOR_PERIODIC_WELLS_ROTATING_LARGE))
#define XYIN_INITIALISED (B_DOMAIN == D_IMAGE)

double light[2] = {0.40824829, 0.816496581};   /* location of light source for SHADE_2D option*/

#include "global_pdes.c"        /* constants and global variables */
#include "sub_maze.c"           /* support for generating mazes */
#include "sub_wave.c"           /* common functions for wave_billiard, heat and schrodinger */
#include "wave_common.c"        /* common functions for wave_billiard, wave_comparison, etc */

FILE *time_series_left, *time_series_right;

/*********************/
/* animation part    */
/*********************/

// void evolve_wave_half(double *phi_in[NX], double *psi_in[NX], double *phi_out[NX], double *psi_out[NX], 
//                       short int *xy_in[NX])
void evolve_wave_half(double *phi_in[NX], double *psi_in[NX], double *phi_out[NX], 
                      short int *xy_in[NX], double *tcc[NX], double *tgamma[NX])
/* time step of field evolution */
/* phi is value of field at time t, psi at time t-1 */
/* this version of the function has been rewritten in order to minimize the number of if-branches */
{
    int i, j, iplus, iminus, jplus, jminus, ij[2];
    double delta, x, y, c, cc, gamma, tb_shift;
    static long time = 0;
    static double tc[NX][NY];
    static short int first = 1;
    
    time++;
    
//     if (OSCILLATE_TOPBOT) tb_shift = (int)((X_SHIFT - XMIN)*(double)NX/(XMAX - XMIN));
    if (OSCILLATE_TOPBOT) tb_shift = (int)((XMAX - XMIN)*(double)NX/(XMAX - XMIN));
    
    /* initialize tables with wave speeds and dissipation */
    if (first)
    {
        for (i=0; i<NX; i++){
            for (j=0; j<NY; j++){
                if (xy_in[i][j] != 0)
                {
                    tc[i][j] = sqrt(tcc[i][j]);
                }
                else if (TWOSPEEDS)
                {
                    tc[i][j] = COURANTB;
                }
            }
        }
        first = 0;
    }
    
    #pragma omp parallel for private(i,j,iplus,iminus,jplus,jminus,delta,x,y)
    /* evolution in the bulk */
    for (i=1; i<NX-1; i++){
        for (j=1; j<NY-1; j++){
            if ((TWOSPEEDS)||(xy_in[i][j] != 0)){
                x = phi_in[i][j];
		y = psi_in[i][j];
                
                /* discretized Laplacian */
                delta = phi_in[i+1][j] + phi_in[i-1][j] + phi_in[i][j+1] + phi_in[i][j-1] - 4.0*x;

                /* evolve phi */
                phi_out[i][j] = -y + 2*x + tcc[i][j]*delta - KAPPA*x - tgamma[i][j]*(x-y);
//                 psi_out[i][j] = x;
            }
        }
    }
    
    /* left boundary */
    if (OSCILLATE_LEFT) for (j=1; j<NY-1; j++) phi_out[0][j] = oscillating_bc(time, j);
    else for (j=1; j<NY-1; j++){
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
//                     delta = phi_in[1][j] + phi_in[0][j+1] + phi_in[0][j-1] - 3.0*x;
                    phi_out[0][j] = x - tc[0][j]*(x - phi_in[1][j]) - KAPPA_SIDES*x - GAMMA_SIDES*(x-y);
                    break;
                }
                case (BC_VPER_HABS):
                {
//                     delta = phi_in[1][j] + phi_in[0][j+1] + phi_in[0][j-1] - 3.0*x;
                    phi_out[0][j] = x - tc[0][j]*(x - phi_in[1][j]) - KAPPA_SIDES*x - GAMMA_SIDES*(x-y);
                    break;
                }
            }
//             psi_out[0][j] = x;
        }
    }
    
    /* right boundary */
    for (j=1; j<NY-1; j++){
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
//                     delta = phi_in[NX-2][j] + phi_in[NX-1][j+1] + phi_in[NX-1][j-1] - 3.0*x;
                    phi_out[NX-1][j] = x - tc[NX-1][j]*(x - phi_in[NX-2][j]) - KAPPA_SIDES*x - GAMMA_SIDES*(x-y);
                    break;
                }
                case (BC_VPER_HABS):
                {
//                     delta = phi_in[NX-2][j] + phi_in[NX-1][j+1] + phi_in[NX-1][j-1] - 3.0*x;
                    phi_out[NX-1][j] = x - tc[NX-1][j]*(x - phi_in[NX-2][j]) - KAPPA_SIDES*x - GAMMA_SIDES*(x-y);
                    break;
                }
            }
//             psi_out[NX-1][j] = x;
        }
    }
    
    /* top boundary */
    for (i=0; i<NX; i++){
        if ((TWOSPEEDS)||(xy_in[i][NY-1] != 0)){
            x = phi_in[i][NY-1];
            y = psi_in[i][NY-1];
                    
            if ((OSCILLATE_TOPBOT)&&(i < tb_shift)&&(i<NX-1)&&(i>0))
            {
                iplus = i+1;
                iminus = i-1;   if (iminus < 0) iminus = 0;
                delta = phi_in[iplus][NY-1] + phi_in[iminus][NY-1] + - 2.0*x;
                phi_out[i][NY-1] = -y + 2*x + tcc[i][NY-1]*delta - KAPPA*x - tgamma[i][NY-1]*(x-y);
            }
            
            else switch (B_COND) {
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
                    iminus = (i-1) % NX;
                    if (iminus < 0) iminus += NX;
                    
                    delta = phi_in[iplus][NY-1] + phi_in[iminus][NY-1] + phi_in[i][NY-2] + phi_in[i][0] - 4.0*x;
                    phi_out[i][NY-1] = -y + 2*x + tcc[i][NY-1]*delta - KAPPA*x - tgamma[i][NY-1]*(x-y);
                    break;
                }
                case (BC_ABSORBING):
                {
//                     iplus = (i+1);   if (iplus == NX) iplus = NX-1;
//                     iminus = (i-1);  if (iminus == -1) iminus = 0;
                    
//                     delta = phi_in[iplus][NY-1] + phi_in[iminus][NY-1] + phi_in[i][NY-2] - 3.0*x;
                    phi_out[i][NY-1] = x - tc[i][NY-1]*(x - phi_in[i][NY-2]) - KAPPA_TOPBOT*x - GAMMA_TOPBOT*(x-y);
                    break;
                }
                case (BC_VPER_HABS):
                {
                    iplus = (i+1);   if (iplus == NX) iplus = NX-1;
                    iminus = (i-1);  if (iminus == -1) iminus = 0;

                    delta = phi_in[iplus][NY-1] + phi_in[iminus][NY-1] + phi_in[i][NY-2] + phi_in[i][0] - 4.0*x;
                    if (i==0) phi_out[0][NY-1] = x - tc[0][NY-1]*(x - phi_in[1][NY-1]) - KAPPA_SIDES*x - GAMMA_SIDES*(x-y);
                    else phi_out[i][NY-1] = -y + 2*x + tcc[i][NY-1]*delta - KAPPA*x - tgamma[i][NY-1]*(x-y);
                    break;
                }
            }
//             psi_out[i][NY-1] = x;
        }
    }
    
    /* bottom boundary */
    for (i=0; i<NX; i++){
        if ((TWOSPEEDS)||(xy_in[i][0] != 0)){
            x = phi_in[i][0];
            y = psi_in[i][0];
                    
            if ((OSCILLATE_TOPBOT)&&(i < tb_shift)&&(i<NX-1)&&(i>0))
            {
                iplus = i+1;
                iminus = i-1;   if (iminus < 0) iminus = 0;
                delta = phi_in[iplus][0] + phi_in[iminus][0] + - 2.0*x;
                phi_out[i][0] = -y + 2*x + tcc[i][0]*delta - KAPPA*x - tgamma[i][0]*(x-y);
            }
            
            else switch (B_COND) {
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
                    iminus = (i-1) % NX;
                    if (iminus < 0) iminus += NX;
                    
                    delta = phi_in[iplus][0] + phi_in[iminus][0] + phi_in[i][1] + phi_in[i][NY-1] - 4.0*x;
                    phi_out[i][0] = -y + 2*x + tcc[i][0]*delta - KAPPA*x - tgamma[i][0]*(x-y);
                    break;
                }
                case (BC_ABSORBING):
                {
//                     iplus = (i+1);   if (iplus == NX) iplus = NX-1;
//                     iminus = (i-1);  if (iminus == -1) iminus = 0;
//                     
//                     delta = phi_in[iplus][0] + phi_in[iminus][0] + phi_in[i][1] - 3.0*x;
                    phi_out[i][0] = x - tc[i][0]*(x - phi_in[i][1]) - KAPPA_TOPBOT*x - GAMMA_TOPBOT*(x-y);
                    break;
                }
                case (BC_VPER_HABS):
                {
                    iplus = (i+1);   if (iplus == NX) iplus = NX-1;
                    iminus = (i-1);  if (iminus == -1) iminus = 0;

                    delta = phi_in[iplus][0] + phi_in[iminus][0] + phi_in[i][1] + phi_in[i][NY-1] - 4.0*x;
                    if (i==0) phi_out[0][0] = x - tc[0][0]*(x - phi_in[1][0]) - KAPPA_SIDES*x - GAMMA_SIDES*(x-y);
                    else phi_out[i][0] = -y + 2*x + tcc[i][0]*delta - KAPPA*x - tgamma[i][0]*(x-y);
                    break;
                }
            }
//             psi_out[i][0] = x;
        }
    }
    
    /* add oscillating boundary condition on the left corners */
    if (OSCILLATE_LEFT)
    {
        phi_out[0][0] = oscillating_bc(time, 0);
        phi_out[0][NY-1] = oscillating_bc(time, NY-1);
    }
    
    /* for debugging purposes/if there is a risk of blow-up */
    if (FLOOR) for (i=0; i<NX; i++){
        for (j=0; j<NY; j++){
            if (xy_in[i][j] != 0) 
            {
                if (phi_out[i][j] > VMAX) phi_out[i][j] = VMAX;
                if (phi_out[i][j] < -VMAX) phi_out[i][j] = -VMAX;
//                 if (psi_out[i][j] > VMAX) psi_out[i][j] = VMAX;
//                 if (psi_out[i][j] < -VMAX) psi_out[i][j] = -VMAX;
            }
        }
    }
}


void evolve_wave(double *phi[NX], double *psi[NX], double *tmp[NX], short int *xy_in[NX], double *tcc_table[NX], double *tgamma_table[NX])
/* time step of field evolution */
/* phi is value of field at time t, psi at time t-1 */
{
   // For the purpose of these comments w[t], w[t-1], w[t+1] are used to refer
    // to phi, psi and the result respectively to avoid confusion with the
    // passed parameter names.
    // At the beginning w[t] is saved in phi, w[t-1] in psi and tmp is space
    // for the next wave state w[t+1]. Take w[t] and w[t-1] to calculate the
    // next wave state. Write this new state in temp
    evolve_wave_half(phi, psi, tmp, xy_in, tcc_table, tgamma_table);
    // now w[t] is saved in tmp, w[t-1] in phi and the result is written to psi
    evolve_wave_half(tmp, phi, psi, xy_in, tcc_table, tgamma_table);
    // now w[t] is saved in psi, w[t-1] in tmp and the result is written to phi
    evolve_wave_half(psi, tmp, phi, xy_in, tcc_table, tgamma_table);
    // now w[t] is saved in phi, w[t-1] in psi and tmp is free again to take
    // the new wave state w[t+1] in the next call to this function, thus
    // matching the given parameter names again
}


void draw_color_bar(int plot, double range)
{
    if (ROTATE_COLOR_SCHEME) draw_color_scheme(-1.0, -0.8, XMAX - 0.1, -1.0, plot, -range, range);
    else draw_color_scheme(XMAX - 0.3, YMIN + 0.1, XMAX - 0.1, YMAX - 0.1, plot, -range, range);
//     else draw_color_scheme(1.7, YMIN + 0.25, 1.9, YMAX - 0.25, plot, -range, range);
}

void draw_color_bar_palette(int plot, double range, int palette, int circular, int fade, double fade_value)
{
    double width = 0.14;
//     double width = 0.2;
    
    if (ROTATE_COLOR_SCHEME) 
        draw_color_scheme_palette_fade(-1.0, -0.8, XMAX - 0.1, -1.0, plot, -range, range, palette, fade, fade_value);
    else if (circular)
        draw_circular_color_scheme_palette_fade(XMAX - 2.0*width, YMIN + 2.0*width, 1.5*width, plot, -range, range, palette, fade, fade_value);
    else 
        draw_color_scheme_palette_fade(XMAX - 1.5*width, YMIN + 0.1, XMAX - 0.5*width, YMAX - 0.1, plot, -range, range, palette, fade, fade_value);
}

void animation()
{
    double time, scale, ratio, startleft[2], startright[2], sign[N_SOURCES], r2, xy[2], fade_value, yshift, speed = 0.0, a, b, c, x, y, angle = 0.0, x1, ior_angle = 0.0, omega, phase_shift, vshift, dsource, finv, source_amp; 
    double *phi[NX], *psi[NX], *tmp[NX], *total_energy[NX], *average_energy[NX], *color_scale[NX], *total_flux, *tcc_table[NX], *tgamma_table[NX], *fade_table;
    short int *xy_in[NX];
    int i, j, k, s, sample_left[2], sample_right[2], period = 0, fade, source_counter = 0, p, q, first_source = 1, imin, imax, ij[2], source, source_period, source_shift[N_SOURCES];
//     static int image_counter = 0;
    int image_counter = 0;
    long int wave_value;
    t_wave_packet *packet;
    t_wave_source wave_source[25];
    
    if (SAVE_TIME_SERIES)
    {
        time_series_left = fopen("wave_left.dat", "w");
        time_series_right = fopen("wave_right.dat", "w");
    }

    /* Since NX and NY are big, it seemed wiser to use some memory allocation here */
    for (i=0; i<NX; i++)
    {
        phi[i] = (double *)malloc(NY*sizeof(double));
        psi[i] = (double *)malloc(NY*sizeof(double));
        tmp[i] = (double *)malloc(NY*sizeof(double));
        total_energy[i] = (double *)malloc(NY*sizeof(double));
        average_energy[i] = (double *)malloc(NY*sizeof(double));
        xy_in[i] = (short int *)malloc(NY*sizeof(short int));
        color_scale[i] = (double *)malloc(NY*sizeof(double));
        tcc_table[i] = (double *)malloc(NX*sizeof(double));
        tgamma_table[i] = (double *)malloc(NX*sizeof(double));
    }
    
    if (FADE_IN_OBSTACLE) fade_table = (double *)malloc(NX*NY*sizeof(double));
    
    if (MEAN_FLUX) total_flux = (double *)malloc(4*NX*NY*sizeof(double));
    
    if (ADD_WAVE_PACKET_SOURCES)  
    {
        packet = (t_wave_packet *)malloc(N_WAVE_PACKETS*sizeof(t_wave_packet));
        init_wave_packets(packet, WAVE_PACKET_RADIUS);
    }
    
    if (ALTERNATE_SOURCE_PHASES) for (i=0; i<N_SOURCES; i++) sign[i] = pow(-1.0,(double)i);
    else for (i=0; i<N_SOURCES; i++) sign[i] = 1.0;
    
    /* initialise positions and radii of circles */
    if ((B_DOMAIN == D_CIRCLES)||(B_DOMAIN == D_CIRCLES_IN_RECT)) ncircles = init_circle_config(circles);
    else if (B_DOMAIN == D_POLYGONS) ncircles = init_polygon_config(polygons);
    printf("Polygons initialized\n");
    
    /* initialise polyline for von Koch and similar domains */
    npolyline = init_polyline(MDEPTH, polyline);
    for (i=0; i<npolyline; i++) printf("vertex %i: (%.3f, %.3f)\n", i, polyline[i].x, polyline[i].y);
    
    npolyrect = init_polyrect(polyrect);
    for (i=0; i<npolyrect; i++) printf("polyrect vertex %i: (%.3f, %.3f) - (%.3f, %.3f)\n", i, polyrect[i].x1, polyrect[i].y1, polyrect[i].x2, polyrect[i].y2);
    printf("Rectangles initialized\n");
    
    init_polyrect_arc(polyrectrot, polyarc, &npolyrect_rot, &npolyarc);
    printf("Rotated rectangles and arcs initialized\n");
    printf("%i rotated rectangles, %i arcs\n", npolyrect_rot, npolyarc);
    
    if ((DRAW_WAVE_TIMESERIES)||(USE_INPUT_TIMESERIES)) init_input_signal();
    
    courant2 = COURANT*COURANT;
    courantb2 = COURANTB*COURANTB;
    c = COURANT*(XMAX - XMIN)/(double)NX;
    
    /* initialize color scale, for option RESCALE_COLOR_IN_CENTER */
    if (RESCALE_COLOR_IN_CENTER)
    {
        for (i=0; i<NX; i++)
            for (j=0; j<NY; j++)
            {
                ij_to_xy(i, j, xy);
                r2 = xy[0]*xy[0] + xy[1]*xy[1];
                color_scale[i][j] = 1.0 - exp(-4.0*r2/LAMBDA*LAMBDA);
            }
    }

    /* initialize wave with a drop at one point, zero elsewhere */
//     init_circular_wave(0.0, -LAMBDA, phi, psi, xy_in);
    
    /* initialize total energy table */
    if ((PLOT == P_MEAN_ENERGY)||(PLOT_B == P_MEAN_ENERGY)||(PLOT == P_LOG_MEAN_ENERGY)||(PLOT_B == P_LOG_MEAN_ENERGY))
        for (i=0; i<NX; i++)
            for (j=0; j<NY; j++) 
                total_energy[i][j] = 0.0;
    
    /* initialize average energy table */
    if ((PLOT == P_AVERAGE_ENERGY)||(PLOT_B == P_AVERAGE_ENERGY)||(PLOT == P_LOG_AVERAGE_ENERGY)||(PLOT_B == P_LOG_AVERAGE_ENERGY))
        for (i=0; i<NX; i++)
            for (j=0; j<NY; j++) 
                average_energy[i][j] = 0.0;
            
    if (MEAN_FLUX) 
        for (i=0; i<4*NX*NY; i++)
            total_flux[i] = 0.0;

    ratio = (XMAX - XMIN)/8.4;  /* for Tokarsky billiard */
    
    if (B_DOMAIN == D_MICHELSON_MOVING)
    {
        xy_to_ij(michelson_schedule(0) - 0.1, YMIN, ij);
        imin = ij[0];
        xy_to_ij(michelson_schedule(NSTEPS) + 0.1, YMIN, ij);
        imax = ij[0];
//         imin = NX/2;
//         imax = NX;
        printf("imin = %i, imax = %i\n", imin, imax);
    }
    
    if (XYIN_INITIALISED) init_xyin_from_image(xy_in);
    
//     isospectral_initial_point(0.2, 0.0, startleft, startright);    /* for isospectral billiards */
//     homophonic_initial_point(0.5, -0.25, 1.5, -0.25, startleft, startright);
//     homophonic_initial_point(0.5, -0.25, 1.5, -0.25, startleft, startright);
//     printf("xleft = (%.3f, %.3f) xright = (%.3f, %.3f)\n", startleft[0], startleft[1], startright[0], startright[1]);    
    
//     xy_to_ij(startleft[0], startleft[1], sample_left);
//     xy_to_ij(startright[0], startright[1], sample_right);
//     printf("xleft = (%.3f, %.3f) xright = (%.3f, %.3f)\n", xin_left, yin_left, xin_right, yin_right);
    
    init_wave_flat(phi, psi, xy_in);
    
    init_ior_2d(xy_in, tcc_table, tgamma_table, ior_angle);
    if (FADE_IN_OBSTACLE) init_fade_table(tcc_table, fade_table); 

//     init_circular_wave(-LAMBDA, 0.0, phi, psi, xy_in);
//     x = XMIN + (XMAX - XMIN)*rand()/RAND_MAX;
//     y = YMIN + (YMAX - YMIN)*rand()/RAND_MAX;
//     init_circular_wave(0.0, -0.8, phi, psi, xy_in);
//     add_circular_wave(-1.0, -1.5, -0.8, phi, psi, xy_in);
//     add_circular_wave(-1.0, 1.5, -0.8, phi, psi, xy_in);
//     sign = -sign;
    
//         init_circular_wave(2.0*LAMBDA*cos(APOLY*PID), 2.0*LAMBDA*sin(APOLY*PID), phi, psi, xy_in);
//         angle = DPI/(double)NPOLY;
//         for (j=1; j<NPOLY; j++)
//             add_circular_wave(1.0, 2.0*LAMBDA*cos((double)j*angle + APOLY*PID), 2.0*LAMBDA*sin((double)j*angle + APOLY*PID), phi, psi, xy_in);

//     init_wave_plus(LAMBDA - 0.3*MU, 0.5*MU, phi, psi, xy_in);
//     init_wave(LAMBDA - 0.3*MU, 0.5*MU, phi, psi, xy_in);
//     init_circular_wave(X_SHOOTER, Y_SHOOTER, phi, psi, xy_in);
//     printf("Initializing wave\n");
//     init_circular_wave(-0.5, 0.0, phi, psi, xy_in);
//     printf("Wave initialized\n");

        
//     init_circular_wave(0.6*cos((double)(period)*DPI/3.0), 0.6*sin((double)(period)*DPI/3.0), phi, psi, xy_in);
//     period++;
//     for (i=0; i<3; i++)
//     {
//         add_circular_wave(-1.0, 0.6*cos(PID + (double)(i)*DPI/3.0), 0.6*sin(PID + (double)(i)*DPI/3.0), phi, psi, xy_in);
//     }
//     add_circular_wave(-1.0, 0.0, LAMBDA, phi, psi, xy_in);
//     add_circular_wave(1.0, -LAMBDA, 0.0, phi, psi, xy_in);
//     add_circular_wave(-1.0, 0.0, -LAMBDA, phi, psi, xy_in);
    
//     init_circular_wave_xplusminus(startleft[0], startleft[1], startright[0], startright[1], phi, psi, xy_in);

//     init_circular_wave_xplusminus(-0.9, 0.0, 0.81, 0.0, phi, psi, xy_in);
//     init_circular_wave(-2.0*ratio, 0.0, phi, psi, xy_in);
//     init_planar_wave(XMIN + 0.015, 0.0, phi, psi, xy_in);
//     init_planar_wave(XMIN + 0.02, 0.0, phi, psi, xy_in);
//     init_planar_travelling_wave(XMIN + 0.1, 0.3, 1.0, phi, psi, xy_in);
//     add_planar_travelling_wave(XMIN + 0.05, 0.0, 1.0, 0.3, phi, psi, xy_in);
    
//     init_wave(-1.5, 0.0, phi, psi, xy_in);
//     init_wave(0.0, 0.0, phi, psi, xy_in);

    /* add a drop at another point */
//     add_drop_to_wave(1.0, 0.7, 0.0, phi, psi);
//     add_drop_to_wave(1.0, -0.7, 0.0, phi, psi);
//     add_drop_to_wave(1.0, 0.0, -0.7, phi, psi);

    blank();
    glColor3f(0.0, 0.0, 0.0);
//     draw_wave(phi, psi, xy_in, 1.0, 0, PLOT);
    if (HIGHRES) draw_wave_highres_palette(2, phi, psi, total_energy, average_energy, fade_table, total_flux, color_scale, xy_in, 1.0, 0, PLOT, COLOR_PALETTE, 1, 0, 1.0);
    else draw_wave_epalette(phi, psi, total_energy, average_energy, total_flux, color_scale, xy_in, 1.0, 0, PLOT, COLOR_PALETTE, 0, 1.0);

//     draw_billiard(0, 1.0);
    
    
    if (DRAW_COLOR_SCHEME) draw_color_bar_palette(PLOT, COLORBAR_RANGE, COLOR_PALETTE, CIRC_COLORBAR, fade, fade_value);
    if (PRINT_SPEED) 
    {   
        a = 0.0075;
        b = 0.00015;
//         speed = a/((double)(NVID)*c);
//         speed = 0.55*a/((double)(NVID*OSCILLATING_SOURCE_PERIOD)*c);
        speed = a/((double)(3*NVID*OSCILLATING_SOURCE_PERIOD)*c);
        /* the factor 3 is due to evolve_wave calling evolve_wave_half 3 times */
        print_speed(speed, 0, 1.0);
    }

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

//         draw_wave(phi, psi, xy_in, scale, i, PLOT);
        if (HIGHRES) draw_wave_highres_palette(2, phi, psi, total_energy, average_energy, fade_table, total_flux, color_scale, xy_in, scale, i, PLOT, COLOR_PALETTE, 0, 0, 1.0);
        else draw_wave_epalette(phi, psi, total_energy, average_energy, total_flux, color_scale, xy_in, scale, i, PLOT, COLOR_PALETTE, 0, 1.0);
        for (j=0; j<NVID; j++) 
        {
            evolve_wave(phi, psi, tmp, xy_in, tcc_table, tgamma_table);
            if (SAVE_TIME_SERIES)
            {
                wave_value = (long int)(phi[sample_left[0]][sample_left[1]]*1.0e16);
                fprintf(time_series_left, "%019ld\n", wave_value);
                wave_value = (long int)(phi[sample_right[0]][sample_right[1]]*1.0e16);
                fprintf(time_series_right, "%019ld\n", wave_value);
                if ((j == 0)&&(i%10 == 0)) printf("Frame %i of %i\n", i, NSTEPS);
//                 fprintf(time_series_right, "%.15f\n", phi[sample_right[0]][sample_right[1]]);
            }
        }
        
//         draw_billiard(0, 1.0);
        
        if (DRAW_COLOR_SCHEME) draw_color_bar_palette(PLOT, COLORBAR_RANGE, COLOR_PALETTE, CIRC_COLORBAR, fade, fade_value);
        
        /* add oscillating waves */
        for (source = 0; source < N_SOURCES; source++)
        {
            angle = APOLY*PID + (double)source*DPI/(double)NPOLY;
            wave_source_x[source] = LAMBDA*cos(angle) + 2.0*MU*cos(angle - 0.5*PID);
            wave_source_y[source] = LAMBDA*sin(angle) + 2.0*MU*sin(angle - 0.5*PID);
            source_shift[source] = source*OSCILLATING_SOURCE_PERIOD/N_SOURCES;
            if ((ADD_OSCILLATING_SOURCE)&&(i%OSCILLATING_SOURCE_PERIOD == source_shift[source]))
            {
                if (ALTERNATE_OSCILLATING_SOURCE) sign[source] = -sign[source];
                add_circular_wave(-sign[source]*INITIAL_AMP, wave_source_x[source], wave_source_y[source], phi, psi, xy_in);
            }
        }

        
        if (ADD_WAVE_PACKET_SOURCES) add_wave_packets(phi, psi, xy_in, packet, i, WAVE_PACKET_RADIUS, 1, 4, 1);
        if (PRINT_SPEED) print_speed(speed, 0, 1.0);
        if (PRINT_FREQUENCY) print_frequency(phase_shift, 0, 1.0);
        if ((VARIABLE_IOR)&&(REFRESH_IOR)&&(i%3 == 0))
        {
            ior_angle = ior_angle_schedule(i);
            printf("IOR angle = %.5lg\n", ior_angle); 
            init_ior_2d(xy_in, tcc_table, tgamma_table, ior_angle);
            printf("speed = %.5lg\n", tcc_table[3*NX/4][NY/2]);
        }
        if (B_DOMAIN == D_MICHELSON_MOVING) 
        {
            michelson_position = michelson_schedule(i);
            printf("Michelson mirror position %.3lg\n", michelson_position);
            init_xy_tcc_in_xrange(xy_in, tcc_table, tgamma_table, imin, imax);
        }

	if (!((NO_EXTRA_BUFFER_SWAP)&&(MOVIE))) glutSwapBuffers();

	if (MOVIE)
        {
            if (i >= INITIAL_TIME) save_frame();
            else printf("Initial phase time %i of %i\n", i, INITIAL_TIME);
            
            if ((i >= INITIAL_TIME)&&(DOUBLE_MOVIE))
            {
//                 draw_wave(phi, psi, xy_in, scale, i, PLOT_B);
                if (HIGHRES) 
                    draw_wave_highres_palette(2, phi, psi, total_energy, average_energy, fade_table, total_flux, color_scale, xy_in, scale, i, PLOT_B, COLOR_PALETTE_B, 1, 0, 1.0);
                else draw_wave_epalette(phi, psi, total_energy, average_energy, total_flux, color_scale, xy_in, scale, i, PLOT_B, COLOR_PALETTE_B, 0, 1.0);
//                 draw_billiard(0, 1.0);
                if (DRAW_COLOR_SCHEME) draw_color_bar_palette(PLOT_B, COLORBAR_RANGE_B, COLOR_PALETTE_B, CIRC_COLORBAR_B, 0, 1.0);  
                if (PRINT_SPEED) print_speed(speed, 0, 1.0);
                if (PRINT_FREQUENCY) print_frequency(phase_shift, 0, 1.0);
                glutSwapBuffers();
//                 printf("image_counter = %i\n", image_counter);
//                 printf("image number = %i\n", NSTEPS + MID_FRAMES + 1 + image_counter);
                save_frame_counter(NSTEPS + MID_FRAMES + 1 + image_counter);
                image_counter++;
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

    }

    if (MOVIE) 
    {
        if (DOUBLE_MOVIE) 
        {
//             draw_wave(phi, psi, xy_in, scale, i, PLOT);
            if (HIGHRES) draw_wave_highres_palette(2, phi, psi, total_energy, average_energy, fade_table, total_flux, color_scale, xy_in, scale, NSTEPS, PLOT, COLOR_PALETTE, 0, 0, 1.0);
            else draw_wave_epalette(phi, psi, total_energy, average_energy, total_flux, color_scale, xy_in, scale, NSTEPS, PLOT, COLOR_PALETTE, 0, 1.0);
//             draw_billiard(0, 1.0);
            if (DRAW_COLOR_SCHEME) draw_color_bar_palette(PLOT, COLORBAR_RANGE, COLOR_PALETTE, CIRC_COLORBAR, 0, 1.0);   
            if (PRINT_SPEED) print_speed(speed, 0, 1.0);
            if (PRINT_FREQUENCY) print_frequency(phase_shift, 0, 1.0);
            glutSwapBuffers();
        }
        if (!FADE) for (i=0; i<MID_FRAMES; i++) save_frame();
        else for (i=0; i<MID_FRAMES; i++) 
        {
            fade_value = 1.0 - (double)i/(double)MID_FRAMES;
            if (HIGHRES) 
                draw_wave_highres_palette(2, phi, psi, total_energy, average_energy, fade_table, total_flux, color_scale, xy_in, scale, NSTEPS, PLOT, COLOR_PALETTE, 0, 1, fade_value);
            else draw_wave_epalette(phi, psi, total_energy, average_energy, total_flux, color_scale, xy_in, scale, NSTEPS, PLOT, COLOR_PALETTE, 1, fade_value);
//             draw_billiard(1, fade_value);
            if (DRAW_COLOR_SCHEME) draw_color_bar_palette(PLOT, COLORBAR_RANGE, COLOR_PALETTE, CIRC_COLORBAR, 1, fade_value); 
            if (PRINT_SPEED) print_speed(speed, 1, fade_value);
            if (PRINT_FREQUENCY) print_frequency(phase_shift, 1, fade_value);
            if (!NO_EXTRA_BUFFER_SWAP) glutSwapBuffers();
            save_frame_counter(NSTEPS + i + 1);
        }
        if (DOUBLE_MOVIE) 
        {
//             draw_wave(phi, psi, xy_in, scale, i, PLOT_B);
            if (HIGHRES) 
                draw_wave_highres_palette(2, phi, psi, total_energy, average_energy, fade_table, total_flux, color_scale, xy_in, scale, NSTEPS, PLOT_B, COLOR_PALETTE_B, 1, 0, 1.0);
            else draw_wave_epalette(phi, psi, total_energy, average_energy, total_flux, color_scale, xy_in, scale, NSTEPS, PLOT_B, COLOR_PALETTE_B, 0, 1.0);
//             draw_billiard(0, 1.0);
            if (DRAW_COLOR_SCHEME) draw_color_bar_palette(PLOT_B, COLORBAR_RANGE_B, COLOR_PALETTE_B, CIRC_COLORBAR_B, 0, 1.0); 
            if (PRINT_SPEED) print_speed(speed, 0, 1.0);
            if (PRINT_FREQUENCY) print_frequency(phase_shift, 0, 1.0);
            glutSwapBuffers();

            if (!FADE) for (i=0; i<END_FRAMES; i++) save_frame_counter(NSTEPS + MID_FRAMES + 1 + image_counter + i);
            else for (i=0; i<END_FRAMES; i++) 
            {
                fade_value = 1.0 - (double)i/(double)END_FRAMES;
                if (HIGHRES) 
                    draw_wave_highres_palette(2, phi, psi, total_energy, average_energy, fade_table, total_flux, color_scale, xy_in, scale, NSTEPS, PLOT_B, COLOR_PALETTE_B, 1, 1, fade_value);
                else draw_wave_epalette(phi, psi, total_energy, average_energy, total_flux, color_scale, xy_in, scale, NSTEPS, PLOT_B, COLOR_PALETTE_B, 1, fade_value);
//                 draw_billiard(1, fade_value);
                if (DRAW_COLOR_SCHEME) draw_color_bar_palette(PLOT_B, COLORBAR_RANGE_B, COLOR_PALETTE_B, CIRC_COLORBAR_B, 1, fade_value); 
                if (PRINT_SPEED) print_speed(speed, 1, fade_value);
                if (PRINT_FREQUENCY) print_frequency(phase_shift, 0, 1.0);
                glutSwapBuffers();
                save_frame_counter(NSTEPS + MID_FRAMES + 1 + image_counter + i);
            }
        }
        
        s = system("mv wave*.tif tif_wave/");
    }
    for (i=0; i<NX; i++)
    {
        free(phi[i]);
        free(psi[i]);
        free(tmp[i]);
        free(total_energy[i]);
        free(average_energy[i]);
        free(xy_in[i]);
        free(color_scale[i]);
        free(tcc_table[i]);
        free(tgamma_table[i]);
    }
    
    if (FADE_IN_OBSTACLE) free(fade_table); 
    
    if (MEAN_FLUX) free(total_flux);
    
    if (ADD_WAVE_PACKET_SOURCES) free(packet);
    
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

    init();

    glutDisplayFunc(display);

    glutMainLoop();

    return 0;
}

