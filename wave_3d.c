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
#define DOUBLE_MOVIE 0  /* set to 1 to produce movies for wave height and energy simultaneously */
#define SAVE_MEMORY 1   /* set to 1 to save memory when writing tiff images */
#define NO_EXTRA_BUFFER_SWAP 1    /* some OS require one less buffer swap when recording images */

/* General geometrical parameters */

#define WINWIDTH 	1920  /* window width */
#define WINHEIGHT 	1150  /* window height */
// // #define NX 1920          /* number of grid points on x axis */
// #define NY 1150          /* number of grid points on y axis */
#define NX 3000          /* number of grid points on x axis */
#define NY 1600          /* number of grid points on y axis */
// #define NX 3840          /* number of grid points on x axis */
// #define NY 2300          /* number of grid points on y axis */

// #define XMIN -2.0
// #define XMAX 2.0	/* x interval */
// #define YMIN -1.197916667
// #define YMAX 1.197916667	/* y interval for 9/16 aspect ratio */
#define XMIN -1.669565217
#define XMAX 1.669565217	/* x interval */
#define YMIN -1.0
#define YMAX 1.0	/* y interval for 9/16 aspect ratio */

#define HIGHRES 0        /* set to 1 if resolution of grid is double that of displayed image */

// #define WINWIDTH 	1280  /* window width */
// #define WINHEIGHT 	720   /* window height */
// 
// // #define NX 1280          /* number of grid points on x axis */
// // #define NX 720         /* number of grid points on x axis */
// // #define NY 720         /* number of grid points on y axis */
// #define NX 2560          /* number of grid points on x axis */
// // #define NX 1440          /* number of grid points on x axis */
// #define NY 1440         /* number of grid points on y axis */
// 
// // #define NX 360         /* number of grid points on x axis */
// // #define NY 360         /* number of grid points on y axis */
//  
// #define XMIN -2.0
// #define XMAX 2.0	/* x interval  */
// #define YMIN -1.125
// #define YMAX 1.125	/* y interval for 9/16 aspect ratio */

#define JULIA_SCALE 0.8 /* scaling for Julia sets */

/* Choice of the billiard table */

#define B_DOMAIN 17        /* choice of domain shape, see list in global_pdes.c */

#define CIRCLE_PATTERN 2   /* pattern of circles or polygons, see list in global_pdes.c */

#define COMPARISON 0        /* set to 1 to compare two different patterns */
#define B_DOMAIN_B 20       /* second domain shape, for comparisons */
#define CIRCLE_PATTERN_B 0  /* second pattern of circles or polygons */

#define VARIABLE_IOR 0      /* set to 1 for a variable index of refraction */
#define IOR 9               /* choice of index of refraction, see list in global_pdes.c */
#define IOR_TOTAL_TURNS 1.0 /* total angle of rotation for IOR_PERIODIC_WELLS_ROTATING */
#define MANDEL_IOR_SCALE -0.05   /* parameter controlling dependence of IoR on Mandelbrot escape speed */

#define P_PERCOL 0.25       /* probability of having a circle in C_RAND_PERCOL arrangement */
#define NPOISSON 1000        /* number of points for Poisson C_RAND_POISSON arrangement */
#define RANDOM_POLY_ANGLE 1 /* set to 1 to randomize angle of polygons */

#define LAMBDA 3.0	    /* parameter controlling the dimensions of domain */
#define MU 0.14             /* parameter controlling the dimensions of domain */
#define NPOLY 6             /* number of sides of polygon */
#define APOLY 0.0           /* angle by which to turn polygon, in units of Pi/2 */ 
#define MDEPTH 2            /* depth of computation of Menger gasket */
#define MRATIO 3            /* ratio defining Menger gasket */
#define MANDELLEVEL 2000    /* iteration level for Mandelbrot set */
#define MANDELLIMIT 20.0    /* limit value for approximation of Mandelbrot set */
#define FOCI 1              /* set to 1 to draw focal points of ellipse */
#define NGRIDX 30           /* number of grid point for grid of disks */
#define NGRIDY 18           /* number of grid point for grid of disks */
#define WALL_WIDTH 0.1      /* width of wall separating lenses */

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
#define OSCILLATION_SCHEDULE 3  /* oscillation schedule, see list in global_pdes.c */
#define OSCIL_YMAX 0.35      /* defines oscillation range */

#define OMEGA 0.015        /* frequency of periodic excitation */
#define AMPLITUDE 1.0     /* amplitude of periodic excitation */ 
#define ACHIRP 0.2        /* acceleration coefficient in chirp */
#define DAMPING 0.0       /* damping of periodic excitation */
#define COURANT 0.1       /* Courant number */
#define COURANTB 0.01      /* Courant number in medium B */
#define GAMMA 0.0          /* damping factor in wave equation */
#define GAMMAB 1.0e-6         /* damping factor in wave equation */
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

#define ADD_OSCILLATING_SOURCE 0        /* set to 1 to add an oscillating wave source */
#define OSCILLATING_SOURCE_PERIOD 30    /* period of oscillating source */
#define ALTERNATE_OSCILLATING_SOURCE 1  /* set to 1 to alternate sign of oscillating source */

#define ADD_WAVE_PACKET_SOURCES 0       /* set to 1 to add several sources emitting wave packets */
#define WAVE_PACKET_SOURCE_TYPE 1       /* type of wave packet sources */
#define N_WAVE_PACKETS 15               /* number of wave packets */
#define WAVE_PACKET_RADIUS 20            /* radius of wave packets */

/* Boundary conditions, see list in global_pdes.c  */

// #define B_COND 1
#define B_COND 3

#define PRECOMPUTE_BC 0     /* set to 1 to compute neighbours for Laplacian in advance */

/* Parameters for length and speed of simulation */

#define NSTEPS 1000       /* number of frames of movie */
// #define NSTEPS 300         /* number of frames of movie */
#define NVID 20            /* number of iterations between images displayed on screen */
// #define NVID 10            /* number of iterations between images displayed on screen */
#define NSEG 1000          /* number of segments of boundary */
#define INITIAL_TIME 0      /* time after which to start saving frames */
#define BOUNDARY_WIDTH 2    /* width of billiard boundary */
#define PRINT_SPEED 0       /* set to 1 to print speed of moving source */

#define PAUSE 100       /* number of frames after which to pause */
#define PSLEEP 3         /* sleep time during pause */
#define SLEEP1  1        /* initial sleeping time */
#define SLEEP2  1        /* final sleeping time */
#define MID_FRAMES 100    /* number of still frames between parts of two-part movie */
#define END_FRAMES 100   /* number of still frames at end of movie */
#define FADE 1           /* set to 1 to fade at end of movie */
#define ROTATE_VIEW_WHILE_FADE 1    /* set to 1 to keep rotating viewpoint during fade */

/* Parameters of initial condition */

#define INITIAL_AMP 0.75            /* amplitude of initial condition */
// #define INITIAL_VARIANCE 0.000025    /* variance of initial condition */
#define INITIAL_VARIANCE 0.00005    /* variance of initial condition */
#define INITIAL_WAVELENGTH  0.05  /* wavelength of initial condition */


/* Plot type, see list in global_pdes.c  */

#define ZPLOT 108     /* wave height */
#define CPLOT 108     /* color scheme */

#define ZPLOT_B 103 
#define CPLOT_B 103        /* plot type for second movie */

#define CHANGE_LUMINOSITY 1     /* set to 1 to let luminosity depend on energy flux intensity */
#define FLUX_WINDOW 30           /* size of averaging window of flux intensity */
#define AMPLITUDE_HIGH_RES 1    /* set to 1 to increase resolution of plot */
#define SHADE_3D 1              /* set to 1 to change luminosity according to normal vector */
#define NON_DIRICHLET_BC 0      /* set to 1 to draw only facets in domain, if field is not zero on boundary */
#define FLOOR_ZCOORD 1          /* set to 1 to draw only facets with z not too negative */
#define DRAW_BILLIARD 0         /* set to 1 to draw boundary */
#define DRAW_BILLIARD_FRONT 0   /* set to 1 to draw front of boundary after drawing wave */
#define DRAW_CONSTRUCTION_LINES 0   /* set to 1 to draw construction lines of certain domains */
#define FADE_IN_OBSTACLE 1      /* set to 1 to fade color inside obstacles */
#define DRAW_OUTSIDE_GRAY 0     /* experimental, draw outside of billiard in gray */

#define PLOT_SCALE_ENERGY 0.4          /* vertical scaling in energy plot */
#define PLOT_SCALE_LOG_ENERGY 0.5       /* vertical scaling in log energy plot */

/* 3D representation */

#define REPRESENTATION_3D 1     /* choice of 3D representation */ 

#define REP_AXO_3D 0        /* linear projection (axonometry) */
#define REP_PROJ_3D 1       /* projection on plane orthogonal to observer line of sight */

#define ROTATE_VIEW 0       /* set to 1 to rotate position of observer */
#define ROTATE_ANGLE 360.0  /* total angle of rotation during simulation */
// #define ROTATE_ANGLE 45.0  /* total angle of rotation during simulation */

/* Color schemes */

#define COLOR_PALETTE 13      /* Color palette, see list in global_pdes.c  */
#define COLOR_PALETTE_B 17     /* Color palette, see list in global_pdes.c  */

#define BLACK 1          /* background */

#define COLOR_SCHEME 3   /* choice of color scheme, see list in global_pdes.c  */

#define SCALE 0          /* set to 1 to adjust color scheme to variance of field */
#define SLOPE 1.0        /* sensitivity of color on wave amplitude */
#define VSCALE_AMPLITUDE 2.0   /* additional scaling factor for color scheme P_3D_AMPLITUDE */
#define VSCALE_ENERGY 1.5      /* additional scaling factor for color scheme P_3D_ENERGY */
#define PHASE_FACTOR 20.0      /* factor in computation of phase in color scheme P_3D_PHASE */
#define PHASE_SHIFT 0.0      /* shift of phase in color scheme P_3D_PHASE */
#define ATTENUATION 0.0    /* exponential attenuation coefficient of contrast with time */
#define E_SCALE 50.0      /* scaling factor for energy representation */
#define LOG_SCALE 0.75     /* scaling factor for energy log representation */
#define LOG_SHIFT 0.5      /* shift of colors on log scale */
#define LOG_ENERGY_FLOOR -10.0    /* floor value for log of (total) energy */
#define LOG_MEAN_ENERGY_SHIFT 1.0   /* additional shift for log of mean energy */
#define FLUX_SCALE 4000.0    /* scaling factor for energy flux representation */
#define FLUX_CSCALE 2.0      /* scaling factor for color in energy flux representation */
#define AVRG_E_FACTOR 0.95   /* controls time window size in P_AVERAGE_ENERGY scheme */
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
#define RAND_SHIFT 5        /* seed of random number generator */
#define MAZE_XSHIFT 0.0     /* horizontal shift of maze */
#define MAZE_WIDTH 0.02     /* half width of maze walls */

#define DRAW_COLOR_SCHEME 1       /* set to 1 to plot the color scheme */
#define COLORBAR_RANGE 1.0      /* scale of color scheme bar */
#define COLORBAR_RANGE_B 1.5    /* scale of color scheme bar for 2nd part */
#define ROTATE_COLOR_SCHEME 0     /* set to 1 to draw color scheme horizontally */

#define SAVE_TIME_SERIES 0      /* set to 1 to save wave time series at a point */

#define ADD_POTENTIAL 1         /* set to 1 to add potential to z coordinate */
// #define POT_MAZE 7
#define POTENTIAL 10
#define POT_FACT 20.0
#define DRAW_WAVE_PROFILE 0     /* set to 1 to draw a profile of the wave */
#define MU_B 1.0           /* parameter controlling the dimensions of domain */
#define HORIZONTAL_WAVE_PROFILE 0 /* set to 1 to draw wave profile vertically */
#define VERTICAL_WAVE_PROFILE 0 /* set to 1 to draw wave profile vertically */
#define DRAW_WAVE_TIMESERIES 0  /* set to 1 to draw a time series of the wave */
#define WAVE_PROFILE_X 2.1      /* value of x to sample wave profile */
#define WAVE_PROFILE_Y -1.0      /* value of y to sample wave profile */
#define PROFILE_AT_BOTTOM 1     /* draw wave profile at bottom instead of top */
#define AVERAGE_WAVE_PROFILE 1  /* set to 1 to draw time-average of wave profile squared*/
#define TIMESERIES_NVALUES 400  /* number of values plotted in time series */
#define DRAW_WAVE_SOURCE 0      /* set to 1 to draw source of wave at (wave_source_x, wave_source_y) */
/* end of constants only used by sub_wave and sub_maze */


/* For debugging purposes only */
#define FLOOR 1         /* set to 1 to limit wave amplitude to VMAX */
#define VMAX 10.0       /* max value of wave amplitude */

/* Parameters controlling 3D projection */

double u_3d[2] = {0.75, -0.45};     /* projections of basis vectors for REP_AXO_3D representation */
double v_3d[2] = {-0.75, -0.45};
double w_3d[2] = {0.0, 0.015};
double light[3] = {0.816496581, -0.40824829, 0.40824829};      /* vector of "light" direction for P_3D_ANGLE color scheme */
double observer[3] = {8.0, 6.0, 6.0};    /* location of observer for REP_PROJ_3D representation */ 
int reset_view = 0;         /* switch to reset 3D view parameters (for option ROTATE_VIEW) */

#define Z_SCALING_FACTOR 0.15    /* overall scaling factor of z axis for REP_PROJ_3D representation */
#define XY_SCALING_FACTOR 1.7   /* overall scaling factor for on-screen (x,y) coordinates after projection */
#define ZMAX_FACTOR 1.0         /* max value of z coordinate for REP_PROJ_3D representation */
#define XSHIFT_3D 0.1           /* overall x shift for REP_PROJ_3D representation */
#define YSHIFT_3D 0.1           /* overall y shift for REP_PROJ_3D representation */


#include "global_pdes.c"        /* constants and global variables */
#include "global_3d.c"          /* constants and global variables */
#include "sub_maze.c"           /* support for generating mazes */
#include "sub_wave.c"           /* common functions for wave_billiard, heat and schrodinger */
#include "wave_common.c"        /* common functions for wave_billiard, wave_comparison, etc */

#include "sub_wave_3d.c"        /* graphical functions specific to wave_3d */

FILE *time_series_left, *time_series_right;

double courant2, courantb2;  /* Courant parameters squared */


void evolve_wave_half_new(double phi_in[NX*NY], double psi_in[NX*NY], double phi_out[NX*NY], 
                      short int xy_in[NX*NY], double tc[NX*NY], double tcc[NX*NY], double tgamma[NX*NY], 
                      t_laplacian laplace[NX*NY])
/* time step of field evolution */
/* phi is value of field at time t, psi at time t-1 */
/* this version of the function has been rewritten in order to minimize the number of if-branches */
{
    int i, j, iplus, iminus, jplus, jminus;
    double x, y, c, cc, gamma, tb_shift;
    static long time = 0;
    double *delta;
    
    delta = (double *)malloc(NX*NY*sizeof(double));
    
    /* compute the Laplacian */
    compute_laplacian(phi_in, laplace, delta, xy_in);
    
    time++;
    
    if (OSCILLATE_TOPBOT) tb_shift = (int)((XMAX - XMIN)*(double)NX/(XMAX - XMIN));
    
    /* evolution in the bulk */
    #pragma omp parallel for private(i,j,iplus,iminus,jplus,jminus,x,y)
    for (i=1; i<NX-1; i++){
        for (j=1; j<NY-1; j++){
            if ((TWOSPEEDS)||(xy_in[i*NY+j] != 0)){
                x = phi_in[i*NY+j];
		y = psi_in[i*NY+j];
                
                /* evolve phi */
                phi_out[i*NY+j] = -y + 2*x + tcc[i*NY+j]*delta[i*NY+j] - KAPPA*x - tgamma[i*NY+j]*(x-y);
            }
        }
    }
    
    /* left boundary */
//     if (OSCILLATE_LEFT) for (j=1; j<NY-1; j++) phi_out[j] = AMPLITUDE*cos((double)time*OMEGA);
    if (OSCILLATE_LEFT) 
    {
        for (j=1; j<NY-1; j++) phi_out[j] = oscillating_bc(time, j);
        printf("Boundary condition %.3lg\n", oscillating_bc(time, 0));
    }
    else for (j=1; j<NY-1; j++){
        if ((TWOSPEEDS)||(xy_in[j] != 0)){
            x = phi_in[j];
            y = psi_in[j];
                    
            switch (B_COND) {
                case (BC_DIRICHLET):
                {
                    phi_out[j] = -y + 2*x + tcc[j]*delta[j] - KAPPA*x - tgamma[j]*(x-y);
                    break;
                }
                case (BC_PERIODIC):
                {
                    phi_out[j] = -y + 2*x + tcc[j]*delta[j] - KAPPA*x - tgamma[j]*(x-y);
                    break;
                }
                case (BC_ABSORBING):
                {
                    phi_out[j] = x - tc[j]*(x - phi_in[NY+j]) - KAPPA_SIDES*x - GAMMA_SIDES*(x-y);
                    break;
                }
                case (BC_VPER_HABS):
                {
                    phi_out[j] = x - tc[j]*(x - phi_in[NY+j]) - KAPPA_SIDES*x - GAMMA_SIDES*(x-y);
                    break;
                }
            }
        }
    }
    
    /* right boundary */
    for (j=1; j<NY-1; j++){
        if ((TWOSPEEDS)||(xy_in[(NX-1)*NY+j] != 0)){
            x = phi_in[(NX-1)*NY+j];
            y = psi_in[(NX-1)*NY+j];
                    
            switch (B_COND) {
                case (BC_DIRICHLET):
                {
                    phi_out[(NX-1)*NY+j] = -y + 2*x + tcc[(NX-1)*NY+j]*delta[(NX-1)*NY+j] - KAPPA*x - tgamma[(NX-1)*NY+j]*(x-y);
                    break;
                }
                case (BC_PERIODIC):
                {
                    phi_out[(NX-1)*NY+j] = -y + 2*x + tcc[(NX-1)*NY+j]*delta[(NX-1)*NY+j] - KAPPA*x - tgamma[(NX-1)*NY+j]*(x-y);
                    break;
                }
                case (BC_ABSORBING):
                {
                    phi_out[(NX-1)*NY+j] = x - tc[(NX-1)*NY+j]*(x - phi_in[(NX-2)*NY+j]) - KAPPA_SIDES*x - GAMMA_SIDES*(x-y);
                    break;
                }
                case (BC_VPER_HABS):
                {
                    phi_out[(NX-1)*NY+j] = x - tc[(NX-1)*NY+j]*(x - phi_in[(NX-2)*NY+j]) - KAPPA_SIDES*x - GAMMA_SIDES*(x-y);
                    break;
                }
            }
        }
    }
    
    /* top boundary */
    for (i=0; i<NX; i++){
        if ((TWOSPEEDS)||(xy_in[i*NY+NY-1] != 0)){
            x = phi_in[i*NY+NY-1];
            y = psi_in[i*NY+NY-1];
                    
            if ((OSCILLATE_TOPBOT)&&(i < tb_shift)&&(i<NX-1)&&(i>0))
            {
                iplus = i+1;
                iminus = i-1;   if (iminus < 0) iminus = 0;
                /* TO ADAPT */
                phi_out[i*NY+NY-1] = -y + 2*x + tcc[i*NY+NY-1]*delta[i*NY+NY-1] - KAPPA*x - tgamma[i*NY+NY-1]*(x-y);
            }
            
            else switch (B_COND) {
                case (BC_DIRICHLET):
                {
                    phi_out[i*NY+NY-1] = -y + 2*x + tcc[i*NY+NY-1]*delta[i*NY+NY-1] - KAPPA*x - tgamma[i*NY+NY-1]*(x-y);
                    break;
                }
                case (BC_PERIODIC):
                {
                    phi_out[i*NY+NY-1] = -y + 2*x + tcc[i*NY+NY-1]*delta[i*NY+NY-1] - KAPPA*x - tgamma[i*NY+NY-1]*(x-y);
                    break;
                }
                case (BC_ABSORBING):
                {
                    phi_out[i*NY+NY-1] = x - tc[i*NY+NY-1]*(x - phi_in[i*NY+NY-2]) - KAPPA_TOPBOT*x - GAMMA_TOPBOT*(x-y);
                    break;
                }
                case (BC_VPER_HABS):
                {
                    if (i==0) phi_out[NY-1] = x - tc[NY-1]*(x - phi_in[1*NY+NY-1]) - KAPPA_SIDES*x - GAMMA_SIDES*(x-y);
                    else phi_out[i*NY+NY-1] = -y + 2*x + tcc[i*NY+NY-1]*delta[i*NY+NY-1] - KAPPA*x - tgamma[i*NY+NY-1]*(x-y);
                    break;
                }
            }
        }
    }
    
    /* bottom boundary */
    for (i=0; i<NX; i++){
        if ((TWOSPEEDS)||(xy_in[i*NY] != 0)){
            x = phi_in[i*NY];
            y = psi_in[i*NY];
                    
            if ((OSCILLATE_TOPBOT)&&(i < tb_shift)&&(i<NX-1)&&(i>0))
            {
                /* TO ADAPT */
                phi_out[i*NY] = -y + 2*x + tcc[i*NY]*delta[i*NY] - KAPPA*x - tgamma[i*NY]*(x-y);
            }
            
            else switch (B_COND) {
                case (BC_DIRICHLET):
                {
                    phi_out[i*NY] = -y + 2*x + tcc[i*NY]*delta[i*NY] - KAPPA*x - tgamma[i*NY]*(x-y);
                    break;
                }
                case (BC_PERIODIC):
                {
                    phi_out[i*NY] = -y + 2*x + tcc[i*NY]*delta[i*NY] - KAPPA*x - tgamma[i*NY]*(x-y);
                    break;
                }
                case (BC_ABSORBING):
                {
                    phi_out[i*NY] = x - tc[i*NY]*(x - phi_in[i*NY+1]) - KAPPA_TOPBOT*x - GAMMA_TOPBOT*(x-y);
                    break;
                }
                case (BC_VPER_HABS):
                {
                    if (i==0) phi_out[0] = x - tc[0]*(x - phi_in[NY]) - KAPPA_SIDES*x - GAMMA_SIDES*(x-y);
                    else phi_out[i*NY] = -y + 2*x + tcc[i*NY]*delta[i*NY] - KAPPA*x - tgamma[i*NY]*(x-y);
                    break;
                }
            }
        }
    }
    
    /* add oscillating boundary condition on the left corners */
    if (OSCILLATE_LEFT)
    {
        phi_out[0] = AMPLITUDE*cos((double)time*OMEGA);
        phi_out[NY-1] = AMPLITUDE*cos((double)time*OMEGA);
    }
    
    /* for debugging purposes/if there is a risk of blow-up */
    if (FLOOR) 
    {
        #pragma omp parallel for private(i,j)
        for (i=0; i<NX; i++){
            for (j=0; j<NY; j++){
                if (xy_in[i*NY+j] != 0) 
                {
                    if (phi_out[i*NY+j] > VMAX) phi_out[i*NY+j] = VMAX;
                    if (phi_out[i*NY+j] < -VMAX) phi_out[i*NY+j] = -VMAX;
                }
            }
        }
    }
    
    free(delta);
}


void evolve_wave_half(double phi_in[NX*NY], double psi_in[NX*NY], double phi_out[NX*NY], 
                      short int xy_in[NX*NY], double tc[NX*NY], double tcc[NX*NY], double tgamma[NX*NY])
/* time step of field evolution */
/* phi is value of field at time t, psi at time t-1 */
/* this version of the function has been rewritten in order to minimize the number of if-branches */
{
    int i, j, iplus, iminus, jplus, jminus, jtop, jbot;
    double delta, x, y, c, cc, gamma;
    static long time = 0;
    static short int first = 1;
    static double tb_shift;
    
    time++;
    
    if ((OSCILLATE_TOPBOT)&&(first)) 
    {
        tb_shift = (int)((XMAX - XMIN)*(double)NX/(XMAX - XMIN));
        if (tb_shift >= NX-1) tb_shift = NY - 2;
        first = 0;
    }
    
    #pragma omp parallel for private(i,j,iplus,iminus,jplus,jminus,delta,x,y)
    /* evolution in the bulk */
    for (i=1; i<NX-1; i++){
        for (j=1; j<NY-1; j++){
            if ((TWOSPEEDS)||(xy_in[i*NY+j] != 0)){
                x = phi_in[i*NY+j];
		y = psi_in[i*NY+j];
                
                /* discretized Laplacian */
                delta = phi_in[(i+1)*NY+j] + phi_in[(i-1)*NY+j] + phi_in[i*NY+j+1] + phi_in[i*NY+j-1] - 4.0*x;

                /* evolve phi */
                phi_out[i*NY+j] = -y + 2*x + tcc[i*NY+j]*delta - KAPPA*x - tgamma[i*NY+j]*(x-y);
            }
        }
    }
    
    /* left boundary */
//     if (OSCILLATE_LEFT) for (j=1; j<NY-1; j++) phi_out[j] = AMPLITUDE*cos((double)time*OMEGA);
    if (OSCILLATE_LEFT) 
    {
        for (j=1; j<NY-1; j++) phi_out[j] = oscillating_bc(time, j);
//         printf("Boundary condition %.3lg\n", oscillating_bc(time));
    }
    else for (j=1; j<NY-1; j++){
        if ((TWOSPEEDS)||(xy_in[j] != 0)){
            x = phi_in[j];
            y = psi_in[j];
                    
            switch (B_COND) {
                case (BC_DIRICHLET):
                {
                    delta = phi_in[NY+j] + phi_in[j+1] + phi_in[j-1] - 3.0*x;
                    phi_out[j] = -y + 2*x + tcc[j]*delta - KAPPA*x - tgamma[j]*(x-y);
                    break;
                }
                case (BC_PERIODIC):
                {
                    delta = phi_in[NY+j] + phi_in[(NX-1)*NY+j] + phi_in[j+1] + phi_in[j-1] - 4.0*x;
                    phi_out[j] = -y + 2*x + tcc[j]*delta - KAPPA*x - tgamma[j]*(x-y);
                    break;
                }
                case (BC_ABSORBING):
                {
                    delta = phi_in[NY+j] + phi_in[j+1] + phi_in[j-1] - 3.0*x;
                    phi_out[j] = x - tc[j]*(x - phi_in[NY+j]) - KAPPA_SIDES*x - GAMMA_SIDES*(x-y);
                    break;
                }
                case (BC_VPER_HABS):
                {
                    delta = phi_in[NY+j] + phi_in[j+1] + phi_in[j-1] - 3.0*x;
                    phi_out[j] = x - tc[j]*(x - phi_in[NY+j]) - KAPPA_SIDES*x - GAMMA_SIDES*(x-y);
                    break;
                }
            }
        }
    }
    
    /* right boundary */
    for (j=1; j<NY-1; j++){
        if ((TWOSPEEDS)||(xy_in[(NX-1)*NY+j] != 0)){
            x = phi_in[(NX-1)*NY+j];
            y = psi_in[(NX-1)*NY+j];
                    
            switch (B_COND) {
                case (BC_DIRICHLET):
                {
                    delta = phi_in[(NX-2)*NY+j] + phi_in[(NX-1)*NY+j+1] + phi_in[(NX-1)*NY+j-1] - 3.0*x;
                    phi_out[(NX-1)*NY+j] = -y + 2*x + tcc[(NX-1)*NY+j]*delta - KAPPA*x - tgamma[(NX-1)*NY+j]*(x-y);
                    break;
                }
                case (BC_PERIODIC):
                {
                    delta = phi_in[(NX-2)*NY+j] + phi_in[j] + phi_in[(NX-1)*NY+j+1] + phi_in[(NX-1)*NY+j-1] - 4.0*x;
                    phi_out[(NX-1)*NY+j] = -y + 2*x + tcc[(NX-1)*NY+j]*delta - KAPPA*x - tgamma[(NX-1)*NY+j]*(x-y);
                    break;
                }
                case (BC_ABSORBING):
                {
                    delta = phi_in[(NX-2)*NY+j] + phi_in[(NX-1)*NY+j+1] + phi_in[(NX-1)*NY+j-1] - 3.0*x;
                    phi_out[(NX-1)*NY+j] = x - tc[(NX-1)*NY+j]*(x - phi_in[(NX-2)*NY+j]) - KAPPA_SIDES*x - GAMMA_SIDES*(x-y);
                    break;
                }
                case (BC_VPER_HABS):
                {
                    delta = phi_in[(NX-2)*NY+j] + phi_in[(NX-1)*NY+j+1] + phi_in[(NX-1)*NY+j-1] - 3.0*x;
                    phi_out[(NX-1)*NY+j] = x - tc[(NX-1)*NY+j]*(x - phi_in[(NX-2)*NY+j]) - KAPPA_SIDES*x - GAMMA_SIDES*(x-y);
                    break;
                }
            }
        }
    }
    
    /* top boundary */
    if (COMPARISON) jbot = NY/2;
    else jbot = 0;
    for (i=0; i<NX; i++){
        if ((TWOSPEEDS)||(xy_in[i*NY+NY-1] != 0)){
            x = phi_in[i*NY+NY-1];
            y = psi_in[i*NY+NY-1];
                    
            if ((OSCILLATE_TOPBOT)&&(i < tb_shift)&&(i<NX-1)&&(i>0))
            {
                iplus = i+1;    
                iminus = i-1;   if (iminus < 0) iminus = 0;
                delta = phi_in[iplus*NY+NY-1] + phi_in[iminus*NY+NY-1] + - 2.0*x;
                phi_out[i*NY+NY-1] = -y + 2*x + tcc[i*NY+NY-1]*delta - KAPPA*x - tgamma[i*NY+NY-1]*(x-y);
            }
            
            else switch (B_COND) {
                case (BC_DIRICHLET):
                {
                    iplus = i+1;   if (iplus == NX) iplus = NX-1;
                    iminus = i-1;  if (iminus == -1) iminus = 0;
                    
                    delta = phi_in[iplus*NY+NY-1] + phi_in[iminus*NY+NY-1] + phi_in[i*NY+NY-2] - 3.0*x;
                    phi_out[i*NY+NY-1] = -y + 2*x + tcc[i*NY+NY-1]*delta - KAPPA*x - tgamma[i*NY+NY-1]*(x-y);
                    break;
                }
                case (BC_PERIODIC):
                {
                    iplus = (i+1) % NX;
                    iminus = (i-1) % NX;
                    if (iminus < 0) iminus += NX;
                    
                    delta = phi_in[iplus*NY+NY-1] + phi_in[iminus*NY+NY-1] + phi_in[i*NY+NY-2] + phi_in[i*NY+jbot] - 4.0*x;
                    phi_out[i*NY+NY-1] = -y + 2*x + tcc[i*NY+NY-1]*delta - KAPPA*x - tgamma[i*NY+NY-1]*(x-y);
                    break;
                }
                case (BC_ABSORBING):
                {
                    iplus = (i+1);   if (iplus == NX) iplus = NX-1;
                    iminus = (i-1);  if (iminus == -1) iminus = 0;
                    
                    delta = phi_in[iplus*NY+NY-1] + phi_in[iminus*NY+NY-1] + phi_in[i*NY+NY-2] - 3.0*x;
                    phi_out[i*NY+NY-1] = x - tc[i*NY+NY-1]*(x - phi_in[i*NY+NY-2]) - KAPPA_TOPBOT*x - GAMMA_TOPBOT*(x-y);
                    break;
                }
                case (BC_VPER_HABS):
                {
                    iplus = (i+1);   if (iplus == NX) iplus = NX-1;
                    iminus = (i-1);  if (iminus == -1) iminus = 0;

                    delta = phi_in[iplus*NY+NY-1] + phi_in[iminus*NY+NY-1] + phi_in[i*NY+NY-2] + phi_in[i*NY+jbot] - 4.0*x;
                    if (i==0) phi_out[NY-1] = x - tc[NY-1]*(x - phi_in[1*NY+NY-1]) - KAPPA_SIDES*x - GAMMA_SIDES*(x-y);
                    else phi_out[i*NY+NY-1] = -y + 2*x + tcc[i*NY+NY-1]*delta - KAPPA*x - tgamma[i*NY+NY-1]*(x-y);
                    break;
                }
            }
        }
    }
    
    /* bottom boundary */
    if (COMPARISON) jtop = NY/2-1;
    else jtop = NY-1;
    for (i=0; i<NX; i++){
        if ((TWOSPEEDS)||(xy_in[i*NY] != 0)){
            x = phi_in[i*NY];
            y = psi_in[i*NY];
                    
            if ((OSCILLATE_TOPBOT)&&(i < tb_shift)&&(i<NX-1)&&(i>0))
            {
                iplus = i+1;
                iminus = i-1;   if (iminus < 0) iminus = 0;
                delta = phi_in[iplus*NY] + phi_in[iminus*NY] + - 2.0*x;
                phi_out[i*NY] = -y + 2*x + tcc[i*NY]*delta - KAPPA*x - tgamma[i*NY]*(x-y);
            }
            
            else switch (B_COND) {
                case (BC_DIRICHLET):
                {
                    iplus = i+1;   if (iplus == NX) iplus = NX-1;
                    iminus = i-1;  if (iminus == -1) iminus = 0;
                    
                    delta = phi_in[iplus*NY] + phi_in[iminus*NY] + phi_in[i*NY+1] - 3.0*x;
                    phi_out[i*NY] = -y + 2*x + tcc[i*NY]*delta - KAPPA*x - tgamma[i*NY]*(x-y);
                    break;
                }
                case (BC_PERIODIC):
                {
                    iplus = (i+1) % NX;
                    iminus = (i-1) % NX;
                    if (iminus < 0) iminus += NX;
                    
                    delta = phi_in[iplus*NY] + phi_in[iminus*NY] + phi_in[i*NY+1] + phi_in[i*NY+jtop] - 4.0*x;
                    phi_out[i*NY] = -y + 2*x + tcc[i*NY]*delta - KAPPA*x - tgamma[i*NY]*(x-y);
                    break;
                }
                case (BC_ABSORBING):
                {
                    iplus = (i+1);   if (iplus == NX) iplus = NX-1;
                    iminus = (i-1);  if (iminus == -1) iminus = 0;
                    
                    delta = phi_in[iplus*NY] + phi_in[iminus*NY] + phi_in[i*NY+1] - 3.0*x;
                    phi_out[i*NY] = x - tc[i*NY]*(x - phi_in[i*NY+1]) - KAPPA_TOPBOT*x - GAMMA_TOPBOT*(x-y);
                    break;
                }
                case (BC_VPER_HABS):
                {
                    iplus = (i+1);   if (iplus == NX) iplus = NX-1;
                    iminus = (i-1);  if (iminus == -1) iminus = 0;

                    delta = phi_in[iplus*NY] + phi_in[iminus*NY] + phi_in[i*NY+1] + phi_in[i*NY+jtop] - 4.0*x;
                    if (i==0) phi_out[0] = x - tc[0]*(x - phi_in[NY]) - KAPPA_SIDES*x - GAMMA_SIDES*(x-y);
                    else phi_out[i*NY] = -y + 2*x + tcc[i*NY]*delta - KAPPA*x - tgamma[i*NY]*(x-y);
                    break;
                }
            }
        }
    }
    
    /* case of comparisons - top of bottom half */
    if (COMPARISON) for (i=0; i<NX; i++){
        if ((TWOSPEEDS)||(xy_in[i*NY+jtop] != 0)){
            x = phi_in[i*NY+jtop];
            y = psi_in[i*NY+jtop];
                    
            switch (B_COND) {
                case (BC_DIRICHLET):
                {
                    iplus = i+1;   if (iplus == NX) iplus = NX-1;
                    iminus = i-1;  if (iminus == -1) iminus = 0;
                    
                    delta = phi_in[iplus*NY+jtop] + phi_in[iminus*NY+jtop] + phi_in[i*NY+jtop-1] - 3.0*x;
                    phi_out[i*NY+jtop] = -y + 2*x + tcc[i*NY+jtop]*delta - KAPPA*x - tgamma[i*NY+jtop]*(x-y);
                    break;
                }
                case (BC_PERIODIC):
                {
                    iplus = (i+1) % NX;
                    iminus = (i-1) % NX;
                    if (iminus < 0) iminus += NX;
                    
                    delta = phi_in[iplus*NY+jtop] + phi_in[iminus*NY+jtop] + phi_in[i*NY+jtop-1] + phi_in[i*NY] - 4.0*x;
                    phi_out[i*NY+jtop] = -y + 2*x + tcc[i*NY+jtop]*delta - KAPPA*x - tgamma[i*NY+jtop]*(x-y);
                    break;
                }
                case (BC_ABSORBING):
                {
                    iplus = (i+1);   if (iplus == NX) iplus = NX-1;
                    iminus = (i-1);  if (iminus == -1) iminus = 0;
                    
                    delta = phi_in[iplus*NY+jtop] + phi_in[iminus*NY+jtop] + phi_in[i*NY+jtop-1] - 3.0*x;
                    phi_out[i*NY+jtop] = x - tc[i*NY+jtop]*(x - phi_in[i*NY+jtop-1]) - KAPPA_TOPBOT*x - GAMMA_TOPBOT*(x-y);
                    break;
                }
                case (BC_VPER_HABS):
                {
                    iplus = (i+1);   if (iplus == NX) iplus = NX-1;
                    iminus = (i-1);  if (iminus == -1) iminus = 0;

                    delta = phi_in[iplus*NY+jtop] + phi_in[iminus*NY+jtop] + phi_in[i*NY+jtop-1] + phi_in[i*NY] - 4.0*x;
                    if (i==0) phi_out[jtop] = x - tc[jtop]*(x - phi_in[1*NY+jtop]) - KAPPA_SIDES*x - GAMMA_SIDES*(x-y);
                    else phi_out[i*NY+jtop] = -y + 2*x + tcc[i*NY+jtop]*delta - KAPPA*x - tgamma[i*NY+jtop]*(x-y);
                    break;
                }
            }
        }
    }
    
    /* case of comparisons - bottom of top half */
    if (COMPARISON) for (i=0; i<NX; i++){
        if ((TWOSPEEDS)||(xy_in[i*NY+jbot] != 0)){
            x = phi_in[i*NY+jbot];
            y = psi_in[i*NY+jbot];
                    
            switch (B_COND) {
                case (BC_DIRICHLET):
                {
                    iplus = i+1;   if (iplus == NX) iplus = NX-1;
                    iminus = i-1;  if (iminus == -1) iminus = 0;
                    
                    delta = phi_in[iplus*NY+jbot] + phi_in[iminus*NY+jbot] + phi_in[i*NY+jbot+1] - 3.0*x;
                    phi_out[i*NY+jbot] = -y + 2*x + tcc[i*NY+jbot]*delta - KAPPA*x - tgamma[i*NY+jbot]*(x-y);
                    break;
                }
                case (BC_PERIODIC):
                {
                    iplus = (i+1) % NX;
                    iminus = (i-1) % NX;
                    if (iminus < 0) iminus += NX;
                    
                    delta = phi_in[iplus*NY+jbot] + phi_in[iminus*NY+jbot] + phi_in[i*NY+jbot+1] + phi_in[i*NY+jbot] - 4.0*x;
                    phi_out[i*NY+jbot] = -y + 2*x + tcc[i*NY+jbot]*delta - KAPPA*x - tgamma[i*NY+jbot]*(x-y);
                    break;
                }
                case (BC_ABSORBING):
                {
                    iplus = (i+1);   if (iplus == NX) iplus = NX-1;
                    iminus = (i-1);  if (iminus == -1) iminus = 0;
                    
                    delta = phi_in[iplus*NY+jbot] + phi_in[iminus*NY+jbot] + phi_in[i*NY+jbot+1] - 3.0*x;
                    phi_out[i*NY+jbot] = x - tc[i*NY+jbot]*(x - phi_in[i*NY+jbot+1]) - KAPPA_TOPBOT*x - GAMMA_TOPBOT*(x-y);
                    break;
                }
                case (BC_VPER_HABS):
                {
                    iplus = (i+1);   if (iplus == NX) iplus = NX-1;
                    iminus = (i-1);  if (iminus == -1) iminus = 0;

                    delta = phi_in[iplus*NY+jbot] + phi_in[iminus*NY+jbot] + phi_in[i*NY+jbot+1] + phi_in[i*NY+NY-1] - 4.0*x;
                    if (i==0) phi_out[jbot] = x - tc[jbot]*(x - phi_in[NY]) - KAPPA_SIDES*x - GAMMA_SIDES*(x-y);
                    else phi_out[i*NY+jbot] = -y + 2*x + tcc[i*NY+jbot]*delta - KAPPA*x - tgamma[i*NY+jbot]*(x-y);
                    break;
                }
            }
        }
    }
    
    /* add oscillating boundary condition on the left corners */
    if (OSCILLATE_LEFT)
    {
        phi_out[0] = oscillating_bc(time, 0);
        phi_out[NY-1] = oscillating_bc(time, NY-1);
    }
    
    
    
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
    double tc[NX*NY], double tcc[NX*NY], double tgamma[NX*NY], t_laplacian laplace[NX*NY], t_laplacian laplace1[NX*NY], 
    t_laplacian laplace2[NX*NY])
/* time step of field evolution */
/* phi is value of field at time t, psi at time t-1 */
{
    if (PRECOMPUTE_BC)
    {
//         evolve_wave_half_new(phi, psi, phi_tmp, psi_tmp, xy_in, tc, tcc, tgamma, laplace);
//         evolve_wave_half_new(phi_tmp, psi_tmp, phi, psi, xy_in, tc, tcc, tgamma, laplace_tmp);
        evolve_wave_half_new(phi, psi, tmp, xy_in, tc, tcc, tgamma, laplace);
        evolve_wave_half_new(tmp, phi, psi, xy_in, tc, tcc, tgamma, laplace1);
        evolve_wave_half_new(psi, tmp, phi, xy_in, tc, tcc, tgamma, laplace2);
    }
    else
    {
//         evolve_wave_half(phi, psi, phi_tmp, psi_tmp, xy_in, tc, tcc, tgamma);
//         evolve_wave_half(phi_tmp, psi_tmp, phi, psi, xy_in, tc, tcc, tgamma);
        evolve_wave_half(phi, psi, tmp, xy_in, tc, tcc, tgamma);
        evolve_wave_half(tmp, phi, psi, xy_in, tc, tcc, tgamma);
        evolve_wave_half(psi, tmp, phi, xy_in, tc, tcc, tgamma);
    }
}


void draw_color_bar_palette(int plot, double range, int palette, int fade, double fade_value)
{
    double width = 0.1;
//     double width = 0.14;
//     double width = 0.2;
    
    width *= (double)NX/(double)WINWIDTH;
    
    if (ROTATE_COLOR_SCHEME) 
        draw_color_scheme_palette_3d(-1.0, -0.8, XMAX - 0.1, -1.0, plot, -range, range, palette, fade, fade_value);
    else 
        draw_color_scheme_palette_3d(XMAX - 1.5*width, YMIN + 0.1, XMAX - 0.5*width, YMAX - 0.1, plot, -range, range, palette, fade, fade_value);
}

void viewpoint_schedule(int i)
/* change position of observer */
{
    int j;
    double angle, ca, sa;
    static double observer_initial[3];
    static int first = 1;
    
    if (first)
    {
        for (j=0; j<3; j++) observer_initial[j] = observer[j];
        first = 0;
    }
    
    angle = (ROTATE_ANGLE*DPI/360.0)*(double)i/(double)NSTEPS;
    ca = cos(angle);
    sa = sin(angle);
    observer[0] = ca*observer_initial[0] - sa*observer_initial[1];
    observer[1] = sa*observer_initial[0] + ca*observer_initial[1];
    printf("Angle %.3lg, Observer position (%.3lg, %.3lg, %.3lg)\n", angle, observer[0], observer[1], observer[2]);
}


void animation()
{
    double time, scale, ratio, startleft[2], startright[2], sign = 1.0, r2, xy[2], fade_value, yshift, speed = 0.0, a, b, c, angle = 0.0, lambda1, y, x1, sign1, omega, phase_shift; 
    double *phi, *psi, *tmp, *color_scale, *tc, *tcc, *tgamma;
//     double *total_energy;
    short int *xy_in;
    int i, j, s, sample_left[2], sample_right[2], period = 0, fade, source_counter = 0, k, p, q;
    static int counter = 0, first_source = 1;
    long int wave_value;
    t_wave *wave;
    t_laplacian *laplace, *laplace1, *laplace2;
    t_wave_source wave_source[25];
    
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
    
    courant2 = COURANT*COURANT;
    courantb2 = COURANTB*COURANTB;
    c = COURANT*(XMAX - XMIN)/(double)NX;
//     a = 0.015;
//     b = 0.0003;
//     a = 0.04;
//     b = 0.0018;
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

    /* initialize coordinates of neighbours for discrete Laplacian */
    if (PRECOMPUTE_BC) 
    {
        init_laplacian_coords(laplace, phi);
        init_laplacian_coords(laplace1, tmp);
        init_laplacian_coords(laplace2, psi);
    }
    
    
    init_wave_fields(wave);
    
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
        
    init_wave_flat_mod(phi, psi, xy_in);
//     init_circular_wave_mod(-0.5, 0.0, phi, psi, xy_in);
//     add_circular_wave_mod(1.0, 1.0, 0.0, phi, psi, xy_in);

//     printf("Wave initialized\n");


    /* initialize table of wave speeds/dissipation */
    init_speed_dissipation(xy_in, tc, tcc, tgamma);
    
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
    draw_wave_3d(0, phi, psi, xy_in, wave, ZPLOT, CPLOT, COLOR_PALETTE, 0, 1.0, 1);
//     draw_billiard();
    
    
    if (DRAW_COLOR_SCHEME) draw_color_bar_palette(CPLOT, COLORBAR_RANGE, COLOR_PALETTE, 0, 1.0);

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
        
        draw_wave_3d(0, phi, psi, xy_in, wave, ZPLOT, CPLOT, COLOR_PALETTE, 0, 1.0, 1);
        for (j=0; j<NVID; j++) 
        {
            evolve_wave(phi, psi, tmp, xy_in, tc, tcc, tgamma, laplace, laplace1, laplace2);
            if (SAVE_TIME_SERIES)
            {
                wave_value = (long int)(phi[sample_left[0]*NY+sample_left[1]]*1.0e16);
                fprintf(time_series_left, "%019ld\n", wave_value);
                wave_value = (long int)(phi[sample_right[0]*NY+sample_right[1]]*1.0e16);
                fprintf(time_series_right, "%019ld\n", wave_value);
                if ((j == 0)&&(i%10 == 0)) printf("Frame %i of %i\n", i, NSTEPS);
//                 fprintf(time_series_right, "%.15f\n", phi[sample_right[0]][sample_right[1]]);
            }
//             if (i % 10 == 9) oscillate_linear_wave(0.2*scale, 0.15*(double)(i*NVID + j), -1.5, YMIN, -1.5, YMAX, phi, psi);
        }
        
//         draw_billiard();
        
        if (DRAW_COLOR_SCHEME) draw_color_bar_palette(CPLOT, COLORBAR_RANGE, COLOR_PALETTE, fade, fade_value); 
        
        /* add oscillating waves */
//         if ((ADD_OSCILLATING_SOURCE)&&(i%OSCILLATING_SOURCE_PERIOD == OSCILLATING_SOURCE_PERIOD - 1))
        if ((ADD_OSCILLATING_SOURCE)&&(i%OSCILLATING_SOURCE_PERIOD == 1))
        {
            if (ALTERNATE_OSCILLATING_SOURCE) sign = -sign;
            add_circular_wave_mod(sign, -0.5, 0.0, phi, psi, xy_in);

//             phase_shift = 1/30.0;
// //             phase_shift = 1/60.0;
// 
//             if (first_source) for (k=0; k<25; k++)
//             {
//                 angle = 0.0;
//                 omega = DPI/25.0;
//                 wave_source[k].xc = 0.05*cos((double)k*omega);;
//                 wave_source[k].yc = 0.05*sin((double)k*omega);
//                 wave_source[k].phase = 0.99 - 1.4*sin(0.7*(1.0 + wave_source[k].xc/0.05));
//                 wave_source[k].amp = 1.0;
// //                 if (wave_source[k].phase > 0.0) wave_source[k].sign = 1;
// //                 else wave_source[k].sign = -1;
//                 wave_source[k].sign = 1;
//             }
//             first_source = 0;
//             
//             for (k=0; k<25; k++)
//                 wave_source[k].phase += 1.4*sin(0.7*(1.0 + wave_source[k].xc*cos(angle)/0.05 + wave_source[k].yc*sin(angle)/0.05));
// 
//             angle = DPI*(double)i/(double)NSTEPS;
// //             phase_shift = 0.02 + 0.08*(double)i/(double)NSTEPS;
// //             printf("Phase shift = %.3lg\n", phase_shift);
//             
//             for (k=0; k<25; k++)
//             {
// //                 wave_source[k].phase += 0.07;
//                 wave_source[k].phase += phase_shift;
//                 wave_source[k].phase -= 1.4*sin(0.7*(1.0 + wave_source[k].xc*cos(angle)/0.05 + wave_source[k].yc*sin(angle)/0.05));
//                 
//                 if (wave_source[k].phase > 1.0)
//                 {
//                     add_circular_wave_mod((double)wave_source[k].sign*wave_source[k].amp, wave_source[k].xc, wave_source[k].yc, phi, psi, xy_in);
//                     printf("Adding wave at (%.2lg, %.2lg)\n", wave_source[k].xc, wave_source[k].yc);
//                     wave_source[k].phase -= 1.0;
//                     wave_source[k].sign *= -1;
//                 }
//             }
            
//             for (j=0; j<NPOLY; j++)
//                 add_circular_wave_mod(sign, lambda1*cos(((double)j+0.5)*angle), lambda1*sin(((double)j+0.5)*angle), phi, psi, xy_in);
            
//             p = phased_array_schedule(i);
//             p = 2;
//             y = -1.0;
//             sign1 = sign;
//             printf("p = %i\n", p);
//             for (k=-8; k<9; k++)
//             {
//                 x1 = 0.05*((double)source_counter/(double)p + (double)k);
//                 if ((x1 > 0.083333333*XMIN)&&(x1 < 0.083333333*XMAX)) 
//                 {
//                     add_circular_wave_mod(sign1, x1, y, phi, psi, xy_in);
//                     printf("Adding wave at (%.2lg, %.2lg)\n", x1, y);
//                 }
//                 sign1 = -sign1;
//             }
//             source_counter++;
//             if (p > 0) q = p;
//             else q = -p;
//             if (source_counter >= q) 
//             {
//                 source_counter = 0;
//                 sign = -sign;
//             }

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
                draw_wave_3d(1, phi, psi, xy_in, wave, ZPLOT_B, CPLOT_B, COLOR_PALETTE_B, 0, 1.0, 1);
                if (DRAW_COLOR_SCHEME) draw_color_bar_palette(CPLOT_B, COLORBAR_RANGE_B, COLOR_PALETTE_B, 0, 1.0);  
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
            draw_wave_3d(0, phi, psi, xy_in, wave, ZPLOT, CPLOT, COLOR_PALETTE, 0, 1.0, 1);
            if (DRAW_COLOR_SCHEME) draw_color_bar_palette(CPLOT, COLORBAR_RANGE, COLOR_PALETTE, 0, 1.0);   
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
                draw_wave_3d(0, phi, psi, xy_in, wave, ZPLOT, CPLOT, COLOR_PALETTE, 1, fade_value, 0);
                if (DRAW_COLOR_SCHEME) draw_color_bar_palette(CPLOT, COLORBAR_RANGE, COLOR_PALETTE, 1, fade_value);   
                if (PRINT_SPEED) print_speed_3d(speed, 1, fade_value);
                if (!NO_EXTRA_BUFFER_SWAP) glutSwapBuffers();
                save_frame_counter(NSTEPS + i + 1);
            }
            
            if ((ROTATE_VIEW)&&(ROTATE_VIEW_WHILE_FADE)) 
            {
                viewpoint_schedule(NSTEPS - INITIAL_TIME);
                reset_view = 1;
            }
            draw_wave_3d(1, phi, psi, xy_in, wave, ZPLOT_B, CPLOT_B, COLOR_PALETTE_B, 0, 1.0, 1);
            if (DRAW_COLOR_SCHEME) draw_color_bar_palette(CPLOT_B, COLORBAR_RANGE_B, COLOR_PALETTE_B, 0, 1.0); 
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
                draw_wave_3d(1, phi, psi, xy_in, wave, ZPLOT_B, CPLOT_B, COLOR_PALETTE_B, 1, fade_value, 0);
                if (DRAW_COLOR_SCHEME) draw_color_bar_palette(CPLOT_B, COLORBAR_RANGE_B, COLOR_PALETTE_B, 1, fade_value);   
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
                draw_wave_3d(0, phi, psi, xy_in, wave, ZPLOT, CPLOT, COLOR_PALETTE, 1, fade_value, 0);
                if (DRAW_COLOR_SCHEME) draw_color_bar_palette(CPLOT, COLORBAR_RANGE, COLOR_PALETTE, 1, fade_value); 
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

    init_3d();

    glutDisplayFunc(display);

    glutMainLoop();

    return 0;
}

