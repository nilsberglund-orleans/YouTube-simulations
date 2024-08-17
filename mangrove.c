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

#define MOVIE 0         /* set to 1 to generate movie */
#define SAVE_MEMORY 0   /* set to 1 to save memory when writing tiff images */

/* General geometrical parameters */

#define WINWIDTH 	1280  /* window width */
#define WINHEIGHT 	720   /* window height */

#define NX 640          /* number of grid points on x axis */
#define NY 360          /* number of grid points on y axis */
// #define NX 1280          /* number of grid points on x axis */
// #define NY 720          /* number of grid points on y axis */

#define XMIN -2.0
#define XMAX 2.0	/* x interval */
#define YMIN -1.125
#define YMAX 1.125	/* y interval for 9/16 aspect ratio */

#define JULIA_SCALE 1.0 /* scaling for Julia sets */

/* Choice of the billiard table */

#define B_DOMAIN 20      /* choice of domain shape, see list in global_pdes.c */

// #define CIRCLE_PATTERN 1   /* pattern of circles, see list in global_pdes.c */
#define CIRCLE_PATTERN 8   /* pattern of circles, see list in global_pdes.c */

#define P_PERCOL 0.25       /* probability of having a circle in C_RAND_PERCOL arrangement */
#define NPOISSON 340        /* number of points for Poisson C_RAND_POISSON arrangement */
#define PDISC_FACTOR 3.25   /* controls density of Poisson disc process (default: 3.25) */
#define RANDOM_POLY_ANGLE 0 /* set to 1 to randomize angle of polygons */

#define LAMBDA 0.85	    /* parameter controlling the dimensions of domain */
#define MU 0.03	            /* parameter controlling the dimensions of domain */
#define NPOLY 3             /* number of sides of polygon */
#define APOLY 1.0           /* angle by which to turn polygon, in units of Pi/2 */ 
#define MDEPTH 4            /* depth of computation of Menger gasket */
#define MRATIO 3            /* ratio defining Menger gasket */
#define MANDELLEVEL 1000    /* iteration level for Mandelbrot set */
#define MANDELLIMIT 10.0    /* limit value for approximation of Mandelbrot set */
#define FOCI 1              /* set to 1 to draw focal points of ellipse */
#define NGRIDX 15           /* number of grid point for grid of disks */
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
#define X_SHIFT -0.9        /* x range on which to apply OSCILLATE_TOPBOT */

#define OMEGA 0.00133333333        /* frequency of periodic excitation */
#define K_BC 3.0             /* spatial period of periodic excitation in y direction */
#define KX_BC 10.0           /* spatial period of periodic excitation in x direction */
#define KY_BC 3.3333           /* spatial period of periodic excitation in y direction */
// #define KX_BC 20.0           /* spatial period of periodic excitation in x direction */
// #define KY_BC 6.66666           /* spatial period of periodic excitation in y direction */
// #define OMEGA 0.002        /* frequency of periodic excitation */
// #define K_BC 3.0           /* spatial period of periodic excitation in y direction */
// #define KX_BC 30.0           /* spatial period of periodic excitation in x direction */
// #define KY_BC 10.0           /* spatial period of periodic excitation in y direction */
#define AMPLITUDE 1.0      /* amplitude of periodic excitation */ 
#define COURANT 0.02       /* Courant number */
#define COURANTB 0.015      /* Courant number in medium B */
// #define COURANTB 0.00666   /* Courant number in medium B */
#define GAMMA 3.0e-6          /* damping factor in wave equation */
#define GAMMAB 5.0e-4           /* damping factor in wave equation */
// #define GAMMA 2.0e-6          /* damping factor in wave equation */
// #define GAMMAB 2.5e-4           /* damping factor in wave equation */
#define GAMMA_SIDES 1.0e-4      /* damping factor on boundary */
#define GAMMA_TOPBOT 1.0e-6     /* damping factor on boundary */
#define KAPPA 0.0               /* "elasticity" term enforcing oscillations */
#define KAPPAB 1.0e-6           /* "elasticity" term enforcing oscillations */
#define KAPPA_SIDES 5.0e-4  /* "elasticity" term on absorbing boundary */
#define KAPPA_TOPBOT 0.0    /* "elasticity" term on absorbing boundary */
/* The Courant number is given by c*DT/DX, where DT is the time step and DX the lattice spacing */
/* The physical damping coefficient is given by GAMMA/(DT)^2 */
/* Increasing COURANT speeds up the simulation, but decreases accuracy */
/* For similar wave forms, COURANT^2*GAMMA should be kept constant */

/* Boundary conditions, see list in global_pdes.c  */

#define B_COND 3

/* Parameters for length and speed of simulation */

#define NSTEPS 2000      /* number of frames of movie */
// #define NSTEPS 5500      /* number of frames of movie */
#define NVID 60          /* number of iterations between images displayed on screen */
#define NSEG 100         /* number of segments of boundary */
#define INITIAL_TIME 100    /* time after which to start saving frames */
#define BOUNDARY_WIDTH 2    /* width of billiard boundary */

#define PAUSE 1000         /* number of frames after which to pause */
#define PSLEEP 1         /* sleep time during pause */
#define SLEEP1  1        /* initial sleeping time */
#define SLEEP2  1   /* final sleeping time */
#define END_FRAMES 100   /* number of still frames at end of movie */

/* Parameters of initial condition */

#define INITIAL_AMP 0.2          /* amplitude of initial condition */
#define INITIAL_VARIANCE 0.002   /* variance of initial condition */
#define INITIAL_WAVELENGTH  0.1 /* wavelength of initial condition */

/* Plot type, see list in global_pdes.c  */

#define PLOT 0

/* Color schemes */

#define COLOR_PALETTE 0     /* Color palette, see list in global_pdes.c  */

#define BLACK 1          /* background */

#define COLOR_SCHEME 1   /* choice of color scheme, see list in global_pdes.c  */

#define SCALE 0          /* set to 1 to adjust color scheme to variance of field */
#define SLOPE 1.0        /* sensitivity of color on wave amplitude */
#define PHASE_FACTOR 1.0       /* factor in computation of phase in color scheme P_3D_PHASE */
#define PHASE_SHIFT 0.0      /* shift of phase in color scheme P_3D_PHASE */
#define ATTENUATION 0.0  /* exponential attenuation coefficient of contrast with time */
#define VSHIFT_AMPLITUDE 0.0   /* additional shift for wave amplitude */
#define VSCALE_AMPLITUDE 0.5   /* additional scaling factor for wave amplitude */
#define E_SCALE 2500.0     /* scaling factor for energy representation */
#define LOG_SCALE 1.0    /* scaling factor for energy log representation */
#define LOG_SHIFT 0.0     /* shift of colors on log scale */
#define FLUX_SCALE 1.0e4    /* scaling factor for enegy flux represtnation */
#define RESCALE_COLOR_IN_CENTER 0   /* set to 1 to decrease color intentiy in the center (for wave escaping ring) */

#define COLORHUE 260     /* initial hue of water color for scheme C_LUM */
#define COLORDRIFT 0.0   /* how much the color hue drifts during the whole simulation */
#define LUMMEAN 0.5      /* amplitude of luminosity variation for scheme C_LUM */
#define LUMAMP 0.3       /* amplitude of luminosity variation for scheme C_LUM */
#define HUEMEAN 220.0    /* mean value of hue for color scheme C_HUE */
#define HUEAMP -50.0      /* amplitude of variation of hue for color scheme C_HUE */

/* mangrove properties */

#define MANGROVE_HUE_MIN 180.0      /* color of original mangrove */
#define MANGROVE_HUE_MAX -50.0        /* color of saturated mangrove */
// #define MANGROVE_EMAX 5.0e-3           /* max energy for mangrove to survive */
#define MANGROVE_EMAX 1.5e-3           /* max energy for mangrove to survive */

#define RANDOM_RADIUS 1     /* set to 1 for random circle radius */
#define ERODE_MANGROVES 1   /* set to 1 for mangroves to be eroded */
#define RECOVER_MANGROVES 1 /* set to 1 to allow mangroves to recover */
#define MOVE_MANGROVES 1    /* set to 1 for mobile mangroves */
#define DETACH_MANGROVES 1  /* set to 1 for mangroves to be able to detach */
#define INERTIA 1           /* set to 1 for taking inertia into account */
#define REPELL_MANGROVES 1  /* set to 1 for mangroves to repell each other */
#define DT_MANGROVE 0.1     /* time step for mangrove displacement */
#define KSPRING 0.05        /* spring constant of mangroves */
#define KWAVE 4.0           /* constant in force due to wave gradient */
#define KREPEL 5.0          /* constant in repelling force between mangroves */
#define REPEL_RADIUS 1.1    /* radius in which repelling force acts (in units of mangrove radius) */
#define DXMAX 0.02          /* max displacement of mangrove in one time step */
#define L_DETACH 0.25       /* spring length beyond which mangroves detach */ 
#define DAMP_MANGROVE 0.1   /* damping coefficient of mangroves */
#define MANGROVE_MASS 1.5   /* mass of mangrove of radius MU */

#define HASHX 25    /* size of hashgrid in x direction */
#define HASHY 15    /* size of hashgrid in y direction */
#define HASHMAX 10  /* maximal number of mangroves per hashgrid cell */
#define HASHGRID_PADDING 0.1    /* padding of hashgrid outside simulation window */

#define DRAW_COLOR_SCHEME 0     /* set to 1 to plot the color scheme */
#define COLORBAR_RANGE 8.0    /* scale of color scheme bar */
#define COLORBAR_RANGE_B 12.0    /* scale of color scheme bar for 2nd part */
#define ROTATE_COLOR_SCHEME 0   /* set to 1 to draw color scheme horizontally */

/* For debugging purposes only */
#define FLOOR 1         /* set to 1 to limit wave amplitude to VMAX */
#define VMAX 10.0       /* max value of wave amplitude */

/* the following constants are only used by wave_billiard and wave_3d so far */
#define COMPARISON 0        /* set to 1 to compare two different patterns */
#define B_DOMAIN_B 20       /* second domain shape, for comparisons */
#define CIRCLE_PATTERN_B 0  /* second pattern of circles or polygons */
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
#define RADIUS_FACTOR 0.3   /* controls inner radius for C_RING arrangements */
#define VARIABLE_IOR 0      /* set to 1 for a variable index of refraction */
#define IOR 7               /* choice of index of refraction, see list in global_pdes.c */
#define IOR_TOTAL_TURNS 1.5 /* total angle of rotation for IOR_PERIODIC_WELLS_ROTATING */
#define MANDEL_IOR_SCALE -0.05   /* parameter controlling dependence of IoR on Mandelbrot escape speed */
#define WAVE_PACKET_SOURCE_TYPE 1       /* type of wave packet sources */
#define N_WAVE_PACKETS 15               /* number of wave packets */
#define OSCIL_LEFT_YSHIFT 0.0   /* y-dependence of left oscillation (for non-horizontal waves) */
#define DRAW_WAVE_PROFILE 0     /* set to 1 to draw a profile of the wave */
#define OSCILLATING_SOURCE_PERIOD 20    /* period of oscillating source */
#define MU_B 1.0           /* parameter controlling the dimensions of domain */
#define DRAW_WAVE_PROFILE 0     /* set to 1 to draw a profile of the wave */
#define VERTICAL_WAVE_PROFILE 0 /* set to 1 to draw wave profile vertically */
#define DRAW_WAVE_TIMESERIES 0  /* set to 1 to draw a time series of the wave */
#define WALL_WIDTH 0.1      /* width of wall separating lenses */
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

#include "global_pdes.c"
#include "sub_maze.c"           /* support for generating mazes */
#include "sub_wave.c"
#include "wave_common.c"

double courant2, courantb2;  /* Courant parameters squared */

typedef struct
{
    double xc, yc, radius;      /* center and radius of circle */
    short int active;           /* circle is active */
    double energy;              /* dissipated energy */
    double yc_wrapped;          /* position of circle centers wrapped vertically */
    double anchorx;             /* points moving circles are attached to */
    double anchory;             /* points moving circles are attached to */
    double vx;                  /* x velocity of circles */
    double vy;                  /* y velocity of circles */
    double radius_initial;      /* initial circle radii */
    double mass_inv;            /* inverse of mangrove mass */
    short int attached;         /* has value 1 if the circle is attached to its anchor */
    int hashx;                  /* hash grid positions of mangroves */
    int hashy;                  /* hash grid positions of mangroves */
} t_mangrove;

typedef struct
{
    int number;                 /* total number of mangroves in cell */
    int mangroves[HASHMAX];     /* numbers of mangroves in cell */
} t_hashgrid;

/*********************/
/* animation part    */
/*********************/


void init_bc_phase(double left_bc[NY], double top_bc[NX], double bot_bc[NX])
/* initialize boundary condition phase KX*x + KY*y */
{
    int i, j;
    double xy[2];
    
    for (j=0; j<NY; j++)
    {
        ij_to_xy(0, j, xy);
        left_bc[j] = KX_BC*XMIN + KY_BC*xy[1];
    }
    for (i=0; i<NX; i++)
    {
        ij_to_xy(i, 0, xy);
        bot_bc[i] = KX_BC*xy[0] + KY_BC*YMIN;
        top_bc[i] = KX_BC*xy[0] + KY_BC*YMAX;        
    }
    
}


// void evolve_wave_half_old(double *phi_in[NX], double *psi_in[NX], double *phi_out[NX], double *psi_out[NX], 
//                       short int *xy_in[NX])
// // void evolve_wave_half(phi_in, psi_in, phi_out, psi_out, xy_in)
// /* time step of field evolution */
// /* phi is value of field at time t, psi at time t-1 */
// {
//     int i, j, iplus, iminus, jplus, jminus, tb_shift;
//     double delta, x, y, c, cc, gamma, kappa, phase, phasemin;
//     static long time = 0;
//     static int init_bc = 1;
//     static double left_bc[NY], top_bc[NX], bot_bc[NX];
//     
//     time++;
//     
//     /* initialize boundary condition phase KX*x + KY*y */
//     if ((OSCILLATE_LEFT)&&(init_bc)) 
//     {
//         init_bc_phase(left_bc, top_bc, bot_bc);
//         tb_shift = (int)((X_SHIFT - XMIN)*(double)NX/(XMAX - XMIN));
//         printf("tb_shift %i\n", tb_shift);
//         init_bc = 0;
//     }
//     
//     #pragma omp parallel for private(i,j,iplus,iminus,jplus,jminus,delta,x,y,c,cc,gamma,kappa)
//     for (i=0; i<NX; i++){
//         for (j=0; j<NY; j++){
//             if (xy_in[i][j])
//             {
//                 c = COURANT;
//                 cc = courant2;
//                 gamma = GAMMA;
//                 kappa = KAPPA;
//             }
//             else if (TWOSPEEDS)
//             {
//                 c = COURANTB;
//                 cc = courantb2;
//                 gamma = GAMMAB;
//                 kappa = KAPPAB;
//             }
// 
//             if ((TWOSPEEDS)||(xy_in[i][j])){
//                 /* discretized Laplacian for various boundary conditions */
//                 if ((B_COND == BC_DIRICHLET)||(B_COND == BC_ABSORBING))
//                 {
//                     iplus = (i+1);   if (iplus == NX) iplus = NX-1;
//                     iminus = (i-1);  if (iminus == -1) iminus = 0;
//                     jplus = (j+1);   if (jplus == NY) jplus = NY-1;
//                     jminus = (j-1);  if (jminus == -1) jminus = 0;
//                 }
//                 else if (B_COND == BC_PERIODIC)
//                 {
//                     iplus = (i+1) % NX;
//                     iminus = (i-1) % NX;
//                     if (iminus < 0) iminus += NX;
//                     jplus = (j+1) % NY;
//                     jminus = (j-1) % NY;
//                     if (jminus < 0) jminus += NY;
//                 }
//                 else if (B_COND == BC_VPER_HABS)
//                 {
//                     iplus = (i+1);   if (iplus == NX) iplus = NX-1;
//                     iminus = (i-1);  if (iminus == -1) iminus = 0;
//                     jplus = (j+1) % NY;
//                     jminus = (j-1) % NY;
//                     if (jminus < 0) jminus += NY;
//                 }
//                 
//                 /* imposing linear wave on top and bottom by making Laplacian 1d */
//                 if ((OSCILLATE_TOPBOT)&&(i < tb_shift))
//                 {
//                     if (j == NY-1) 
//                     {
//                         jminus = NY-1;
//                         jplus = NY-1;
//                     }
//                     else if (j == 0) 
//                     {   
//                         jminus = 0;
//                         jplus = 0;
//                     }
//                 }
//                 
//                 delta = phi_in[iplus][j] + phi_in[iminus][j] + phi_in[i][jplus] + phi_in[i][jminus] - 4.0*phi_in[i][j];
// 
//                 x = phi_in[i][j];
// 		y = psi_in[i][j];
// 
//                 /* evolve phi */
//                 if ((B_COND == BC_PERIODIC)||(B_COND == BC_DIRICHLET)) 
//                     phi_out[i][j] = -y + 2*x + cc*delta - kappa*x - gamma*(x-y);
//                 else if (B_COND == BC_ABSORBING)
//                 {
//                     if ((i>0)&&(i<NX-1)&&(j>0)&&(j<NY-1))
//                         phi_out[i][j] = -y + 2*x + cc*delta - kappa*x - gamma*(x-y);
//                 
//                     /* upper border */
//                     else if (j==NY-1) 
//                         phi_out[i][j] = x - c*(x - phi_in[i][NY-2]) - KAPPA_TOPBOT*x - GAMMA_TOPBOT*(x-y);
//                     
//                     /* lower border */
//                     else if (j==0) 
//                         phi_out[i][j] = x - c*(x - phi_in[i][1]) - KAPPA_TOPBOT*x - GAMMA_TOPBOT*(x-y);
//                 
//                     /* right border */
//                     if (i==NX-1) 
//                         phi_out[i][j] = x - c*(x - phi_in[NX-2][j]) - KAPPA_SIDES*x - GAMMA_SIDES*(x-y);
//                     
//                     /* left border */
//                     else if (i==0) 
//                         phi_out[i][j] = x - c*(x - phi_in[1][j]) - KAPPA_SIDES*x - GAMMA_SIDES*(x-y);
//                 }
//                 else if (B_COND == BC_VPER_HABS)
//                 {
//                     if ((i>0)&&(i<NX-1))
//                         phi_out[i][j] = -y + 2*x + cc*delta - kappa*x - gamma*(x-y);
//                 
//                     /* right border */
//                     else if (i==NX-1) 
//                         phi_out[i][j] = x - c*(x - phi_in[NX-2][j]) - KAPPA_SIDES*x - GAMMA_SIDES*(x-y);
//                     
//                     /* left border */
//                     else if (i==0) 
//                         phi_out[i][j] = x - c*(x - phi_in[1][j]) - KAPPA_SIDES*x - GAMMA_SIDES*(x-y);
//                 }
//                 psi_out[i][j] = x;
//                 
//                 /* add oscillating boundary condition on the left */
// //                 if ((i == 0)&&(OSCILLATE_LEFT)) 
// //                 {
// //                     phase =  (double)time*OMEGA - DPI*K_BC*(double)j/(double)NY;
// //                     if (phase < 0.0) phase = 0.0;
// //                     phi_out[i][j] = AMPLITUDE*sin(phase);
// //                 }
// 
//                 /* add oscillating boundary condition on the left */
//                 if (OSCILLATE_LEFT)
//                 {
//                     phasemin = left_bc[0]; 
//                     if (i == 0)
//                     {
//                         phase =  (double)time*OMEGA - left_bc[j] + phasemin;
//                         if (phase < 0.0) phase = 0.0;
//                         phi_out[i][j] = AMPLITUDE*sin(phase);
//                     }
//                     if ((j == 0)&&(i < tb_shift))
//                     {
//                         phase =  (double)time*OMEGA - bot_bc[i] + phasemin;
//                         if (phase < 0.0) phase = 0.0;
//                         phi_out[i][j] = AMPLITUDE*sin(phase);
//                     }
//                     else if ((j == NY-1)&&(i < tb_shift))
//                     {
//                         phase =  (double)time*OMEGA - top_bc[i] + phasemin;
//                         if (phase < 0.0) phase = 0.0;
//                         phi_out[i][j] = AMPLITUDE*sin(phase);
//                     }                    
//                 }
//                 
//                 if (FLOOR)
//                 {
//                     if (phi_out[i][j] > VMAX) phi_out[i][j] = VMAX;
//                     if (phi_out[i][j] < -VMAX) phi_out[i][j] = -VMAX;
//                     if (psi_out[i][j] > VMAX) psi_out[i][j] = VMAX;
//                     if (psi_out[i][j] < -VMAX) psi_out[i][j] = -VMAX;
//                 }
//             }
//         }
//     }
// //     printf("phi(0,0) = %.3lg, psi(0,0) = %.3lg\n", phi[NX/2][NY/2], psi[NX/2][NY/2]);
// }

void evolve_wave_half(double *phi_in[NX], double *psi_in[NX], double *phi_out[NX],
                      short int *xy_in[NX])
/* time step of field evolution */
/* phi is value of field at time t, psi at time t-1 */
/* this version of the function has been rewritten in order to minimize the number of if-branches */
{
    int i, j, iplus, iminus, jplus, jminus, tb_shift;
    double delta, x, y, c, cc, gamma, kappa, phase, phasemin;
    static long time = 0;
    static double tc[NX][NY], tcc[NX][NY], tgamma[NX][NY], left_bc[NY], top_bc[NX], bot_bc[NX];
    static short int first = 1, init_bc = 1;
    
    time++;
    
    /* initialize boundary condition phase KX*x + KY*y */
    if ((OSCILLATE_LEFT)&&(init_bc)) 
    {
        init_bc_phase(left_bc, top_bc, bot_bc);
        tb_shift = (int)((X_SHIFT - XMIN)*(double)NX/(XMAX - XMIN));
        printf("tb_shift %i\n", tb_shift);
        init_bc = 0;
    }
    
    /* initialize tables with wave speeds and dissipation */
//     if (first)
    {
        for (i=0; i<NX; i++){
            for (j=0; j<NY; j++){
                if (xy_in[i][j] != 0)
                {
                    tc[i][j] = COURANT;
                    tcc[i][j] = courant2;
                    if (xy_in[i][j] == 1) tgamma[i][j] = GAMMA;
                    else tgamma[i][j] = GAMMAB;
                }
                else if (TWOSPEEDS)
                {
                    tc[i][j] = COURANTB;
                    tcc[i][j] = courantb2;
                    tgamma[i][j] = GAMMAB;
                }
            }
        }
//         first = 0;
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
            }
        }
    }
    
    /* left boundary */
//     if (OSCILLATE_LEFT) for (j=1; j<NY-1; j++) phi_out[0][j] = AMPLITUDE*cos((double)time*OMEGA);
    if (OSCILLATE_LEFT) for (j=1; j<NY-1; j++) 
    {
        phasemin = left_bc[0]; 
        phase =  (double)time*OMEGA - left_bc[j] + phasemin;
        if (phase < 0.0) phase = 0.0;
        phi_out[0][j] = AMPLITUDE*sin(phase);
    }
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
            }
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
            }
        }
    }
    
    /* top boundary */
    for (i=0; i<NX; i++){
        if ((TWOSPEEDS)||(xy_in[i][NY-1] != 0)){
            x = phi_in[i][NY-1];
            y = psi_in[i][NY-1];
            
            if ((OSCILLATE_TOPBOT)&&(i < tb_shift))
            {
                iplus = i+1;
                iminus = i-1;   if (iminus < 0) iminus = 0;
                delta = phi_in[iplus][NY-1] + phi_in[iminus][NY-1] + - 2.0*x;
                phi_out[i][NY-1] = -y + 2*x + tcc[i][NY-1]*delta - KAPPA*x - tgamma[i][NY-1]*(x-y);
            }
            else if ((OSCILLATE_LEFT)&&(i < tb_shift))
            {
                phasemin = left_bc[0]; 
                phase =  (double)time*OMEGA - top_bc[i] + phasemin;
                if (phase < 0.0) phase = 0.0;
                phi_out[i][NY-1] = AMPLITUDE*sin(phase);
            }
            else switch (B_COND) {
                case (BC_DIRICHLET):
                {
                    iplus = (i+1);   if (iplus == NX) iplus = NX-1;
                    iminus = (i-1);  if (iminus == -1) iminus = 0;
                    
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

                    delta = phi_in[iplus][NY-1] + phi_in[iminus][NY-1] + phi_in[i][NY-2] + phi_in[i][0] - 4.0*x;
                    phi_out[i][NY-1] = -y + 2*x + tcc[i][NY-1]*delta - KAPPA*x - tgamma[i][NY-1]*(x-y);
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
                    
            if ((OSCILLATE_TOPBOT)&&(i < tb_shift))
            {
                iplus = i+1;
                iminus = i-1;   if (iminus < 0) iminus = 0;
                delta = phi_in[iplus][0] + phi_in[iminus][0] + - 2.0*x;
                phi_out[i][0] = -y + 2*x + tcc[i][0]*delta - KAPPA*x - tgamma[i][0]*(x-y);
            }
            else if ((OSCILLATE_LEFT)&&(i < tb_shift))
            {
                phasemin = left_bc[0]; 
                phase =  (double)time*OMEGA - bot_bc[i] + phasemin;
                if (phase < 0.0) phase = 0.0;
                phi_out[i][0] = AMPLITUDE*sin(phase);
            }
            else switch (B_COND) {
                case (BC_DIRICHLET):
                {
                    iplus = (i+1);   if (iplus == NX) iplus = NX-1;
                    iminus = (i-1);  if (iminus == -1) iminus = 0;
                    
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

                    delta = phi_in[iplus][0] + phi_in[iminus][0] + phi_in[i][1] + phi_in[i][NY-1] - 4.0*x;
                    phi_out[i][0] = -y + 2*x + tcc[i][0]*delta - KAPPA*x - tgamma[i][0]*(x-y);
                    break;
                }
            }
        }
    }
    
    /* add oscillating boundary condition on the left corners - NEEDED ? */
    if ((i == 0)&&(OSCILLATE_LEFT))
    {
        phi_out[i][0] = AMPLITUDE*cos((double)time*OMEGA);
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
}


void evolve_wave(double *phi[NX], double *psi[NX], double *tmp[NX], short int *xy_in[NX])
/* time step of field evolution */
/* phi is value of field at time t, psi at time t-1 */
{
    evolve_wave_half(phi, psi, tmp, xy_in);
    evolve_wave_half(tmp, phi, psi, xy_in);
    evolve_wave_half(psi, tmp, phi, xy_in);
}


void hash_xy_to_ij(double x, double y, int ij[2])
{
    static int first = 1;
    static double lx, ly;
    int i, j;
    
    if (first)
    {
        lx = XMAX - XMIN + 2.0*HASHGRID_PADDING;
        ly = YMAX - YMIN + 2.0*HASHGRID_PADDING;
        first = 0;
    }
    
    i = (int)((double)HASHX*(x - XMIN + HASHGRID_PADDING)/lx);
    j = (int)((double)HASHY*(y - YMIN + HASHGRID_PADDING)/ly);
    
    if (i<0) i = 0;
    if (i>=HASHX) i = HASHX-1;
    if (j<0) j = 0;
    if (j>=HASHY) j = HASHY-1;
    
    ij[0] = i;
    ij[1] = j;

//     printf("Mapped (%.3f,%.3f) to (%i, %i)\n", x, y, ij[0], ij[1]);
}


void compute_repelling_force(int i, int j, double force[2], t_mangrove* mangrove)
/* compute repelling force of mangrove j on mangrove i */
{
    double x1, y1, x2, y2, distance, r, f;
    
    x1 = mangrove[i].xc;
    y1 = mangrove[i].yc;
    x2 = mangrove[j].xc;
    y2 = mangrove[j].yc;
    
    distance = module2(x2 - x1, y2 - y1);
    r = mangrove[i].radius + mangrove[j].radius;
    if (r <= 0.0) r = 0.001*MU;
    f = KREPEL/(0.001 + distance*distance);
    
    if ((distance > 0.0)&&(distance < REPEL_RADIUS*r))
    {
        force[0] = f*(x1 - x2)/distance;
        force[1] = f*(y1 - y2)/distance;
    }
    else
    {
        force[0] = 0.0;
        force[1] = 0.0;
    }
}


void update_hashgrid(t_mangrove* mangrove, int* hashgrid_number, int* hashgrid_mangroves)
{
    int i, j, k, n, m, ij[2], max = 0;
    
    printf("Updating hashgrid_number\n");
    for (i=0; i<HASHX*HASHY; i++) hashgrid_number[i] = 0;
    printf("Updated hashgrid_number\n");
        
    /* place each mangrove in hash grid */
    for (k=1; k<ncircles; k++)
//         if (circleactive[k])
        {
//             printf("placing circle %i\t", k);
            hash_xy_to_ij(mangrove[k].xc, mangrove[k].yc, ij);
            i = ij[0];  j = ij[1];
//             printf("ij = (%i, %i)\t", i, j);
            n = hashgrid_number[i*HASHY + j];
            m = i*HASHY*HASHMAX + j*HASHMAX + n;
//             printf("n = %i, m = %i\n", n, m);
            if (m < HASHX*HASHY*HASHMAX) hashgrid_mangroves[m] = k;
            else printf("Too many mangroves in hash cell, try increasing HASHMAX\n");
            hashgrid_number[i*HASHY + j]++;
            mangrove[k].hashx = i;
            mangrove[k].hashy = j;
            
            if (n > max) max = n;
//                 printf("Placed mangrove %i at (%i,%i) in hashgrid\n", k, ij[0], ij[1]);
//                 printf("%i mangroves at (%i,%i)\n", hashgrid_number[ij[0]][ij[1]], ij[0], ij[1]);
        }
    
    printf("Maximal number of mangroves per hash cell: %i\n", max);
}


void animation()
{
    double time, scale, diss, rgb[3], hue, y, dissip, ej, gradient[2], dx, dy, dt, xleft, xright, 
            length, fx, fy, force[2];
    double *phi[NX], *psi[NX], *tmp[NX];
    short int *xy_in[NX], redraw = 0;
    int i, j, k, n, s, ij[2], i0, iplus, iminus, j0, jplus, jminus, p, q;
    static int imin, imax;
    static short int first = 1;
    t_mangrove *mangrove;
    int *hashgrid_number, *hashgrid_mangroves;  
    t_hashgrid *hashgrid;

    /* Since NX and NY are big, it seemed wiser to use some memory allocation here */
    for (i=0; i<NX; i++)
    {
        phi[i] = (double *)malloc(NY*sizeof(double));
        psi[i] = (double *)malloc(NY*sizeof(double));
        tmp[i] = (double *)malloc(NY*sizeof(double));
        xy_in[i] = (short int *)malloc(NY*sizeof(short int));
    }
    mangrove = (t_mangrove *)malloc(NMAXCIRCLES*sizeof(t_mangrove));    /* mangroves */  
    
    hashgrid = (t_hashgrid *)malloc(HASHX*HASHY*sizeof(t_hashgrid));    /* hashgrid */      
        
    hashgrid_number = (int *)malloc(HASHX*HASHY*sizeof(int));   /* total number of mangroves in each hash grid cell */
    hashgrid_mangroves = (int *)malloc(HASHX*HASHY*HASHMAX*sizeof(int)); /* numbers of mangroves in each hash grid cell */
    
    /* initialise positions and radii of circles */
    if ((B_DOMAIN == D_CIRCLES)||(B_DOMAIN == D_CIRCLES_IN_RECT)) ncircles = init_circle_config(circles);
    else if (B_DOMAIN == D_POLYGONS) ncircles = init_polygon_config(polygons);

    

    courant2 = COURANT*COURANT;
    courantb2 = COURANTB*COURANTB;
//     dt = 0.01;

    /* initialize wave with a drop at one point, zero elsewhere */
    init_wave_flat(phi, psi, xy_in);
    
    /* initialise mangroves */
    for (i=0; i < ncircles; i++) 
    {
        /* to avoid having to recode init_circle_config, would be more elegant in C++ */
        mangrove[i].xc = circles[i].xc;
        mangrove[i].yc = circles[i].yc;
        mangrove[i].radius = circles[i].radius;
        mangrove[i].active = circles[i].active;

        mangrove[i].energy = 0.0;
        y = mangrove[i].yc;
        if (y >= YMAX) y -= mangrove[i].radius;
        if (y <= YMIN) y += mangrove[i].radius;
//         if (y >= YMAX) y -= (YMAX - YMIN);
//         if (y <= YMIN) y += (YMAX - YMIN);
        mangrove[i].yc_wrapped = y;
//         mangrove[i].active = 1;
        
        if (RANDOM_RADIUS) mangrove[i].radius = mangrove[i].radius*(0.75 + 0.5*((double)rand()/RAND_MAX));
        
        mangrove[i].radius_initial = mangrove[i].radius;
        mangrove[i].attached = 1; 
        mangrove[i].mass_inv = MU*MU/(MANGROVE_MASS*mangrove[i].radius*mangrove[i].radius);
        
        if (MOVE_MANGROVES)
        {
            mangrove[i].anchorx = mangrove[i].xc;
            mangrove[i].anchory = mangrove[i].yc_wrapped;
//             mangrove[i].anchory = mangrove[i].yc;
        }
        
        if (INERTIA)
        {
            mangrove[i].vx = 0.0;
            mangrove[i].vy = 0.0;
        }
    }
    
    /* initialise hash table for interacting mangroves */
    if (REPELL_MANGROVES) update_hashgrid(mangrove, hashgrid_number, hashgrid_mangroves);
    
    if (first) /* compute box limits where circles are reset */
    {
        /* find leftmost and rightmost circle */
        for (i=0; i<ncircles; i++) 
            if ((mangrove[i].active)&&(mangrove[i].xc - mangrove[i].radius < xleft)) xleft = mangrove[i].xc - mangrove[i].radius; 
        for (i=0; i<ncircles; i++) 
            if ((mangrove[i].active)&&(mangrove[i].xc + mangrove[i].radius > xright)) xright = mangrove[i].xc + mangrove[i].radius; 
        
        xy_to_ij(xleft, 0.0, ij);
        imin = ij[0] - 10;
        if (imin < 0) imin = 0;
        xy_to_ij(xright, 0.0, ij);
        imax = ij[0];
        if (imax >= NX) imax = NX-1;
        first = 0;
        
        printf("xleft = %.3lg, xright = %.3lg, imin = %i, imax = %i\n", xleft, xright, imin, imax);
    }

    blank();
    glColor3f(0.0, 0.0, 0.0);
    draw_wave(phi, psi, xy_in, 1.0, 0, PLOT);
    draw_billiard(0, 1.0);

    glutSwapBuffers();



    sleep(SLEEP1);

    for (i=0; i<=INITIAL_TIME + NSTEPS; i++)
    {
	printf("Computing frame %d\n",i);
        /* compute the variance of the field to adjust color scheme */
        /* the color depends on the field divided by sqrt(1 + variance) */
        if (SCALE)
        {
            scale = sqrt(1.0 + compute_variance(phi,psi, xy_in));
//             printf("Scaling factor: %5lg\n", scale);
        }
        else scale = 1.0;

        printf("Drawing wave\n");
        draw_wave(phi, psi, xy_in, scale, i, PLOT);


        printf("Evolving wave\n");
        for (j=0; j<NVID; j++) 
        {
//             printf("%i ", j);
            evolve_wave(phi, psi, tmp, xy_in);
//             if (i % 10 == 9) oscillate_linear_wave(0.2*scale, 0.15*(double)(i*NVID + j), -1.5, YMIN, -1.5, YMAX, phi, psi);
        }
        
        
        /* move mangroves */
        if (MOVE_MANGROVES) for (j=0; j<ncircles; j++) if (mangrove[j].active)
        {
            compute_gradient(phi, psi, mangrove[j].xc, mangrove[j].yc_wrapped, gradient);
//             printf("gradient = (%.3lg, %.3lg)\t", gradient[0], gradient[1]);
            
//             if (j%NGRIDY == 0) printf("gradient (%.3lg, %.3lg)\n", gradient[0], gradient[1]);
//             if (j%NGRIDY == 0) printf("circle %i (%.3lg, %.3lg) -> ", j, mangrove[j].xc, mangrove[j].yc);
            
            /* compute force of wave */
            dx = DT_MANGROVE*KWAVE*gradient[0];
            dy = DT_MANGROVE*KWAVE*gradient[1];

            /* compute force of spring */
            if (mangrove[j].attached)
            {
                dx += DT_MANGROVE*(-KSPRING*(mangrove[j].xc - mangrove[j].anchorx));
                dy += DT_MANGROVE*(-KSPRING*(mangrove[j].yc_wrapped - mangrove[j].anchory));
            }
            
            /* compute repelling force from other mangroves */
            if (REPELL_MANGROVES)
            {
                /* determine neighboring grid points */
                i0 = mangrove[j].hashx;
                iminus = i0 - 1;    if (iminus < 0) iminus = 0;
                iplus = i0 + 1;     if (iplus >= HASHX) iplus = HASHX-1;

                j0 = mangrove[j].hashy;
                jminus = j0 - 1;    if (jminus < 0) jminus = 0;
                jplus = j0 + 1;     if (jplus >= HASHY) jplus = HASHY-1;
                
                fx = 0.0; 
                fy = 0.0;
                for (p=iminus; p<= iplus; p++)
                    for (q=jminus; q<= jplus; q++)
                        for (k=0; k<hashgrid_number[p*HASHY+q]; k++) 
                            if (mangrove[hashgrid_mangroves[p*HASHY*HASHMAX + q*HASHMAX + k]].active)
                            {
                                compute_repelling_force(j, hashgrid_mangroves[p*HASHY*HASHMAX + q*HASHMAX + k], force, mangrove);
                                fx += force[0];
                                fy += force[1];
                            }
                        
//                 if (fx*fx + fy*fy > 0.001) printf("Force on mangrove %i: (%.3f, %.3f)\n", j, fx, fy);
                        
                dx += DT_MANGROVE*fx;
                dy += DT_MANGROVE*fy;
            }
            
            /* detach mangrove if spring is too long */
            if (DETACH_MANGROVES)
            {
                length = module2(mangrove[j].xc - mangrove[j].anchorx, mangrove[j].yc_wrapped - mangrove[j].anchory);
//                 if (j%NGRIDY == 0) printf("spring length %.i:  %.3lg\n", j, length);
//                 if (length > L_DETACH) mangrove[j].attached = 0;
                if (length*mangrove[j].mass_inv > L_DETACH) mangrove[j].attached = 0;
            }
            
            if (dx > DXMAX) dx = DXMAX;
            if (dx < -DXMAX) dx = -DXMAX;
            if (dy > DXMAX) dy = DXMAX;
            if (dy < -DXMAX) dy = -DXMAX;
            
            if (INERTIA)
            {
                mangrove[j].vx += (dx - DAMP_MANGROVE*mangrove[j].vx)*mangrove[j].mass_inv;
                mangrove[j].vy += (dy - DAMP_MANGROVE*mangrove[j].vy)*mangrove[j].mass_inv;
                mangrove[j].xc += mangrove[j].vx*DT_MANGROVE;
                mangrove[j].yc += mangrove[j].vy*DT_MANGROVE;
                mangrove[j].yc_wrapped += mangrove[j].vy*DT_MANGROVE;
//                 if (j%NGRIDY == 0) 
//                     printf("circle %.i: (dx,dy) = (%.3lg,%.3lg), (vx,vy) = (%.3lg,%.3lg)\n", 
//                            j, mangrove[j].xc-mangrove[j].anchorx, mangrove[j].yc-mangrove[j].anchory, mangrove[j].vx, mangrove[j].vy);
            }
            else
            {
                mangrove[j].xc += dx*mangrove[j].mass_inv*DT_MANGROVE;
                mangrove[j].yc += dy*mangrove[j].mass_inv*DT_MANGROVE;
                mangrove[j].yc_wrapped += dy*mangrove[j].mass_inv*DT_MANGROVE;
            }
            
            if (mangrove[j].xc <= XMIN) mangrove[j].xc = XMIN;
            if (mangrove[j].xc >= XMAX) mangrove[j].xc = XMAX;
            if (mangrove[j].yc_wrapped <= YMIN) mangrove[j].yc_wrapped = YMIN;
            if (mangrove[j].yc_wrapped >= YMAX) mangrove[j].yc_wrapped = YMAX;
        
//             if (j%NGRIDY == 0) printf("(%.3lg, %.3lg)\n", mangrove[j].xc, mangrove[j].yc);
                        
            redraw = 1;
        }
        
        /* test for debugging */
        if (1) for (j=0; j<ncircles; j++) 
        {
            dissip = compute_dissipation(phi, psi, xy_in, mangrove[j].xc, mangrove[j].yc_wrapped);
            
            /* make sure the dissipation does not grow too fast because of round-off/blow-up */
            if (dissip > 0.1*MANGROVE_EMAX) 
            {
                dissip = 0.1*MANGROVE_EMAX;
                printf("Flooring dissipation!\n");
            }
            
            if (mangrove[j].active) 
            {
                mangrove[j].energy += dissip;
                ej = mangrove[j].energy;
//                 printf("ej = %.3f\n", ej);
                if (ej <= MANGROVE_EMAX)
                {
                    if (ej > 0.0) 
                    {
                        hue = MANGROVE_HUE_MIN + (MANGROVE_HUE_MAX - MANGROVE_HUE_MIN)*ej/MANGROVE_EMAX;
                        if (hue < 0.0) hue += 360.0;
                    }
                    else hue = MANGROVE_HUE_MIN;
                    hsl_to_rgb(hue, 0.9, 0.5, rgb);
//                     if (j%NGRIDY == 0) printf("Circle %i, energy %.5lg, hue %.5lg\n", j, ej, hue);
                    draw_colored_circle(mangrove[j].xc, mangrove[j].yc, mangrove[j].radius, NSEG, rgb);
                    
                    /* shrink mangrove */
                    if ((ERODE_MANGROVES)&&(ej > 0.0))
                    {
                        mangrove[j].radius = mangrove[j].radius_initial*(1.0 - ej*ej/(MANGROVE_EMAX*MANGROVE_EMAX));
                        redraw = 1;
                    }
                    else mangrove[j].radius = mangrove[j].radius_initial;

                }
                else    /* remove mangrove */
                {
                    mangrove[j].active = 0;
                    /* reinitialize table xy_in */
                    redraw = 1;
                } 
            }
            else if (RECOVER_MANGROVES)   /* allow disabled mangroves to recover */
            {
                mangrove[j].energy -= 0.15*dissip;
                printf("Circle %i energy %.3lg\n", j, mangrove[j].energy);
                if (mangrove[j].energy < 0.0)
                {
                    printf("Reactivating circle %i?\n", j);
                    /* THE PROBLEM occurs when circleactive[0] is set to 1 again */
                    if (j>0) mangrove[j].active = 1;
                    mangrove[j].radius = mangrove[j].radius_initial;
                    mangrove[j].energy = -MANGROVE_EMAX;
                    /* reinitialize table xy_in */
                    redraw = 1;
                }   
                
            }
        }
        
    /* for compatibility with draw_billiard, may be improvable */
    for (j=0; j<ncircles; j++)
    {
            circles[j].xc = mangrove[j].xc;
            circles[j].yc = mangrove[j].yc;
            circles[j].radius = mangrove[j].radius;
    }
        
        /* compute energy dissipated in obstacles */
/*        if (ERODE_MANGROVES) for (j=0; j<ncircles; j++) 
        {
//             printf("j = %i\t", j);
            dissip = compute_dissipation(phi, psi, xy_in, mangrove[j].xc, mangrove[j].yc_wrapped);
            printf("dissip = %.3f\t", dissip);
            
            /* make sure the dissipation does not grow too fast because of round-off/blow-up */
//             if (dissip > 0.1*MANGROVE_EMAX) 
//             {
//                 dissip = 0.1*MANGROVE_EMAX;
//                 printf("Flooring dissipation!\n");
//             }
//             
//             if (mangrove[j].active) 
//             {
//                 mangrove[j].energy += dissip;
//                 ej = mangrove[j].energy;
//                 printf("ej = %.3f\n", ej);
//                 if (ej <= MANGROVE_EMAX)
//                 {
//                     if (ej > 0.0) 
//                     {
//                         hue = MANGROVE_HUE_MIN + (MANGROVE_HUE_MAX - MANGROVE_HUE_MIN)*ej/MANGROVE_EMAX;
//                         if (hue < 0.0) hue += 360.0;
//                     }
//                     else hue = MANGROVE_HUE_MIN;
//                     hsl_to_rgb(hue, 0.9, 0.5, rgb);
//                     if (j%NGRIDY == 0) printf("Circle %i, energy %.5lg, hue %.5lg\n", j, ej, hue);
//                     draw_colored_circle(mangrove[j].xc, mangrove[j].yc, mangrove[j].radius, NSEG, rgb);
//                     
//                     /* shrink mangrove */
//                     if (ej > 0.0)
//                     {
//                         mangrove[j].radius -= MU*ej*ej/(MANGROVE_EMAX*MANGROVE_EMAX);
//                         if (mangrove[j].radius < 0.0) mangrove[j].radius = 0.0;
//                         mangrove[j].radius = mangrove[j].radius_initial*(1.0 - ej*ej/(MANGROVE_EMAX*MANGROVE_EMAX));
//                         redraw = 1;
//                     }
//                     else mangrove[j].radius = mangrove[j].radius_initial;
//                 }
//                 else    /* remove mangrove */
//                 {
//                     mangrove[j].active = 0;
                    /* reinitialize table xy_in */
//                     redraw = 1;
//                 }   
//             }
//             else    /* allow disabled mangroves to recover */
//             {
//                 mangrove[j].energy -= 0.15*dissip;
//                 printf("ej = %.3f\n", mangrove[j].energy);
//                 mangrove[j].radius += 0.005*MU;
//                 if (mangrove[j].radius > MU) mangrove[j].radius = MU;
//                 if ((mangrove[j].energy < 0.0)&&(mangrove[j].radius > 0.0))
//                 if (mangrove[j].energy < 0.0)
//                 {
//                     mangrove[j].active = 1;
//                     mangrove[j].radius = mangrove[j].radius*(0.75 + 0.5*((double)rand()/RAND_MAX));
//                     mangrove[j].radius = mangrove[j].radius_initial;
//                     mangrove[j].energy = -MANGROVE_EMAX;
                    /* reinitialize table xy_in */
//                     redraw = 1;
//                 }   
                
//             }
            
//             printf("Circle %i, energy %.5lg\n", j, mangrove[j].energy);
//         }

        printf("Updating hashgrid\n");
        if (REPELL_MANGROVES) update_hashgrid(mangrove, hashgrid_number, hashgrid_mangroves);

        
        printf("Drawing billiard\n");
        draw_billiard(0, 1.0);

	glutSwapBuffers();

        if (redraw) 
        {
            printf("Reinitializing xy_in\n");
            init_xyin_xrange(xy_in, imin, NX);
//             init_xyin_xrange(xy_in, imin, imax);
        }
        redraw = 0;

	if (MOVIE)
        {
            if (i >= INITIAL_TIME) save_frame();
            else printf("Initial phase time %i of %i\n", i, INITIAL_TIME);

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
        for (i=0; i<END_FRAMES; i++) save_frame();
        s = system("mv wave*.tif tif_wave/");
    }
    for (i=0; i<NX; i++)
    {
        free(phi[i]);
        free(psi[i]);
        free(tmp[i]);
        free(xy_in[i]);
    }
    free(mangrove);
        
    free(hashgrid_number);
    free(hashgrid_mangroves);
}


void display(void)
{
    glPushMatrix();

    blank();
    glutSwapBuffers();
    blank();
    glutSwapBuffers();

    animation();
    sleep(SLEEP2);

    glPopMatrix();

    glutDestroyWindow(glutGetWindow());

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
