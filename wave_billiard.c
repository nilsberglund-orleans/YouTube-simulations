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

#define MOVIE 1         /* set to 1 to generate movie */
#define DOUBLE_MOVIE 1  /* set to 1 to produce movies for wave height and energy simultaneously */

/* General geometrical parameters */

#define WINWIDTH 	1920  /* window width */
#define WINHEIGHT 	1000  /* window height */
#define NX 1920          /* number of grid points on x axis */
#define NY 1000          /* number of grid points on y axis */

#define XMIN -2.0
#define XMAX 2.0	/* x interval */
#define YMIN -1.0329
#define YMAX 1.1129	/* y interval for 9/16 aspect ratio */

// #define WINWIDTH 	1280  /* window width */
// #define WINHEIGHT 	720   /* window height */

// #define NX 1280          /* number of grid points on x axis */
// #define NY 720          /* number of grid points on y axis */
// #define NX 2560          /* number of grid points on x axis */
// #define NY 1440          /* number of grid points on y axis */

// #define YMIN -1.125
// #define YMAX 1.125	/* y interval for 9/16 aspect ratio */

#define JULIA_SCALE 1.0 /* scaling for Julia sets */

/* Choice of the billiard table */

#define B_DOMAIN 42       /* choice of domain shape, see list in global_pdes.c */

#define CIRCLE_PATTERN 2   /* pattern of circles or polygons, see list in global_pdes.c */

#define P_PERCOL 0.25       /* probability of having a circle in C_RAND_PERCOL arrangement */
#define NPOISSON 300        /* number of points for Poisson C_RAND_POISSON arrangement */
#define RANDOM_POLY_ANGLE 1 /* set to 1 to randomize angle of polygons */

#define LAMBDA 1.1	    /* parameter controlling the dimensions of domain */
#define MU 0.4              /* parameter controlling the dimensions of domain */
#define NPOLY 14             /* number of sides of polygon */
#define APOLY -1.0           /* angle by which to turn polygon, in units of Pi/2 */ 
#define MDEPTH 5            /* depth of computation of Menger gasket */
#define MRATIO 3            /* ratio defining Menger gasket */
#define MANDELLEVEL 1000    /* iteration level for Mandelbrot set */
#define MANDELLIMIT 10.0    /* limit value for approximation of Mandelbrot set */
#define FOCI 1              /* set to 1 to draw focal points of ellipse */
#define NGRIDX 10           /* number of grid point for grid of disks */
#define NGRIDY 12           /* number of grid point for grid of disks */

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

#define OMEGA 0.002        /* frequency of periodic excitation */
#define AMPLITUDE 1.0      /* amplitude of periodic excitation */ 
#define COURANT 0.04       /* Courant number */
#define COURANTB 0.01      /* Courant number in medium B */
#define GAMMA 0.0          /* damping factor in wave equation */
#define GAMMAB 1.0e-6           /* damping factor in wave equation */
#define GAMMA_SIDES 1.0e-4      /* damping factor on boundary */
#define GAMMA_TOPBOT 1.0e-7     /* damping factor on boundary */
#define KAPPA 0.0           /* "elasticity" term enforcing oscillations */
#define KAPPA_SIDES 5.0e-4  /* "elasticity" term on absorbing boundary */
#define KAPPA_TOPBOT 0.0    /* "elasticity" term on absorbing boundary */
/* The Courant number is given by c*DT/DX, where DT is the time step and DX the lattice spacing */
/* The physical damping coefficient is given by GAMMA/(DT)^2 */
/* Increasing COURANT speeds up the simulation, but decreases accuracy */
/* For similar wave forms, COURANT^2*GAMMA should be kept constant */

/* Boundary conditions, see list in global_pdes.c  */

#define B_COND 3

/* Parameters for length and speed of simulation */

#define NSTEPS 2100      /* number of frames of movie */
#define NVID 25          /* number of iterations between images displayed on screen */
#define NSEG 100         /* number of segments of boundary */
#define INITIAL_TIME 0      /* time after which to start saving frames */
#define BOUNDARY_WIDTH 3    /* width of billiard boundary */

#define PAUSE 1000       /* number of frames after which to pause */
#define PSLEEP 1         /* sleep time during pause */
#define SLEEP1  1        /* initial sleeping time */
#define SLEEP2  1        /* final sleeping time */
#define MID_FRAMES 20    /* number of still frames between parts of two-part movie */
#define END_FRAMES 50    /* number of still frames at end of movie */

/* Parameters of initial condition */

#define INITIAL_AMP 0.75         /* amplitude of initial condition */
#define INITIAL_VARIANCE 0.0003  /* variance of initial condition */
#define INITIAL_WAVELENGTH  0.015  /* wavelength of initial condition */

/* Plot type, see list in global_pdes.c  */

#define PLOT 1

#define PLOT_B 0        /* plot type for second movie */

/* Color schemes */

#define COLOR_PALETTE 12     /* Color palette, see list in global_pdes.c  */

#define BLACK 1          /* background */

#define COLOR_SCHEME 3   /* choice of color scheme, see list in global_pdes.c  */

#define SCALE 0          /* set to 1 to adjust color scheme to variance of field */
#define SLOPE 0.5        /* sensitivity of color on wave amplitude */
#define ATTENUATION 0.0  /* exponential attenuation coefficient of contrast with time */
#define E_SCALE 150.0     /* scaling factor for energy representation */
#define LOG_SCALE 2.5     /* scaling factor for energy log representation */

#define COLORHUE 260     /* initial hue of water color for scheme C_LUM */
#define COLORDRIFT 0.0   /* how much the color hue drifts during the whole simulation */
#define LUMMEAN 0.5      /* amplitude of luminosity variation for scheme C_LUM */
#define LUMAMP 0.3       /* amplitude of luminosity variation for scheme C_LUM */
#define HUEMEAN 180.0    /* mean value of hue for color scheme C_HUE */
#define HUEAMP -180.0      /* amplitude of variation of hue for color scheme C_HUE */

#define DRAW_COLOR_SCHEME 1     /* set to 1 to plot the color scheme */
#define COLORBAR_RANGE 6.0    /* scale of color scheme bar */
#define COLORBAR_RANGE_B 10.0    /* scale of color scheme bar for 2nd part */
#define ROTATE_COLOR_SCHEME 0   /* set to 1 to draw color scheme horizontally */

#define SAVE_TIME_SERIES 0      /* set to 1 to save wave time series at a point */

/* For debugging purposes only */
#define FLOOR 0         /* set to 1 to limit wave amplitude to VMAX */
#define VMAX 10.0       /* max value of wave amplitude */

#include "global_pdes.c"        /* constants and global variables */
#include "sub_wave.c"           /* common functions for wave_billiard, heat and schrodinger */
#include "wave_common.c"        /* common functions for wave_billiard, wave_comparison, etc */

FILE *time_series_left, *time_series_right;

double courant2, courantb2;  /* Courant parameters squared */

/*********************/
/* animation part    */
/*********************/

void evolve_wave_half_old(double *phi_in[NX], double *psi_in[NX], double *phi_out[NX], double *psi_out[NX], 
                      short int *xy_in[NX])
/* time step of field evolution */
/* phi is value of field at time t, psi at time t-1 */
{
    int i, j, iplus, iminus, jplus, jminus;
    double delta, x, y, c, cc, gamma;
    static long time = 0;
    
    time++;
    
//     c = COURANT;
//     cc = courant2;

    #pragma omp parallel for private(i,j,iplus,iminus,jplus,jminus,delta,x,y,c,cc,gamma)
    for (i=0; i<NX; i++){
        for (j=0; j<NY; j++){
//             if (xy_in[i][j])
//             {
//                 c = COURANT;
//                 cc = courant2;
//                 gamma = GAMMA;
//             }
            if (xy_in[i][j] != 0)
            {
                c = COURANT;
                cc = courant2;
                if (xy_in[i][j] == 1) gamma = GAMMA;
                else gamma = GAMMAB;
            }
            else if (TWOSPEEDS)
            {
                c = COURANTB;
                cc = courantb2;
                gamma = GAMMAB;
            }

            if ((TWOSPEEDS)||(xy_in[i][j] != 0)){
                /* discretized Laplacian for various boundary conditions */
                if ((B_COND == BC_DIRICHLET)||(B_COND == BC_ABSORBING))
                {
                    iplus = (i+1);   if (iplus == NX) iplus = NX-1;
                    iminus = (i-1);  if (iminus == -1) iminus = 0;
                    jplus = (j+1);   if (jplus == NY) jplus = NY-1;
                    jminus = (j-1);  if (jminus == -1) jminus = 0;
                }
                else if (B_COND == BC_PERIODIC)
                {
                    iplus = (i+1) % NX;
                    iminus = (i-1) % NX;
                    if (iminus < 0) iminus += NX;
                    jplus = (j+1) % NY;
                    jminus = (j-1) % NY;
                    if (jminus < 0) jminus += NY;
                }
                else if (B_COND == BC_VPER_HABS)
                {
                    iplus = (i+1);   if (iplus == NX) iplus = NX-1;
                    iminus = (i-1);  if (iminus == -1) iminus = 0;
                    jplus = (j+1) % NY;
                    jminus = (j-1) % NY;
                    if (jminus < 0) jminus += NY;
                }
                
                /* imposing linear wave on top and bottom by making Laplacian 1d */
                if (OSCILLATE_TOPBOT)
                {
                    if (j == NY-1) jminus = NY-1;
                    else if (j == 0) jplus = 0;
                }
                
                delta = phi_in[iplus][j] + phi_in[iminus][j] + phi_in[i][jplus] + phi_in[i][jminus] - 4.0*phi_in[i][j];

                x = phi_in[i][j];
		y = psi_in[i][j];

                /* evolve phi */
                if ((B_COND == BC_PERIODIC)||(B_COND == BC_DIRICHLET)) 
                    phi_out[i][j] = -y + 2*x + cc*delta - KAPPA*x - gamma*(x-y);
                else if (B_COND == BC_ABSORBING)
                {
                    if ((i>0)&&(i<NX-1)&&(j>0)&&(j<NY-1))
                        phi_out[i][j] = -y + 2*x + cc*delta - KAPPA*x - gamma*(x-y);
                
                    /* upper border */
                    else if (j==NY-1) 
                        phi_out[i][j] = x - c*(x - phi_in[i][NY-2]) - KAPPA_TOPBOT*x - GAMMA_TOPBOT*(x-y);
                    
                    /* lower border */
                    else if (j==0) 
                        phi_out[i][j] = x - c*(x - phi_in[i][1]) - KAPPA_TOPBOT*x - GAMMA_TOPBOT*(x-y);
                
                    /* right border */
                    if (i==NX-1) 
                        phi_out[i][j] = x - c*(x - phi_in[NX-2][j]) - KAPPA_SIDES*x - GAMMA_SIDES*(x-y);
                    
                    /* left border */
                    else if (i==0) 
                        phi_out[i][j] = x - c*(x - phi_in[1][j]) - KAPPA_SIDES*x - GAMMA_SIDES*(x-y);
                }
                else if (B_COND == BC_VPER_HABS)
                {
                    if ((i>0)&&(i<NX-1))
                        phi_out[i][j] = -y + 2*x + cc*delta - KAPPA*x - gamma*(x-y);
                
                    /* right border */
                    else if (i==NX-1) 
                        phi_out[i][j] = x - c*(x - phi_in[NX-2][j]) - KAPPA_SIDES*x - GAMMA_SIDES*(x-y);
                    
                    /* left border */
                    else if (i==0) 
                        phi_out[i][j] = x - c*(x - phi_in[1][j]) - KAPPA_SIDES*x - GAMMA_SIDES*(x-y);
                }
                psi_out[i][j] = x;
                
                /* add oscillating boundary condition on the left */
                if ((i == 0)&&(OSCILLATE_LEFT)) phi_out[i][j] = AMPLITUDE*cos((double)time*OMEGA);
//                 psi_out[i][j] = x;

                if (FLOOR)
                {
                    if (phi_out[i][j] > VMAX) phi_out[i][j] = VMAX;
                    if (phi_out[i][j] < -VMAX) phi_out[i][j] = -VMAX;
                    if (psi_out[i][j] > VMAX) psi_out[i][j] = VMAX;
                    if (psi_out[i][j] < -VMAX) psi_out[i][j] = -VMAX;
                }
            }
        }
    }
//     printf("phi(0,0) = %.3lg, psi(0,0) = %.3lg\n", phi[NX/2][NY/2], psi[NX/2][NY/2]);
}


void evolve_wave_half(double *phi_in[NX], double *psi_in[NX], double *phi_out[NX], double *psi_out[NX], 
                      short int *xy_in[NX])
/* time step of field evolution */
/* phi is value of field at time t, psi at time t-1 */
/* this version of the function has been rewritten in order to minimize the number of if-branches */
{
    int i, j, iplus, iminus, jplus, jminus;
    double delta, x, y, c, cc, gamma;
    static long time = 0;
    static double tc[NX][NY], tcc[NX][NY], tgamma[NX][NY];
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
                psi_out[i][j] = x;
            }
        }
    }
    
    /* left boundary */
    if (OSCILLATE_LEFT) for (j=1; j<NY-1; j++) phi_out[0][j] = AMPLITUDE*cos((double)time*OMEGA);
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
            psi_out[0][j] = x;
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
            psi_out[NX-1][j] = x;
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
                    if (i==0) phi_out[0][NY-1] = x - tc[0][NY-1]*(x - phi_in[1][NY-1]) - KAPPA_SIDES*x - GAMMA_SIDES*(x-y);
                    else phi_out[i][NY-1] = -y + 2*x + tcc[i][NY-1]*delta - KAPPA*x - tgamma[i][NY-1]*(x-y);
                    break;
                }
            }
            psi_out[i][NY-1] = x;
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
                    if (i==0) phi_out[0][0] = x - tc[0][0]*(x - phi_in[1][0]) - KAPPA_SIDES*x - GAMMA_SIDES*(x-y);
                    else phi_out[i][0] = -y + 2*x + tcc[i][0]*delta - KAPPA*x - tgamma[i][0]*(x-y);
                    break;
                }
            }
            psi_out[i][0] = x;
        }
    }
    
    /* add oscillating boundary condition on the left corners */
    if (OSCILLATE_LEFT)
    {
        phi_out[0][0] = AMPLITUDE*cos((double)time*OMEGA);
        phi_out[0][NY-1] = AMPLITUDE*cos((double)time*OMEGA);
    }
    
    /* for debugging purposes/if there is a risk of blow-up */
    if (FLOOR) for (i=0; i<NX; i++){
        for (j=0; j<NY; j++){
            if (xy_in[i][j] != 0) 
            {
                if (phi_out[i][j] > VMAX) phi_out[i][j] = VMAX;
                if (phi_out[i][j] < -VMAX) phi_out[i][j] = -VMAX;
                if (psi_out[i][j] > VMAX) psi_out[i][j] = VMAX;
                if (psi_out[i][j] < -VMAX) psi_out[i][j] = -VMAX;
            }
        }
    }
}


void evolve_wave(double *phi[NX], double *psi[NX], double *phi_tmp[NX], double *psi_tmp[NX], short int *xy_in[NX])
/* time step of field evolution */
/* phi is value of field at time t, psi at time t-1 */
{
    evolve_wave_half(phi, psi, phi_tmp, psi_tmp, xy_in);
    evolve_wave_half(phi_tmp, psi_tmp, phi, psi, xy_in);
//     evolve_wave_half_old(phi, psi, phi_tmp, psi_tmp, xy_in);
//     evolve_wave_half_old(phi_tmp, psi_tmp, phi, psi, xy_in);
}


void draw_color_bar(int plot, double range)
{
    if (ROTATE_COLOR_SCHEME) draw_color_scheme(-1.0, -0.8, XMAX - 0.1, -1.0, plot, -range, range);
    else draw_color_scheme(1.7, YMIN + 0.1, 1.9, YMAX - 0.1, plot, -range, range);
}


void animation()
{
    double time, scale, ratio, startleft[2], startright[2], sign; 
    double *phi[NX], *psi[NX], *phi_tmp[NX], *psi_tmp[NX], *total_energy[NX];
    short int *xy_in[NX];
    int i, j, s, sample_left[2], sample_right[2];
    static int counter = 0;
    long int wave_value;
    
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
        phi_tmp[i] = (double *)malloc(NY*sizeof(double));
        psi_tmp[i] = (double *)malloc(NY*sizeof(double));
        total_energy[i] = (double *)malloc(NY*sizeof(double));
        xy_in[i] = (short int *)malloc(NY*sizeof(short int));
    }
    
    /* initialise positions and radii of circles */
    if ((B_DOMAIN == D_CIRCLES)||(B_DOMAIN == D_CIRCLES_IN_RECT)) init_circle_config(circles);
    else if (B_DOMAIN == D_POLYGONS) init_polygon_config(polygons);
    printf("Polygons initialized\n");
    
    /* initialise polyline for von Koch and similar domains */
    npolyline = init_polyline(MDEPTH, polyline);
    for (i=0; i<npolyline; i++) printf("vertex %i: (%.3f, %.3f)\n", i, polyline[i].x, polyline[i].y);

    courant2 = COURANT*COURANT;
    courantb2 = COURANTB*COURANTB;

    /* initialize wave with a drop at one point, zero elsewhere */
//     init_circular_wave(0.0, -LAMBDA, phi, psi, xy_in);
    
    /* initialize total energy table */
    if ((PLOT == P_MEAN_ENERGY)||(PLOT_B == P_MEAN_ENERGY)||(PLOT == P_LOG_MEAN_ENERGY)||(PLOT_B == P_LOG_MEAN_ENERGY))
        for (i=0; i<NX; i++)
            for (j=0; j<NY; j++) 
                total_energy[i][j] = 0.0;
    
    ratio = (XMAX - XMIN)/8.4;  /* for Tokarsky billiard */
    
//     isospectral_initial_point(0.2, 0.0, startleft, startright);    /* for isospectral billiards */
//     homophonic_initial_point(0.5, -0.25, 1.5, -0.25, startleft, startright);
//     homophonic_initial_point(0.5, -0.25, 1.5, -0.25, startleft, startright);
//     printf("xleft = (%.3f, %.3f) xright = (%.3f, %.3f)\n", startleft[0], startleft[1], startright[0], startright[1]);    
    
//     xy_to_ij(startleft[0], startleft[1], sample_left);
//     xy_to_ij(startright[0], startright[1], sample_right);
//     printf("xleft = (%.3f, %.3f) xright = (%.3f, %.3f)\n", xin_left, yin_left, xin_right, yin_right);
    
//     init_wave_flat(phi, psi, xy_in);
    
//     init_wave_plus(LAMBDA - 0.3*MU, 0.5*MU, phi, psi, xy_in);
//     init_wave(LAMBDA - 0.3*MU, 0.5*MU, phi, psi, xy_in);
//     init_circular_wave(X_SHOOTER, Y_SHOOTER, phi, psi, xy_in);
//     printf("Initializing wave\n");
//     init_circular_wave(-1.0, 0.0, phi, psi, xy_in);
//     printf("Wave initialized\n");
    init_circular_wave(0.0, 0.0, phi, psi, xy_in);
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
//     init_planar_wave(XMIN + 0.5, 0.0, phi, psi, xy_in);
//     init_wave(-1.5, 0.0, phi, psi, xy_in);
//     init_wave(0.0, 0.0, phi, psi, xy_in);

    /* add a drop at another point */
//     add_drop_to_wave(1.0, 0.7, 0.0, phi, psi);
//     add_drop_to_wave(1.0, -0.7, 0.0, phi, psi);
//     add_drop_to_wave(1.0, 0.0, -0.7, phi, psi);

    blank();
    glColor3f(0.0, 0.0, 0.0);
//     draw_wave(phi, psi, xy_in, 1.0, 0, PLOT);
    draw_wave_e(phi, psi, total_energy, xy_in, 1.0, 0, PLOT);
    draw_billiard();
    
    
    if (DRAW_COLOR_SCHEME) draw_color_bar(PLOT, COLORBAR_RANGE);

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
        draw_wave_e(phi, psi, total_energy, xy_in, scale, i, PLOT);
        for (j=0; j<NVID; j++) 
        {
            evolve_wave(phi, psi, phi_tmp, psi_tmp, xy_in);
            if (SAVE_TIME_SERIES)
            {
                wave_value = (long int)(phi[sample_left[0]][sample_left[1]]*1.0e16);
                fprintf(time_series_left, "%019ld\n", wave_value);
                wave_value = (long int)(phi[sample_right[0]][sample_right[1]]*1.0e16);
                fprintf(time_series_right, "%019ld\n", wave_value);
                if ((j == 0)&&(i%10 == 0)) printf("Frame %i of %i\n", i, NSTEPS);
//                 fprintf(time_series_right, "%.15f\n", phi[sample_right[0]][sample_right[1]]);
            }
//             if (i % 10 == 9) oscillate_linear_wave(0.2*scale, 0.15*(double)(i*NVID + j), -1.5, YMIN, -1.5, YMAX, phi, psi);
        }
        
        draw_billiard();
        
        if (DRAW_COLOR_SCHEME) draw_color_bar(PLOT, COLORBAR_RANGE); 
        
        /* add oscillating waves */
//         if (i%345 == 344)
//         {
//                 add_circular_wave(1.0, -LAMBDA, 0.0, phi, psi, xy_in);
//         }

	glutSwapBuffers();

	if (MOVIE)
        {
            if (i >= INITIAL_TIME) save_frame();
            else printf("Initial phase time %i of %i\n", i, INITIAL_TIME);
            
            if ((i >= INITIAL_TIME)&&(DOUBLE_MOVIE))
            {
//                 draw_wave(phi, psi, xy_in, scale, i, PLOT_B);
                draw_wave_e(phi, psi, total_energy, xy_in, scale, i, PLOT_B);
                if (DRAW_COLOR_SCHEME) draw_color_bar(PLOT_B, COLORBAR_RANGE_B);  
                draw_billiard();
                glutSwapBuffers();
                save_frame_counter(NSTEPS + 21 + counter);
                counter++;
            }

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
            draw_wave_e(phi, psi, total_energy, xy_in, scale, NSTEPS, PLOT);
            draw_billiard();
            if (DRAW_COLOR_SCHEME) draw_color_bar(PLOT, COLORBAR_RANGE);   
            glutSwapBuffers();
        }
        for (i=0; i<MID_FRAMES; i++) save_frame();
        if (DOUBLE_MOVIE) 
        {
//             draw_wave(phi, psi, xy_in, scale, i, PLOT_B);
            draw_wave_e(phi, psi, total_energy, xy_in, scale, NSTEPS, PLOT_B);
            draw_billiard();
            if (DRAW_COLOR_SCHEME) draw_color_bar(PLOT_B, COLORBAR_RANGE_B); 
            glutSwapBuffers();
        }
        for (i=0; i<END_FRAMES; i++) save_frame_counter(NSTEPS + MID_FRAMES + 1 + counter + i);
        
        s = system("mv wave*.tif tif_wave/");
    }
    for (i=0; i<NX; i++)
    {
        free(phi[i]);
        free(psi[i]);
        free(phi_tmp[i]);
        free(psi_tmp[i]);
        free(total_energy[i]);
        free(xy_in[i]);
    }
    
    if (SAVE_TIME_SERIES)
    {
        fclose(time_series_left);
        fclose(time_series_right);
    }


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

