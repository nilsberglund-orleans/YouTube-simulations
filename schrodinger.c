/*********************************************************************************/
/*                                                                               */
/*  Animation of Schrödinger equation in a planar domain                         */
/*                                                                               */
/*  N. Berglund, May 2021                                                        */
/*                                                                               */
/*  Feel free to reuse, but if doing so it would be nice to drop a               */
/*  line to nils.berglund@univ-orleans.fr - Thanks!                              */
/*                                                                               */
/*  compile with                                                                 */
/*  gcc -o schrodinger schrodinger.c                                             */
/* -L/usr/X11R6/lib -ltiff -lm -lGL -lGLU -lX11 -lXmu -lglut -O3 -fopenmp        */
/*                                                                               */
/*  To make a video, set MOVIE to 1 and create subfolder tif_schrod              */
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

// #define NX 1280          /* number of grid points on x axis */
// #define NX 720          /* number of grid points on x axis */
#define NX 640          /* number of grid points on x axis */
#define NY 360          /* number of grid points on y axis */

/* setting NX to WINWIDTH and NY to WINHEIGHT increases resolution */
/* but will multiply run time by 4                                 */

#define XMIN -2.0
#define XMAX 2.0	/* x interval */
#define YMIN -1.125
#define YMAX 1.125	/* y interval for 9/16 aspect ratio */

#define JULIA_SCALE 1.0 /* scaling for Julia sets */

/* Choice of the billiard table, see list in global_pdes.c  */

#define B_DOMAIN 10      /* choice of domain shape */

#define CIRCLE_PATTERN 0    /* pattern of circles, see list in global_pdes.c */

#define P_PERCOL 0.25       /* probability of having a circle in C_RAND_PERCOL arrangement */
#define NPOISSON 300        /* number of points for Poisson C_RAND_POISSON arrangement */
#define RANDOM_POLY_ANGLE 1 /* set to 1 to randomize angle of polygons */

#define LAMBDA 0.1	    /* parameter controlling the dimensions of domain */
#define MU 0.03	            /* parameter controlling the dimensions of domain */
#define NPOLY 6             /* number of sides of polygon */
#define APOLY 1.0           /* angle by which to turn polygon, in units of Pi/2 */ 
#define MDEPTH 5            /* depth of computation of Menger gasket */
#define MRATIO 3            /* ratio defining Menger gasket */
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

#define DT 0.00000001
// #define DT 0.00000001
// #define DT 0.000000005
// #define DT 0.000000005
#define HBAR 1.0

/* Boundary conditions, see list in global_pdes.c  */

#define B_COND 1

/* Parameters for length and speed of simulation */

#define NSTEPS 2500      /* number of frames of movie */
// #define NVID 2000         /* number of iterations between images displayed on screen */
#define NVID 1200         /* number of iterations between images displayed on screen */
#define NSEG 100         /* number of segments of boundary */
#define BOUNDARY_WIDTH 2    /* width of billiard boundary */

#define PAUSE 1000       /* number of frames after which to pause */
#define PSLEEP 1         /* sleep time during pause */
#define SLEEP1  1        /* initial sleeping time */
#define SLEEP2  1        /* final sleeping time */
#define END_FRAMES 100   /* still frames at end of movie */

/* For debugging purposes only */
#define FLOOR 0         /* set to 1 to limit wave amplitude to VMAX */
#define VMAX 10.0       /* max value of wave amplitude */


/* Plot type, see list in global_pdes.c  */

#define PLOT 11


/* Color schemes, see list in global_pdes.c  */

#define COLOR_PALETTE 10     /* Color palette, see list in global_pdes.c  */

#define BLACK 1          /* black background */

#define COLOR_SCHEME 3   /* choice of color scheme */

#define SCALE 1          /* set to 1 to adjust color scheme to variance of field */
#define SLOPE 1.0        /* sensitivity of color on wave amplitude */
#define ATTENUATION 0.0  /* exponential attenuation coefficient of contrast with time */
#define E_SCALE 150.0     /* scaling factor for energy representation */
#define LOG_SCALE 1.0    /* scaling factor for energy log representation */
#define LOG_SHIFT 0.0     /* shift of colors on log scale */
#define RESCALE_COLOR_IN_CENTER 0   /* set to 1 to decrease color intentiy in the center (for wave escaping ring) */

#define COLORHUE 260     /* initial hue of water color for scheme C_LUM */
#define COLORDRIFT 0.0   /* how much the color hue drifts during the whole simulation */
#define LUMMEAN 0.5      /* amplitude of luminosity variation for scheme C_LUM */
#define LUMAMP 0.3       /* amplitude of luminosity variation for scheme C_LUM */
#define HUEMEAN 180.0    /* mean value of hue for color scheme C_HUE */
#define HUEAMP 180.0      /* amplitude of variation of hue for color scheme C_HUE */

#define DRAW_COLOR_SCHEME 1     /* set to 1 to plot the color scheme */
#define COLORBAR_RANGE 2.0    /* scale of color scheme bar */
#define COLORBAR_RANGE_B 12.0    /* scale of color scheme bar for 2nd part */
#define ROTATE_COLOR_SCHEME 0   /* set to 1 to draw color scheme horizontally */

/* the following constants are only used by wave_billiard and wave_3d so far */
#define COMPARISON 0        /* set to 1 to compare two different patterns */
#define B_DOMAIN_B 20       /* second domain shape, for comparisons */
#define CIRCLE_PATTERN_B 0  /* second pattern of circles or polygons */
#define OSCILLATION_SCHEDULE 3  /* oscillation schedule, see list in global_pdes.c */
#define ACHIRP 0.2        /* acceleration coefficient in chirp */
#define DAMPING 0.0        /* damping of periodic excitation */
#define OMEGA 0.001       /* frequency of periodic excitation */
#define AMPLITUDE 0.8      /* amplitude of periodic excitation */ 
/* end of constants only used by wave_billiard and wave_3d */

/* for compatibility with sub_wave and sub_maze */
#define NXMAZE 7      /* width of maze */
#define NYMAZE 7      /* height of maze */
#define MAZE_MAX_NGBH 4     /* max number of neighbours of maze cell */
#define RAND_SHIFT 24        /* seed of random number generator */
#define MAZE_XSHIFT 0.0     /* horizontal shift of maze */
#define MAZE_WIDTH 0.02     /* half width of maze walls */
#define ADD_POTENTIAL 0
#define POT_MAZE 7
#define POTENTIAL 0
#define VARIABLE_IOR 1      /* set to 1 for a variable index of refraction */
#define IOR 7               /* choice of index of refraction, see list in global_pdes.c */
#define IOR_TOTAL_TURNS 1.5 /* total angle of rotation for IOR_PERIODIC_WELLS_ROTATING */
#define MANDEL_IOR_SCALE -0.05   /* parameter controlling dependence of IoR on Mandelbrot escape speed */
#define COURANT 0.04       /* Courant number */
#define COURANTB 0.0       /* Courant number in medium B */
#define INITIAL_AMP 0.5            /* amplitude of initial condition */
#define INITIAL_VARIANCE 0.0003    /* variance of initial condition */
#define INITIAL_WAVELENGTH  0.015  /* wavelength of initial condition */
#define TWOSPEEDS 0          /* set to 1 to replace hardcore boundary by medium with different speed */
#define WAVE_PACKET_SOURCE_TYPE 1       /* type of wave packet sources */
#define N_WAVE_PACKETS 15               /* number of wave packets */
#define OSCIL_LEFT_YSHIFT 0.0   /* y-dependence of left oscillation (for non-horizontal waves) */
#define OSCILLATING_SOURCE_PERIOD 20    /* period of oscillating source */
#define MU_B 1.0           /* parameter controlling the dimensions of domain */
#define GAMMA 0.0          /* damping factor in wave equation */
#define GAMMAB 0.0         /* damping factor in wave equation */
#define DRAW_WAVE_PROFILE 0     /* set to 1 to draw a profile of the wave */
#define VERTICAL_WAVE_PROFILE 0 /* set to 1 to draw wave profile vertically */
#define WALL_WIDTH 0.1      /* width of wall separating lenses */
#define INITIAL_TIME 50      /* time after which to start saving frames */
#define OSCIL_YMAX 0.35      /* defines oscillation range */
#define MESSAGE_LDASH 14         /* length of dash for Morse code message */
#define MESSAGE_LDOT 8           /* length of dot for Morse code message */
#define MESSAGE_LINTERVAL 54     /* length of interval between dashes/dots for Morse code message */
#define MESSAGE_LINTERLETTER 60  /* length of interval between letters for Morse code message */
#define MESSAGE_LSPACE 48        /* length of space for Morse code message */
#define MESSAGE_INITIAL_TIME 100 /* initial time before starting message for Morse code message */    
#define WAVE_PROFILE_X 2.1      /* value of x to sample wave profile */
#define HRES 1          /* dummy, only used by rde.c */
#define INITIAL_SHIFT 20.0          /* time shift of initial wave packet (in oscillation periods) */
#define WAVE_PACKET_SHIFT 200.0     /* time shift between wave packets (in oscillation periods) */
#define FADE_IN_OBSTACLE 0      /* set to 1 to fade color inside obstacles */
#define SHADE_2D 0       /* set to 1 to add pseudo-3d shading effect */ 
#define SHADE_SCALE_2D 0.05  /* lower value increases sensitivity of shading */
#define N_SOURCES 1                     /* number of sources, for option draw_sources */
double light[2] = {0.40824829, 0.816496581};   /* location of light source for SHADE_2D option*/
/* end of constants only used by sub_wave and sub_maze */


#include "global_pdes.c"
#include "sub_maze.c"
#include "sub_wave.c"

double courant2;  /* Courant parameter squared */
double dx2;       /* spatial step size squared */
double intstep;   /* integration step */
double intstep1;  /* integration step used in absorbing boundary conditions */



void init_coherent_state(double x, double y, double px, double py, double scalex, double *phi[NX], 
                         double *psi[NX], short int *xy_in[NX])
/* initialise field with coherent state of position (x,y) and momentum (px, py) */
/* phi is real part, psi is imaginary part */
{
    int i, j;
    double xy[2], dist2, module, phase, scale2;    

    scale2 = scalex*scalex;
    for (i=0; i<NX; i++)
        for (j=0; j<NY; j++)
        {
            ij_to_xy(i, j, xy);
	    xy_in[i][j] = xy_in_billiard(xy[0],xy[1]);

            if (xy_in[i][j])
            {
                dist2 = (xy[0]-x)*(xy[0]-x) + (xy[1]-y)*(xy[1]-y);
                module = exp(-dist2/scale2);
                if (module < 1.0e-15) module = 1.0e-15;
                phase = (px*(xy[0]-x) + py*(xy[1]-y))/scalex;

                phi[i][j] = module*cos(phase);
                psi[i][j] = module*sin(phase);
            }
            else
            {
                phi[i][j] = 0.0;
                psi[i][j] = 0.0;
            }
        }
}



/*********************/
/* animation part    */
/*********************/

void schrodinger_color_scheme(double phi, double psi, double scale, int time, double rgb[3])
// double phi, psi, scale, rgb[3];
// int time;
{
    double phase, amp, lum;
    
    if (PLOT == P_MODULE)
        color_scheme(COLOR_SCHEME, 2.0*module2(phi, psi)-1.0, scale, time, rgb);
    else if (PLOT == P_PHASE)
    {
        amp = module2(phi,psi);
//         if (amp < 1.0e-10) amp = 1.0e-10;
        phase = argument(phi/amp, psi/amp);
        if (phase < 0.0) phase += DPI;
        lum = (color_amplitude(amp, scale, time))*0.5;
        if (lum < 0.0) lum = 0.0;
        hsl_to_rgb(phase*360.0/DPI, 0.9, lum, rgb);
    }
    else if (PLOT == P_REAL) color_scheme(COLOR_SCHEME, phi, scale, time, rgb);
    else if (PLOT == P_IMAGINARY) color_scheme(COLOR_SCHEME, psi, scale, time, rgb);
}


void draw_wave(double *phi[NX], double *psi[NX], short int *xy_in[NX], double scale, int time)
/* draw the field */
{
    int i, j;
    double rgb[3], xy[2], x1, y1, x2, y2, amp, phase;

    glBegin(GL_QUADS);

    for (i=0; i<NX; i++)
        for (j=0; j<NY; j++)
        {
            if (xy_in[i][j])
            {
                schrodinger_color_scheme(phi[i][j],psi[i][j], scale, time, rgb);
                    
                glColor3f(rgb[0], rgb[1], rgb[2]);

                glVertex2i(i, j);
                glVertex2i(i+1, j);
                glVertex2i(i+1, j+1);
                glVertex2i(i, j+1);
            }
        }

    glEnd ();
}

void evolve_wave_half_old(double *phi_in[NX], double *psi_in[NX], double *phi_out[NX], double *psi_out[NX], 
                      short int *xy_in[NX])
// void evolve_wave_half(phi_in, psi_in, phi_out, psi_out, xy_in)
// /* time step of field evolution */
// /* phi is real part, psi is imaginary part */
//     double *phi_in[NX], *psi_in[NX], *phi_out[NX], *psi_out[NX]; short int *xy_in[NX];
{
    int i, j, iplus, iminus, jplus, jminus;
    double delta1, delta2, x, y;
    
    #pragma omp parallel for private(i,j,iplus,iminus,jplus,jminus,delta1,delta2,x,y)
    for (i=0; i<NX; i++){
        for (j=0; j<NY; j++){
            if (xy_in[i][j]){
                /* discretized Laplacian depending on boundary conditions */
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
                
                delta1 = phi_in[iplus][j] + phi_in[iminus][j] + phi_in[i][jplus] + phi_in[i][jminus] - 4.0*phi_in[i][j];
                delta2 = psi_in[iplus][j] + psi_in[iminus][j] + psi_in[i][jplus] + psi_in[i][jminus] - 4.0*psi_in[i][j];

                x = phi_in[i][j];
		y = psi_in[i][j];

                /* evolve phi and psi */
                if (B_COND != BC_ABSORBING)
                {
                    phi_out[i][j] = x - intstep*delta2;
                    psi_out[i][j] = y + intstep*delta1;
                }
                else        /* case of absorbing b.c. - this is only an approximation of correct way of implementing */
                {
                    /* in the bulk */
                    if ((i>0)&&(i<NX-1)&&(j>0)&&(j<NY-1))
                    {
                        phi_out[i][j] = x - intstep*delta2;
                        psi_out[i][j] = y + intstep*delta1;
                    }
                     /* right border */
                    else if (i==NX-1) 
                    {
                        phi_out[i][j] = x - intstep1*(y - psi_in[i-1][j]);
                        psi_out[i][j] = y + intstep1*(x - phi_in[i-1][j]);
                    }
                    /* upper border */
                    else if (j==NY-1) 
                    {
                        phi_out[i][j] = x - intstep1*(y - psi_in[i][j-1]);
                        psi_out[i][j] = y + intstep1*(x - phi_in[i][j-1]);
                    }
                    /* left border */
                    else if (i==0) 
                    {
                        phi_out[i][j] = x - intstep1*(y - psi_in[1][j]);
                        psi_out[i][j] = y + intstep1*(x - phi_in[1][j]);
                    }
                   /* lower border */
                    else if (j==0) 
                    {
                        phi_out[i][j] = x - intstep1*(y - psi_in[i][1]);
                        psi_out[i][j] = y + intstep1*(x - phi_in[i][1]);
                    }
                }


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
// void evolve_wave_half(phi_in, psi_in, phi_out, psi_out, xy_in)
// /* time step of field evolution */
// /* phi is real part, psi is imaginary part */
{
    int i, j, iplus, iminus, jplus, jminus;
    double delta1, delta2, x, y;
    
    #pragma omp parallel for private(i,j,iplus,iminus,jplus,jminus,delta1,delta2,x,y)
    for (i=1; i<NX-1; i++){
        for (j=1; j<NY-1; j++){
            if (xy_in[i][j]){
                x = phi_in[i][j];
		y = psi_in[i][j];
                
                delta1 = phi_in[i+1][j] + phi_in[i-1][j] + phi_in[i][j+1] + phi_in[i][j-1] - 4.0*x;
                delta2 = psi_in[i+1][j] + psi_in[i-1][j] + psi_in[i][j+1] + psi_in[i][j-1] - 4.0*y;

                /* evolve phi and psi */
                phi_out[i][j] = x - intstep*delta2;
                psi_out[i][j] = y + intstep*delta1;
            }
        }
    }
    
    /* left boundary */
    for (j=1; j<NY-1; j++){
        if (xy_in[0][j]){
            x = phi_in[0][j];
            y = psi_in[0][j];
                    
            switch (B_COND) {
                case (BC_DIRICHLET):
                {
                    delta1 = phi_in[1][j] + phi_in[0][j+1] + phi_in[0][j-1] - 3.0*x;
                    delta2 = psi_in[1][j] + psi_in[0][j+1] + psi_in[0][j-1] - 3.0*y;
                    phi_out[0][j] = x - intstep*delta2;
                    psi_out[0][j] = y + intstep*delta1;
                    break;
                }
                case (BC_PERIODIC):
                {
                    delta1 = phi_in[1][j] + phi_in[NX-1][j] + phi_in[0][j+1] + phi_in[0][j-1] - 4.0*x;
                    delta2 = psi_in[1][j] + psi_in[NX-1][j] + psi_in[0][j+1] + psi_in[0][j-1] - 4.0*y;
                    phi_out[0][j] = x - intstep*delta2;
                    psi_out[0][j] = y + intstep*delta1;
                    break;
                }
            }
        }
    }
    
    /* right boundary */
    for (j=1; j<NY-1; j++){
        if (xy_in[0][j]){
            x = phi_in[NX-1][j];
            y = psi_in[NX-1][j];
                    
            switch (B_COND) {
                case (BC_DIRICHLET):
                {
                    delta1 = phi_in[NX-2][j] + phi_in[NX-1][j+1] + phi_in[NX-1][j-1] - 3.0*x;
                    delta2 = psi_in[NX-2][j] + psi_in[NX-1][j+1] + psi_in[NX-1][j-1] - 3.0*y;
                    phi_out[NX-1][j] = x - intstep*delta2;
                    psi_out[NX-1][j] = y + intstep*delta1;
                    break;
                }
                case (BC_PERIODIC):
                {
                    delta1 = phi_in[NX-2][j] + phi_in[0][j] + phi_in[NX-1][j+1] + phi_in[NX-1][j-1] - 4.0*x;
                    delta2 = psi_in[NX-2][j] + psi_in[0][j] + psi_in[NX-1][j+1] + psi_in[NX-1][j-1] - 4.0*y;
                    phi_out[NX-1][j] = x - intstep*delta2;
                    psi_out[NX-1][j] = y + intstep*delta1;
                    break;
                }
            }
        }
    }

    /* top boundary */
    for (i=0; i<NX; i++){
        if (xy_in[i][NY-1]){
            x = phi_in[i][NY-1];
            y = psi_in[i][NY-1];
                    
            switch (B_COND) {
                case (BC_DIRICHLET):
                {
                    iplus = i+1;   if (iplus == NX) iplus = NX-1;
                    iminus = i-1;  if (iminus == -1) iminus = 0;
                    
                    delta1 = phi_in[iplus][NY-1] + phi_in[iminus][NY-1] + phi_in[i][NY-2] - 3.0*x;
                    delta2 = psi_in[iplus][NY-1] + psi_in[iminus][NY-1] + psi_in[i][NY-2] - 3.0*x;
                    phi_out[i][NY-1] = x - intstep*delta2;
                    psi_out[i][NY-1] = y + intstep*delta1;
                    break;
                }
                case (BC_PERIODIC):
                {
                    iplus = (i+1) % NX;
                    iminus = (i-1) % NX;
                    if (iminus < 0) iminus += NX;
                    
                    delta1 = phi_in[iplus][NY-1] + phi_in[iminus][NY-1] + phi_in[i][NY-2] + phi_in[i][0] - 4.0*x;
                    delta2 = psi_in[iplus][NY-1] + psi_in[iminus][NY-1] + psi_in[i][NY-2] + psi_in[i][0] - 4.0*y;
                    phi_out[i][NY-1] = x - intstep*delta2;
                    psi_out[i][NY-1] = y + intstep*delta1;
                    break;
                }
            }
        }
    }
    
    /* bottom boundary */
    for (i=0; i<NX; i++){
        if (xy_in[i][0]){
            x = phi_in[i][0];
            y = psi_in[i][0];
                    
            switch (B_COND) {
                case (BC_DIRICHLET):
                {
                    iplus = i+1;   if (iplus == NX) iplus = NX-1;
                    iminus = i-1;  if (iminus == -1) iminus = 0;
                    
                    delta1 = phi_in[iplus][0] + phi_in[iminus][0] + phi_in[i][1] - 3.0*x;
                    delta2 = psi_in[iplus][0] + psi_in[iminus][0] + psi_in[i][1] - 3.0*x;
                    phi_out[i][0] = x - intstep*delta2;
                    psi_out[i][0] = y + intstep*delta1;
                    break;
                }
                case (BC_PERIODIC):
                {
                    iplus = (i+1) % NX;
                    iminus = (i-1) % NX;
                    if (iminus < 0) iminus += NX;
                    
                    delta1 = phi_in[iplus][0] + phi_in[iminus][0] + phi_in[i][1] + phi_in[i][NY-1] - 4.0*x;
                    delta2 = psi_in[iplus][0] + psi_in[iminus][0] + psi_in[i][1] + psi_in[i][NY-1] - 4.0*y;
                    phi_out[i][0] = x - intstep*delta2;
                    psi_out[i][0] = y + intstep*delta1;
                    break;
                }
            }
        }
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
/* phi is real part, psi is imaginary part */
{
    evolve_wave_half(phi, psi, phi_tmp, psi_tmp, xy_in);
    evolve_wave_half(phi_tmp, psi_tmp, phi, psi, xy_in);
}


double compute_variance(double *phi[NX], double *psi[NX], short int *xy_in[NX])
// double compute_variance(phi, psi, xy_in)
/* compute the variance (total probability) of the field */
// double *phi[NX], *psi[NX]; short int * xy_in[NX];
{
    int i, j, n = 0;
    double variance = 0.0;

    for (i=1; i<NX; i++)
        for (j=1; j<NY; j++)
        {
            if (xy_in[i][j])
            {
                n++;
                variance += phi[i][j]*phi[i][j] + psi[i][j]*psi[i][j];
            }
        }
    if (n==0) n=1;
    return(variance/(double)n);
}

void renormalise_field(double *phi[NX], double *psi[NX], short int *xy_in[NX], double variance)
/* renormalise variance of field */
{
    int i, j;
    double stdv;
    
    stdv = sqrt(variance);

    for (i=1; i<NX; i++)
        for (j=1; j<NY; j++)
        {
            if (xy_in[i][j])
            {
                phi[i][j] = phi[i][j]/stdv;
                psi[i][j] = psi[i][j]/stdv;
            }
        }
}


void draw_color_bar(int plot, double range)
{
    if (ROTATE_COLOR_SCHEME) draw_color_scheme(-1.0, -0.8, XMAX - 0.1, -1.0, plot, -range, range);
    else draw_color_scheme(1.7, YMIN + 0.1, 1.9, YMAX - 0.1, plot, -range, range);
}


void animation()
{
    double time, scale, dx, var;
    double *phi[NX], *psi[NX], *phi_tmp[NX], *psi_tmp[NX];
    short int *xy_in[NX];
    int i, j, s;

    /* Since NX and NY are big, it seemed wiser to use some memory allocation here */
    for (i=0; i<NX; i++)
    {
        phi[i] = (double *)malloc(NY*sizeof(double));
        psi[i] = (double *)malloc(NY*sizeof(double));
        phi_tmp[i] = (double *)malloc(NY*sizeof(double));
        psi_tmp[i] = (double *)malloc(NY*sizeof(double));
        xy_in[i] = (short int *)malloc(NY*sizeof(short int));
    }
    
    /* initialise polyline for von Koch and simular domains */
    npolyline = init_polyline(MDEPTH, polyline);
//     for (i=0; i<npolyline; i++) printf("vertex %i: (%.3f, %.3f)\n", i, polyline[i].x, polyline[i].y);

    dx = (XMAX-XMIN)/((double)NX);
    intstep = DT/(dx*dx*HBAR);
    intstep1 = DT/(dx*HBAR);
    
    printf("Integration step %.3lg\n", intstep);

    /* initialize wave wave function */
    init_coherent_state(-0.5, 0.0, 15.0, 0.0, 0.15, phi, psi, xy_in);
//     init_coherent_state(0.0, 0.0, 0.0, 5.0, 0.03, phi, psi, xy_in);
//     init_coherent_state(-0.5, 0.0, 1.0, 1.0, 0.05, phi, psi, xy_in);
    
    
    
    if (SCALE)
    {
        var = compute_variance(phi,psi, xy_in);
        scale = sqrt(1.0 + var);
        renormalise_field(phi, psi, xy_in, var);
    }
    
    blank();
    
    if (DRAW_COLOR_SCHEME) draw_color_bar(PLOT, COLORBAR_RANGE);

    glColor3f(0.0, 0.0, 0.0);

    glutSwapBuffers();

    sleep(SLEEP1);

    for (i=0; i<=NSTEPS; i++)
    {
        /* compute the variance of the field to adjust color scheme */
        /* the color depends on the field divided by sqrt(1 + variance) */
        if (SCALE)
        {
            var = compute_variance(phi,psi, xy_in);
            scale = sqrt(1.0 + var);
//             printf("Norm: %5lg\t Scaling factor: %5lg\n", var, scale);
            renormalise_field(phi, psi, xy_in, var);
        }
        else scale = 1.0;

        draw_wave(phi, psi, xy_in, scale, i);
        
//         printf("Wave drawn\n");
        
        for (j=0; j<NVID; j++) evolve_wave(phi, psi, phi_tmp, psi_tmp, xy_in);
        
        draw_billiard(0, 1.0);
        
        if (DRAW_COLOR_SCHEME) draw_color_bar(PLOT, COLORBAR_RANGE); 
        
        glutSwapBuffers();

	if (MOVIE)
        {
            save_frame();

            /* it seems that saving too many files too fast can cause trouble with the file system */
            /* so this is to make a pause from time to time - parameter PAUSE may need adjusting   */
            if (i % PAUSE == PAUSE - 1)
            {
                printf("Making a short pause\n");
                sleep(PSLEEP);
                s = system("mv wave*.tif tif_schrod/");
            }
        }
    }

    if (MOVIE)
    {
        for (i=0; i<END_FRAMES; i++) save_frame();
        s = system("mv wave*.tif tif_schrod/");
    }
    
    for (i=0; i<NX; i++)
    {
        free(phi[i]);
        free(psi[i]);
        free(phi_tmp[i]);
        free(psi_tmp[i]);
        free(xy_in[i]);
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
    glutCreateWindow("Schrodinger equation in a planar domain");

    init();

    glutDisplayFunc(display);

    glutMainLoop();

    return 0;
}

