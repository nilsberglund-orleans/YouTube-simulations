/*********************************************************************************/
/*                                                                               */
/*  Animation of heat equation in a planar domain                                */
/*                                                                               */
/*  N. Berglund, May 2021                                                        */
/*                                                                               */
/*  Feel free to reuse, but if doing so it would be nice to drop a               */
/*  line to nils.berglund@univ-orleans.fr - Thanks!                              */
/*                                                                               */
/*  compile with                                                                 */
/*  gcc -o heat heat.c                                                           */
/* -L/usr/X11R6/lib -ltiff -lm -lGL -lGLU -lX11 -lXmu -lglut -O3 -fopenmp        */
/*                                                                               */
/*  To make a video, set MOVIE to 1 and create subfolder tif_heat                */
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

// #define WINWIDTH 	1280  /* window width */
#define WINWIDTH 	720  /* window width */
#define WINHEIGHT 	720   /* window height */

// #define NX 1280          /* number of grid points on x axis */
#define NX 720          /* number of grid points on x axis */
#define NY 720          /* number of grid points on y axis */
// #define NX 640          /* number of grid points on x axis */
// #define NY 360          /* number of grid points on y axis */

/* setting NX to WINWIDTH and NY to WINHEIGHT increases resolution */
/* but will multiply run time by 4                                 */

// #define XMIN -2.0
// #define XMAX 2.0	/* x interval */
#define XMIN -1.125
#define XMAX 1.125	/* x interval */
#define YMIN -1.125
#define YMAX 1.125	/* y interval for 9/16 aspect ratio */

#define JULIA_SCALE 0.5 /* scaling for Julia sets */

/* Choice of the billiard table */

#define B_DOMAIN 27      /* choice of domain shape, see list in global_pdes.c  */

#define CIRCLE_PATTERN 0    /* pattern of circles, see list in global_pdes.c */

#define P_PERCOL 0.25       /* probability of having a circle in C_RAND_PERCOL arrangement */
#define NPOISSON 300        /* number of points for Poisson C_RAND_POISSON arrangement */
#define PDISC_FACTOR 3.25   /* controls density of Poisson disc process (default: 3.25) */
#define RANDOM_POLY_ANGLE 0 /* set to 1 to randomize angle of polygons */

#define LAMBDA 1.1	    /* parameter controlling the dimensions of domain */
#define MU 0.8	            /* parameter controlling the dimensions of domain */
#define NPOLY 6             /* number of sides of polygon */
#define APOLY 1.0           /* angle by which to turn polygon, in units of Pi/2 */
#define MDEPTH 1           /* depth of computation of Menger gasket */
#define MRATIO 5            /* ratio defining Menger gasket */
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

// #define DT 0.00001
#define DT 0.000004
// #define DT 0.000002
// #define DT 0.00000002
// #define DT 0.000000005
#define VISCOSITY 10.0
#define T_OUT 2.0       /* outside temperature */
#define T_IN 0.0        /* inside temperature */
// #define T_OUT 0.0       /* outside temperature */
// #define T_IN 2.0        /* inside temperature */
#define SPEED 0.0       /* speed of drift to the right */

/* Boundary conditions, see list in global_pdes.c  */

#define B_COND 1

/* Parameters for length and speed of simulation */

#define NSTEPS 1200      /* number of frames of movie */
#define NVID 30          /* number of iterations between images displayed on screen */
// #define NVID 100          /* number of iterations between images displayed on screen */
#define NSEG 100         /* number of segments of boundary */
#define BOUNDARY_WIDTH 1    /* width of billiard boundary */
#define DRAW_BILLIARD 0     /* set to 1 to draw billiard boundary */

#define PAUSE 100       /* number of frames after which to pause */
#define PSLEEP 1         /* sleep time during pause */
#define SLEEP1  2        /* initial sleeping time */
#define SLEEP2  1        /* final sleeping time */

/* For debugging purposes only */
#define FLOOR 0         /* set to 1 to limit wave amplitude to VMAX */
#define VMAX 10.0       /* max value of wave amplitude */

/* Field representation */

#define FIELD_REP 0

#define F_INTENSITY 0   /* color represents intensity */
#define F_GRADIENT 1    /* color represents norm of gradient */ 

#define DRAW_FIELD_LINES 0  /* set to 1 to draw field lines */
#define FIELD_LINE_WIDTH 1  /* width of field lines */
#define N_FIELD_LINES 120   /* number of field lines */
#define FIELD_LINE_FACTOR 120 /* factor controlling precision when computing origin of field lines */

/* Color schemes, see list in global_pdes.c  */

#define COLOR_PALETTE 10     /* Color palette, see list in global_pdes.c  */

#define BLACK 1          /* black background */

#define COLOR_SCHEME 1   /* choice of color scheme */

#define SCALE 0          /* set to 1 to adjust color scheme to variance of field */
// #define SLOPE 0.1        /* sensitivity of color on wave amplitude */
#define SLOPE 0.2        /* sensitivity of color on wave amplitude */
#define ATTENUATION 0.0  /* exponential attenuation coefficient of contrast with time */
#define PHASE_SHIFT 0.0      /* shift of phase in color scheme P_3D_PHASE */

#define COLORHUE 260     /* initial hue of water color for scheme C_LUM */
#define COLORDRIFT 0.0   /* how much the color hue drifts during the whole simulation */
#define LUMMEAN 0.5      /* amplitude of luminosity variation for scheme C_LUM */
#define LUMAMP 0.3       /* amplitude of luminosity variation for scheme C_LUM */
// #define HUEMEAN 180.0    /* mean value of hue for color scheme C_HUE */
// #define HUEAMP -180.0      /* amplitude of variation of hue for color scheme C_HUE */
#define HUEMEAN 359.0    /* mean value of hue for color scheme C_HUE */
#define HUEAMP -359.0      /* amplitude of variation of hue for color scheme C_HUE */
// #define HUEMEAN 270.0    /* mean value of hue for color scheme C_HUE */
// #define HUEAMP -130.0      /* amplitude of variation of hue for color scheme C_HUE */
#define E_SCALE 100.0     /* scaling factor for energy representation */
#define LOG_SCALE 1.0     /* scaling factor for energy log representation */
#define LOG_SHIFT 0.0   

#define DRAW_COLOR_SCHEME 1     /* set to 1 to plot the color scheme */
#define COLORBAR_RANGE 2.0    /* scale of color scheme bar */
#define COLORBAR_RANGE_B 12.0    /* scale of color scheme bar for 2nd part */
#define ROTATE_COLOR_SCHEME 1   /* set to 1 to draw color scheme horizontally */

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
#define RADIUS_FACTOR 0.3   /* controls inner radius for C_RING arrangements */
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



void init_gaussian(double x, double y, double mean, double amplitude, double scalex, 
                   double *phi[NX], short int * xy_in[NX])
/* initialise field with gaussian at position (x,y) */
{
    int i, j, in;
    double xy[2], dist2, module, phase, scale2;    

    scale2 = scalex*scalex;
    printf("Initialising field\n");
    for (i=0; i<NX; i++)
        for (j=0; j<NY; j++)
        {
            ij_to_xy(i, j, xy);
	    xy_in[i][j] = xy_in_billiard(xy[0],xy[1]);

            in = xy_in[i][j];
            if (in == 1)
            {
                dist2 = (xy[0]-x)*(xy[0]-x) + (xy[1]-y)*(xy[1]-y);
                module = amplitude*exp(-dist2/scale2);
                if (module < 1.0e-15) module = 1.0e-15;

                phi[i][j] = mean + module/scalex;
            }   /* boundary temperatures */
            else if (in >= 2) phi[i][j] = T_IN*pow(0.75, (double)(in-2));
//             else if (in >= 2) phi[i][j] = T_IN*pow(1.0 - 0.5*(double)(in-2), (double)(in-2));
//             else if (in >= 2) phi[i][j] = T_IN*(1.0 - (double)(in-2)/((double)MDEPTH))*(1.0 - (double)(in-2)/((double)MDEPTH));
            else phi[i][j] = T_OUT;
        }
}

void init_julia_set(double *phi[NX], short int * xy_in[NX])
/* change Julia set boundary condition */
{
    int i, j, in;
    double xy[2], dist2, module, phase, scale2;    

//     printf("Changing Julia set\n");
    for (i=0; i<NX; i++)
        for (j=0; j<NY; j++)
        {
            ij_to_xy(i, j, xy);
	    xy_in[i][j] = xy_in_billiard(xy[0],xy[1]);

            in = xy_in[i][j];
            if (in >= 2) phi[i][j] = T_IN;
        }
}


/*********************/
/* animation part    */
/*********************/


void compute_gradient(double *phi[NX], double *nablax[NX], double *nablay[NX])
/* compute the gradient of the field */
{
    int i, j, iplus, iminus, jplus, jminus; 
    double dx;
    
    dx = (XMAX-XMIN)/((double)NX);
    
    for (i=0; i<NX; i++)
        for (j=0; j<NY; j++)
        {
            iplus = i+1;  if (iplus == NX) iplus = NX-1;
            iminus = i-1; if (iminus == -1) iminus = 0;
            jplus = j+1;  if (jplus == NX) jplus = NY-1;
            jminus = j-1; if (jminus == -1) jminus = 0;
            nablax[i][j] = (phi[iplus][j] - phi[iminus][j])/dx;
            nablay[i][j] = (phi[i][jplus] - phi[i][jminus])/dx;
        }
}

void draw_field_line(double x, double y, short int *xy_in[NX], double *nablax[NX], 
                     double *nablay[NX], double delta, int nsteps)
/* draw a field line of the gradient, starting in (x,y) */
{
    double x1, y1, x2, y2, pos[2], nabx, naby, norm2, norm;
    int i = 0, ij[2], cont = 1;
    
    glColor3f(1.0, 1.0, 1.0);
//     glColor3f(0.0, 0.0, 0.0);
    glLineWidth(FIELD_LINE_WIDTH);
    x1 = x;
    y1 = y;
    
//     printf("Drawing field line \n");

    glEnable(GL_LINE_SMOOTH);
    glBegin(GL_LINE_STRIP);
    xy_to_pos(x1, y1, pos);
    glVertex2d(pos[0], pos[1]);
    
    i = 0;
    while ((cont)&&(i < nsteps))
    {
        xy_to_ij(x1, y1, ij);
        
        if (ij[0] < 0) ij[0] = 0;
        if (ij[0] > NX-1) ij[0] = NX-1;
        if (ij[1] < 0) ij[1] = 0;
        if (ij[1] > NY-1) ij[1] = NY-1;
        
        nabx = nablax[ij[0]][ij[1]];
        naby = nablay[ij[0]][ij[1]];
        
        norm2 = nabx*nabx + naby*naby;
        
        if (norm2 > 1.0e-14)
        {
            /* avoid too large step size */
            if (norm2 < 1.0e-9) norm2 = 1.0e-9;
            norm = sqrt(norm2);
            x1 = x1 + delta*nabx/norm;
            y1 = y1 + delta*naby/norm;
        }
        else cont = 0;
        
        if (!xy_in[ij[0]][ij[1]]) cont = 0;
        
        /* stop if the boundary is hit */
//         if (xy_in[ij[0]][ij[1]] != 1) cont = 0;
        
//         printf("x1 = %.3lg \t y1 = %.3lg \n", x1, y1);
                
        xy_to_pos(x1, y1, pos);
        glVertex2d(pos[0], pos[1]);
        
        i++;
    }
    glEnd();
}

void draw_wave(double *phi[NX], short int *xy_in[NX], double scale, int time)
/* draw the field */
{
    int i, j, iplus, iminus, jplus, jminus, ij[2], counter = 0;
    static int first = 1;
    double rgb[3], xy[2], x1, y1, x2, y2, dx, value, angle, dangle, intens, deltaintens, sum = 0.0;
    double *nablax[NX], *nablay[NX];
    static double linex[N_FIELD_LINES*FIELD_LINE_FACTOR], liney[N_FIELD_LINES*FIELD_LINE_FACTOR], distance[N_FIELD_LINES*FIELD_LINE_FACTOR], integral[N_FIELD_LINES*FIELD_LINE_FACTOR + 1];

    for (i=0; i<NX; i++) 
    {
        nablax[i] = (double *)malloc(NY*sizeof(double));
        nablay[i] = (double *)malloc(NY*sizeof(double));
    }
    
    /* compute the gradient */
    compute_gradient(phi, nablax, nablay);
    
    /* compute the position of origins of field lines */
    if ((first)&&(DRAW_FIELD_LINES))
    {
        first = 0;
        
        printf("computing linex\n");
        
//         x1 = LAMBDA + MU*1.01;
//         y1 = 1.0;
        x1 = 0.99*LAMBDA;
        y1 = 0.0;
        linex[0] = x1;
        liney[0] = y1;
        dangle = DPI/((double)(N_FIELD_LINES*FIELD_LINE_FACTOR));
            
        for (i = 1; i < N_FIELD_LINES*FIELD_LINE_FACTOR; i++)
        {
            angle = (double)i*dangle;
//             x2 = LAMBDA + MU*1.01*cos(angle);
//             y2 = 0.5 + MU*1.01*sin(angle);
            x2 = 0.99*LAMBDA*cos(angle);
            y2 = 0.99*LAMBDA*sin(angle);
            linex[i] = x2;
            liney[i] = y2;
            distance[i-1] = module2(x2-x1,y2-y1);
            x1 = x2;
            y1 = y2;
        }
        distance[N_FIELD_LINES*FIELD_LINE_FACTOR - 1] = module2(x2- 0.99*LAMBDA,y2);
//         distance[N_FIELD_LINES*FIELD_LINE_FACTOR - 1] = module2(x2-LAMBDA,y2-0.5);
    }

    dx = (XMAX-XMIN)/((double)NX);
    glBegin(GL_QUADS);

    for (i=0; i<NX; i++)
        for (j=0; j<NY; j++)
        {
            if (FIELD_REP == F_INTENSITY) value = phi[i][j];
            else if (FIELD_REP == F_GRADIENT)
            {
                value = module2(nablax[i][j], nablay[i][j]);
            }
            
            if (xy_in[i][j] == 1) 
            {
                color_scheme(COLOR_SCHEME, value, scale, time, rgb);
                glColor3f(rgb[0], rgb[1], rgb[2]);
            }
            else glColor3f(0.0, 0.0, 0.0);

            glVertex2i(i, j);
            glVertex2i(i+1, j);
            glVertex2i(i+1, j+1);
            glVertex2i(i, j+1);
        }
    glEnd ();
        
    /* draw a field line */
    if (DRAW_FIELD_LINES)
    {
        /* compute gradient norm along boundary and its integral */
        for (i = 0; i < N_FIELD_LINES*FIELD_LINE_FACTOR; i++)
        {
            xy_to_ij(linex[i], liney[i], ij);
            intens = module2(nablax[ij[0]][ij[1]], nablay[ij[0]][ij[1]])*distance[i];
            if (i > 0) integral[i] = integral[i-1] + intens;
            else integral[i] = intens;
        }
        deltaintens = integral[N_FIELD_LINES*FIELD_LINE_FACTOR-1]/((double)N_FIELD_LINES);
        
//         printf("delta = %.5lg\n", deltaintens);
        
        i = 0;
        draw_field_line(linex[0], liney[0], xy_in, nablax, nablay, 0.00002, 100000);
        for (j = 1; j < N_FIELD_LINES+1; j++)
        {
            while ((integral[i] <= j*deltaintens)&&(i < N_FIELD_LINES*FIELD_LINE_FACTOR)) i++; 
            draw_field_line(linex[i], liney[i], xy_in, nablax, nablay, 0.00002, 100000);
            counter++;
        }
        printf("%i lines\n", counter);
    }

    
    for (i=0; i<NX; i++)
    {
        free(nablax[i]);
        free(nablay[i]);
    }
}



void evolve_wave_half(double *phi_in[NX], double *phi_out[NX], short int *xy_in[NX])
/* time step of field evolution */
{
    int i, j, iplus, iminus, jplus, jminus;
    double delta1, delta2, x, y;
    
    #pragma omp parallel for private(i,j,iplus,iminus,jplus,jminus,delta1,delta2,x,y)
    for (i=0; i<NX; i++){
        for (j=0; j<NY; j++){
            if (xy_in[i][j] == 1){
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

                x = phi_in[i][j];

                /* evolve phi */
                if (B_COND != BC_ABSORBING)
                {
                    phi_out[i][j] = x + intstep*(delta1 - SPEED*(phi_in[iplus][j] - phi_in[i][j]));
                }
                else        /* case of absorbing b.c. - this is only an approximation of correct way of implementing */
                {
                    /* in the bulk */
                    if ((i>0)&&(i<NX-1)&&(j>0)&&(j<NY-1))
                    {
                        phi_out[i][j] = x - intstep*delta2;
                    }
                     /* right border */
                    else if (i==NX-1) 
                    {
                        phi_out[i][j] = x - intstep1*(x - phi_in[i-1][j]);
                    }
                    /* upper border */
                    else if (j==NY-1) 
                    {
                        phi_out[i][j] = x - intstep1*(x - phi_in[i][j-1]);
                    }
                    /* left border */
                    else if (i==0) 
                    {
                        phi_out[i][j] = x - intstep1*(x - phi_in[1][j]);
                    }
                   /* lower border */
                    else if (j==0) 
                    {
                        phi_out[i][j] = x - intstep1*(x - phi_in[i][1]);
                    }
                }


                if (FLOOR)
                {
                    if (phi_out[i][j] > VMAX) phi_out[i][j] = VMAX;
                    if (phi_out[i][j] < -VMAX) phi_out[i][j] = -VMAX;
                }
            }
        }
    }
    
//     printf("phi(0,0) = %.3lg, psi(0,0) = %.3lg\n", phi[NX/2][NY/2], psi[NX/2][NY/2]);
}

void evolve_wave(double *phi[NX], double *phi_tmp[NX], short int *xy_in[NX])
/* time step of field evolution */
{
    evolve_wave_half(phi, phi_tmp, xy_in);
    evolve_wave_half(phi_tmp, phi, xy_in);
}




double compute_variance(double *phi[NX], short int * xy_in[NX])
/* compute the variance (total probability) of the field */
{
    int i, j, n = 0;
    double variance = 0.0;

    for (i=1; i<NX; i++)
        for (j=1; j<NY; j++)
        {
            if (xy_in[i][j])
            {
                n++;
                variance += phi[i][j]*phi[i][j];
            }
        }
    if (n==0) n=1;
    return(variance/(double)n);
}

void renormalise_field(double *phi[NX], short int * xy_in[NX], double variance)
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
            }
        }
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


void print_Julia_parameters()
{
    double pos[2];
    char message[50];
    
    glColor3f(1.0, 1.0, 1.0);
    if (julia_y >= 0.0) sprintf(message, "c = %.5f + %.5f i", julia_x, julia_y);
    else sprintf(message, "c = %.5f %.5f i", julia_x, julia_y);
    xy_to_pos(XMIN + 0.1, YMAX - 0.2, pos);
    write_text(pos[0], pos[1], message);
}

void set_Julia_parameters(int time, double *phi[NX], short int *xy_in[NX])
{
    double jangle, cosj, sinj, radius = 0.15;

    jangle = (double)time*DPI/(double)NSTEPS;
//     jangle = (double)time*0.001;
//     jangle = (double)time*0.0001;

    cosj = cos(jangle);
    sinj = sin(jangle);
    julia_x = -0.9 + radius*cosj;
    julia_y = radius*sinj;
    init_julia_set(phi, xy_in);
    
    printf("Julia set parameters : i = %i, angle = %.5lg, cx = %.5lg, cy = %.5lg \n", time, jangle, julia_x, julia_y);
}

void set_Julia_parameters_cardioid(int time, double *phi[NX], short int *xy_in[NX])
{
    double jangle, cosj, sinj, yshift;

    jangle = pow(1.05 + (double)time*0.00003, 0.333);
    yshift = 0.02*sin((double)time*PID*0.002);
//     jangle = pow(1.0 + (double)time*0.00003, 0.333);
//     jangle = pow(0.05 + (double)time*0.00003, 0.333);
//     jangle = pow(0.1 + (double)time*0.00001, 0.333);
//     yshift = 0.04*sin((double)time*PID*0.002);

    cosj = cos(jangle);
    sinj = sin(jangle);
    julia_x = 0.5*(cosj*(1.0 - 0.5*cosj) + 0.5*sinj*sinj);
    julia_y = 0.5*sinj*(1.0-cosj) + yshift;
//     julia_x = 0.5*(cosj*(1.0 - 0.5*cosj) + 0.5*sinj*sinj);
//     julia_y = 0.5*sinj*(1.0-cosj);
    init_julia_set(phi, xy_in);
    
    printf("Julia set parameters : i = %i, angle = %.5lg, cx = %.5lg, cy = %.5lg \n", time, jangle, julia_x, julia_y);
}

void animation()
{
    double time, scale, dx, var, jangle, cosj, sinj;
    double *phi[NX], *phi_tmp[NX];
    short int *xy_in[NX];
    int i, j, s;

    /* Since NX and NY are big, it seemed wiser to use some memory allocation here */
    for (i=0; i<NX; i++)
    {
        phi[i] = (double *)malloc(NY*sizeof(double));
        phi_tmp[i] = (double *)malloc(NY*sizeof(double));
        xy_in[i] = (short int *)malloc(NY*sizeof(short int));
    }

    npolyline = init_polyline(MDEPTH, polyline);
    for (i=0; i<npolyline; i++) printf("vertex %i: (%.3f, %.3f)\n", i, polyline[i].x, polyline[i].y);

    dx = (XMAX-XMIN)/((double)NX);
    intstep = DT/(dx*dx*VISCOSITY);
    intstep1 = DT/(dx*VISCOSITY);
    
//     julia_x = 0.1; 
//     julia_y = 0.6; 
    
//     set_Julia_parameters(0, phi, xy_in);
    
    printf("Integration step %.3lg\n", intstep);

    /* initialize wave wave function */
    init_gaussian(-1.0, 0.0, 0.1, 0.0, 0.01, phi, xy_in);
//     init_gaussian(x, y, mean, amplitude, scalex, phi, xy_in)
    
    if (SCALE)
    {
        var = compute_variance(phi, xy_in);
        scale = sqrt(1.0 + var);
        renormalise_field(phi, xy_in, var);
    }

    blank();
    glColor3f(0.0, 0.0, 0.0);
    

    glutSwapBuffers();
    
    draw_wave(phi, xy_in, 1.0, 0);
    if (DRAW_BILLIARD) draw_billiard(0, 1.0);
//     print_Julia_parameters(i);
    
//     print_level(MDEPTH);

    glutSwapBuffers();

    sleep(SLEEP1);
    if (MOVIE) for (i=0; i<SLEEP1*25; i++) save_frame();

    for (i=0; i<=NSTEPS; i++)
    {
        /* compute the variance of the field to adjust color scheme */
        /* the color depends on the field divided by sqrt(1 + variance) */
        if (SCALE)
        {
            var = compute_variance(phi, xy_in);
            scale = sqrt(1.0 + var);
//             printf("Norm: %5lg\t Scaling factor: %5lg\n", var, scale);
            renormalise_field(phi, xy_in, var);
        }
        else scale = 1.0;
        
        draw_wave(phi, xy_in, scale, i);
        
        for (j=0; j<NVID; j++) evolve_wave(phi, phi_tmp, xy_in);

        if (DRAW_BILLIARD) draw_billiard(0, 1.0);
        
//         print_level(MDEPTH);
//         print_Julia_parameters(i);

	glutSwapBuffers();
        
        /* modify Julia set */
//         set_Julia_parameters(i, phi, xy_in);

	if (MOVIE)
        {
            save_frame();

            /* it seems that saving too many files too fast can cause trouble with the file system */
            /* so this is to make a pause from time to time - parameter PAUSE may need adjusting   */
            if (i % PAUSE == PAUSE - 1)
            {
                printf("Making a short pause\n");
                sleep(PSLEEP);
                s = system("mv wave*.tif tif_heat/");
            }
        }

    }

    if (MOVIE)
    {
        for (i=0; i<20; i++) save_frame();
        s = system("mv wave*.tif tif_heat/");
    }
    for (i=0; i<NX; i++)
    {
        free(phi[i]);
        free(phi_tmp[i]);
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
    glutCreateWindow("Heat equation in a planar domain");

    init();

    glutDisplayFunc(display);

    glutMainLoop();

    return 0;
}

