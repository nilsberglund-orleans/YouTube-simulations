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
#define DOUBLE_MOVIE 0  /* set to 1 to produce movies for wave height and energy simultaneously */

/* General geometrical parameters */

#define WINWIDTH 	1920  /* window width */
#define WINHEIGHT 	1000  /* window height */
#define NX 480          /* number of grid points on x axis */
#define NY 240          /* number of grid points on y axis */

#define XMIN -2.0
#define XMAX 2.0	/* x interval */
#define YMIN -1.041666667
#define YMAX 1.041666667	/* y interval for 9/16 aspect ratio */

// #define WINWIDTH 	1280  /* window width */
// #define WINHEIGHT 	720   /* window height */
// 
// // #define NX 160          /* number of grid points on x axis */
// // #define NY 90          /* number of grid points on y axis */
// #define NX 320          /* number of grid points on x axis */
// #define NY 180          /* number of grid points on y axis */
// 
// // #define NX 640          /* number of grid points on x axis */
// // #define NY 360          /* number of grid points on y axis */
// 
// // #define NX 1280          /* number of grid points on x axis */
// // #define NY 720          /* number of grid points on y axis */
// 
// #define XMIN -2.0
// #define XMAX 2.0	/* x interval */
// #define YMIN -1.125
// #define YMAX 1.125	/* y interval for 9/16 aspect ratio */

/* Choice of simulated equation */

#define RDE_EQUATION 5  /* choice of reaction term, see list in global_3d.c */
#define NFIELDS 2       /* number of fields in reaction-diffusion equation */
#define NLAPLACIANS 2   /* number of fields for which to compute Laplacian */

#define ADD_POTENTIAL 1 /* set to 1 to add a potential (for Schrodiner equation) */
#define POTENTIAL 2     /* type of potential, see list in global_3d.c  */

#define JULIA_SCALE 0.5 /* scaling for Julia sets */

/* Choice of the billiard table */

#define B_DOMAIN 999      /* choice of domain shape, see list in global_pdes.c  */

#define CIRCLE_PATTERN 99    /* pattern of circles, see list in global_pdes.c */

#define P_PERCOL 0.25       /* probability of having a circle in C_RAND_PERCOL arrangement */
#define NPOISSON 300        /* number of points for Poisson C_RAND_POISSON arrangement */
#define RANDOM_POLY_ANGLE 0 /* set to 1 to randomize angle of polygons */

#define LAMBDA 1.0	    /* parameter controlling the dimensions of domain */
#define MU 1.0	            /* parameter controlling the dimensions of domain */
#define NPOLY 6             /* number of sides of polygon */
#define APOLY 0.333333333   /* angle by which to turn polygon, in units of Pi/2 */
#define MDEPTH 7            /* depth of computation of Menger gasket */
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

#define DT 0.00000002

#define VISCOSITY 2.0

#define RPSA 0.75         /* parameter in Rock-Paper-Scissors-type interaction */
#define RPSLZB 0.75       /* second parameter in Rock-Paper-Scissors-Lizard-Spock type interaction */

#define EPSILON 0.8     /* time scale separation */
#define DELTA 0.1       /* time scale separation */
#define FHNA 1.0        /* parameter in FHN equation */
#define FHNC -0.01      /* parameter in FHN equation */
#define K_HARMONIC 0.2  /* spring constant of harmonic potential */
#define K_COULOMB 0.5   /* constant in Coulomb potential */
#define BZQ 0.0008      /* parameter in BZ equation */
#define BZF 1.2         /* parameter in BZ equation */

#define T_OUT 2.0       /* outside temperature */
#define T_IN 0.0        /* inside temperature */
#define SPEED 0.0       /* speed of drift to the right */

#define ADD_NOISE 0     /* set to 1 to add noise, set to 2 to add noise in right half */
#define NOISE_INTENSITY 0.005      /* noise intensity */
#define CHANGE_NOISE 1      /* set to 1 to increase noise intensity */
#define NOISE_FACTOR 40.0   /* factor by which to increase noise intensity */
#define NOISE_INITIAL_TIME 100  /* initial time during which noise remains constant */

#define CHANGE_VISCOSITY 0      /* set to 1 to change the viscosity in the course of the simulation */
#define ADJUST_INTSTEP 0       /* set to 1 to decrease integration step when viscosity increases */
#define VISCOSITY_INITIAL_TIME 10  /* initial time during which viscosity remains constant */
#define VISCOSITY_FACTOR 100.0   /* factor by which to change viscosity */
#define VISCOSITY_MAX 2.0        /* max value of viscosity beyond which NVID is increased and integration step is decrase, 
                                    for numerical stability */
                                        
#define CHANGE_RPSLZB 0         /* set to 1 to change second parameter in Rock-Paper-Scissors-Lizard-Spock equation */
#define RPSLZB_CHANGE 0.75      /* factor by which to rpslzb parameter */
#define RPSLZB_INITIAL_TIME 0   /* initial time during which rpslzb remains constant */
#define RPSLZB_FINAL_TIME 500   /* final time during which rpslzb remains constant */
                                      

/* Boundary conditions, see list in global_pdes.c  */

#define B_COND 1

/* Parameters for length and speed of simulation */

#define NSTEPS 1150      /* number of frames of movie */
#define NVID 850          /* number of iterations between images displayed on screen */
#define ACCELERATION_FACTOR 1.0 /* factor by which to increase NVID in course of simulation */
#define DT_ACCELERATION_FACTOR 1.0 /* factor by which to increase time step in course of simulation  */
#define MAX_DT 0.024     /* maximal value of integration step */
#define NSEG 100         /* number of segments of boundary */
#define BOUNDARY_WIDTH 4    /* width of billiard boundary */

#define PAUSE 100       /* number of frames after which to pause */
#define PSLEEP 2         /* sleep time during pause */
#define SLEEP1  2        /* initial sleeping time */
#define SLEEP2  1        /* final sleeping time */
#define INITIAL_TIME 0  /* initial still time */
#define MID_FRAMES 50    /* number of still frames between parts of two-part movie */
#define END_FRAMES 50    /* number of still frames at end of movie */
#define FADE 1           /* set to 1 to fade at end of movie */

/* Visualisation */

#define PLOT_3D 1    /* controls whether plot is 2D or 3D */

#define ROTATE_VIEW 1       /* set to 1 to rotate position of observer */
#define ROTATE_ANGLE 360.0  /* total angle of rotation during simulation */

/* Plot type - color scheme */

#define CPLOT 30
#define CPLOT_B 31

/* Plot type - height of 3D plot */

// #define ZPLOT 30     /* z coordinate in 3D plot */
// #define ZPLOT_B 32    /* z coordinate in second 3D plot */
#define ZPLOT 30     /* z coordinate in 3D plot */
#define ZPLOT_B 30    /* z coordinate in second 3D plot */

#define AMPLITUDE_HIGH_RES 1    /* set to 1 to increase resolution of P_3D_AMPLITUDE plot */
#define SHADE_3D 1              /* set to 1 to change luminosity according to normal vector */
#define NON_DIRICHLET_BC 0      /* set to 1 to draw only facets in domain, if field is not zero on boundary */
#define WRAP_ANGLE 1            /* experimental: wrap angle to [0, 2Pi) for interpolation in angle schemes */
#define FADE_IN_OBSTACLE 0      /* set to 1 to fade color inside obstacles */
#define DRAW_OUTSIDE_GRAY 0     /* experimental - draw outside of billiard in gray */
#define ADD_POTENTIAL_TO_Z 1    /* set to 1 to add the external potential to z-coordinate of plot */
#define ADD_POT_CONSTANT 0.5    /* constant in front of added potential */

#define PLOT_SCALE_ENERGY 0.05      /* vertical scaling in energy plot */

#define PRINT_TIME 0        /* set to 1 to print running time */
#define PRINT_VISCOSITY 0   /* set to 1 to print viscosity */
#define PRINT_RPSLZB 0      /* set to 1 to print rpslzb parameter */
#define PRINT_PROBABILITIES 0   /* set to 1 to print probabilities (for Ehrenfest urn configuration) */
#define PRINT_NOISE 0       /* set to 1 to print noise intensity */

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

#define COLOR_PALETTE 14     /* Color palette, see list in global_pdes.c  */
#define COLOR_PALETTE_B 10     /* Color palette, see list in global_pdes.c  */

#define BLACK 1          /* black background */

#define COLOR_SCHEME 3   /* choice of color scheme */

#define COLOR_PHASE_SHIFT 0.0   /* phase shift of color scheme, in units of Pi */

#define SCALE 0          /* set to 1 to adjust color scheme to variance of field */
#define SLOPE 1.0        /* sensitivity of color on wave amplitude */
#define VSCALE_AMPLITUDE 7.5      /* additional scaling factor for color scheme P_3D_AMPLITUDE */
#define ATTENUATION 0.0  /* exponential attenuation coefficient of contrast with time */
#define CURL_SCALE 0.000015   /* scaling factor for curl representation */
#define RESCALE_COLOR_IN_CENTER 0   /* set to 1 to decrease color intentiy in the center (for wave escaping ring) */
#define SLOPE_SCHROD_LUM 15.0       /* sensitivity of luminosity on module, for color scheme Z_ARGUMENT */
#define MIN_SCHROD_LUM 0.075      /* minimal luminosity in color scheme Z_ARGUMENT*/

#define COLORHUE 260     /* initial hue of water color for scheme C_LUM */
#define COLORDRIFT 0.0   /* how much the color hue drifts during the whole simulation */
#define LUMMEAN 0.5      /* amplitude of luminosity variation for scheme C_LUM */
#define LUMAMP 0.3       /* amplitude of luminosity variation for scheme C_LUM */
#define HUEMEAN 359.0    /* mean value of hue for color scheme C_HUE */
#define HUEAMP -359.0      /* amplitude of variation of hue for color scheme C_HUE */
#define E_SCALE 100.0     /* scaling factor for energy representation */
#define LOG_SCALE 1.0     /* scaling factor for energy log representation */
#define LOG_SHIFT 0.0 

#define DRAW_COLOR_SCHEME 1     /* set to 1 to plot the color scheme */
#define COLORBAR_RANGE 3.0      /* scale of color scheme bar */
#define COLORBAR_RANGE_B 3.0   /* scale of color scheme bar for 2nd part */
#define ROTATE_COLOR_SCHEME 0   /* set to 1 to draw color scheme horizontally */

/* only for compatibility with wave_common.c */
#define TWOSPEEDS 0          /* set to 1 to replace hardcore boundary by medium with different speed */
#define OMEGA 0.005        /* frequency of periodic excitation */
#define COURANT 0.08       /* Courant number */
#define COURANTB 0.03      /* Courant number in medium B */
#define INITIAL_AMP 0.5         /* amplitude of initial condition */
#define INITIAL_VARIANCE 0.0002  /* variance of initial condition */
#define INITIAL_WAVELENGTH  0.1  /* wavelength of initial condition */
#define VSCALE_ENERGY 200.0       /* additional scaling factor for color scheme P_3D_ENERGY */
#define PHASE_FACTOR 20.0       /* factor in computation of phase in color scheme P_3D_PHASE */
#define PHASE_SHIFT 0.0      /* shift of phase in color scheme P_3D_PHASE */
/* end of constants added only for compatibility with wave_common.c */


double u_3d[2] = {0.75, -0.45};     /* projections of basis vectors for REP_AXO_3D representation */
double v_3d[2] = {-0.75, -0.45};
double w_3d[2] = {0.0, 0.015};
double light[3] = {0.816496581, -0.40824829, 0.40824829};      /* vector of "light" direction for P_3D_ANGLE color scheme */
double observer[3] = {8.0, 8.0, 6.0};    /* location of observer for REP_PROJ_3D representation */ 
int reset_view = 0;         /* switch to reset 3D view parameters (for option ROTATE_VIEW) */

#define Z_SCALING_FACTOR 0.75   /* overall scaling factor of z axis for REP_PROJ_3D representation */
#define XY_SCALING_FACTOR 2.2   /* overall scaling factor for on-screen (x,y) coordinates after projection */
#define ZMAX_FACTOR 1.0         /* max value of z coordinate for REP_PROJ_3D representation */
#define XSHIFT_3D -0.1           /* overall x shift for REP_PROJ_3D representation */
#define YSHIFT_3D 0.2          /* overall y shift for REP_PROJ_3D representation */


/* For debugging purposes only */
#define FLOOR 1         /* set to 1 to limit wave amplitude to VMAX */
#define VMAX 2.0        /* max value of wave amplitude */

#define REFRESH_B (ZPLOT_B != ZPLOT)||(CPLOT_B != CPLOT)    /* to save computing time, to be improved */
#define COMPUTE_WRAP_ANGLE ((WRAP_ANGLE)&&((cplot == Z_ANGLE_GRADIENT)||(cplot == Z_ANGLE_GRADIENTX)||(cplot == Z_ARGUMENT)||(cplot == Z_ANGLE_GRADIENTX)))
#define PRINT_PARAMETERS ((PRINT_TIME)||(PRINT_VISCOSITY)||(PRINT_RPSLZB)||(PRINT_PROBABILITIES)||(PRINT_NOISE))

#include "global_pdes.c"
#include "sub_wave.c"
#include "wave_common.c"        /* common functions for wave_billiard, wave_comparison, etc */

#include "global_3d.c"          /* constants and global variables */
#include "sub_wave_3d_rde.c"    /* should be later replaced by sub_wave_rde.c */

#include "sub_rde.c"    


double potential(int i, int j)
/* compute potential (e.g. for Schrödinger equation) */
{
    double x, y, xy[2], r, small = 2.0e-1, kx, ky;
    
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
        default:
        {
            return(0.0);
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

void evolve_wave_half(double *phi_in[NFIELDS], double *phi_out[NFIELDS], short int xy_in[NX*NY], double potential_field[NX*NY])
/* time step of field evolution */
{
    int i, j, k, iplus, iminus, jplus, jminus;
    double x, y, z, deltax, deltay, deltaz, rho, pot;
    double *delta_phi[NLAPLACIANS];
    static double invsqr3 = 0.577350269;    /* 1/sqrt(3) */
    
    for (i=0; i<NLAPLACIANS; i++) delta_phi[i] = (double *)malloc(NX*NY*sizeof(double));
    
    /* compute the Laplacian of phi */
    for (i=0; i<NLAPLACIANS; i++) compute_laplacian(phi_in[i], delta_phi[i], xy_in);
    
    #pragma omp parallel for private(i,j,k,x,y,z,deltax,deltay,deltaz,rho)
    for (i=0; i<NX; i++){
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
                    phi_out[0][i*NY+j] = phi_in[0][i*NY+j] + intstep*(deltax + x*(1.0-x*x));
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
                    if (ADD_POTENTIAL)
                    {
                        pot = potential_field[i*NY+j];
                        phi_out[0][i*NY+j] += intstep*pot*phi_in[1][i*NY+j];
                        phi_out[1][i*NY+j] -= intstep*pot*phi_in[0][i*NY+j];
                    }
                }
            }
        }
    }
                
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
}

void evolve_wave(double *phi[NFIELDS], double *phi_tmp[NFIELDS], short int xy_in[NX*NY], double potential_field[NX*NY])
/* time step of field evolution */
{
    evolve_wave_half(phi, phi_tmp, xy_in, potential_field);
    evolve_wave_half(phi_tmp, phi, xy_in, potential_field);
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



void print_parameters(t_rde rde[NX*NY], short int xy_in[NX*NY], double time, short int left, double viscosity, double noise)
{
    char message[100];
    double density, hue, rgb[3], logratio, x, y, pos[2], probas[2];
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
                xbox = XMAX - 0.39;
                xtext = XMAX - 0.55;
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
                xbox = XMAX - 0.39;
                xtext = XMAX - 0.61;
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
        if (PLOT_3D) write_text(xtext, y, message);
        else
        {
            xy_to_pos(xtext, y, pos);
            write_text(pos[0], pos[1], message);
        }
    }
}

void draw_color_bar_palette(int plot, double range, int palette, int fade, double fade_value)
{
    double width = 0.14;
//     double width = 0.2;
    
    if (ROTATE_COLOR_SCHEME) 
        draw_color_scheme_palette_3d(-1.0, -0.8, XMAX - 0.1, -1.0, plot, -range, range, palette, fade, fade_value);
    else 
        draw_color_scheme_palette_3d(XMAX - 1.5*width, YMIN + 0.1, XMAX - 0.5*width, YMAX - 0.1, plot, -range, range, palette, fade, fade_value);
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
    double time = 0.0, scale, dx, var, jangle, cosj, sinj, sqrintstep, 
        intstep0, viscosity_printed, fade_value, noise = NOISE_INTENSITY;
    double *phi[NFIELDS], *phi_tmp[NFIELDS], *potential_field;
    short int *xy_in;
    int i, j, k, s, nvid, field;
    static int counter = 0;
    t_rde *rde;

    /* Since NX and NY are big, it seemed wiser to use some memory allocation here */
    for (i=0; i<NFIELDS; i++)
    {
        phi[i] = (double *)malloc(NX*NY*sizeof(double));
        phi_tmp[i] = (double *)malloc(NX*NY*sizeof(double));
    }

    xy_in = (short int *)malloc(NX*NY*sizeof(short int));
    rde = (t_rde *)malloc(NX*NY*sizeof(t_rde));
    
    if (ADD_POTENTIAL) 
    {
        potential_field = (double *)malloc(NX*NY*sizeof(double));
        initialize_potential(potential_field);
    }

    npolyline = init_polyline(MDEPTH, polyline);
    for (i=0; i<npolyline; i++) printf("vertex %i: (%.3f, %.3f)\n", i, polyline[i].x, polyline[i].y);

    dx = (XMAX-XMIN)/((double)NX);
    intstep = DT/(dx*dx);
    
    intstep0 = intstep;
    intstep1 = DT/dx;
    
    viscosity = VISCOSITY;
    
    sqrintstep = sqrt(intstep*(double)NVID);
        
    printf("Integration step %.3lg\n", intstep);

    /* initialize field */
//     init_random(0.5, 0.4, phi, xy_in);
//     init_random(0.0, 0.4, phi, xy_in);
//     init_gaussian(x, y, mean, amplitude, scalex, phi, xy_in)
    init_coherent_state(-0.7, 0.0, 3.5, 0.0, 0.15, phi, xy_in);
    
    init_cfield_rde(phi, xy_in, CPLOT, rde, 0);
    if (PLOT_3D) init_zfield_rde(phi, xy_in, ZPLOT, rde, 0);
    
    if (DOUBLE_MOVIE)
    {
        init_cfield_rde(phi, xy_in, CPLOT_B, rde, 1);
        if (PLOT_3D) init_zfield_rde(phi, xy_in, ZPLOT_B, rde, 1);
    }
    
    blank();
    glColor3f(0.0, 0.0, 0.0);
    

    glutSwapBuffers();
    
    printf("Drawing wave\n");
    draw_wave_rde(0, phi, xy_in, rde, potential_field, ZPLOT, CPLOT, COLOR_PALETTE, 0, 1.0, 1);
//     draw_billiard();
    if (PRINT_PARAMETERS) print_parameters(rde, xy_in, time, 0, VISCOSITY, noise);
    if (DRAW_COLOR_SCHEME) draw_color_bar_palette(CPLOT, COLORBAR_RANGE, COLOR_PALETTE, 0, 1.0);

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
        
        if (ROTATE_VIEW) 
        {
            viewpoint_schedule(i - INITIAL_TIME);
            reset_view = 1;
        }
        
        printf("Drawing wave %i\n", i);
        draw_wave_rde(0, phi, xy_in, rde, potential_field, ZPLOT, CPLOT, COLOR_PALETTE, 0, 1.0, 1);
        
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
        for (j=0; j<nvid; j++) evolve_wave(phi, phi_tmp, xy_in, potential_field);
        
        for (j=0; j<NFIELDS; j++) printf("field[%i] = %.3lg\t", j, phi[j][0]);
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
        if (PRINT_PARAMETERS) print_parameters(rde, xy_in, time, 0, viscosity_printed, noise);
        if (DRAW_COLOR_SCHEME) draw_color_bar_palette(CPLOT, COLORBAR_RANGE, COLOR_PALETTE, 0, 1.0); 
        
//         print_level(MDEPTH);
//         print_Julia_parameters(i);

	glutSwapBuffers();
        
        /* modify Julia set */
//         set_Julia_parameters(i, phi, xy_in);

	if (MOVIE)
        {
            printf("Saving frame %i\n", i);
            save_frame();
            
            if ((i >= INITIAL_TIME)&&(DOUBLE_MOVIE))
            {
                draw_wave_rde(1, phi, xy_in, rde, potential_field, ZPLOT_B, CPLOT_B, COLOR_PALETTE_B, 0, 1.0, REFRESH_B);
//                 draw_billiard();
                if (PRINT_PARAMETERS) print_parameters(rde, xy_in, time, 0, viscosity_printed, noise);
                if (DRAW_COLOR_SCHEME) draw_color_bar_palette(CPLOT_B, COLORBAR_RANGE_B, COLOR_PALETTE_B, 0, 1.0);  
                glutSwapBuffers();
                save_frame_counter(NSTEPS + MID_FRAMES + 1 + counter);
                counter++;
            }

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
            draw_wave_rde(0, phi, xy_in, rde, potential_field, ZPLOT, CPLOT, COLOR_PALETTE, 0, 1.0, 1);
//             draw_billiard();
            if (PRINT_PARAMETERS) print_parameters(rde, xy_in, time, 0, viscosity_printed, noise);
            if (DRAW_COLOR_SCHEME) draw_color_bar_palette(CPLOT, COLORBAR_RANGE, COLOR_PALETTE, 0, 1.0);   
            glutSwapBuffers();
            
            if (!FADE) for (i=0; i<MID_FRAMES; i++) save_frame();
            else for (i=0; i<MID_FRAMES; i++) 
            {
                fade_value = 1.0 - (double)i/(double)MID_FRAMES;
                draw_wave_rde(0, phi, xy_in, rde, potential_field, ZPLOT, CPLOT, COLOR_PALETTE, 1, fade_value, 0);
//                 draw_billiard();
                if (PRINT_PARAMETERS) print_parameters(rde, xy_in, time, 0, viscosity_printed, noise);
                if (DRAW_COLOR_SCHEME) draw_color_bar_palette(CPLOT, COLORBAR_RANGE, COLOR_PALETTE, 1, fade_value);   
                glutSwapBuffers();
                save_frame_counter(NSTEPS + i + 1);
            }
            draw_wave_rde(1, phi, xy_in, rde, potential_field, ZPLOT_B, CPLOT_B, COLOR_PALETTE_B, 0, 1.0, REFRESH_B);
            if (DRAW_COLOR_SCHEME) draw_color_bar_palette(CPLOT_B, COLORBAR_RANGE_B, COLOR_PALETTE_B, 0, 1.0); 
            glutSwapBuffers();
            
            if (!FADE) for (i=0; i<END_FRAMES; i++) save_frame_counter(NSTEPS + MID_FRAMES + 1 + counter + i);
            else for (i=0; i<END_FRAMES; i++) 
            {
                fade_value = 1.0 - (double)i/(double)END_FRAMES;
                draw_wave_rde(1, phi, xy_in, rde, potential_field, ZPLOT_B, CPLOT_B, COLOR_PALETTE_B, 1, fade_value, 0);
                if (DRAW_COLOR_SCHEME) draw_color_bar_palette(CPLOT_B, COLORBAR_RANGE_B, COLOR_PALETTE_B, 1, fade_value);   
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
                draw_wave_rde(0, phi, xy_in, rde, potential_field, ZPLOT, CPLOT, COLOR_PALETTE, 1, fade_value, 0);
                if (DRAW_COLOR_SCHEME) draw_color_bar_palette(CPLOT, COLORBAR_RANGE, COLOR_PALETTE, 1, fade_value); 
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
    if (ADD_POTENTIAL) free(potential_field);
    
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
    else init();

    glutDisplayFunc(display);

    glutMainLoop();

    return 0;
}

