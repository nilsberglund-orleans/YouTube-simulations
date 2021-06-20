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

/* General geometrical parameters */

#define WINWIDTH 	1280  /* window width */
#define WINHEIGHT 	720   /* window height */

#define NX 1280          /* number of grid points on x axis */
#define NY 720          /* number of grid points on y axis */
// #define NX 640          /* number of grid points on x axis */
// #define NY 360          /* number of grid points on y axis */

/* setting NX to WINWIDTH and NY to WINHEIGHT increases resolution */
/* but will multiply run time by 4                                 */

#define XMIN -2.0
#define XMAX 2.0	/* x interval */
#define YMIN -1.125
#define YMAX 1.125	/* y interval for 9/16 aspect ratio */

#define JULIA_SCALE 1.1 /* scaling for Julia sets */

/* Choice of the billiard table */

#define B_DOMAIN 16      /* choice of domain shape */

#define D_RECTANGLE 0   /* rectangular domain */
#define D_ELLIPSE 1     /* elliptical domain */
#define D_STADIUM 2     /* stadium-shaped domain */
#define D_SINAI 3       /* Sinai billiard */
#define D_DIAMOND 4     /* diamond-shaped billiard */
#define D_TRIANGLE 5    /* triangular billiard */
#define D_FLAT 6        /* flat interface */
#define D_ANNULUS 7     /* annulus */
#define D_POLYGON 8     /* polygon */
#define D_YOUNG 9       /* Young diffraction slits */
#define D_GRATING 10    /* diffraction grating */
#define D_EHRENFEST 11  /* Ehrenfest urn type geometry */

#define D_MENGER 15     /* Menger-Sierpinski carpet */ 
#define D_JULIA_INT 16  /* interior of Julia set */ 

/* Billiard tables for heat equation */

#define D_ANNULUS_HEATED 21 /* annulus with different temperatures */
#define D_MENGER_HEATED 22  /* Menger gasket with different temperatures */
#define D_MENGER_H_OPEN 23  /* Menger gasket with different temperatures and larger domain */
#define D_MANDELBROT 24     /* Mandelbrot set */
#define D_JULIA 25          /* Julia set */
#define D_MANDELBROT_CIRCLE 26     /* Mandelbrot set with circular conductor */

#define LAMBDA 1.0	    /* parameter controlling the dimensions of domain */
#define MU 0.05	            /* parameter controlling the dimensions of domain */
#define NPOLY 8             /* number of sides of polygon */
#define APOLY 1.0           /* angle by which to turn polygon, in units of Pi/2 */ 
#define MDEPTH 4            /* depth of computation of Menger gasket */
#define MRATIO 3            /* ratio defining Menger gasket */
#define MANDELLEVEL 1000      /* iteration level for Mandelbrot set */
#define MANDELLIMIT 10.0     /* limit value for approximation of Mandelbrot set */
#define FOCI 1              /* set to 1 to draw focal points of ellipse */

/* You can add more billiard tables by adapting the functions */
/* xy_in_billiard and draw_billiard below */

/* Physical parameters of wave equation */

#define OMEGA 0.9           /* frequency of periodic excitation */
#define COURANT 0.01       /* Courant number */
#define GAMMA 0.0      /* damping factor in wave equation */
// #define GAMMA 5.0e-10      /* damping factor in wave equation */
#define KAPPA 0.0       /* "elasticity" term enforcing oscillations */
// #define KAPPA 5.0e-6       /* "elasticity" term enforcing oscillations */
// #define KAPPA 5.0e-9       /* "elasticity" term enforcing oscillations */
// #define KAPPA 5.0e-8       /* "elasticity" term enforcing oscillations */
/* The Courant number is given by c*DT/DX, where DT is the time step and DX the lattice spacing */
/* The physical damping coefficient is given by GAMMA/(DT)^2 */
/* Increasing COURANT speeds up the simulation, but decreases accuracy */
/* For similar wave forms, COURANT^2*GAMMA should be kept constant */

/* Boundary conditions */

#define B_COND 2

#define BC_DIRICHLET 0   /* Dirichlet boundary conditions */
#define BC_PERIODIC 1    /* periodic boundary conditions */
#define BC_ABSORBING 2   /* absorbing boundary conditions (beta version) */


/* For debugging purposes only */
#define FLOOR 0         /* set to 1 to limit wave amplitude to VMAX */
#define VMAX 10.0       /* max value of wave amplitude */

/* Parameters for length and speed of simulation */

#define NSTEPS 4000      /* number of frames of movie */
#define NVID 25          /* number of iterations between images displayed on screen */
#define NSEG 100         /* number of segments of boundary */

#define PAUSE 1000         /* number of frames after which to pause */
#define PSLEEP 1         /* sleep time during pause */
#define SLEEP1  1        /* initial sleeping time */
#define SLEEP2  1   /* final sleeping time */

/* Color schemes */

#define BLACK 1          /* background */

#define COLOR_SCHEME 1   /* choice of color scheme */

#define C_LUM 0          /* color scheme modifies luminosity (with slow drift of hue) */
#define C_HUE 1          /* color scheme modifies hue */

#define SCALE 0          /* set to 1 to adjust color scheme to variance of field */
#define SLOPE 1.0        /* sensitivity of color on wave amplitude */
#define ATTENUATION 0.0  /* exponential attenuation coefficient of contrast with time */

#define COLORHUE 260     /* initial hue of water color for scheme C_LUM */
#define COLORDRIFT 0.0   /* how much the color hue drifts during the whole simulation */
#define LUMMEAN 0.5      /* amplitude of luminosity variation for scheme C_LUM */
#define LUMAMP 0.3       /* amplitude of luminosity variation for scheme C_LUM */
#define HUEMEAN 230.0    /* mean value of hue for color scheme C_HUE */
#define HUEAMP 50.0      /* amplitude of variation of hue for color scheme C_HUE */
// #define HUEMEAN 320.0    /* mean value of hue for color scheme C_HUE */
// #define HUEAMP 100.0      /* amplitude of variation of hue for color scheme C_HUE */

/* Basic math */

#define PI 	3.141592654
#define DPI 	6.283185307
#define PID 	1.570796327

double julia_x = -0.5, julia_y = 0.5;    /* parameters for Julia sets */
// double julia_x = 0.33267, julia_y = 0.06395;    /* parameters for Julia sets */
// double julia_x = 0.37468, julia_y = 0.21115;    /* parameters for Julia sets */

#include "sub_wave.c"

double courant2;  /* Courant parameter squared */

void init_wave(x, y, phi, psi, xy_in)
/* initialise field with drop at (x,y) - phi is wave height, psi is phi at time t-1 */
    double x, y, *phi[NX], *psi[NX]; short int * xy_in[NX];

{
    int i, j;
    double xy[2], dist2;

    for (i=0; i<NX; i++)
        for (j=0; j<NY; j++)
        {
            ij_to_xy(i, j, xy);
            dist2 = (xy[0]-x)*(xy[0]-x) + (xy[1]-y)*(xy[1]-y);
	    xy_in[i][j] = xy_in_billiard(xy[0],xy[1]);
	    phi[i][j] = 0.2*exp(-dist2/0.001)*cos(-sqrt(dist2)/0.01);
            psi[i][j] = 0.0;
        }
}

void init_wave_flat(phi, psi, xy_in)
/* initialise flat field - phi is wave height, psi is phi at time t-1 */
    double *phi[NX], *psi[NX]; short int * xy_in[NX];

{
    int i, j;
    double xy[2], dist2;

    for (i=0; i<NX; i++)
        for (j=0; j<NY; j++)
        {
            ij_to_xy(i, j, xy);
	    xy_in[i][j] = xy_in_billiard(xy[0],xy[1]);
	    phi[i][j] = 0.0;
            psi[i][j] = 0.0;
        }
}

void add_drop_to_wave(factor, x, y, phi, psi)
/* add drop at (x,y) to the field with given prefactor */
double factor, x, y, *phi[NX], *psi[NX];
{
    int i, j;
    double xy[2], dist2;

    for (i=0; i<NX; i++)
        for (j=0; j<NY; j++)
        {
            ij_to_xy(i, j, xy);
            dist2 = (xy[0]-x)*(xy[0]-x) + (xy[1]-y)*(xy[1]-y);
            phi[i][j] += 0.2*factor*exp(-dist2/0.001)*cos(-sqrt(dist2)/0.01);
        }
}

void oscillate_linear_wave(amplitude, t, x1, y1, x2, y2, phi, psi)
/* oscillating boundary condition at (x,y) */
double amplitude, t, x1, y1, x2, y2, *phi[NX], *psi[NX];
{
    int i, j, ij1[2], ij2[2], imin, imax, jmin, jmax, d = 5;
    double xy[2], dist2;
    
    xy_to_ij(x1, y1, ij1);
    xy_to_ij(x2, y2, ij2);
    imin = ij1[0] - d;     if (imin < 0) imin = 0;
    imax = ij2[0] + d;     if (imax >= NX) imax = NX-1;
    jmin = ij1[1] - d;     if (jmin < 0) jmin = 0;
    jmax = ij2[1] + d;     if (jmax >= NY) jmax = NY-1;

    for (i = imin; i < imax; i++)
        for (j = jmin; j < jmax; j++)
        {
            ij_to_xy(i, j, xy);
            dist2 = (xy[0]-x1)*(xy[0]-x1);      /* to be improved */
//             dist2 = (xy[0]-x)*(xy[0]-x) + (xy[1]-y)*(xy[1]-y);
//             if (dist2 < 0.01)
            if (dist2 < 0.001)
                phi[i][j] = amplitude*exp(-dist2/0.001)*cos(-sqrt(dist2)/0.01)*cos(t*OMEGA);
//                 phi[i][j] += 0.2*exp(-dist2/0.001)*cos(-sqrt(dist2)/0.01)*cos(t*OMEGA);
        }
}




/*********************/
/* animation part    */
/*********************/


void draw_wave(phi, psi, xy_in, scale, time)
/* draw the field */
double *phi[NX], *psi[NX], scale;
short int *xy_in[NX];
int time;
{
    int i, j;
    double rgb[3], xy[2], x1, y1, x2, y2;

    glBegin(GL_QUADS);

    for (i=0; i<NX; i++)
        for (j=0; j<NY; j++)
        {
            if (xy_in[i][j])
            {
                color_scheme(COLOR_SCHEME, phi[i][j], scale, time, rgb);
                glColor3f(rgb[0], rgb[1], rgb[2]);

                glVertex2i(i, j);
                glVertex2i(i+1, j);
                glVertex2i(i+1, j+1);
                glVertex2i(i, j+1);
            }
        }

    glEnd ();
}

void evolve_wave_half(phi_in, psi_in, phi_out, psi_out, xy_in)
/* time step of field evolution */
/* phi is value of field at time t, psi at time t-1 */
    double *phi_in[NX], *psi_in[NX], *phi_out[NX], *psi_out[NX]; short int *xy_in[NX];
{
    int i, j, iplus, iminus, jplus, jminus;
    double delta, x, y, c, cc;
    
    c = COURANT;
    cc = courant2;

    #pragma omp parallel for private(i,j,iplus,iminus,jplus,jminus,delta,x,y)
    for (i=0; i<NX; i++){
        for (j=0; j<NY; j++){
            if (xy_in[i][j]){
                /* discretized Laplacian */
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
                delta = phi_in[iplus][j] + phi_in[iminus][j] + phi_in[i][jplus] + phi_in[i][jminus] - 4.0*phi_in[i][j];

                x = phi_in[i][j];
		y = psi_in[i][j];

                /* evolve phi */
                if (B_COND != BC_ABSORBING) phi_out[i][j] = -y + 2*x + cc*delta - KAPPA*x - GAMMA*(x-y);
                else
                {
                    if ((i>0)&&(i<NX-1)&&(j>0)&&(j<NY-1))
                        phi_out[i][j] = -y + 2*x + cc*delta - KAPPA*x - GAMMA*(x-y);
                
                    /* upper border */
                    else if (j==NY-1) 
                        phi_out[i][j] = x - c*(x - phi_in[i][NY-2]) - KAPPA*x - GAMMA*(x-y);
                    
                    /* lower border */
                    else if (j==0) 
                        phi_out[i][j] = x - c*(x - phi_in[i][1]) - KAPPA*x - GAMMA*(x-y);
                
                    /* right border */
                    if (i==NX-1) 
                        phi_out[i][j] = x - c*(x - phi_in[NX-2][j]) - KAPPA*x - GAMMA*(x-y);
                    
                    /* left border */
                    else if (i==0) 
                        phi_out[i][j] = x - c*(x - phi_in[1][j]) - KAPPA*x - GAMMA*(x-y);
                }
                psi_out[i][j] = x;

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


void evolve_wave(phi, psi, phi_tmp, psi_tmp, xy_in)
/* time step of field evolution */
/* phi is value of field at time t, psi at time t-1 */
    double *phi[NX], *psi[NX], *phi_tmp[NX], *psi_tmp[NX]; short int *xy_in[NX];
{
    evolve_wave_half(phi, psi, phi_tmp, psi_tmp, xy_in);
    evolve_wave_half(phi_tmp, psi_tmp, phi, psi, xy_in);
}

void old_evolve_wave(phi, psi, xy_in)
/* time step of field evolution */
/* phi is value of field at time t, psi at time t-1 */
    double *phi[NX], *psi[NX]; short int *xy_in[NX];
{
    int i, j, iplus, iminus, jplus, jminus;
    double delta, x, y;

    #pragma omp parallel for private(i,j,iplus,iminus,jplus,jminus,delta,x,y)
    for (i=0; i<NX; i++){
        for (j=0; j<NY; j++){
            if (xy_in[i][j]){
                /* discretized Laplacian */
		iplus = (i+1) % NX;
		iminus = (i-1) % NX;
		if (iminus < 0) iminus += NX;
                jplus = (j+1) % NY;
                jminus = (j-1) % NY;
                if (jminus < 0) jminus += NY;
                delta = phi[iplus][j] + phi[iminus][j] + phi[i][jplus] + phi[i][jminus] - 4.0*phi[i][j];

                x = phi[i][j];
		y = psi[i][j];

                /* evolve phi */
                phi[i][j] = -y + 2*x + courant2*delta - KAPPA*x - GAMMA*(x-y);

                /* Old versions of the simulation used this: */
//                 phi[i][j] = (-psi[i][j] + 2*phi[i][j] + courant2*delta)*damping;
//                 where damping = 1.0 - 0.0001;

                psi[i][j] = x;

                if (FLOOR)
                {
                    if (phi[i][j] > VMAX) phi[i][j] = VMAX;
                    if (phi[i][j] < -VMAX) phi[i][j] = -VMAX;
                    if (psi[i][j] > VMAX) psi[i][j] = VMAX;
                    if (psi[i][j] < -VMAX) psi[i][j] = -VMAX;
                }
            }
        }
    }
//     printf("phi(0,0) = %.3lg, psi(0,0) = %.3lg\n", phi[NX/2][NY/2], psi[NX/2][NY/2]);
}


double compute_variance(phi, psi, xy_in)
/* compute the variance of the field, to adjust color scheme */
    double *phi[NX], *psi[NX]; short int * xy_in[NX];
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


void animation()
{
    double time, scale;
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

    courant2 = COURANT*COURANT;

    /* initialize wave with a drop at one point, zero elsewhere */
//     init_wave_flat(phi, psi, xy_in);
    init_wave(0.0, 0.0, phi, psi, xy_in);

    /* add a drop at another point */
//     add_drop_to_wave(1.0, 0.7, 0.0, phi, psi);
//     add_drop_to_wave(1.0, -0.7, 0.0, phi, psi);
//     add_drop_to_wave(1.0, 0.0, -0.7, phi, psi);

    blank();
    glColor3f(0.0, 0.0, 0.0);
    draw_wave(phi, psi, xy_in, 1.0, 0);
    draw_billiard();

    glutSwapBuffers();



    sleep(SLEEP1);

    for (i=0; i<=NSTEPS; i++)
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

//         /* TO BE ADAPTED */
//         scale = 1.0;

        draw_wave(phi, psi, xy_in, scale, i);
        for (j=0; j<NVID; j++) 
        {
            evolve_wave(phi, psi, phi_tmp, psi_tmp, xy_in);
//             if (i % 10 == 9) oscillate_linear_wave(0.2*scale, 0.15*(double)(i*NVID + j), -1.5, YMIN, -1.5, YMAX, phi, psi);
        }
//         for (j=0; j<NVID; j++) evolve_wave(phi, psi, phi_tmp, psi_tmp, xy_in);
        
        draw_billiard();


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
                s = system("mv wave*.tif tif_wave/");
            }
        }

    }

    if (MOVIE) 
    {
        for (i=0; i<20; i++) save_frame();
        s = system("mv wave*.tif tif_wave/");
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
    glutCreateWindow("Wave equation in a planar domain");

    init();

    glutDisplayFunc(display);

    glutMainLoop();

    return 0;
}

