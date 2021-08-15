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

#define WINWIDTH 	1280  /* window width */
#define WINHEIGHT 	720   /* window height */

#define NX 1280          /* number of grid points on x axis */
#define NY 720          /* number of grid points on y axis */

#define XMIN -2.0
#define XMAX 2.0	/* x interval */
#define YMIN -1.125
#define YMAX 1.125	/* y interval for 9/16 aspect ratio */

#define JULIA_SCALE 1.0 /* scaling for Julia sets */

/* Choice of the billiard table */

#define B_DOMAIN 20        /* choice of domain shape, see list in global_pdes.c */
#define B_DOMAIN_B 20      /* choice of domain shape, see list in global_pdes.c */

#define CIRCLE_PATTERN 1       /* pattern of circles, see list in global_pdes.c */
#define CIRCLE_PATTERN_B 10    /* pattern of circles, see list in global_pdes.c */

#define P_PERCOL 0.25       /* probability of having a circle in C_RAND_PERCOL arrangement */
#define NPOISSON 300        /* number of points for Poisson C_RAND_POISSON arrangement */

#define LAMBDA 0.8	    /* parameter controlling the dimensions of domain */
#define MU 0.03 	    /* parameter controlling the dimensions of domain */
#define MUB 0.03	    /* parameter controlling the dimensions of domain */
#define NPOLY 3             /* number of sides of polygon */
#define APOLY 1.0           /* angle by which to turn polygon, in units of Pi/2 */ 
#define MDEPTH 4            /* depth of computation of Menger gasket */
#define MRATIO 3            /* ratio defining Menger gasket */
#define MANDELLEVEL 1000      /* iteration level for Mandelbrot set */
#define MANDELLIMIT 10.0     /* limit value for approximation of Mandelbrot set */
#define FOCI 1              /* set to 1 to draw focal points of ellipse */
#define NGRIDX 15           /* number of grid point for grid of disks */
#define NGRIDY 20           /* number of grid point for grid of disks */

/* You can add more billiard tables by adapting the functions */
/* xy_in_billiard and draw_billiard below */

/* Physical parameters of wave equation */

#define TWOSPEEDS 0         /* set to 1 to replace hardcore boundary by medium with different speed */
#define OSCILLATE_LEFT 0     /* set to 1 to add oscilating boundary condition on the left */
#define OSCILLATE_TOPBOT 0   /* set to 1 to enforce a planar wave on top and bottom boundary */

#define OMEGA 0.0035        /* frequency of periodic excitation */
#define AMPLITUDE 1.0      /* amplitude of periodic excitation */ 
#define COURANT 0.02       /* Courant number */
#define COURANTB 0.0075      /* Courant number in medium B */
#define GAMMA 0.0      /* damping factor in wave equation */
#define GAMMAB 1.0e-5           /* damping factor in wave equation */
#define GAMMA_SIDES 1.0e-4      /* damping factor on boundary */
#define GAMMA_TOPBOT 1.0e-6      /* damping factor on boundary */
#define KAPPA 0.0       /* "elasticity" term enforcing oscillations */
#define KAPPA_SIDES 5.0e-4       /* "elasticity" term on absorbing boundary */
#define KAPPA_TOPBOT 0.0       /* "elasticity" term on absorbing boundary */
/* The Courant number is given by c*DT/DX, where DT is the time step and DX the lattice spacing */
/* The physical damping coefficient is given by GAMMA/(DT)^2 */
/* Increasing COURANT speeds up the simulation, but decreases accuracy */
/* For similar wave forms, COURANT^2*GAMMA should be kept constant */

/* Boundary conditions, see list in global_pdes.c  */

#define B_COND 3

/* Parameters for length and speed of simulation */

#define NSTEPS 6000      /* number of frames of movie */
#define NVID 25          /* number of iterations between images displayed on screen */
#define NSEG 100         /* number of segments of boundary */
#define INITIAL_TIME 200    /* time after which to start saving frames */
#define COMPUTE_ENERGIES 1  /* set to 1 to compute and print energies */
#define BOUNDARY_WIDTH 2    /* width of billiard boundary */

#define PAUSE 1000         /* number of frames after which to pause */
#define PSLEEP 1         /* sleep time during pause */
#define SLEEP1  1        /* initial sleeping time */
#define SLEEP2  1   /* final sleeping time */

/* Plot type, see list in global_pdes.c  */

#define PLOT 1

/* Color schemes */

#define BLACK 1          /* background */

#define COLOR_SCHEME 1   /* choice of color scheme, see list in global_pdes.c  */

#define SCALE 0          /* set to 1 to adjust color scheme to variance of field */
#define SLOPE 1.5        /* sensitivity of color on wave amplitude */
#define ATTENUATION 0.0  /* exponential attenuation coefficient of contrast with time */
#define E_SCALE 2000.0     /* scaling factor for energy representation */

#define COLORHUE 260     /* initial hue of water color for scheme C_LUM */
#define COLORDRIFT 0.0   /* how much the color hue drifts during the whole simulation */
#define LUMMEAN 0.5      /* amplitude of luminosity variation for scheme C_LUM */
#define LUMAMP 0.3       /* amplitude of luminosity variation for scheme C_LUM */
#define HUEMEAN 220.0    /* mean value of hue for color scheme C_HUE */
#define HUEAMP -220.0      /* amplitude of variation of hue for color scheme C_HUE */

/* For debugging purposes only */
#define FLOOR 0         /* set to 1 to limit wave amplitude to VMAX */
#define VMAX 10.0       /* max value of wave amplitude */


#include "global_pdes.c"        /* constants and global variables */
#include "sub_wave.c"           /* common functions for wave_billiard, heat and schrodinger */
#include "wave_common.c"        /* common functions for wave_billiard, wave_comparison, etc */
#include "sub_wave_comp.c"      /* some functions specific to wave_comparison */

double courant2, courantb2;  /* Courant parameters squared */


/*********************/
/* animation part    */
/*********************/

void evolve_wave_half(double *phi_in[NX], double *psi_in[NX], double *phi_out[NX], double *psi_out[NX], 
                      short int *xy_in[NX])
/* time step of field evolution */
/* phi is value of field at time t, psi at time t-1 */
{
    int i, j, iplus, iminus, jplus, jminus, jmid = NY/2;
    double delta, x, y, c, cc, gamma;
    static long time = 0;
    
    time++;

    #pragma omp parallel for private(i,j,iplus,iminus,jplus,jminus,delta,x,y)
    for (i=0; i<NX; i++){
        for (j=0; j<NY; j++){
            if (xy_in[i][j])
            {
                c = COURANT;
                cc = courant2;
                gamma = GAMMA;
            }
            else if (TWOSPEEDS)
            {
                c = COURANTB;
                cc = courantb2;
                gamma = GAMMAB;
            }

            if (((TWOSPEEDS)&&(xy_in[i][j] != 2))||(xy_in[i][j] == 1)){
                /* discretized Laplacian for various boundary conditions */
                if ((B_COND == BC_DIRICHLET)||(B_COND == BC_ABSORBING)||(B_COND == BC_ABS_REFLECT))
                {
                    iplus = (i+1);   if (iplus == NX) iplus = NX-1;
                    iminus = (i-1);  if (iminus == -1) iminus = 0;
                    jplus = (j+1);   
                    if (jplus == NY) jplus = NY-1;
                    else if (jplus == jmid) jplus = jmid-1;
                    jminus = (j-1);  
                    if (jminus == -1) jminus = 0;
                    else if (jminus == jmid-1) jminus = jmid;
                }
                else if (B_COND == BC_PERIODIC)
                {
                    iplus = (i+1) % NX;
                    iminus = (i-1) % NX;
                    if (iminus < 0) iminus += NX;
                    if (j < jmid)   /* lower half */
                    {
                        jplus = (j+1) % jmid;
                        jminus = (j-1) % jmid;
                        if (jminus < 0) jminus += jmid;
                    }
                    else    /* upper half */
                    {
                        jplus = j+1;
                        if (jplus >= NY) jplus -= jmid;
                        jminus = j-1;
                        if (jminus < jmid) jminus += jmid;
                    }
                }
                else if (B_COND == BC_VPER_HABS)
                {
                    iplus = (i+1);   if (iplus == NX) iplus = NX-1;
                    iminus = (i-1);  if (iminus == -1) iminus = 0;
                    if (j < jmid)   /* lower half */
                    {
                        jplus = (j+1) % jmid;
                        jminus = (j-1) % jmid;
                        if (jminus < 0) jminus += jmid;
                    }
                    else    /* upper half */
                    {
                        jplus = j+1;
                        if (jplus >= NY) jplus -= jmid;
                        jminus = j-1;
                        if (jminus < jmid) jminus += jmid;
                    }
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
                else if ((B_COND == BC_ABSORBING)||(B_COND == BC_ABS_REFLECT))
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
                
                /* add oscillating boundary condition on the left */
                if ((i == 0)&&(OSCILLATE_LEFT)) phi_out[i][j] = AMPLITUDE*cos((double)time*OMEGA);
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


void evolve_wave(double *phi[NX], double *psi[NX], double *phi_tmp[NX], double *psi_tmp[NX], short int *xy_in[NX])
/* time step of field evolution */
/* phi is value of field at time t, psi at time t-1 */
{
    evolve_wave_half(phi, psi, phi_tmp, psi_tmp, xy_in);
    evolve_wave_half(phi_tmp, psi_tmp, phi, psi, xy_in);
}



void animation()
{
    double time, scale, energies[6], top_energy, bottom_energy;
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
    
    /* initialise positions and radii of circles */
    printf("initializing circle configuration\n");
    if ((B_DOMAIN == D_CIRCLES)||(B_DOMAIN_B == D_CIRCLES)) init_circle_config_comp();

    courant2 = COURANT*COURANT;
    courantb2 = COURANTB*COURANTB;

    /* initialize wave with a drop at one point, zero elsewhere */
    int_planar_wave_comp(XMIN + 0.015, 0.0, phi, psi, xy_in);
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

    /* initialize energies */
    if (COMPUTE_ENERGIES) 
    {
        compute_energy_tblr(phi, psi, xy_in, energies);
        top_energy = energies[0] + energies[1] + energies[2]; 
        bottom_energy = energies[3] + energies[4] + energies[5];
    }

    blank();
    glColor3f(0.0, 0.0, 0.0);
    printf("drawing wave\n");
    draw_wave_comp(phi, psi, xy_in, 1.0, 0);

    printf("drawing billiard\n");
    draw_billiard_comp();

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

        draw_wave_comp(phi, psi, xy_in, scale, i);
        
        draw_billiard_comp();


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
        
        for (j=0; j<NVID; j++) 
        {
            evolve_wave(phi, psi, phi_tmp, psi_tmp, xy_in);
//             if (i % 10 == 9) oscillate_linear_wave(0.2*scale, 0.15*(double)(i*NVID + j), -1.5, YMIN, -1.5, YMAX, phi, psi);
        }
        
        glutSwapBuffers();

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

