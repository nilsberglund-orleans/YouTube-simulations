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

/* General geometrical parameters */

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

#define B_DOMAIN 20      /* choice of domain shape, see list in global_pdes.c */

#define CIRCLE_PATTERN 8   /* pattern of circles, see list in global_pdes.c */

#define P_PERCOL 0.25       /* probability of having a circle in C_RAND_PERCOL arrangement */
#define NPOISSON 340        /* number of points for Poisson C_RAND_POISSON arrangement */

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

/* You can add more billiard tables by adapting the functions */
/* xy_in_billiard and draw_billiard below */

/* Physical parameters of wave equation */

#define TWOSPEEDS 1         /* set to 1 to replace hardcore boundary by medium with different speed */
#define OSCILLATE_LEFT 1    /* set to 1 to add oscilating boundary condition on the left */
#define OSCILLATE_TOPBOT 1  /* set to 1 to enforce a planar wave on top and bottom boundary */
#define X_SHIFT -0.9        /* x range on which to apply OSCILLATE_TOPBOT */

#define OMEGA 0.002        /* frequency of periodic excitation */
#define K_BC 3.0           /* spatial period of periodic excitation in y direction */
#define KX_BC 30.0           /* spatial period of periodic excitation in x direction */
#define KY_BC 10.0           /* spatial period of periodic excitation in y direction */
#define AMPLITUDE 1.0      /* amplitude of periodic excitation */ 
#define COURANT 0.02       /* Courant number */
#define COURANTB 0.01      /* Courant number in medium B */
// #define COURANTB 0.00666   /* Courant number in medium B */
#define GAMMA 2.0e-6          /* damping factor in wave equation */
#define GAMMAB 2.5e-4           /* damping factor in wave equation */
// #define GAMMAB 5.0e-4           /* damping factor in wave equation */
// #define GAMMAB 1.0e-4           /* damping factor in wave equation */
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

// #define NSTEPS 1000      /* number of frames of movie */
#define NSTEPS 4500      /* number of frames of movie */
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

#define BLACK 1          /* background */

#define COLOR_SCHEME 1   /* choice of color scheme, see list in global_pdes.c  */

#define SCALE 0          /* set to 1 to adjust color scheme to variance of field */
#define SLOPE 1.0        /* sensitivity of color on wave amplitude */
#define ATTENUATION 0.0  /* exponential attenuation coefficient of contrast with time */
#define E_SCALE 2500.0     /* scaling factor for energy representation */

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
#define MANGROVE_EMAX 1.1e-3           /* max energy for mangrove to survive */

#define RANDOM_RADIUS 1     /* set to 1 for random circle radius */
#define ERODE_MANGROVES 0   /* set to 1 for mangroves to be eroded */
#define MOVE_MANGROVES 1    /* set to 1 for mobile mangroves */
#define DETACH_MANGROVES 1  /* set to 1 for mangroves to be able to detach */
#define INERTIA 1           /* set to 1 for taking inertia into account */
#define DT_MANGROVE 0.1     /* time step for mangrove displacement */
#define KSPRING 0.25        /* spring constant of mangroves */
#define KWAVE 2.0           /* constant in force due to wave gradient */
#define DXMAX 0.02          /* max displacement of mangrove in one time step */
#define L_DETACH 0.2        /* spring length beyond which mangroves detach */ 
#define DAMP_MANGROVE 0.2   /* damping coefficient of mangroves */
#define MANGROVE_MASS 1.5   /* mass of mangrove of radius MU */

/* For debugging purposes only */
#define FLOOR 1         /* set to 1 to limit wave amplitude to VMAX */
#define VMAX 10.0       /* max value of wave amplitude */


#include "global_pdes.c"
#include "sub_wave.c"
#include "wave_common.c"

double courant2, courantb2;  /* Courant parameters squared */
double circle_energy[NMAXCIRCLES];      /* energy dissipated by the circles */
double circley_wrapped[NMAXCIRCLES];  /* position of circle centers wrapped vertically */
double anchor_x[NMAXCIRCLES];         /* points moving circles are attached to */
double anchor_y[NMAXCIRCLES];         /* points moving circles are attached to */
double vx[NMAXCIRCLES];               /* x velocity of circles */
double vy[NMAXCIRCLES];               /* y velocity of circles */
double circlerad_initial[NMAXCIRCLES];  /* initial circle radii */
double mass_inverse[NMAXCIRCLES];       /* inverse of mangrove mass */
short int circle_attached[NMAXCIRCLES]; /* has value 1 if the circle is attached to its anchor */

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


void evolve_wave_half(double *phi_in[NX], double *psi_in[NX], double *phi_out[NX], double *psi_out[NX], 
                      short int *xy_in[NX])
// void evolve_wave_half(phi_in, psi_in, phi_out, psi_out, xy_in)
/* time step of field evolution */
/* phi is value of field at time t, psi at time t-1 */
{
    int i, j, iplus, iminus, jplus, jminus, tb_shift;
    double delta, x, y, c, cc, gamma, kappa, phase, phasemin;
    static long time = 0;
    static int init_bc = 1;
    static double left_bc[NY], top_bc[NX], bot_bc[NX];
    
    time++;
    
    /* initialize boundary condition phase KX*x + KY*y */
    if ((OSCILLATE_LEFT)&&(init_bc)) 
    {
        init_bc_phase(left_bc, top_bc, bot_bc);
        tb_shift = (int)((X_SHIFT - XMIN)*(double)NX/(XMAX - XMIN));
        printf("tb_shift %i\n", tb_shift);
        init_bc = 0;
    }
    
    #pragma omp parallel for private(i,j,iplus,iminus,jplus,jminus,delta,x,y,c,cc,gamma,kappa)
    for (i=0; i<NX; i++){
        for (j=0; j<NY; j++){
            if (xy_in[i][j])
            {
                c = COURANT;
                cc = courant2;
                gamma = GAMMA;
                kappa = KAPPA;
            }
            else if (TWOSPEEDS)
            {
                c = COURANTB;
                cc = courantb2;
                gamma = GAMMAB;
                kappa = KAPPAB;
            }

            if ((TWOSPEEDS)||(xy_in[i][j])){
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
                if ((OSCILLATE_TOPBOT)&&(i < tb_shift))
                {
                    if (j == NY-1) 
                    {
                        jminus = NY-1;
                        jplus = NY-1;
                    }
                    else if (j == 0) 
                    {   
                        jminus = 0;
                        jplus = 0;
                    }
                }
                
                delta = phi_in[iplus][j] + phi_in[iminus][j] + phi_in[i][jplus] + phi_in[i][jminus] - 4.0*phi_in[i][j];

                x = phi_in[i][j];
		y = psi_in[i][j];

                /* evolve phi */
                if ((B_COND == BC_PERIODIC)||(B_COND == BC_DIRICHLET)) 
                    phi_out[i][j] = -y + 2*x + cc*delta - kappa*x - gamma*(x-y);
                else if (B_COND == BC_ABSORBING)
                {
                    if ((i>0)&&(i<NX-1)&&(j>0)&&(j<NY-1))
                        phi_out[i][j] = -y + 2*x + cc*delta - kappa*x - gamma*(x-y);
                
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
                        phi_out[i][j] = -y + 2*x + cc*delta - kappa*x - gamma*(x-y);
                
                    /* right border */
                    else if (i==NX-1) 
                        phi_out[i][j] = x - c*(x - phi_in[NX-2][j]) - KAPPA_SIDES*x - GAMMA_SIDES*(x-y);
                    
                    /* left border */
                    else if (i==0) 
                        phi_out[i][j] = x - c*(x - phi_in[1][j]) - KAPPA_SIDES*x - GAMMA_SIDES*(x-y);
                }
                psi_out[i][j] = x;
                
                /* add oscillating boundary condition on the left */
//                 if ((i == 0)&&(OSCILLATE_LEFT)) 
//                 {
//                     phase =  (double)time*OMEGA - DPI*K_BC*(double)j/(double)NY;
//                     if (phase < 0.0) phase = 0.0;
//                     phi_out[i][j] = AMPLITUDE*sin(phase);
//                 }

                /* add oscillating boundary condition on the left */
                if (OSCILLATE_LEFT)
                {
                    phasemin = left_bc[0]; 
                    if (i == 0)
                    {
                        phase =  (double)time*OMEGA - left_bc[j] + phasemin;
                        if (phase < 0.0) phase = 0.0;
                        phi_out[i][j] = AMPLITUDE*sin(phase);
                    }
                    if ((j == 0)&&(i < tb_shift))
                    {
                        phase =  (double)time*OMEGA - bot_bc[i] + phasemin;
                        if (phase < 0.0) phase = 0.0;
                        phi_out[i][j] = AMPLITUDE*sin(phase);
                    }
                    else if ((j == NY-1)&&(i < tb_shift))
                    {
                        phase =  (double)time*OMEGA - top_bc[i] + phasemin;
                        if (phase < 0.0) phase = 0.0;
                        phi_out[i][j] = AMPLITUDE*sin(phase);
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


void evolve_wave(double *phi[NX], double *psi[NX], double *phi_tmp[NX], double *psi_tmp[NX], short int *xy_in[NX])
/* time step of field evolution */
/* phi is value of field at time t, psi at time t-1 */
{
    evolve_wave_half(phi, psi, phi_tmp, psi_tmp, xy_in);
    evolve_wave_half(phi_tmp, psi_tmp, phi, psi, xy_in);
}


void animation()
{
    double time, scale, diss, rgb[3], hue, y, dissip, ej, gradient[2], dx, dy, dt, xleft, xright, length;
    double *phi[NX], *psi[NX], *phi_tmp[NX], *psi_tmp[NX];
    short int *xy_in[NX], redraw = 0;
    int i, j, s, ij[2];
    static int imin, imax;
    static short int first = 1;

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
    if (B_DOMAIN == D_CIRCLES) init_circle_config();

    courant2 = COURANT*COURANT;
    courantb2 = COURANTB*COURANTB;
//     dt = 0.01;

    /* initialize wave with a drop at one point, zero elsewhere */
    init_wave_flat(phi, psi, xy_in);
    
//     init_planar_wave(XMIN + 0.01, 0.0, phi, psi, xy_in);
//     init_planar_wave(XMIN + 0.02, 0.0, phi, psi, xy_in);
//     init_planar_wave(XMIN + 1.0, 0.0, phi, psi, xy_in);
//     init_wave(-1.5, 0.0, phi, psi, xy_in);
//     init_wave(0.0, 0.0, phi, psi, xy_in);

    /* add a drop at another point */
//     add_drop_to_wave(1.0, 0.7, 0.0, phi, psi);
//     add_drop_to_wave(1.0, -0.7, 0.0, phi, psi);
//     add_drop_to_wave(1.0, 0.0, -0.7, phi, psi);

    /* initialise mangroves */
    for (i=0; i < ncircles; i++) 
    {
        circle_energy[i] = 0.0;
        y = circley[i];
        if (y >= YMAX) y -= circlerad[i];
        if (y <= YMIN) y += circlerad[i];
//         if (y >= YMAX) y -= (YMAX - YMIN);
//         if (y <= YMIN) y += (YMAX - YMIN);
        circley_wrapped[i] = y;
        
        if (RANDOM_RADIUS) circlerad[i] = circlerad[i]*(0.75 + 0.5*((double)rand()/RAND_MAX));
        
        circlerad_initial[i] = circlerad[i];
        circle_attached[i] = 1; 
        mass_inverse[i] = MU*MU/(MANGROVE_MASS*circlerad[i]*circlerad[i]);
        
        if (MOVE_MANGROVES)
        {
            anchor_x[i] = circlex[i];
            anchor_y[i] = circley_wrapped[i];
//             anchor_y[i] = circley[i];
        }
        
        if (INERTIA)
        {
            vx[i] = 0.0;
            vy[i] = 0.0;
        }
    }
    
    if (first) /* compute box limits where circles are reset */
    {
        /* find leftmost and rightmost circle */
        for (i=0; i<ncircles; i++) 
            if ((circleactive[i])&&(circlex[i] - circlerad[i] < xleft)) xleft = circlex[i] - circlerad[i]; 
        for (i=0; i<ncircles; i++) 
            if ((circleactive[i])&&(circlex[i] + circlerad[i] > xright)) xright = circlex[i] + circlerad[i]; 
        
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
    draw_billiard();

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

        draw_wave(phi, psi, xy_in, scale, i, PLOT);
        for (j=0; j<NVID; j++) 
        {
            evolve_wave(phi, psi, phi_tmp, psi_tmp, xy_in);
//             if (i % 10 == 9) oscillate_linear_wave(0.2*scale, 0.15*(double)(i*NVID + j), -1.5, YMIN, -1.5, YMAX, phi, psi);
        }
        
        /* compute energy dissipated in obstacles */
        if (ERODE_MANGROVES) for (j=0; j<ncircles; j++) 
        {
            dissip = compute_dissipation(phi, psi, xy_in, circlex[j], circley_wrapped[j]);
            
            /* make sure the dissipation does not grow too fast because of round-off/blow-up */
            if (dissip > 0.1*MANGROVE_EMAX) 
            {
                dissip = 0.1*MANGROVE_EMAX;
                printf("Flooring dissipation!\n");
            }
            
            if (circleactive[j]) 
            {
                circle_energy[j] += dissip;
                ej = circle_energy[j];
                if (ej <= MANGROVE_EMAX)
                {
                    if (ej > 0.0) 
                    {
                        hue = MANGROVE_HUE_MIN + (MANGROVE_HUE_MAX - MANGROVE_HUE_MIN)*ej/MANGROVE_EMAX;
                        if (hue < 0.0) hue += 360.0;
                    }
                    else hue = MANGROVE_HUE_MIN;
                    hsl_to_rgb(hue, 0.9, 0.5, rgb);
                    if (j%NGRIDY == 0) printf("Circle %i, energy %.5lg, hue %.5lg\n", j, ej, hue);
                    draw_colored_circle(circlex[j], circley[j], circlerad[j], NSEG, rgb);
                    
                    /* shrink mangrove */
                    if (ej > 0.0)
                    {
//                         circlerad[j] -= MU*ej*ej/(MANGROVE_EMAX*MANGROVE_EMAX);
//                         if (circlerad[j] < 0.0) circlerad[j] = 0.0;
                        circlerad[j] = circlerad_initial[j]*(1.0 - ej*ej/(MANGROVE_EMAX*MANGROVE_EMAX));
                        redraw = 1;
                    }
                    else circlerad[j] = circlerad_initial[j];
                }
                else    /* remove mangrove */
                {
                    circleactive[j] = 0;
                    /* reinitialize table xy_in */
                    redraw = 1;
                }   
            }
            else    /* allow disabled mangroves to recover */
            {
                circle_energy[j] -= 0.15*dissip;
//                 circlerad[j] += 0.005*MU;
//                 if (circlerad[j] > MU) circlerad[j] = MU;
//                 if ((circle_energy[j] < 0.0)&&(circlerad[j] > 0.0))
                if (circle_energy[j] < 0.0)
                {
                    circleactive[j] = 1;
//                     circlerad[j] = circlerad[j]*(0.75 + 0.5*((double)rand()/RAND_MAX));
                    circlerad[j] = circlerad_initial[j];
                    circle_energy[j] = -MANGROVE_EMAX;
                    /* reinitialize table xy_in */
                    redraw = 1;
                }   
                
            }
            
//             printf("Circle %i, energy %.5lg\n", j, circle_energy[j]);
        }
        
        /* move mangroves */
        if (MOVE_MANGROVES) for (j=0; j<ncircles; j++) if (circleactive[j])
        {
            compute_gradient(phi, psi, circlex[j], circley_wrapped[j], gradient);
            
//             if (j%NGRIDY == 0) printf("gradient (%.3lg, %.3lg)\n", gradient[0], gradient[1]);
//             if (j%NGRIDY == 0) printf("circle %i (%.3lg, %.3lg) -> ", j, circlex[j], circley[j]);
            
            /* compute force of wave */
            dx = DT_MANGROVE*KWAVE*gradient[0];
            dy = DT_MANGROVE*KWAVE*gradient[1];

            /* compute force of spring */
            if (circle_attached[j])
            {
                dx += DT_MANGROVE*(-KSPRING*(circlex[j] - anchor_x[j]));
                dy += DT_MANGROVE*(-KSPRING*(circley_wrapped[j] - anchor_y[j]));
            }
            
            /* detach mangrove if spring is too long */
            if (DETACH_MANGROVES)
            {
                length = module2(circlex[j] - anchor_x[j], circley_wrapped[j] - anchor_y[j]);
                if (j%NGRIDY == 0) printf("spring length %.i:  %.3lg\n", j, length);
//                 if (length > L_DETACH) circle_attached[j] = 0;
                if (length*mass_inverse[j] > L_DETACH) circle_attached[j] = 0;
            }
            
            if (dx > DXMAX) dx = DXMAX;
            if (dx < -DXMAX) dx = -DXMAX;
            if (dy > DXMAX) dy = DXMAX;
            if (dy < -DXMAX) dy = -DXMAX;
            
            if (INERTIA)
            {
                vx[j] += (dx - DAMP_MANGROVE*vx[j])*mass_inverse[j];
                vy[j] += (dy - DAMP_MANGROVE*vy[j])*mass_inverse[j];
                circlex[j] += vx[j]*DT_MANGROVE;
                circley[j] += vy[j]*DT_MANGROVE;
                circley_wrapped[j] += vy[j]*DT_MANGROVE;
                if (j%NGRIDY == 0) 
                    printf("circle %.i: (dx,dy) = (%.3lg,%.3lg), (vx,vy) = (%.3lg,%.3lg)\n", 
                           j, circlex[j]-anchor_x[j], circley[j]-anchor_y[j], vx[j], vy[j]);
            }
            else
            {
                circlex[j] += dx*mass_inverse[j]*DT_MANGROVE;
                circley[j] += dy*mass_inverse[j]*DT_MANGROVE;
                circley_wrapped[j] += dy*mass_inverse[j]*DT_MANGROVE;
            }
            
            if (circlex[j] <= XMIN) circlex[j] = XMIN;
            if (circlex[j] >= XMAX) circlex[j] = XMAX;
            if (circley_wrapped[j] <= YMIN) circley_wrapped[j] = YMIN;
            if (circley_wrapped[j] >= YMAX) circley_wrapped[j] = YMAX;
        
//             if (j%NGRIDY == 0) printf("(%.3lg, %.3lg)\n", circlex[j], circley[j]);

            redraw = 1;
        }
        
        draw_billiard();

	glutSwapBuffers();

        if (redraw) 
        {
            printf("Reinitializing xy_in\n");
            init_xyin_xrange(xy_in, imin, NX-1);
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

