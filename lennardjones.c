/*********************************************************************************/
/*                                                                               */
/*  Animation of interacting particles in a planar domain                        */
/*                                                                               */
/*  N. Berglund, november 2021                                                   */
/*                                                                               */
/*  UPDATE 24/04: distinction between damping and "elasticity" parameters        */
/*  UPDATE 27/04: new billiard shapes, bug in color scheme fixed                 */
/*  UPDATE 28/04: code made more efficient, with help of Marco Mancini           */
/*                                                                               */
/*  Feel free to reuse, but if doing so it would be nice to drop a               */
/*  line to nils.berglund@univ-orleans.fr - Thanks!                              */
/*                                                                               */
/*  compile with                                                                 */
/*  gcc -o lennardjones lennardjones.c                                           */
/* -L/usr/X11R6/lib -ltiff -lm -lGL -lGLU -lX11 -lXmu -lglut -O3 -fopenmp        */
/*                                                                               */
/*  OMP acceleration may be more effective after executing                       */
/*  export OMP_NUM_THREADS=2 in the shell before running the program             */
/*                                                                               */
/*  To make a video, set MOVIE to 1 and create subfolder tif_ljones              */
/*  It may be possible to increase parameter PAUSE                               */
/*                                                                               */
/*  create movie using                                                           */
/*  ffmpeg -i lj.%05d.tif -vcodec libx264 lj.mp4                             */
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
#define DOUBLE_MOVIE 0  /* set to 1 to produce movies for wave height and energy simultaneously */

/* General geometrical parameters */

#define WINWIDTH 	1280  /* window width */
#define WINHEIGHT 	720   /* window height */

#define XMIN -2.0
#define XMAX 2.0	/* x interval */
#define YMIN -1.125
#define YMAX 1.125	/* y interval for 9/16 aspect ratio */

#define INITXMIN -2.0
#define INITXMAX 2.0	/* x interval for initial condition */
#define INITYMIN -0.8
#define INITYMAX 1.125	/* y interval for initial condition */

/* Choice of the billiard table */

// #define B_DOMAIN 20      /* choice of domain shape, see list in global_ljones.c */

#define CIRCLE_PATTERN 8   /* pattern of circles, see list in global_ljones.c */

#define INTERACTION 1       /* particle interaction, see list in global_ljones.c */

#define P_PERCOL 0.25       /* probability of having a circle in C_RAND_PERCOL arrangement */
#define NPOISSON 100        /* number of points for Poisson C_RAND_POISSON arrangement */
#define PDISC_DISTANCE 5.0  /* minimal distance in Poisson disc process, controls density of particles */
#define PDISC_CANDIDATES 100 /* number of candidates in construction of Poisson disc process */
#define RANDOM_POLY_ANGLE 0 /* set to 1 to randomize angle of polygons */

#define LAMBDA 2.0	    /* parameter controlling the dimensions of domain */
// #define MU 0.02 	    /* parameter controlling radius of particles */
#define MU 0.015 	    /* parameter controlling radius of particles */
#define NPOLY 3             /* number of sides of polygon */
#define APOLY 1.0           /* angle by which to turn polygon, in units of Pi/2 */ 
#define MDEPTH 4            /* depth of computation of Menger gasket */
#define MRATIO 3            /* ratio defining Menger gasket */
#define MANDELLEVEL 1000    /* iteration level for Mandelbrot set */
#define MANDELLIMIT 10.0    /* limit value for approximation of Mandelbrot set */
#define FOCI 1              /* set to 1 to draw focal points of ellipse */
// #define NGRIDX 50           /* number of grid point for grid of disks */
#define NGRIDX 32           /* number of grid point for grid of disks */
#define NGRIDY 26           /* number of grid point for grid of disks */

#define X_SHOOTER -0.2
#define Y_SHOOTER -0.6
#define X_TARGET 0.4
#define Y_TARGET 0.7        /* shooter and target positions in laser fight */

/* Boundary conditions, see list in global_ljones.c  */

#define B_COND 3

/* Parameters for length and speed of simulation */

#define NSTEPS 3600      /* number of frames of movie */
// #define NSTEPS 100      /* number of frames of movie */
#define NVID 250          /* number of iterations between images displayed on screen */
#define NSEG 100         /* number of segments of boundary */
#define INITIAL_TIME 0    /* time after which to start saving frames */
#define BOUNDARY_WIDTH 2    /* width of billiard boundary */

#define PAUSE 1000         /* number of frames after which to pause */
#define PSLEEP 1         /* sleep time during pause */
#define SLEEP1  1        /* initial sleeping time */
#define SLEEP2  1   /* final sleeping time */
#define MID_FRAMES 20    /* number of still frames between parts of two-part movie */
#define END_FRAMES 100   /* number of still frames at end of movie */

/* Parameters of initial condition */


/* Plot type, see list in global_ljones.c  */

#define PLOT 1
#define PLOT_B 3        /* plot type for second movie */


/* Color schemes */

#define COLOR_PALETTE 10     /* Color palette, see list in global_ljones.c  */

#define BLACK 1          /* background */

#define COLOR_SCHEME 1   /* choice of color scheme, see list in global_ljones.c  */

#define SCALE 0          /* set to 1 to adjust color scheme to variance of field */
#define SLOPE 0.5        /* sensitivity of color on wave amplitude */
#define ATTENUATION 0.0  /* exponential attenuation coefficient of contrast with time */

#define COLORHUE 260     /* initial hue of water color for scheme C_LUM */
#define COLORDRIFT 0.0   /* how much the color hue drifts during the whole simulation */
#define LUMMEAN 0.5      /* amplitude of luminosity variation for scheme C_LUM */
#define LUMAMP 0.3       /* amplitude of luminosity variation for scheme C_LUM */
#define HUEMEAN 220.0    /* mean value of hue for color scheme C_HUE */
#define HUEAMP -50.0      /* amplitude of variation of hue for color scheme C_HUE */

/* particle properties */

#define PARTICLE_HUE_MIN 330.0      /* color of original particle */
#define PARTICLE_HUE_MAX 30.0        /* color of saturated particle */
#define PARTICLE_EMAX 2.0e1           /* max energy for particle to survive */

#define RANDOM_RADIUS 0     /* set to 1 for random circle radius */
#define MOVE_PARTICLES 1    /* set to 1 for mobile particles */
#define INERTIA 1           /* set to 1 for taking inertia into account */
#define DT_PARTICLE 2.0e-6     /* time step for particle displacement */
#define KREPEL 10.0          /* constant in repelling force between particles */
#define EQUILIBRIUM_DIST 4.0  /* Lennard-Jones equilibrium distance */
// #define EQUILIBRIUM_DIST 15.0  /* Lennard-Jones equilibrium distance */
#define REPEL_RADIUS 20.0    /* radius in which repelling force acts (in units of particle radius) */
#define DAMPING 1.5e5   /* damping coefficient of particles */
// #define DAMPING 1.0e-10      /* damping coefficient of particles */
#define PARTICLE_MASS 1.0   /* mass of particle of radius MU */
// #define V_INITIAL 0.0       /* initial velocity range */
#define V_INITIAL 5.0        /* initial velocity range */
#define SIGMA 5.0           /* noise intensity in thermostat */
#define BETA 1.0e-2          /* initial inverse temperature */
#define MU_XI 0.1            /* friction constant in thermostat */
#define KSPRING_BOUNDARY 1.0e5    /* confining harmonic potential outside simulation region */
#define NBH_DIST_FACTOR 4.5        /* radius in which to count neighbours */
#define GRAVITY 300.0         /* gravity acting on all particles */

#define INCREASE_BETA 1   /* set to 1 to increase BETA during simulation */
#define BETA_FACTOR 1.0e1   /* factor by which to change BETA during simulation */
#define N_TOSCILLATIONS 0.25   /* number of temperature oscillations in BETA schedule */
// #define BETA_FACTOR 2.0e3   /* factor by which to change BETA during simulation */

#define INCREASE_KREPEL 0   /* set to 1 to increase KREPEL during simulation */
#define KREPEL_FACTOR 1000.0   /* factor by which to change KREPEL during simulation */

#define PART_AT_BOTTOM 1     /* set to 1 to include "seed" particles at bottom */
#define MASS_PART_BOTTOM 10000.0 /* mass of particles at bottom */
#define NPART_BOTTOM 100     /* number of particles at the bottom */

#define ADD_PARTICLES 1    /* set to 1 to add particles */
#define ADD_TIME 50       /* time at which to add first particle */
#define ADD_PERIOD 7      /* time interval between adding further particles */
#define SAFETY_FACTOR 3.0  /* no particles are added at distance less than MU*SAFETY_FACTOR of other particles */

#define FLOOR_FORCE 0      /* set to 1 to limit force on particle to FMAX */
#define FMAX 2.0e10         /* maximal force */

#define HASHX 32    /* size of hashgrid in x direction */
#define HASHY 18   /* size of hashgrid in y direction */
#define HASHMAX 100  /* maximal number of particles per hashgrid cell */
#define HASHGRID_PADDING 0.1    /* padding of hashgrid outside simulation window */

#define DRAW_COLOR_SCHEME 0     /* set to 1 to plot the color scheme */
#define COLORBAR_RANGE 8.0    /* scale of color scheme bar */
#define COLORBAR_RANGE_B 12.0    /* scale of color scheme bar for 2nd part */
#define ROTATE_COLOR_SCHEME 0   /* set to 1 to draw color scheme horizontally */

/* For debugging purposes only */
#define FLOOR 1         /* set to 1 to limit wave amplitude to VMAX */
#define VMAX 10.0       /* max value of wave amplitude */

#include "global_ljones.c"
#include "sub_lj.c"


double gaussian()
/* returns standard normal random variable, using Box-Mueller algorithm */
{
    static double V1, V2, S;
    static int phase = 0;
    double X;

    if (phase == 0) 
    {
        do 
        {
        double U1 = (double)rand() / RAND_MAX;
        double U2 = (double)rand() / RAND_MAX;
        V1 = 2 * U1 - 1;
        V2 = 2 * U2 - 1;
        S = V1 * V1 + V2 * V2;
        } 
        while(S >= 1 || S == 0);
        X = V1 * sqrt(-2 * log(S) / S);
    } 
    else X = V2 * sqrt(-2 * log(S) / S);

    phase = 1 - phase;

    return X;
}


/*********************/
/* animation part    */
/*********************/



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

double lennard_jones_force_aniso(double r, double req)
{
    int i;
    double rmin = 0.01, rplus, ratio = 1.0;
    
    if (r > REPEL_RADIUS*MU) return(0.0);
    else
    {
        if (r > rmin) rplus = r;
        else rplus = rmin;
    
        for (i=0; i<6; i++) ratio *= req*MU/rplus;
    
        return((ratio - 2.0*ratio*ratio)/rplus);
    }
}

double lennard_jones_force(double r)
{
    int i;
    double rmin = 0.01, rplus, ratio = 1.0;
    
    if (r > REPEL_RADIUS*MU) return(0.0);
    else
    {
        if (r > rmin) rplus = r;
        else rplus = rmin;
    
//     ratio = pow(EQUILIBRIUM_DIST*MU/rplus, 6.0);
        for (i=0; i<6; i++) ratio *= EQUILIBRIUM_DIST*MU/rplus;
    
        return((ratio - 2.0*ratio*ratio)/rplus);
    }
}

void aniso_lj_force(double r, double ca, double sa, double force[2])
{
    int i;
    double rmin = 0.01, rplus, ratio = 1.0, c2, s2, c4, s4, a, aprime, f1, f2;
    
    if (r > REPEL_RADIUS*MU)
    {
        force[0] = 0.0;
        force[1] = 0.0;
    }
    else
    {
        if (r > rmin) rplus = r;
        else rplus = rmin;
    
        for (i=0; i<6; i++) ratio *= EQUILIBRIUM_DIST*MU/rplus;
    
        /* cos(2phi) and sin(2phi) */
        c2 = ca*ca - sa*sa;
        s2 = 2.0*ca*sa;

        /* cos(4phi) and sin(4phi) */
        c4 = c2*c2 - s2*s2;
        s4 = 2.0*c2*s2;
        
        a = 0.5*(9.0 - 7.0*c4);
        aprime = 14.0*s4;
        
        f1 = ratio*(a - ratio)/rplus;
        f2 = ratio*aprime/rplus;
        
        force[0] = f1*ca - f2*sa;
        force[1] = f1*sa + f2*ca;
    }
}

void penta_lj_force(double r, double ca, double sa, double force[2])
{
    int i;
    double rmin = 0.01, rplus, ratio = 1.0, c2, s2, c4, s4, c5, s5, a, aprime, f1, f2;
    static double a0, b0;
    static int first = 1;
    
    if (first)
    {
        a0 = cos(0.1*PI) + 0.5;
        b0 = a0 - 1.0;
        first = 0;
    }
    
    if (r > REPEL_RADIUS*MU)
    {
        force[0] = 0.0;
        force[1] = 0.0;
    }
    else
    {
        if (r > rmin) rplus = r;
        else rplus = rmin;
    
        for (i=0; i<6; i++) ratio *= EQUILIBRIUM_DIST*MU/rplus;
    
        /* cos(2phi) and sin(2phi) */
        c2 = ca*ca - sa*sa;
        s2 = 2.0*ca*sa;

        /* cos(4phi) and sin(4phi) */
        c4 = c2*c2 - s2*s2;
        s4 = 2.0*c2*s2;
        
        /* cos(5phi) and sin(5phi) */
        c5 = ca*c4 - sa*s4;
        s5 = sa*c4 + ca*s4;
        
        a = a0 - b0*c5;
        aprime = 5.0*b0*s5;
        
        f1 = ratio*(a - ratio)/rplus;
        f2 = ratio*aprime/rplus;
        
        force[0] = f1*ca - f2*sa;
        force[1] = f1*sa + f2*ca;
    }
}


int compute_repelling_force(int i, int j, double force[2], t_particle* particle, double krepel)
/* compute repelling force of particle j on particle i */
/* returns 1 if distance between particles is smaller than NBH_DIST_FACTOR*MU */
{
    double x1, y1, x2, y2, distance, r, f, angle, ca, sa, aniso, fx, fy, ff[2];
    
    x1 = particle[i].xc;
    y1 = particle[i].yc;
    x2 = particle[j].xc;
    y2 = particle[j].yc;
    
    distance = module2(x2 - x1, y2 - y1);
    
    if (distance == 0.0)
    {
        force[0] = 0.0;
        force[1] = 0.0;
        return(1);
    }
    else
    {
        ca = (x2 - x1)/distance;
        sa = (y2 - y1)/distance;
    
        switch (INTERACTION) {
            case (I_COULOMB):
            {
                f = krepel/(1.0e-8 + distance*distance);
                force[0] = f*ca;
                force[1] = f*sa;   
                break;
            }
            case (I_LENNARD_JONES):
            {
                f = krepel*lennard_jones_force(distance);
                force[0] = f*ca;
                force[1] = f*sa;   
                break;
            }
            case (I_LJ_DIRECTIONAL):
            {
                aniso_lj_force(distance, ca, sa, ff);
                force[0] = krepel*ff[0];
                force[1] = krepel*ff[1];
                break;
            }
            case (I_LJ_PENTA):
            {
                penta_lj_force(distance, ca, sa, ff);
                force[0] = krepel*ff[0];
                force[1] = krepel*ff[1];
                break;
            }
        }
    }

    if ((distance < NBH_DIST_FACTOR*MU)&&(j != i)) return(1);
    else return(0);
}


void update_hashgrid(t_particle* particle, int* hashgrid_number, int* hashgrid_particles)
{
    int i, j, k, n, m, ij[2], max = 0;
    
//     printf("Updating hashgrid_number\n");
    for (i=0; i<HASHX*HASHY; i++) hashgrid_number[i] = 0;
//     printf("Updated hashgrid_number\n");
        
    /* place each particle in hash grid */
    for (k=1; k<ncircles; k++)
//         if (circleactive[k])
        {
//             printf("placing circle %i\t", k);
            hash_xy_to_ij(particle[k].xc, particle[k].yc, ij);
            i = ij[0];  j = ij[1];
//             printf("ij = (%i, %i)\t", i, j);
            n = hashgrid_number[i*HASHY + j];
            m = i*HASHY*HASHMAX + j*HASHMAX + n;
//             printf("n = %i, m = %i\n", n, m);
            if (m < HASHX*HASHY*HASHMAX) hashgrid_particles[m] = k;
            else printf("Too many particles in hash cell, try increasing HASHMAX\n");
            hashgrid_number[i*HASHY + j]++;
            particle[k].hashx = i;
            particle[k].hashy = j;
            
            if (n > max) max = n;
//                 printf("Placed particle %i at (%i,%i) in hashgrid\n", k, ij[0], ij[1]);
//                 printf("%i particles at (%i,%i)\n", hashgrid_number[ij[0]][ij[1]], ij[0], ij[1]);
        }
    
    printf("Maximal number of particles per hash cell: %i\n", max);
}

int add_particle(double x, double y, double vx, double vy, t_particle particle[NMAXCIRCLES])
{
    int i, closeby = 0;
    double dist;
    
    /* test distance to other particles */
    for (i=0; i<ncircles; i++)
    {
        dist = module2(x - particle[i].xc, y - particle[i].yc);
        if (dist < SAFETY_FACTOR*MU) closeby = 1;
    }
    
    if ((closeby)||(ncircles >= NMAXCIRCLES)) 
    {
        printf("Cannot add particle at (%.3lg, %.3lg)\n", x, y);
        return(0);
    }
    else
    {
        i = ncircles;
        
        particle[i].xc = x;
        particle[i].yc = y;
        particle[i].radius = MU;
        particle[i].active = 1;
        particle[i].neighb = 0;

        particle[i].energy = 0.0;

        if (RANDOM_RADIUS) particle[i].radius = particle[i].radius*(0.75 + 0.5*((double)rand()/RAND_MAX));
        
        particle[i].mass_inv = 1.0;
                
        particle[i].vx = vx;
        particle[i].vy = vy;
        particle[i].energy = (particle[i].vx*particle[i].vx + particle[i].vy*particle[i].vy)*particle[i].mass_inv;
    
        ncircles++;
        
        printf("Added particle at (%.3lg, %.3lg)\n", x, y);
        printf("Number of particles: %i\n", ncircles);

        return(1);
    }
}

double neighbour_color(int nnbg)
{
    if (nnbg > 7) nnbg = 7;
    switch(nnbg){
        case (7): return(340.0);
        case (6): return(310.0);
        case (5): return(260.0);
        case (4): return(200.0);
        case (3): return(140.0);
        case (2): return(100.0);
        case (1): return(70.0);
        default:  return(30.0);
    }   
}

void draw_particles(t_particle particle[NMAXCIRCLES], int plot)
{
    int j, k, m, width, nnbg;
    double ej, hue, rgb[3], radius, x1, y1, x2, y2, angle;
    
    blank();
                
    for (j=0; j<ncircles; j++) if (particle[j].active)
    {
        glLineWidth(BOUNDARY_WIDTH);
        glColor3f(1.0, 1.0, 1.0);
        
        switch (plot) {
            case (P_KINETIC): 
            {
                ej = particle[j].energy;
                if (ej > 0.0) 
                {
                    hue = PARTICLE_HUE_MIN + (PARTICLE_HUE_MAX - PARTICLE_HUE_MIN)*ej/PARTICLE_EMAX;
                    if (hue > PARTICLE_HUE_MIN) hue = PARTICLE_HUE_MIN;
                    if (hue < PARTICLE_HUE_MAX) hue = PARTICLE_HUE_MAX;
                }
                radius = particle[j].radius;
                width = BOUNDARY_WIDTH;
                break;
            }
            case (P_NEIGHBOURS): 
            {
                hue = neighbour_color(particle[j].neighb);
                radius = particle[j].radius;
                width = BOUNDARY_WIDTH;
                break;
            }
            case (P_BONDS):
            {
                glLineWidth(BOUNDARY_WIDTH);
                hue = neighbour_color(particle[j].neighb);
                radius = 0.015;
//                 radius = particle[j].radius;
//                 radius = 0.8*particle[j].radius;
                width = 1;
                for (k = 0; k < particle[j].neighb; k++)
                {
                    m = particle[j].neighbours[k];
                    angle = particle[j].nghangle[k];
                    x1 = particle[j].xc + radius*cos(angle);
                    y1 = particle[j].yc + radius*sin(angle);
                    x2 = particle[m].xc - radius*cos(angle);
                    y2 = particle[m].yc - radius*sin(angle);
                    draw_line(x1, y1, x2, y2);
//                     draw_line(particle[j].xc, particle[j].yc, particle[m].xc, particle[m].yc);
                }
                break;
            }
        }
            
        hsl_to_rgb(hue, 0.9, 0.5, rgb);
        draw_colored_circle(particle[j].xc, particle[j].yc, radius, NSEG, rgb);
        
        glLineWidth(width);
        glColor3f(1.0, 1.0, 1.0);
        draw_circle(particle[j].xc, particle[j].yc, radius, NSEG);
    }
}


void print_parameters(double beta, double krepel)
{
    char message[100];
    
    if (INCREASE_BETA)  /* print force constant */
    {
        erase_area_hsl(XMAX - 0.39, YMAX - 0.1 + 0.025, 0.37, 0.05, 0.0, 0.9, 0.0);
        glColor3f(1.0, 1.0, 1.0);
        sprintf(message, "Temperature %.3f", 1.0/beta);
        write_text(XMAX - 0.7, YMAX - 0.1, message);
    }
    else if (INCREASE_KREPEL)  /* print force constant */
    {
        erase_area_hsl(XMAX - 0.24, YMAX - 0.1 + 0.025, 0.22, 0.05, 0.0, 0.9, 0.0);
        glColor3f(1.0, 1.0, 1.0);
        sprintf(message, "Force %.0f", krepel);
        write_text(XMAX - 0.42, YMAX - 0.1, message);
    }   
}


double repel_schedule(int i)
{
    static double kexponent;
    static int first = 1;
    double krepel;
    
    if (first) 
    {
        kexponent = log(KREPEL_FACTOR)/(double)(INITIAL_TIME + NSTEPS);
        first = 0;
    }
    krepel = KREPEL*exp(kexponent*(double)i);
    printf("krepel = %.3lg\n", krepel);
    return(krepel);
}


double temperature_schedule(int i)
{
    static double bexponent, omega;
    static int first = 1;
    double beta;
    
    if (first)
    {
        bexponent = log(BETA_FACTOR)/(double)(INITIAL_TIME + NSTEPS);
        omega = N_TOSCILLATIONS*DPI/(double)(INITIAL_TIME + NSTEPS);
        first = 0;
    }
    beta = BETA*exp(bexponent*(double)i);
    beta = beta*2.0/(1.0 + cos(omega*(double)i));
    printf("beta = %.3lg\n", beta);
    return(beta);
}

void animation()
{
    double time, scale, diss, rgb[3], dissip, gradient[2], x, y, dx, dy, dt, xleft, xright, a, b, 
            length, fx, fy, force[2], totalenergy = 0.0, krepel = KREPEL, pos[2], 
            beta = BETA, xi = 0.0;
    double *qx, *qy, *px, *py;
    int i, j, k, n, m, s, ij[2], i0, iplus, iminus, j0, jplus, jminus, p, q, total_neighbours = 0, min_nb, max_nb, close;
    static int imin, imax;
    static short int first = 1;
    t_particle *particle;
    int *hashgrid_number, *hashgrid_particles;  
    t_hashgrid *hashgrid;
    char message[100];


    particle = (t_particle *)malloc(NMAXCIRCLES*sizeof(t_particle));    /* particles */  
    
    hashgrid = (t_hashgrid *)malloc(HASHX*HASHY*sizeof(t_hashgrid));    /* hashgrid */      
        
    hashgrid_number = (int *)malloc(HASHX*HASHY*sizeof(int));   /* total number of particles in each hash grid cell */
    hashgrid_particles = (int *)malloc(HASHX*HASHY*HASHMAX*sizeof(int)); /* numbers of particles in each hash grid cell */
    
    qx = (double *)malloc(NMAXCIRCLES*sizeof(double));
    qy = (double *)malloc(NMAXCIRCLES*sizeof(double));
    px = (double *)malloc(NMAXCIRCLES*sizeof(double));
    py = (double *)malloc(NMAXCIRCLES*sizeof(double));
    
    /* initialise positions and radii of circles */
    init_particle_config(particle);


    /* initialise particles */
    for (i=0; i < NMAXCIRCLES; i++) 
    {
        particle[i].neighb = 0;

        particle[i].energy = 0.0;
        y = particle[i].yc;
        if (y >= YMAX) y -= particle[i].radius;
        if (y <= YMIN) y += particle[i].radius;

        if (RANDOM_RADIUS) particle[i].radius = particle[i].radius*(0.75 + 0.5*((double)rand()/RAND_MAX));
        
        particle[i].mass_inv = 1.0;
                
        particle[i].vx = V_INITIAL*gaussian();
        particle[i].vy = V_INITIAL*gaussian();
        particle[i].energy = (particle[i].vx*particle[i].vx + particle[i].vy*particle[i].vy)*particle[i].mass_inv;
        
        px[i] = particle[i].vx;
        py[i] = particle[i].vy;
    }
    for (i=ncircles; i < NMAXCIRCLES; i++) 
    {
        particle[i].active = 0;
        particle[i].neighb = 0;
        particle[i].energy = 0.0;
        particle[i].mass_inv = 1.0;
        particle[i].vx = 0.0;
        particle[i].vy = 0.0;
        px[i] = 0.0;
        py[i] = 0.0;        
    }
    
    /* add particles at the bottom as seed */
    if (PART_AT_BOTTOM) for (i=0; i<=NPART_BOTTOM; i++)
    {   
        x = XMIN + (double)i*(XMAX - XMIN)/(double)NPART_BOTTOM;
        y = YMIN + 2.0*MU;
        add_particle(x, y, 0.0, 0.0, particle);
        particle[ncircles-1].mass_inv = 1.0/MASS_PART_BOTTOM;
//         particle[ncircles-1].radius *= 1.2;
    }
    for (i=0; i<=NPART_BOTTOM; i++)
    {   
        x = XMIN + (double)i*(XMAX - XMIN)/(double)NPART_BOTTOM;
        y = YMIN + 4.0*MU;
        add_particle(x, y, 0.0, 0.0, particle);
        particle[ncircles-1].mass_inv = 1.0/MASS_PART_BOTTOM;
//         particle[ncircles-1].radius *= 1.2;
    }
    
    xi = 0.0;
    
    for (i=0; i<ncircles; i++)
    {
        printf("Particle %i at (%.3f, %.3f) of energy %.3f\n", i, particle[i].xc, particle[i].yc, particle[i].energy);
    }
    sleep(1);
    
    b = 0.25*SIGMA*SIGMA*DT_PARTICLE/MU_XI;
    
    update_hashgrid(particle, hashgrid_number, hashgrid_particles);
    
    blank();
//     glColor3f(0.0, 0.0, 0.0);

    glutSwapBuffers();
    
    sleep(SLEEP1);

    for (i=0; i<=INITIAL_TIME + NSTEPS; i++)
    {
	printf("Computing frame %d\n",i);
        
        if (INCREASE_KREPEL) krepel = repel_schedule(i);
        if (INCREASE_BETA) beta = temperature_schedule(i);

        blank();
        
        for(n=0; n<NVID; n++)
        {
            /* compute forces on particles */
            for (j=0; j<ncircles; j++) 
            {
                particle[j].neighb = 0;
                
                /* compute repelling force from other particles */
                /* determine neighboring grid points */
                i0 = particle[j].hashx;
                iminus = i0 - 1;    if (iminus < 0) iminus = 0;
                iplus = i0 + 1;     if (iplus >= HASHX) iplus = HASHX-1;

                j0 = particle[j].hashy;
                jminus = j0 - 1;    if (jminus < 0) jminus = 0;
                jplus = j0 + 1;     if (jplus >= HASHY) jplus = HASHY-1;
                
                fx = 0.0; 
                fy = 0.0;
                for (p=iminus; p<= iplus; p++)
                    for (q=jminus; q<= jplus; q++)
                        for (k=0; k<hashgrid_number[p*HASHY+q]; k++) 
                        {
                            m = p*HASHY*HASHMAX + q*HASHMAX + k;
                            if ((hashgrid_particles[m]!=j)&&(particle[hashgrid_particles[m]].active))
                            {
                                close = compute_repelling_force(j, hashgrid_particles[m], force, particle, krepel);
                                fx += force[0];
                                fy += force[1];
                                if ((close)&&(particle[j].neighb < MAXNEIGH)) 
                                {
                                    particle[j].neighbours[particle[j].neighb] = hashgrid_particles[m];
                                    x = particle[hashgrid_particles[m]].xc;
                                    y = particle[hashgrid_particles[m]].yc;
                                    particle[j].nghangle[particle[j].neighb] = argument(x - particle[j].xc, y - particle[j].yc);
                                }
                                if (close) particle[j].neighb++;
                            }
                        }
                            
                /* add harmonic force outside rectangle */
                if (particle[j].xc > XMAX) fx -= KSPRING_BOUNDARY*(particle[j].xc - XMAX);
                else if (particle[j].xc < XMIN) fx += KSPRING_BOUNDARY*(XMIN - particle[j].xc);
                if (particle[j].yc > YMAX) fy -= KSPRING_BOUNDARY*(particle[j].yc - YMAX);
                else if (particle[j].yc < YMIN) fy += KSPRING_BOUNDARY*(YMIN - particle[j].yc);
                
                /* add gravity */
                fy -= GRAVITY;
                                                
                if (FLOOR_FORCE)
                {
                    if (fx > FMAX) fx = FMAX;
                    if (fx < -FMAX) fx = -FMAX;
                    if (fy > FMAX) fy = FMAX;
                    if (fy < -FMAX) fy = -FMAX;
                }
                
                particle[j].fx = fx;
                particle[j].fy = fy;
            }
            
            /* timestep of thermostat algorithm */
            for (j=0; j<ncircles; j++) 
            {
                particle[j].vx = px[j] + 0.5*DT_PARTICLE*particle[j].fx;
                particle[j].vy = py[j] + 0.5*DT_PARTICLE*particle[j].fy;
                
                px[j] = particle[j].vx + 0.5*DT_PARTICLE*particle[j].fx;
                py[j] = particle[j].vy + 0.5*DT_PARTICLE*particle[j].fy;
                
                particle[j].energy = (px[j]*px[j] + py[j]*py[j])*particle[j].mass_inv;
                
                qx[j] = particle[j].xc + 0.5*DT_PARTICLE*px[j]*particle[j].mass_inv;
                qy[j] = particle[j].yc + 0.5*DT_PARTICLE*py[j]*particle[j].mass_inv;
                
                px[j] *= exp(- 0.5*DT_PARTICLE*xi);
                py[j] *= exp(- 0.5*DT_PARTICLE*xi);
            }
                
            /* compute kinetic energy */
            totalenergy = 0.0;
            for (j=0; j<ncircles; j++) totalenergy += particle[j].energy;
            a = DT_PARTICLE*(totalenergy - (double)ncircles/beta)/MU_XI;
            a += SIGMA*sqrt(DT_PARTICLE)*gaussian();
            xi = (xi + a - b*xi)/(1.0 + b);
            
            for (j=0; j<ncircles; j++) 
            {
                px[j] *= exp(- 0.5*DT_PARTICLE*xi);
                py[j] *= exp(- 0.5*DT_PARTICLE*xi);
                
                particle[j].xc = qx[j] + 0.5*DT_PARTICLE*px[j]*particle[j].mass_inv;
                particle[j].yc = qy[j] + 0.5*DT_PARTICLE*py[j]*particle[j].mass_inv;
            }
        }
                
        printf("Mean kinetic energy: %.3f\n", totalenergy/(double)ncircles); 
        
        total_neighbours = 0;
        min_nb = 100;
        max_nb = 0;
        for (j=0; j<ncircles; j++) if (particle[j].active)
        {
            total_neighbours += particle[j].neighb;
            if (particle[j].neighb > max_nb) max_nb = particle[j].neighb; 
            if (particle[j].neighb < min_nb) min_nb = particle[j].neighb; 
        }
        printf("Mean number of neighbours: %.3f\n", (double)total_neighbours/(double)ncircles); 
        printf("Min number of neighbours: %i\n", min_nb); 
        printf("Max number of neighbours: %i\n", max_nb); 
        
        draw_particles(particle, PLOT);

        /* add a particle */
        if ((ADD_PARTICLES)&&((i - ADD_TIME + 1)%ADD_PERIOD == 0))
        {
//             j = 0;
//             while (module2(particle[j].xc,particle[j].yc) > 0.7) j = rand()%ncircles;
//             x =  particle[j].xc + 2.5*MU;
//             y =  particle[j].yc;
            
            x =  XMIN + (XMAX - XMIN)*rand()/RAND_MAX;
            y =  YMAX + 0.1*rand()/RAND_MAX;
            add_particle(x, y, 0.0, 0.0, particle);
        }
        
        update_hashgrid(particle, hashgrid_number, hashgrid_particles);
        
        print_parameters(beta, krepel);

	glutSwapBuffers();


	if (MOVIE)
        {
            if (i >= INITIAL_TIME) save_frame_lj();
            else printf("Initial phase time %i of %i\n", i, INITIAL_TIME);
            
            if ((i >= INITIAL_TIME)&&(DOUBLE_MOVIE))
            {
                draw_particles(particle, PLOT_B);
                print_parameters(beta, krepel);
                glutSwapBuffers();
                save_frame_lj_counter(NSTEPS + 21 + counter);
                counter++;
            }

            /* it seems that saving too many files too fast can cause trouble with the file system */
            /* so this is to make a pause from time to time - parameter PAUSE may need adjusting   */
            if (i % PAUSE == PAUSE - 1)
            {
                printf("Making a short pause\n");
                sleep(PSLEEP);
                s = system("mv lj*.tif tif_ljones/");
            }
        }

    }

    if (MOVIE) 
    {
        if (DOUBLE_MOVIE) 
        {
            draw_particles(particle, PLOT);   
            print_parameters(beta, krepel);
            glutSwapBuffers();
        }
        for (i=0; i<MID_FRAMES; i++) save_frame_lj();
        if (DOUBLE_MOVIE) 
        {
            draw_particles(particle, PLOT_B);
            print_parameters(beta, krepel);
            glutSwapBuffers();
        }
        for (i=0; i<END_FRAMES; i++) save_frame_lj_counter(NSTEPS + MID_FRAMES + 1 + counter + i);

        s = system("mv lj*.tif tif_ljones/");
    }
    free(particle);
    free(hashgrid);
    free(hashgrid_number);
    free(hashgrid_particles);
    free(qx);
    free(qy);
    free(px);
    free(py);
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
    glutCreateWindow("Particles with Lennard-Jones interaction in a planar domain");

    init();

    glutDisplayFunc(display);

    glutMainLoop();

    return 0;
}

