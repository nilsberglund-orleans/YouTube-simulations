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

#define INITXMIN -1.95
#define INITXMAX 1.95	/* x interval for initial condition */
#define INITYMIN -1.05
#define INITYMAX 0.95	/* y interval for initial condition */

#define CIRCLE_PATTERN 8   /* pattern of circles, see list in global_ljones.c */

#define TWO_TYPES 0         /* set to 1 to have two types of particles */
#define TPYE_PROPORTION 0.8 /* proportion of particles of first type */
#define SYMMETRIZE_FORCE 1  /* set to 1 to symmetrize two-particle interaction, only needed if particles are not all the same */

#define INTERACTION 1       /* particle interaction, see list in global_ljones.c */
#define INTERACTION_B 1     /* particle interaction for second type of particle, see list in global_ljones.c */
#define SPIN_INTER_FREQUENCY 1.0 /* angular frequency of spin-spin interaction */
#define SPIN_INTER_FREQUENCY_B 2.0 /* angular frequency of spin-spin interaction for second particle type */

#define P_PERCOL 0.25       /* probability of having a circle in C_RAND_PERCOL arrangement */
#define NPOISSON 100        /* number of points for Poisson C_RAND_POISSON arrangement */
#define PDISC_DISTANCE 5.0  /* minimal distance in Poisson disc process, controls density of particles */
#define PDISC_CANDIDATES 100 /* number of candidates in construction of Poisson disc process */
#define RANDOM_POLY_ANGLE 0 /* set to 1 to randomize angle of polygons */

#define LAMBDA 2.0	    /* parameter controlling the dimensions of domain */
#define MU 0.015	    /* parameter controlling radius of particles */
#define MU_B 0.02427051     /* parameter controlling radius of particles of second type */
#define NPOLY 3             /* number of sides of polygon */
#define APOLY 0.5           /* angle by which to turn polygon, in units of Pi/2 */ 
#define MDEPTH 4            /* depth of computation of Menger gasket */
#define MRATIO 3            /* ratio defining Menger gasket */
#define MANDELLEVEL 1000    /* iteration level for Mandelbrot set */
#define MANDELLIMIT 10.0    /* limit value for approximation of Mandelbrot set */
#define FOCI 1              /* set to 1 to draw focal points of ellipse */
#define NGRIDX 80           /* number of grid point for grid of disks */
#define NGRIDY 25           /* number of grid point for grid of disks */

#define X_SHOOTER -0.2
#define Y_SHOOTER -0.6
#define X_TARGET 0.4
#define Y_TARGET 0.7        /* shooter and target positions in laser fight */

/* Parameters for length and speed of simulation */

#define NSTEPS 3100       /* number of frames of movie */
#define NVID 150          /* number of iterations between images displayed on screen */
#define NSEG 150         /* number of segments of boundary */
#define INITIAL_TIME 200    /* time after which to start saving frames */
#define BOUNDARY_WIDTH 1    /* width of particle boundary */
#define LINK_WIDTH 2        /* width of links between particles */
#define CONTAINER_WIDTH 4    /* width of container boundary */

#define PAUSE 1000         /* number of frames after which to pause */
#define PSLEEP 1         /* sleep time during pause */
#define SLEEP1  1        /* initial sleeping time */
#define SLEEP2  1   /* final sleeping time */
#define MID_FRAMES 20    /* number of still frames between parts of two-part movie */
#define END_FRAMES 100   /* number of still frames at end of movie */

/* Boundary conditions, see list in global_ljones.c */

#define BOUNDARY_COND 1

/* Plot type, see list in global_ljones.c  */

#define PLOT 0
#define PLOT_B 6        /* plot type for second movie */

#define COLOR_BONDS 0   /* set to 1 to color bonds according to length */

/* Color schemes */

#define COLOR_PALETTE 10     /* Color palette, see list in global_ljones.c  */

#define BLACK 1          /* background */

#define COLOR_SCHEME 3   /* choice of color scheme, see list in global_ljones.c  */

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
#define PARTICLE_HUE_MAX 50.0        /* color of saturated particle */
#define PARTICLE_EMAX 1.0e3           /* max energy for particle to survive */

#define RANDOM_RADIUS 0     /* set to 1 for random circle radius */
#define DT_PARTICLE 2.0e-6     /* time step for particle displacement */
#define KREPEL 10.0          /* constant in repelling force between particles */
#define EQUILIBRIUM_DIST 4.0  /* Lennard-Jones equilibrium distance for second type of particle */
#define EQUILIBRIUM_DIST_B 5.0  /* Lennard-Jones equilibrium distance for second type of particle */
#define REPEL_RADIUS 20.0    /* radius in which repelling force acts (in units of particle radius) */
#define DAMPING 0.0     /* damping coefficient of particles */
#define PARTICLE_MASS 1.0   /* mass of particle of radius MU */
#define PARTICLE_MASS_B 0.1   /* mass of particle of radius MU */
#define PARTICLE_INERTIA_MOMENT 0.05     /* moment of inertia of particle */
#define PARTICLE_INERTIA_MOMENT_B 0.02     /* moment of inertia of second type of particle */
#define V_INITIAL 0.5        /* initial velocity range */
#define OMEGA_INITIAL 2.0        /* initial angular velocity range */

#define THERMOSTAT 0        /* set to 1 to switch on thermostat */
#define SIGMA 5.0           /* noise intensity in thermostat */
#define BETA 0.3            /* initial inverse temperature */
#define MU_XI 0.05           /* friction constant in thermostat */
#define KSPRING_BOUNDARY 2.0e7    /* confining harmonic potential outside simulation region */
#define NBH_DIST_FACTOR 5.0       /* radius in which to count neighbours */
#define GRAVITY 0.0         /* gravity acting on all particles */

#define ROTATION 0            /* set to 1 to include rotation of particles */
#define COUPLE_ANGLE_TO_THERMOSTAT 1    /* set to 1 to couple angular degrees of freedom to thermostat */
#define DIMENSION_FACTOR 1.0  /* scaling factor taking into account number of degrees of freedom */  
#define KTORQUE 7.0        /* force constant in angular dynamics */
#define KTORQUE_B 10.0        /* force constant in angular dynamics */
#define KTORQUE_DIFF 150.0    /* force constant in angular dynamics for different particles */
#define DRAW_SPIN 0           /* set to 1 to draw spin vectors of particles */
#define DRAW_SPIN_B 0         /* set to 1 to draw spin vectors of particles */
#define DRAW_CROSS 1          /* set to 1 to draw cross on particles of second type */
#define SPIN_RANGE 5.0       /* range of spin-spin interaction */
#define SPIN_RANGE_B 5.0     /* range of spin-spin interaction for second type of particle */
#define QUADRUPOLE_RATIO 0.6  /* anisotropy in quadrupole potential */ 

#define INCREASE_BETA 0  /* set to 1 to increase BETA during simulation */
#define BETA_FACTOR 0.001   /* factor by which to change BETA during simulation */
#define N_TOSCILLATIONS 0.5   /* number of temperature oscillations in BETA schedule */
#define NO_OSCILLATION 1        /* set to 1 to have exponential BETA change only */
#define FINAL_CONSTANT_PHASE 200  /* final phase in which temperature is constant */

#define DECREASE_CONTAINER_SIZE 1   /* set to 1 to decrease size of container */
#define SYMMETRIC_DECREASE 0        /* set tp 1 to decrease container symmetrically */
#define COMPRESSION_RATIO 0.3       /* final size of container */
#define RESTORE_CONTAINER_SIZE 1    /* set to 1 to restore container to initial size at end of simulation */
#define RESTORE_TIME 700            /* time before end of sim at which to restore size */

#define MOVE_OBSTACLE 0     /* set to 1 to have a moving obstacle */
#define CENTER_VIEW_ON_OBSTACLE 0   /* set to 1 to center display on moving obstacle */
#define OBSTACLE_RADIUS 0.27 /* radius of obstacle for circle boundary conditions */
#define OBSTACLE_XMIN 0.0   /* initial position of obstacle */
#define OBSTACLE_XMAX 3.0   /* final position of obstacle */

#define INCREASE_KREPEL 0   /* set to 1 to increase KREPEL during simulation */
#define KREPEL_FACTOR 1000.0   /* factor by which to change KREPEL during simulation */

#define PART_AT_BOTTOM 0     /* set to 1 to include "seed" particles at bottom */
#define MASS_PART_BOTTOM 10000.0 /* mass of particles at bottom */
#define NPART_BOTTOM 100     /* number of particles at the bottom */

#define ADD_PARTICLES 0    /* set to 1 to add particles */
#define ADD_TIME 100        /* time at which to add first particle */
#define ADD_PERIOD 500       /* time interval between adding further particles */
#define SAFETY_FACTOR 3.0  /* no particles are added at distance less than MU*SAFETY_FACTOR of other particles */

#define FLOOR_FORCE 1      /* set to 1 to limit force on particle to FMAX */
#define FMAX 1.0e10         /* maximal force */

#define HASHX 23   /* size of hashgrid in x direction */
#define HASHY 14   /* size of hashgrid in y direction */
#define HASHMAX 100  /* maximal number of particles per hashgrid cell */
#define HASHGRID_PADDING 0.1    /* padding of hashgrid outside simulation window */

#define DRAW_COLOR_SCHEME 0     /* set to 1 to plot the color scheme */
#define COLORBAR_RANGE 8.0    /* scale of color scheme bar */
#define COLORBAR_RANGE_B 12.0    /* scale of color scheme bar for 2nd part */
#define ROTATE_COLOR_SCHEME 0   /* set to 1 to draw color scheme horizontally */

#include "global_ljones.c"
#include "sub_lj.c"

double xshift = 0.0;      /* x shift of shown window */
double xspeed = 0.0;      /* x speed of obstacle */

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
    static double lx, ly, padding;
    int i, j;
    
    if (first)
    {
        if ((BOUNDARY_COND == BC_PERIODIC)||(BOUNDARY_COND == BC_PERIODIC_CIRCLE)) padding = 0.0;
        else padding = HASHGRID_PADDING;
        lx = XMAX - XMIN + 2.0*padding;
        ly = YMAX - YMIN + 2.0*padding;
        first = 0;
    }
    
    i = (int)((double)HASHX*(x - XMIN + padding)/lx);
    j = (int)((double)HASHY*(y - YMIN + padding)/ly);
    
    if (i<0) i = 0;
    else if (i>=HASHX) i = HASHX-1;
    if (j<0) j = 0;
    else if (j>=HASHY) j = HASHY-1;
    
    ij[0] = i;
    ij[1] = j;

//     printf("Mapped (%.3f,%.3f) to (%i, %i)\n", x, y, ij[0], ij[1]);
}


double lennard_jones_force(double r, t_particle particle)
{
    int i;
    double rmin = 0.01, rplus, ratio = 1.0;
    
    if (r > REPEL_RADIUS*particle.radius) return(0.0);
    else
    {
        if (r > rmin) rplus = r;
        else rplus = rmin;
    
//         ratio = ipow(particle.eq_dist*particle.radius/rplus, 6);
        ratio = ipow(particle.eq_dist*particle.radius/rplus, 3);
        ratio = ratio*ratio;
        
        return((ratio - 2.0*ratio*ratio)/rplus);
    }
}

void aniso_lj_force(double r, double ca, double sa, double ca_rel, double sa_rel, double force[2], t_particle particle)
{
    int i;
    double rmin = 0.01, rplus, ratio = 1.0, c2, s2, c4, s4, a, aprime, f1, f2;
    
    if (r > REPEL_RADIUS*particle.radius)
    {
        force[0] = 0.0;
        force[1] = 0.0;
    }
    else
    {
        if (r > rmin) rplus = r;
        else rplus = rmin;
    
//         ratio = ipow(particle.eq_dist*particle.radius/rplus, 6);
        ratio = ipow(particle.eq_dist*particle.radius/rplus, 3);
        ratio = ratio*ratio;

        /* cos(2phi) and sin(2phi) */
        c2 = ca_rel*ca_rel - sa_rel*sa_rel;
        s2 = 2.0*ca_rel*sa_rel;

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

void penta_lj_force(double r, double ca, double sa, double ca_rel, double sa_rel, double force[2], t_particle particle)
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
    
    if (r > REPEL_RADIUS*particle.radius)
    {
        force[0] = 0.0;
        force[1] = 0.0;
    }
    else
    {
        if (r > rmin) rplus = r;
        else rplus = rmin;
    
//         ratio = ipow(particle.eq_dist*particle.radius/rplus, 6);
        ratio = ipow(particle.eq_dist*particle.radius/rplus, 3);
        ratio = ratio*ratio;

    
        /* cos(2phi) and sin(2phi) */
        c2 = ca_rel*ca_rel - sa_rel*sa_rel;
        s2 = 2.0*ca_rel*sa_rel;

        /* cos(4phi) and sin(4phi) */
        c4 = c2*c2 - s2*s2;
        s4 = 2.0*c2*s2;
        
        /* cos(5phi) and sin(5phi) */
        c5 = ca_rel*c4 - sa_rel*s4;
        s5 = sa_rel*c4 + ca_rel*s4;
        
        a = a0 - b0*c5;
        aprime = 5.0*b0*s5;
        
        f1 = ratio*(a - ratio)/rplus;
        f2 = ratio*aprime/rplus;
        
        force[0] = f1*ca - f2*sa;
        force[1] = f1*sa + f2*ca;
    }
}

double old_golden_ratio_force(double r, t_particle particle)
/* potential with two minima at distances whose ratio is the golden ratio Phi */
/* old version that does not work very well */
{
    int i;
    double x, y, z, rplus, ratio = 1.0, phi, a, phi3;
    static int first = 1;
    static double rmin, b, c, d;
    
    if (first) 
    {
        rmin = 0.5*particle.radius;
        phi = 0.5*(1.0 + sqrt(5.0));
        phi3 = 1.0/(phi*phi*phi);
        a = 0.66;
        b = 1.0 + phi3 + a;
        d = phi3*a;
        c = phi3 + a + d;
//         b = 7.04;
//         c = 13.66;
//         d = 6.7;
        first = 0;
        printf("a = %.4lg, b = %.4lg, c = %.4lg, d = %.4lg\n", a, b, c, d);
    }
    
    if (r > REPEL_RADIUS*particle.radius) return(0.0);
    else
    {
        if (r > rmin) rplus = r;
        else rplus = rmin;
    
        x = particle.eq_dist*particle.radius/rplus;
        y = x*x*x;
        z = d - c*y + b*y*y - y*y*y;
        return(x*z/rplus);
    }
}

double golden_ratio_force(double r, t_particle particle)
/* potential with two minima at distances whose ratio is the golden ratio Phi */
/* piecewise polynomial/LJ version */
{
    int i;
    double x, rplus, xm6, y1;
    static int first = 1;
    static double rmin, phi, a, h1, h2, phi6;
    
    if (first) 
    {
        rmin = 0.5*particle.radius;
        phi = 0.5*(1.0 + sqrt(5.0));
        a = 1.2;
        
        h1 = 1.0;       /* inner potential well depth */
        h2 = 10.0;      /* outer potential well depth */
        phi6 = ipow(phi, 6);
        
        first = 0;
    }
    
    if (r > REPEL_RADIUS*particle.radius) return(0.0);
    else
    {
        if (r > rmin) rplus = r;
        else rplus = rmin;
        
        x = rplus/(particle.eq_dist*particle.radius);
//         xm6 = 1.0/ipow(x, 6);
        xm6 = 1.0/ipow(x, 3);
        xm6 = xm6*xm6;
    
        if (x <= 1.0) return(12.0*h1*xm6*(1.0 - xm6)/x);
        else if (x <= a)
        {
            y1 = ipow(a - 1.0, 3);
            return(6.0*h1*(x - 1.0)*(a - x)/y1);
        }
        else if (x <= phi)
        {
            y1 = ipow(phi - a, 3);
            return(6.0*h2*(x - a)*(x - phi)/y1);
        }
        else return(12.0*h2*phi6*(1.0 - phi6*xm6)*xm6/x);
    }
}

void dipole_lj_force(double r, double ca, double sa, double ca_rel, double sa_rel, double force[2], t_particle particle)
{
    int i;
    double rmin = 0.01, rplus, ratio = 1.0, a, aprime, f1, f2;
    
    if (r > REPEL_RADIUS*MU)
    {
        force[0] = 0.0;
        force[1] = 0.0;
    }
    else
    {
        if (r > rmin) rplus = r;
        else rplus = rmin;
    
//         ratio = ipow(particle.eq_dist*particle.radius/rplus, 6);
        ratio = ipow(particle.eq_dist*particle.radius/rplus, 3);
        ratio = ratio*ratio;
            
        a = 1.0 + 0.25*ca_rel;
        aprime = -0.25*sa_rel;
        
        f1 = ratio*(a - ratio)/rplus;
        f2 = ratio*aprime/rplus;
        
        force[0] = f1*ca - f2*sa;
        force[1] = f1*sa + f2*ca;
    }
}

void quadrupole_lj_force(double r, double ca, double sa, double ca_rel, double sa_rel, double force[2], t_particle particle)
{
    int i;
    double rmin = 0.01, rplus, ratio = 1.0, a, aprime, f1, f2, ca2, sa2, x, y, dplus, dminus;
    static int first = 1;
    static double a0, b0, aplus, aminus;

    if (first)
    {
        dplus = cos(0.2*PI)*cos(0.1*PI);
//         dminus = 0.8*dplus;
        dminus = QUADRUPOLE_RATIO*dplus;
        aplus = ipow(1.0/dplus, 6);
        aminus = ipow(1.0/dminus, 6);
//         aminus = ipow(cos(0.2*PI)*(0.25 + 0.5*sin(0.1*PI)), 6);
        a0 = 0.5*(aplus + aminus);
        b0 = 0.5*(aplus - aminus);
        first = 0;
    }
    
    if (r > REPEL_RADIUS*particle.radius)
    {
        force[0] = 0.0;
        force[1] = 0.0;
    }
    else
    {
        if (r > rmin) rplus = r;
        else rplus = rmin;
    
//         ratio = ipow(particle.eq_dist*particle.radius/rplus, 6);
        ratio = ipow(particle.eq_dist*particle.radius/rplus, 3);
        ratio = ratio*ratio;
        
        /* cos(2*phi) and sin(2*phi) */
        ca2 = ca_rel*ca_rel - sa_rel*sa_rel;
        sa2 = 2.0*ca_rel*sa_rel;
            
        a = a0 + b0*ca2;
//         if (a == 0.0) a = 1.0e-10; 
        aprime = -2.0*b0*sa2;
                
        f1 = ratio*(a - ratio)/rplus;
        f2 = ratio*aprime/rplus;
        
        force[0] = f1*ca - f2*sa;
        force[1] = f1*sa + f2*ca;
    }
}


void quadrupole_lj_force2(double r, double ca, double sa, double ca_rel, double sa_rel, double force[2], t_particle particle)
{
    int i;
    double rmin = 0.01, rplus, ratio = 1.0, a, aprime, f1, f2, ca2, sa2, x, y, eqdist;
    static int first = 1;
    static double aplus, aminus, a0, b0;

    if (first)
    {
        aplus = ipow(cos(0.2*PI)*cos(0.1*PI), 6);
        aminus = 0.1*aplus;
//         aminus = 0.0;
//         aminus = -2.0*ipow(cos(0.2*PI)*(0.5*sin(0.1*PI)), 6);
//         aminus = ipow(cos(0.2*PI)*(0.25 + 0.5*sin(0.1*PI)), 6);
        a0 = 0.5*(aplus + aminus);
        b0 = 0.5*(aplus - aminus);
        first = 0;
    }
    
    if (r > REPEL_RADIUS*particle.radius)
    {
        force[0] = 0.0;
        force[1] = 0.0;
    }
    else
    {        
        if (r > rmin) rplus = r;
        else rplus = rmin;
    
        /* correct distance */
//         ratio = ipow(particle.eq_dist*particle.radius/rplus, 6);
        ratio = ipow(particle.eq_dist*particle.radius/rplus, 3);
        ratio = ratio*ratio;
        
        /* cos(2*phi) and sin(2*phi) */
        
        ca2 = ca_rel*ca_rel - sa_rel*sa_rel;
        sa2 = 2.0*ca_rel*sa_rel;
            
        a = a0 + b0*ca2;
        if (a == 0.0) a = 1.0e-10; 
        aprime = -2.0*b0*sa2;
        
//         f1 = ratio*(a - ratio)/rplus;
//         f2 = ratio*aprime/rplus;
        
        f1 = ratio*(aplus - ratio)/(rplus);
        f2 = ratio*(aminus - ratio)/(rplus);
        
//         force[0] = f1*ca_rel - f2*sa_rel;
//         force[1] = f1*sa_rel + f2*ca_rel;

        force[0] = f1*ca - f2*sa;
        force[1] = f1*sa + f2*ca;
    }
}



int compute_repelling_force(int i, int j, double force[2], double *torque, t_particle* particle, 
                            double krepel, int wrapx, int wrapy)
/* compute repelling force of particle j on particle i */
/* returns 1 if distance between particles is smaller than NBH_DIST_FACTOR*MU */
{
    double x1, y1, x2, y2, x3, y3, distance, r, f, angle, ca, sa, aniso, fx, fy, ff[2], cj, sj, ca_rel, sa_rel, dist_scaled, spin_f;
    double xtemp, ytemp;
    static double deltax = XMAX - XMIN, deltay = YMAX - YMIN, dxhalf = 0.5*(XMAX - XMIN), dyhalf = 0.5*(YMAX - YMIN);
    int periodx, periody, wwrapx;
//     static double factor = 0.5*(XMAX - XMIN);
//     int interaction;
    
    x1 = particle[i].xc;
    y1 = particle[i].yc;
    x2 = particle[j].xc;
    y2 = particle[j].yc;
    
     /* for monitoring purposes only */
    xtemp = x2;
    ytemp = y2;
    wwrapx = (BOUNDARY_COND == BC_KLEIN)&&(vabs(x2 - x1) > dxhalf);

    /* case of periodic boundary conditions */
    if ((BOUNDARY_COND == BC_PERIODIC)||(BOUNDARY_COND == BC_PERIODIC_CIRCLE))
    {
        if (wrapy != 0)
        {
            if (y2 > 0.0) y2 -= deltay;
            else y2 += deltay;
        }
        if (wrapx != 0)
        {
            if (x2 > 0.0) x2 -= deltax;
            else x2 += deltax;
        }        
    }
    distance = module2(x2 - x1, y2 - y1);
    
    if ((distance == 0.0)||(i == j))
    {
        force[0] = 0.0;
        force[1] = 0.0;
        *torque = 0.0;
        return(1);
    }
    else if (distance > REPEL_RADIUS*particle[i].radius) 
    {
        force[0] = 0.0;
        force[1] = 0.0;
        *torque = 0.0;
        return(0);
    }
    else
    {
        ca = (x2 - x1)/distance;
        sa = (y2 - y1)/distance;
        
        /* compute relative angle in case particles can rotate */
        if (ROTATION)
        {
            cj = cos(particle[j].angle);
            sj = sin(particle[j].angle);
            ca_rel = ca*cj + sa*sj;
            sa_rel = sa*cj - ca*sj;
        }
        else
        {
            ca_rel = ca;
            sa_rel = sa;
        }
        
        switch (particle[j].interaction) {
            case (I_COULOMB):
            {
                f = -krepel/(1.0e-8 + distance*distance);
                force[0] = f*ca;
                force[1] = f*sa; 
                break;
            }
            case (I_LENNARD_JONES):
            {
                f = krepel*lennard_jones_force(distance, particle[j]);
                force[0] = f*ca;
                force[1] = f*sa;   
                break;
            }
            case (I_LJ_DIRECTIONAL):
            {
                aniso_lj_force(distance, ca, sa, ca_rel, sa_rel, ff, particle[j]);
                force[0] = krepel*ff[0];
                force[1] = krepel*ff[1];
                break;
            }
            case (I_LJ_PENTA):
            {
                penta_lj_force(distance, ca, sa, ca_rel, sa_rel, ff, particle[j]);
                force[0] = krepel*ff[0];
                force[1] = krepel*ff[1];
                break;
            }
            case (I_GOLDENRATIO):
            {
                f = krepel*golden_ratio_force(distance, particle[j]);
                force[0] = f*ca;
                force[1] = f*sa;   
                break;
            }
            case (I_LJ_DIPOLE):
            {
                dipole_lj_force(distance, ca, sa, ca_rel, sa_rel, ff, particle[j]);
                force[0] = krepel*ff[0];
                force[1] = krepel*ff[1];
                break;
            }
            case (I_LJ_QUADRUPOLE):
            {
                quadrupole_lj_force(distance, ca, sa, ca_rel, sa_rel, ff, particle[j]);
                force[0] = krepel*ff[0];
                force[1] = krepel*ff[1];
                break;
            }
        }
    }
    
    if (ROTATION) 
    {   
        dist_scaled = distance/(particle[i].spin_range*particle[i].radius);
        spin_f = particle[i].spin_freq;
        *torque = sin(spin_f*(particle[j].angle - particle[i].angle))/(1.0e-8 + dist_scaled*dist_scaled);
        
        if (particle[i].type == particle[j].type) 
        {
            if (particle[i].type == 0) *torque *= KTORQUE;
            else *torque *= KTORQUE_B;
        }
        else *torque *= KTORQUE_DIFF;
    }
//     if (ROTATION) *torque = -sin(particle[j].angle - particle[i].angle)/(1.0e-8 + distance*distance);
//     if (ROTATION) *torque = ff[0]*sin(particle[i].angle) - ff[1]*cos(particle[i].angle);
//     if (ROTATION) *torque = ff[0]*sin(particle[j].angle) - ff[1]*cos(particle[j].angle);
    else *torque = 0.0;

    if ((distance < NBH_DIST_FACTOR*particle[i].radius)&&(j != i)) return(1);
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

int add_particle(double x, double y, double vx, double vy, double mass, short int type, t_particle particle[NMAXCIRCLES])
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
        
        particle[i].type = type;
        
        particle[i].xc = x;
        particle[i].yc = y;
        particle[i].radius = MU*sqrt(mass);
        particle[i].active = 1;
        particle[i].neighb = 0;
        particle[i].thermostat = 1;

        particle[i].energy = 0.0;

        if (RANDOM_RADIUS) particle[i].radius = particle[i].radius*(0.75 + 0.5*((double)rand()/RAND_MAX));
        
        particle[i].mass_inv = 1.0/mass;
        if (particle[i].type == 0) particle[i].inertia_moment_inv = 1.0/PARTICLE_INERTIA_MOMENT;
        else particle[i].inertia_moment_inv = 1.0/PARTICLE_INERTIA_MOMENT;
                
        particle[i].vx = vx;
        particle[i].vy = vy;
        particle[i].energy = (particle[i].vx*particle[i].vx + particle[i].vy*particle[i].vy)*particle[i].mass_inv;
        
        particle[i].angle = DPI*(double)rand()/RAND_MAX;
        particle[i].omega = 0.0;
        
        if (particle[i].type == 1)
        {
            particle[i].interaction = INTERACTION_B;
            particle[i].eq_dist = EQUILIBRIUM_DIST_B;
            particle[i].spin_range = SPIN_RANGE_B;
            particle[i].spin_freq = SPIN_INTER_FREQUENCY_B;            
        }
    
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

void draw_one_particle(t_particle particle, double xc, double yc, double radius, double angle, int nsides, double width, double rgb[3])
/* draw one of the particles */ 
{
    double ca, sa, x1, x2, y1, y2, xc1;
    
    if (CENTER_VIEW_ON_OBSTACLE) 
    {
        xc1 = xc - xshift;
        if (xc1 < XMIN) xc1 += XMAX - XMIN;
        else if (xc1 > XMAX) xc1 -= XMAX - XMIN;
    }
    else xc1 = xc;
    glColor3f(rgb[0], rgb[1], rgb[2]);
    if ((particle.interaction == I_LJ_QUADRUPOLE)||(particle.interaction == I_LJ_DIPOLE)) 
        draw_colored_rhombus(xc1, yc, radius, angle + APOLY*PID, rgb);
    else draw_colored_polygon(xc1, yc, radius, nsides, angle + APOLY*PID, rgb);
        
    /* draw crosses on particles of second type */
    if ((TWO_TYPES)&&(DRAW_CROSS))
        if (particle.type == 1)
        {
            if (ROTATION) angle = angle + APOLY*PID;
            else angle = APOLY*PID;
            ca = cos(angle);
            sa = sin(angle);
            glLineWidth(3);
            glColor3f(0.0, 0.0, 0.0);
            x1 = xc1 - MU_B*ca;
            y1 = yc - MU_B*sa;
            x2 = xc1 + MU_B*ca;
            y2 = yc + MU_B*sa;
            draw_line(x1, y1, x2, y2);
            x1 = xc1 - MU_B*sa;
            y1 = yc + MU_B*ca;
            x2 = xc1 + MU_B*sa;
            y2 = yc - MU_B*ca;
            draw_line(x1, y1, x2, y2);
        }
        
    glLineWidth(width);
    glColor3f(1.0, 1.0, 1.0);
    if ((particle.interaction == I_LJ_QUADRUPOLE)||(particle.interaction == I_LJ_DIPOLE)) 
        draw_rhombus(xc1, yc, radius, angle + APOLY*PID);
    else draw_polygon(xc1, yc, radius, nsides, angle + APOLY*PID);   
}

void draw_particles(t_particle particle[NMAXCIRCLES], int plot)
{
    int i, j, k, m, width, nnbg, nsides;
    double ej, hue, rgb[3], radius, x1, y1, x2, y2, angle, ca, sa, length, linkcolor;
    
    blank();
    glColor3f(1.0, 1.0, 1.0);
    
    /* draw the bonds first */
    if (plot == P_BONDS)
    {
        glLineWidth(LINK_WIDTH);
        for (j=0; j<ncircles; j++) if (particle[j].active)
        {
//             radius = particle[j].radius;
            for (k = 0; k < particle[j].neighb; k++)
            {
                m = particle[j].neighbours[k];
//                 angle = particle[j].nghangle[k];
                x1 = particle[j].xc;
                if (CENTER_VIEW_ON_OBSTACLE) x1 -= xshift;
                y1 = particle[j].yc;
                x2 = particle[m].xc;
                if (CENTER_VIEW_ON_OBSTACLE) x2 -= xshift;
                y2 = particle[m].yc;
                
                /* case of periodic boundary conditions */
                if ((BOUNDARY_COND == BC_PERIODIC)||(BOUNDARY_COND == BC_PERIODIC_CIRCLE))
                {
                    if (x2 - x1 > 0.5*(XMAX - XMIN)) x1 += XMAX - XMIN;
                    else if (x2 - x1 < -0.5*(XMAX - XMIN)) x1 -= XMAX - XMIN;
                    if (y2 - y1 > 0.5*(YMAX - YMIN)) y1 += YMAX - YMIN;
                    else if (y2 - y1 < -0.5*(YMAX - YMIN)) y1 -= YMAX - YMIN;        
                }
                
                if (COLOR_BONDS)
                {
                    length = module2(x1 - x2, y1 - y2)/particle[j].radius;
                    if (length < 1.5) linkcolor = 1.0;
                    else linkcolor = 1.0 - 0.75*(length - 1.5)/(NBH_DIST_FACTOR - 1.5);
//                     printf("length = %.3lg\t color = %.3lg\n", length, linkcolor);
                    glColor3f(linkcolor, linkcolor, linkcolor);
                }
                
                draw_line(x1, y1, x2, y2);
            }
        }
    }
                
    /* determine particle color and size */
    for (j=0; j<ncircles; j++) if (particle[j].active)
    {
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
//                 if (particle[j].type == 1) hue = 70.0;        /* to make second particle type more visible */
//                 if (particle[j].type == 1) hue = neighbour_color(7 - particle[j].neighb);
//                 else 
                hue = neighbour_color(particle[j].neighb);
                radius = particle[j].radius;
                width = 1;
                break;
            }
            case (P_ANGLE):
            {
                hue = particle[j].angle*particle[j].spin_freq/DPI;
                hue -= (double)((int)hue);
                hue = PARTICLE_HUE_MIN + (PARTICLE_HUE_MAX - PARTICLE_HUE_MIN)*hue;
                radius = particle[j].radius;
                width = BOUNDARY_WIDTH;
                break;
            }
            case (P_TYPE):
            {
                if (particle[j].type == 0) hue = 310.0;
                else hue = 70.0;
                radius = particle[j].radius;
                width = BOUNDARY_WIDTH;
                break;
            }
            case (P_DIRECTION): 
            {
                hue = argument(particle[j].vx, particle[j].vy);
                if (hue > DPI) hue -= DPI;
                if (hue < 0.0) hue += DPI;
                hue = PARTICLE_HUE_MIN + (PARTICLE_HUE_MAX - PARTICLE_HUE_MIN)*hue/DPI;
                radius = particle[j].radius;
                width = BOUNDARY_WIDTH;
                break;
            }
        }

        switch (particle[j].interaction) {
            case (I_LJ_DIRECTIONAL): 
            {   
                nsides = 4;
                break;
            }
            case (I_LJ_PENTA): 
            {
                nsides = 5;
                break;
            }
            case (I_LJ_QUADRUPOLE):
            {
                nsides = 4;
                break;
            }
            default: nsides = NSEG;
        }

        if (plot == P_BONDS) hsl_to_rgb_turbo(hue, 0.9, 0.5, rgb);
        else if (plot == P_DIRECTION) hsl_to_rgb_twilight(hue, 0.9, 0.5, rgb);
        else hsl_to_rgb(hue, 0.9, 0.5, rgb);
        angle = particle[j].angle + APOLY*DPI;
        
        draw_one_particle(particle[j], particle[j].xc, particle[j].yc, radius, angle, nsides, width, rgb);
                
        /* in case of periodic b.c., draw translates of particles */
        if ((BOUNDARY_COND == BC_PERIODIC)||(BOUNDARY_COND == BC_PERIODIC_CIRCLE))
        {
            x1 = particle[j].xc;
            y1 = particle[j].yc;
            
            for (i=-2; i<3; i++)
                for (k=-1; k<2; k++)
                    draw_one_particle(particle[j], x1 + (double)i*(XMAX - XMIN), y1 + (double)k*(YMAX - YMIN), radius, angle,                   nsides, width, rgb);
        }
    }
    
    /* draw spin vectors */
    if ((DRAW_SPIN)||(DRAW_SPIN_B))
    {
        glLineWidth(width);
        for (j=0; j<ncircles; j++) 
            if ((particle[j].active)&&(((DRAW_SPIN)&&(particle[j].type == 0))||((DRAW_SPIN_B)&&(particle[j].type == 1))))
        {
//             x1 = particle[j].xc - 2.0*MU*cos(particle[j].angle);
//             y1 = particle[j].yc - 2.0*MU*sin(particle[j].angle);
            x1 = particle[j].xc;
//             if (CENTER_VIEW_ON_OBSTACLE) x1 -= xshift;
            y1 = particle[j].yc;
            x2 = particle[j].xc + 2.0*MU*cos(particle[j].angle);
//             if (CENTER_VIEW_ON_OBSTACLE) x2 -= xshift;
            y2 = particle[j].yc + 2.0*MU*sin(particle[j].angle);
            draw_line(x1, y1, x2, y2);
        }
    }
}


void draw_container(double xmin, double xmax)
/* draw the container, for certain boundary conditions */
{
    int i;
    double rgb[3], x;
    char message[100];
    
    switch (BOUNDARY_COND) {
        case (BC_SCREEN): 
        {
            /* do nothing */
            break;
        }
        case (BC_RECTANGLE): 
        {
            glColor3f(1.0, 1.0, 1.0);
            glLineWidth(CONTAINER_WIDTH);
            
            draw_line(INITXMIN, INITYMIN, INITXMAX, INITYMIN);
            draw_line(INITXMIN, INITYMAX, INITXMAX, INITYMAX);
            
            if (!SYMMETRIC_DECREASE) draw_line(INITXMAX, INITYMIN,  INITXMAX, INITYMAX);

            draw_line(xmin, INITYMIN, xmin, INITYMAX);
            draw_line(XMIN, 0.5*(INITYMIN + INITYMAX), xmin, 0.5*(INITYMIN + INITYMAX));
            
            if (SYMMETRIC_DECREASE)
            {
                draw_line(xmax, INITYMIN, xmax, INITYMAX);
                draw_line(XMAX, 0.5*(INITYMIN + INITYMAX), xmax, 0.5*(INITYMIN + INITYMAX));
            }
            
            break;
        }
        case (BC_CIRCLE):
        {
            glLineWidth(CONTAINER_WIDTH);
            hsl_to_rgb(300.0, 0.1, 0.5, rgb);
            for (i=-1; i<2; i++)
            {
                x = xmin + (double)i*(XMAX - XMIN);
                draw_colored_circle(x, 0.0, OBSTACLE_RADIUS, NSEG, rgb);
                glColor3f(1.0, 1.0, 1.0);
                draw_circle(x, 0.0, OBSTACLE_RADIUS, NSEG);
            }
        }
        case (BC_PERIODIC_CIRCLE):
        {
            glLineWidth(CONTAINER_WIDTH);
            hsl_to_rgb(300.0, 0.1, 0.5, rgb);
            for (i=-1; i<2; i++)
            {
                if (CENTER_VIEW_ON_OBSTACLE) draw_colored_circle(0.0, 0.0, OBSTACLE_RADIUS, NSEG, rgb);
//                 x = xmin + (double)i*(XMAX - XMIN);
                else draw_colored_circle(xmin + (double)i*(XMAX - XMIN), 0.0, OBSTACLE_RADIUS, NSEG, rgb);
                
                glColor3f(1.0, 1.0, 1.0);
                draw_circle(x, 0.0, OBSTACLE_RADIUS, NSEG);
                
                glColor3f(0.0, 0.0, 0.0);
                sprintf(message, "Speed %.2f", xspeed);
                write_text(-0.17, -0.025, message); 
            }
        }
    }
}

void print_parameters(double beta, double temperature, double krepel, double lengthcontainer, double boundary_force, 
                      short int left)
{
    char message[100];
    int i;
    double density, hue, rgb[3], logratio, y;
    static double xbox, xtext, xmid, xmidtext, xxbox, xxtext, pressures[100], meanpressure = 0.0;
    static int first = 1, i_pressure, naverage = 100;
    
    if (first)
    {
        if (left)
        {
            xbox = XMIN + 0.4;
            xtext = XMIN + 0.08;
            xxbox = XMAX - 0.39;
            xxtext = XMAX - 0.73;
       }
        else
        {
            xbox = XMAX - 0.39;
            xtext = XMAX - 0.73;
            xxbox = XMIN + 0.4;
            xxtext = XMIN + 0.08;
        }
        xmid = 0.5*(XMIN + XMAX) - 0.1;
        xmidtext = xmid - 0.24;
        for (i=0; i<naverage; i++) pressures[i] = 0.0;
        i_pressure = 0;
        
        first = 0;
    }
    
    /* table of pressures */
    pressures[i_pressure] = boundary_force/(lengthcontainer + INITYMAX - INITYMIN);
    i_pressure++;
    if (i_pressure == naverage) i_pressure = 0;
    
    for (i=0; i<naverage; i++) meanpressure += pressures[i];
    meanpressure = meanpressure/(double)naverage;
    
    y = YMAX - 0.1;
    if (INCREASE_BETA)  /* print temperature */
    {
        logratio = log(beta/BETA)/log(0.5*BETA_FACTOR);
        if (logratio > 1.0) logratio = 1.0;
        if (BETA_FACTOR > 1.0) hue = PARTICLE_HUE_MAX - (PARTICLE_HUE_MAX - PARTICLE_HUE_MIN)*logratio;
        else hue = PARTICLE_HUE_MIN - (PARTICLE_HUE_MIN - PARTICLE_HUE_MAX)*logratio;
        erase_area_hsl_turbo(xbox, y + 0.025, 0.37, 0.05, hue, 0.9, 0.5);
        if ((hue < 90)||(hue > 270)) glColor3f(1.0, 1.0, 1.0);
        else glColor3f(0.0, 0.0, 0.0);
        sprintf(message, "Temperature %.3f", 1.0/beta);
        write_text(xtext, y, message);
        y -= 0.1;
    }
    if (DECREASE_CONTAINER_SIZE)  /* print density */
    {
        density = (double)ncircles/((lengthcontainer)*(INITYMAX - INITYMIN));
        erase_area_hsl(xbox, y + 0.025, 0.37, 0.05, 0.0, 0.9, 0.0);
        glColor3f(1.0, 1.0, 1.0);
        sprintf(message, "Density %.3f", density);
        write_text(xtext, y, message);
        
        erase_area_hsl(xmid, y + 0.025, 0.37, 0.05, 0.0, 0.9, 0.0);
        glColor3f(1.0, 1.0, 1.0);
        sprintf(message, "Temperature %.3f", temperature);
        write_text(xmidtext, y, message);

        erase_area_hsl(xxbox, y + 0.025, 0.37, 0.05, 0.0, 0.9, 0.0);
        glColor3f(1.0, 1.0, 1.0);
        sprintf(message, "Pressure %.3f", meanpressure);
        write_text(xxtext, y, message);

    }   
    else if (INCREASE_KREPEL)  /* print force constant */
    {
        erase_area_hsl(xbox, y + 0.025, 0.22, 0.05, 0.0, 0.9, 0.0);
        glColor3f(1.0, 1.0, 1.0);
        sprintf(message, "Force %.0f", krepel);
        write_text(xtext + 0.28, y, message);
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
        bexponent = log(BETA_FACTOR)/(double)(INITIAL_TIME + NSTEPS - FINAL_CONSTANT_PHASE);
        omega = N_TOSCILLATIONS*DPI/(double)(INITIAL_TIME + NSTEPS - FINAL_CONSTANT_PHASE);
        first = 0;
    }
    if (i < INITIAL_TIME) beta = BETA;
    else if (i > INITIAL_TIME + NSTEPS - FINAL_CONSTANT_PHASE) beta = BETA*BETA_FACTOR;
    else
    {
        beta = BETA*exp(bexponent*(double)(i - INITIAL_TIME));
        if (!NO_OSCILLATION) beta = beta*2.0/(1.0 + cos(omega*(double)(i - INITIAL_TIME)));
    }
    printf("beta = %.3lg\n", beta);
    return(beta);
}

double container_size_schedule(int i)
{
    if ((i < INITIAL_TIME)||(i > INITIAL_TIME + NSTEPS - RESTORE_TIME)) return(INITXMIN);
    else 
        return(INITXMIN + (1.0-COMPRESSION_RATIO)*(INITXMAX-INITXMIN)*(double)(i-INITIAL_TIME)/(double)(NSTEPS-RESTORE_TIME));
}

double obstacle_schedule_old(int i)
{
    double time;
    static double t1 = 0.5, t2 = 0.75, t3 = 0.875;
    
    if (i < INITIAL_TIME) return(OBSTACLE_XMIN);
    else 
    {
        time = (double)(i-INITIAL_TIME)/(double)(NSTEPS);
        
        if (time < t1) return(OBSTACLE_XMIN + (OBSTACLE_XMAX - OBSTACLE_XMIN)*time/t1);
        else if (time < t2) return(OBSTACLE_XMIN + (OBSTACLE_XMAX - OBSTACLE_XMIN)*(time - t1)/(t2 - t1));
        else if (time < t3) return(OBSTACLE_XMIN + (OBSTACLE_XMAX - OBSTACLE_XMIN)*(time - t2)/(t3 - t2));
        else return(OBSTACLE_XMAX);
    }
}

double obstacle_schedule(int i)
{
    double time;
    double x;
//     static double t1 = 0.5, t2 = 0.75, t3 = 0.875;
    
    if (i < INITIAL_TIME) return(OBSTACLE_XMIN);
    else 
    {
        time = (double)(i-INITIAL_TIME)/(double)(NSTEPS);
        
        x = OBSTACLE_XMIN + 50.0*time*time;
        xspeed = 100.0*time;
        while (x > XMAX) x += XMIN - XMAX;
        return(x);
    }
}

double compute_boundary_force(t_particle particle, double *fx, double *fy, double xleft, double xright)
{
    int i;
    double xmin, xmax, ymin, ymax, padding, r, cphi, sphi, f, fperp = 0.0, x;
    
    switch(BOUNDARY_COND){
        case (BC_SCREEN):
        {
            /* add harmonic force outside screen */
            if (particle.xc > XMAX) *fx -= KSPRING_BOUNDARY*(particle.xc - XMAX);
            else if (particle.xc < XMIN) *fx += KSPRING_BOUNDARY*(XMIN - particle.xc);
            if (particle.yc > YMAX) *fy -= KSPRING_BOUNDARY*(particle.yc - YMAX);
            else if (particle.yc < YMIN) *fy += KSPRING_BOUNDARY*(YMIN - particle.yc);
            return(fperp);
        }
        case (BC_RECTANGLE):
        {
            /* add harmonic force outside rectangular box */
            padding = MU + 0.01;
            xmin = xleft + padding;
            xmax = xright - padding;
            ymin = INITYMIN + padding;
            ymax = INITYMAX - padding;
            
            if (particle.xc > xmax) 
            {
                fperp = KSPRING_BOUNDARY*(particle.xc - xmax);
                *fx -= fperp;
            }
            else if (particle.xc < xmin) 
            {
                fperp = KSPRING_BOUNDARY*(xmin - particle.xc);
                *fx += fperp;
            }
            if (particle.yc > ymax) 
            {
                fperp = KSPRING_BOUNDARY*(particle.yc - ymax);
                *fy -= fperp;
            }
            else if (particle.yc < ymin) 
            {
                fperp = KSPRING_BOUNDARY*(ymin - particle.yc);
                *fy += fperp;
            }
//             if (particle.xc > xmax) *fx -= KSPRING_BOUNDARY*(particle.xc - xmax);
//             else if (particle.xc < xmin) *fx += KSPRING_BOUNDARY*(xmin - particle.xc);
//             if (particle.yc > ymax) *fy -= KSPRING_BOUNDARY*(particle.yc - ymax);
//             else if (particle.yc < ymin) *fy += KSPRING_BOUNDARY*(ymin - particle.yc);
            
            return(fperp);
        }
        case (BC_CIRCLE):
        {
            /* add harmonic force outside screen */
            if (particle.xc > XMAX) *fx -= KSPRING_BOUNDARY*(particle.xc - XMAX);
            else if (particle.xc < XMIN) *fx += KSPRING_BOUNDARY*(XMIN - particle.xc);
            if (particle.yc > YMAX) *fy -= KSPRING_BOUNDARY*(particle.yc - YMAX);
            else if (particle.yc < YMIN) *fy += KSPRING_BOUNDARY*(YMIN - particle.yc);
            
            /* add harmonic force from obstacle */
            for (i=-1; i<2; i++) 
            {
                x = xleft + (double)i*(XMAX - XMIN);
                if (vabs(particle.xc - x) < 1.1*OBSTACLE_RADIUS)
                {
                    r = module2(particle.xc - x, particle.yc);
                    if (r < 1.0e-5) r = 1.0e-05;
                    cphi = (particle.xc - x)/r;
                    sphi = particle.yc/r;
                    padding = MU + 0.03;
                    if (r < OBSTACLE_RADIUS + padding)
                    {
                        f = KSPRING_BOUNDARY*(OBSTACLE_RADIUS + padding - r);
                        *fx += f*cphi;
                        *fy += f*sphi;
                    }
                }
            }
            return(fperp);
        }
        case (BC_PERIODIC_CIRCLE):
        {
            /* add harmonic force from obstacle */
            for (i=-1; i<2; i++) 
            {
                x = xleft + (double)i*(XMAX - XMIN);
                if (vabs(particle.xc - x) < 1.1*OBSTACLE_RADIUS)
                {
                    r = module2(particle.xc - x, particle.yc);
                    if (r < 1.0e-5) r = 1.0e-05;
                    cphi = (particle.xc - x)/r;
                    sphi = particle.yc/r;
                    padding = MU + 0.03;
                    if (r < OBSTACLE_RADIUS + padding)
                    {
                        f = KSPRING_BOUNDARY*(OBSTACLE_RADIUS + padding - r);
                        *fx += f*cphi;
                        *fy += f*sphi;
                    }
                }
            }
            return(fperp);
        }
    }
}

void animation()
{
    double time, scale, diss, rgb[3], dissip, gradient[2], x, y, dx, dy, dt, xleft, xright, a, b, 
            length, fx, fy, force[2], totalenergy = 0.0, krepel = KREPEL, pos[2], prop, vx, 
            beta = BETA, xi = 0.0, xmincontainer = INITXMIN, xmaxcontainer = INITXMAX, torque, torque_ij, 
            fboundary = 0.0;
    double *qx, *qy, *px, *py, *qangle, *pangle;
    int i, j, k, n, m, s, ij[2], i0, iplus, iminus, j0, jplus, jminus, p, q, p1, q1, total_neighbours = 0, 
        min_nb, max_nb, close, wrapx = 0, wrapy = 0, nactive = 0, nadd_particle = 0;
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
    qangle = (double *)malloc(NMAXCIRCLES*sizeof(double));
    pangle = (double *)malloc(NMAXCIRCLES*sizeof(double));
    
    /* initialise positions and radii of circles */
    init_particle_config(particle);


    /* initialise particles */
    for (i=0; i < NMAXCIRCLES; i++) 
    {
        /* set particle type */
        particle[i].type = 0;
        if ((TWO_TYPES)&&((double)rand()/RAND_MAX > TPYE_PROPORTION)) 
        {
            particle[i].type = 1;
            particle[i].radius = MU_B;
        }
        
        particle[i].neighb = 0;
        particle[i].thermostat = 1;

//         particle[i].energy = 0.0;
        y = particle[i].yc;
        if (y >= YMAX) y -= particle[i].radius;
        if (y <= YMIN) y += particle[i].radius;

        if (RANDOM_RADIUS) particle[i].radius = particle[i].radius*(0.75 + 0.5*((double)rand()/RAND_MAX));
        
        if (particle[i].type == 0)
        {
            particle[i].interaction = INTERACTION;
            particle[i].eq_dist = EQUILIBRIUM_DIST;
            particle[i].spin_range = SPIN_RANGE;
            particle[i].spin_freq = SPIN_INTER_FREQUENCY;
            particle[i].mass_inv = 1.0/PARTICLE_MASS;
            particle[i].inertia_moment_inv = 1.0/PARTICLE_INERTIA_MOMENT;
        }
        else 
        {
            particle[i].interaction = INTERACTION_B;
            particle[i].eq_dist = EQUILIBRIUM_DIST_B;
            particle[i].spin_range = SPIN_RANGE_B;
            particle[i].spin_freq = SPIN_INTER_FREQUENCY_B;
            particle[i].mass_inv = 1.0/PARTICLE_MASS_B;
            particle[i].inertia_moment_inv = 1.0/PARTICLE_INERTIA_MOMENT_B;
        }
                
        particle[i].vx = V_INITIAL*gaussian();
        particle[i].vy = V_INITIAL*gaussian();
        particle[i].energy = (particle[i].vx*particle[i].vx + particle[i].vy*particle[i].vy)*particle[i].mass_inv;
        
        px[i] = particle[i].vx;
        py[i] = particle[i].vy;
        
        if (ROTATION) 
        {
            particle[i].angle = DPI*(double)rand()/RAND_MAX;
            particle[i].omega = OMEGA_INITIAL*gaussian();
            if (COUPLE_ANGLE_TO_THERMOSTAT) 
                particle[i].energy += particle[i].omega*particle[i].omega*particle[i].inertia_moment_inv;
        }
        else 
        {
            particle[i].angle = 0.0;
            particle[i].omega = 0.0;
        }
        pangle[i] = particle[i].omega;
    }
    /* initialize dummy values in case particles are added */
    for (i=ncircles; i < NMAXCIRCLES; i++) 
    {
        particle[i].type = 0;
        particle[i].active = 0;
        particle[i].neighb = 0;
        particle[i].thermostat = 0;
        particle[i].energy = 0.0;
        particle[i].mass_inv = 1.0/PARTICLE_MASS;
        particle[i].inertia_moment_inv = 1.0/PARTICLE_INERTIA_MOMENT;
        particle[i].vx = 0.0;
        particle[i].vy = 0.0;
        px[i] = 0.0;
        py[i] = 0.0;        
        particle[i].angle = DPI*(double)rand()/RAND_MAX;
        particle[i].omega = 0.0;
        pangle[i] = 0.0;
        particle[i].interaction = INTERACTION;
        particle[i].eq_dist = EQUILIBRIUM_DIST;
        particle[i].spin_range = SPIN_RANGE;
        particle[i].spin_freq = SPIN_INTER_FREQUENCY;
   }
    
    /* add particles at the bottom as seed */
    if (PART_AT_BOTTOM) for (i=0; i<=NPART_BOTTOM; i++)
    {   
        x = XMIN + (double)i*(XMAX - XMIN)/(double)NPART_BOTTOM;
        y = YMIN + 2.0*MU;
        add_particle(x, y, 0.0, 0.0, MASS_PART_BOTTOM, 0, particle);
    }
    if (PART_AT_BOTTOM) for (i=0; i<=NPART_BOTTOM; i++)
    {   
        x = XMIN + (double)i*(XMAX - XMIN)/(double)NPART_BOTTOM;
        y = YMIN + 4.0*MU;
        add_particle(x, y, 0.0, 0.0, MASS_PART_BOTTOM, 0, particle);
    }
    
    /* inactivate particles in obstacle */
    if ((BOUNDARY_COND == BC_CIRCLE)||(BOUNDARY_COND == BC_PERIODIC_CIRCLE))
        for (i=0; i< ncircles; i++)
            if ((module2(particle[i].xc - OBSTACLE_XMIN, particle[i].yc) < 1.2*OBSTACLE_RADIUS)) 
                particle[i].active = 0;
            
    /* count number of active particles */
    for (i=0; i< ncircles; i++) nactive += particle[i].active;
    
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
        if (DECREASE_CONTAINER_SIZE) 
        {
            xmincontainer = container_size_schedule(i);
            if (SYMMETRIC_DECREASE) xmaxcontainer = -container_size_schedule(i);
        }
        if (MOVE_OBSTACLE) 
        {
            xmincontainer = obstacle_schedule(i);
            xshift = xmincontainer;
            while (xshift > XMAX) xshift -= XMAX - XMIN;
        }
        
        blank();
        
        fboundary = 0.0;
        
        for(n=0; n<NVID; n++)
        {
            /* compute forces on particles */
            for (j=0; j<ncircles; j++) if (particle[j].active)
            {
                particle[j].neighb = 0;
                
                /* compute repelling force from other particles */
                /* determine neighboring grid points */
                i0 = particle[j].hashx;
                iminus = i0 - 1;  
                iplus = i0 + 1;   

                j0 = particle[j].hashy;
                jminus = j0 - 1; 
                jplus = j0 + 1;     
                
                if ((BOUNDARY_COND != BC_PERIODIC)&&(BOUNDARY_COND != BC_PERIODIC_CIRCLE))
                {
                    if (iminus < 0) iminus = 0;
                    else if (iplus >= HASHX) iplus = HASHX-1;

                    if (jminus < 0) jminus = 0;
                    else if (jplus >= HASHY) jplus = HASHY-1;
                }
                    
                
                fx = 0.0; 
                fy = 0.0;
                torque = 0.0;
                for (p=iminus; p<= iplus; p++)
                    for (q=jminus; q<= jplus; q++)
                    {
                        p1 = p;
                        q1 = q;
                        wrapx = 0;
                        wrapy = 0;
                        /* case of periodic boundary conditions */
                        if ((BOUNDARY_COND == BC_PERIODIC)||(BOUNDARY_COND == BC_PERIODIC_CIRCLE))
                        {
                            if (p1 < 0) 
                            {
                                p1 = HASHX-1;
                                wrapx = -1;
                            }
                            else if (p1 >= HASHX) 
                            {
                                p1 = 0;
                                wrapx = 1;
                            }
                            if (q1 < 0) 
                            {
                                q1 = HASHY-1;
                                wrapy = -1;
                            }
                            else if (q1 >= HASHY) 
                            {
                                q1 = 0;
                                wrapy = 1;
                            }
                        }
                        for (k=0; k<hashgrid_number[p1*HASHY+q1]; k++) 
                        {
                            m = p1*HASHY*HASHMAX + q1*HASHMAX + k;
                            if ((hashgrid_particles[m]!=j)&&(particle[hashgrid_particles[m]].active))
                            {
                                close = compute_repelling_force(j, hashgrid_particles[m], force, &torque_ij, particle, krepel, wrapx, wrapy);
                                if (SYMMETRIZE_FORCE)
                                {
                                    fx += 0.5*force[0];
                                    fy += 0.5*force[1];
                                    close = compute_repelling_force(j, hashgrid_particles[m], force, &torque_ij, particle, krepel, wrapx, wrapy);
                                    fx += 0.5*force[0];
                                    fy += 0.5*force[1];
                                }
                                else 
                                {
                                    fx += force[0];
                                    fy += force[1];
                                }
                                torque += torque_ij;
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
                    }      
                /* take care of boundary conditions */
                fboundary += compute_boundary_force(particle[j], &fx, &fy, xmincontainer, xmaxcontainer);
                
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
                particle[j].torque = torque;
                
//                 if (j%500 == 1) printf("Particle %i angle = %.5lg omega = %.5lg torque = %.5lg\n", j, particle[j].angle, particle[j].omega, particle[j].torque);
//                 if (j%500 == 1) printf("Particle %i x = %.5lg vx = %.5lg fx = %.5lg\n", j, particle[j].xc, particle[j].vx, particle[j].fx);
            }
            
            /* timestep of thermostat algorithm */
            for (j=0; j<ncircles; j++) if (particle[j].active)
            {
                particle[j].vx = px[j] + 0.5*DT_PARTICLE*particle[j].fx;
                particle[j].vy = py[j] + 0.5*DT_PARTICLE*particle[j].fy;
                particle[j].omega = pangle[j] + 0.5*DT_PARTICLE*particle[j].torque;
                
                px[j] = particle[j].vx + 0.5*DT_PARTICLE*particle[j].fx;
                py[j] = particle[j].vy + 0.5*DT_PARTICLE*particle[j].fy;
                pangle[j] = particle[j].omega + 0.5*DT_PARTICLE*particle[j].torque;
                
                particle[j].energy = (px[j]*px[j] + py[j]*py[j])*particle[j].mass_inv;
                if ((COUPLE_ANGLE_TO_THERMOSTAT)&&(particle[j].thermostat))
                    particle[j].energy += pangle[j]*pangle[j]*particle[j].inertia_moment_inv;
                
                qx[j] = particle[j].xc + 0.5*DT_PARTICLE*px[j]*particle[j].mass_inv;
                qy[j] = particle[j].yc + 0.5*DT_PARTICLE*py[j]*particle[j].mass_inv;
                qangle[j] = particle[j].angle + 0.5*DT_PARTICLE*pangle[j]*particle[j].inertia_moment_inv;
                
                if ((THERMOSTAT)&&(particle[j].thermostat))
                {
                    px[j] *= exp(- 0.5*DT_PARTICLE*xi);
                    py[j] *= exp(- 0.5*DT_PARTICLE*xi);
                }
                if ((COUPLE_ANGLE_TO_THERMOSTAT)&&(particle[j].thermostat)) 
                    pangle[j] *= exp(- 0.5*DT_PARTICLE*xi);
            }
                
            /* compute kinetic energy */
            totalenergy = 0.0;
            for (j=0; j<ncircles; j++)
                if ((particle[j].active)&&(particle[j].thermostat))
                    totalenergy += particle[j].energy;
            totalenergy *= DIMENSION_FACTOR;    /* normalize energy to take number of degrees of freedom into account */
            if (THERMOSTAT)
            {
                a = DT_PARTICLE*(totalenergy - (double)nactive/beta)/MU_XI;
                a += SIGMA*sqrt(DT_PARTICLE)*gaussian();
                xi = (xi + a - b*xi)/(1.0 + b);
            }
            
            for (j=0; j<ncircles; j++) if (particle[j].active) 
            {
                if ((THERMOSTAT)&&(particle[j].thermostat))
                {
                    px[j] *= exp(- 0.5*DT_PARTICLE*xi);
                    py[j] *= exp(- 0.5*DT_PARTICLE*xi);
                }
                else 
                {
                    px[j] *= exp(- DT_PARTICLE*DAMPING);
                    py[j] *= exp(- DT_PARTICLE*DAMPING);
                }
                if ((COUPLE_ANGLE_TO_THERMOSTAT)&&(particle[j].thermostat))
                    pangle[j] *= exp(- 0.5*DT_PARTICLE*xi);
                
                particle[j].xc = qx[j] + 0.5*DT_PARTICLE*px[j]*particle[j].mass_inv;
                particle[j].yc = qy[j] + 0.5*DT_PARTICLE*py[j]*particle[j].mass_inv;
                particle[j].angle = qangle[j] + 0.5*DT_PARTICLE*pangle[j]*particle[j].inertia_moment_inv;
            }
            
        }
        
        /* reset angular values to [0, 2 Pi) */
        if (ROTATION) for (j=0; j<ncircles; j++) 
        {
            while (particle[j].angle > DPI) particle[j].angle -= DPI;
            while (particle[j].angle < 0.0) particle[j].angle += DPI;
        }
        
        /* case of periodic boundary conditions */
        if ((BOUNDARY_COND == BC_PERIODIC)||(BOUNDARY_COND == BC_PERIODIC_CIRCLE)) for (j=0; j<ncircles; j++) 
            if (particle[j].active) 
            {
                if (particle[j].xc > XMAX) particle[j].xc += XMIN - XMAX;
                else if (particle[j].xc < XMIN) particle[j].xc += XMAX - XMIN;
                if (particle[j].yc > YMAX) particle[j].yc += YMIN - YMAX;
                else if (particle[j].yc < YMIN) particle[j].yc += YMAX - YMIN;
            }
                
        printf("Mean kinetic energy: %.3f\n", totalenergy/(double)ncircles); 
        printf("Boundary force: %.3f\n", fboundary/(double)(ncircles*NVID)); 
        
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
        draw_container(xmincontainer, xmaxcontainer);

        /* add a particle */
        if ((ADD_PARTICLES)&&((i - INITIAL_TIME - ADD_TIME + 1)%ADD_PERIOD == 0))
        {
//             add_particle(XMIN + 0.1, 0.0, 50.0, 0.0, 3.0, 0, particle);
//             px[ncircles - 1] = particle[ncircles - 1].vx;
//             py[ncircles - 1] = particle[ncircles - 1].vy;
//             particle[ncircles - 1].radius = 1.5*MU;
//             j = 0;
//             while (module2(particle[j].xc,particle[j].yc) > 0.7) j = rand()%ncircles;
//             x =  particle[j].xc + 2.5*MU;
//             y =  particle[j].yc;
            
//             x =  XMIN + (XMAX - XMIN)*rand()/RAND_MAX;
//             y =  YMAX + 0.01*rand()/RAND_MAX;
//             add_particle(x, y, 0.0, 0.0, 1.0, 0, particle);

            x =  XMIN + 0.25*(XMAX - XMIN);
            y =  YMAX + 0.01;
            prop = 1.0 - (double)nadd_particle/5.0;
            vx = 100.0*prop;
            add_particle(x, y, vx, -10.0, 5.0*prop, 0, particle);
            particle[ncircles - 1].radius = 10.0*MU*prop;
            particle[ncircles - 1].eq_dist = 2.0;
            particle[ncircles - 1].thermostat = 0;
            px[ncircles - 1] = particle[ncircles - 1].vx;
            py[ncircles - 1] = particle[ncircles - 1].vy;
            nadd_particle++;
        }
        
        update_hashgrid(particle, hashgrid_number, hashgrid_particles);
        
        print_parameters(beta, totalenergy/(double)ncircles, krepel, xmaxcontainer - xmincontainer, 
                         fboundary/(double)(ncircles*NVID), 0);

	glutSwapBuffers();


	if (MOVIE)
        {
            if (i >= INITIAL_TIME) save_frame_lj();
            else printf("Initial phase time %i of %i\n", i, INITIAL_TIME);
            
            if ((i >= INITIAL_TIME)&&(DOUBLE_MOVIE))
            {
                draw_particles(particle, PLOT_B);
                draw_container(xmincontainer, xmaxcontainer);
                print_parameters(beta, totalenergy/(double)ncircles, krepel, xmaxcontainer - xmincontainer, 
                         fboundary/(double)(ncircles*NVID), 0);
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
            draw_container(xmincontainer, xmaxcontainer);
            print_parameters(beta, totalenergy/(double)ncircles, krepel, xmaxcontainer - xmincontainer, 
                         fboundary/(double)(ncircles*NVID), 0);
            glutSwapBuffers();
        }
        for (i=0; i<MID_FRAMES; i++) save_frame_lj();
        if (DOUBLE_MOVIE) 
        {
            draw_particles(particle, PLOT_B);
            draw_container(xmincontainer, xmaxcontainer);
            print_parameters(beta, totalenergy/(double)ncircles, krepel, xmaxcontainer - xmincontainer, 
                         fboundary/(double)(ncircles*NVID), 0);
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
    free(qangle);
    free(pangle);
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

