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
#include <time.h>

#define MOVIE 0         /* set to 1 to generate movie */
#define DOUBLE_MOVIE 0  /* set to 1 to produce movies for wave height and energy simultaneously */

#define TIME_LAPSE 0     /* set to 1 to add a time-lapse movie at the end */
                         /* so far incompatible with double movie */
#define TIME_LAPSE_FACTOR 3    /* factor of time-lapse movie */


/* General geometrical parameters */

#define WINWIDTH 	1280  /* window width */
#define WINHEIGHT 	720   /* window height */

#define XMIN -2.0
#define XMAX 2.0	/* x interval */
#define YMIN -1.125
#define YMAX 1.125	/* y interval for 9/16 aspect ratio */

#define INITXMIN -0.7
#define INITXMAX 0.7	/* x interval for initial condition */
#define INITYMIN 0.1
#define INITYMAX 0.6	/* y interval for initial condition */

#define BCXMIN -0.7
#define BCXMAX 0.7	/* x interval for boundary condition */
#define BCYMIN -0.85
#define BCYMAX 0.85	/* y interval for boundary condition */

#define OBSXMIN -2.0
#define OBSXMAX 2.0     /* x interval for motion of obstacle */

#define CIRCLE_PATTERN 8   /* pattern of circles, see list in global_ljones.c */

#define ADD_FIXED_OBSTACLES 0   /* set to 1 do add fixed circular obstacles */
#define OBSTACLE_PATTERN 2  /* pattern of obstacles, see list in global_ljones.c */

#define ADD_FIXED_SEGMENTS 1    /* set to 1 to add fixed segments as obstacles */
#define SEGMENT_PATTERN 2   /* pattern of repelling segments, see list in global_ljones.c */

#define TWO_TYPES 0         /* set to 1 to have two types of particles */
#define TPYE_PROPORTION 0.8 /* proportion of particles of first type */
#define SYMMETRIZE_FORCE 1  /* set to 1 to symmetrize two-particle interaction, only needed if particles are not all the same */
#define CENTER_PX 0         /* set to 1 to center horizontal momentum */
#define CENTER_PY 0         /* set to 1 to center vertical momentum */
#define CENTER_PANGLE 0     /* set to 1 to center angular momentum */

#define INTERACTION 1       /* particle interaction, see list in global_ljones.c */
#define INTERACTION_B 1     /* particle interaction for second type of particle, see list in global_ljones.c */
#define SPIN_INTER_FREQUENCY 5.0 /* angular frequency of spin-spin interaction */
#define SPIN_INTER_FREQUENCY_B 2.0 /* angular frequency of spin-spin interaction for second particle type */

#define P_PERCOL 0.25       /* probability of having a circle in C_RAND_PERCOL arrangement */
#define NPOISSON 100        /* number of points for Poisson C_RAND_POISSON arrangement */
#define PDISC_DISTANCE 2.75  /* minimal distance in Poisson disc process, controls density of particles */
#define PDISC_CANDIDATES 100 /* number of candidates in construction of Poisson disc process */
#define RANDOM_POLY_ANGLE 0 /* set to 1 to randomize angle of polygons */

#define LAMBDA 2.0	    /* parameter controlling the dimensions of domain */
#define MU 0.013  	    /* parameter controlling radius of particles */
#define MU_B 0.0254         /* parameter controlling radius of particles of second type */
#define NPOLY 3             /* number of sides of polygon */
#define APOLY 0.666666666   /* angle by which to turn polygon, in units of Pi/2 */ 
#define MDEPTH 4            /* depth of computation of Menger gasket */
#define MRATIO 3            /* ratio defining Menger gasket */
#define MANDELLEVEL 1000    /* iteration level for Mandelbrot set */
#define MANDELLIMIT 10.0    /* limit value for approximation of Mandelbrot set */
#define FOCI 1              /* set to 1 to draw focal points of ellipse */
#define NGRIDX 46           /* number of grid point for grid of disks */
#define NGRIDY 24           /* number of grid point for grid of disks */
#define EHRENFEST_RADIUS 0.9    /* radius of container for Ehrenfest urn configuration */
#define EHRENFEST_WIDTH 0.035     /* width of tube for Ehrenfest urn configuration */

#define X_SHOOTER -0.2
#define Y_SHOOTER -0.6
#define X_TARGET 0.4
#define Y_TARGET 0.7        /* shooter and target positions in laser fight */

/* Parameters for length and speed of simulation */

#define NSTEPS 5500      /* number of frames of movie */
#define NVID 650         /* number of iterations between images displayed on screen */
#define NSEG 250         /* number of segments of boundary */
#define INITIAL_TIME 10    /* time after which to start saving frames */
#define BOUNDARY_WIDTH 1    /* width of particle boundary */
#define LINK_WIDTH 2        /* width of links between particles */
#define CONTAINER_WIDTH 4   /* width of container boundary */

#define PAUSE 1000         /* number of frames after which to pause */
#define PSLEEP 1         /* sleep time during pause */
#define SLEEP1  1        /* initial sleeping time */
#define SLEEP2  1   /* final sleeping time */
#define MID_FRAMES 20    /* number of still frames between parts of two-part movie */
#define END_FRAMES 100   /* number of still frames at end of movie */

/* Boundary conditions, see list in global_ljones.c */

#define BOUNDARY_COND 0

/* Plot type, see list in global_ljones.c  */

#define PLOT 0
#define PLOT_B 8        /* plot type for second movie */

#define COLOR_BONDS 1   /* set to 1 to color bonds according to length */

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

#define ENERGY_HUE_MIN 330.0      /* color of original particle */
#define ENERGY_HUE_MAX 50.0        /* color of saturated particle */
#define PARTICLE_HUE_MIN 359.0      /* color of original particle */
#define PARTICLE_HUE_MAX 0.0        /* color of saturated particle */
#define PARTICLE_EMAX 2.0e2           /* energy of particle with hottest color */
#define HUE_TYPE0 280.0     /* hue of particles of type 0 */
#define HUE_TYPE1 135.0      /* hue of particles of type 1 */
#define HUE_TYPE2 70.0      /* hue of particles of type 1 */
#define HUE_TYPE3 210.0      /* hue of particles of type 1 */

#define RANDOM_RADIUS 0     /* set to 1 for random circle radius */
#define DT_PARTICLE 5.0e-7    /* time step for particle displacement */
#define KREPEL 12.0          /* constant in repelling force between particles */
#define EQUILIBRIUM_DIST 4.0    /* Lennard-Jones equilibrium distance */
#define EQUILIBRIUM_DIST_B 5.0  /* Lennard-Jones equilibrium distance for second type of particle */
#define REPEL_RADIUS 20.0    /* radius in which repelling force acts (in units of particle radius) */
#define DAMPING 0.0     /* damping coefficient of particles */
#define PARTICLE_MASS 1.0   /* mass of particle of radius MU */
#define PARTICLE_MASS_B 1.0   /* mass of particle of radius MU */
#define PARTICLE_INERTIA_MOMENT 0.2     /* moment of inertia of particle */
#define PARTICLE_INERTIA_MOMENT_B 0.02     /* moment of inertia of second type of particle */
#define V_INITIAL 10.0        /* initial velocity range */
#define OMEGA_INITIAL 10.0        /* initial angular velocity range */

#define THERMOSTAT 1        /* set to 1 to switch on thermostat */
#define VARY_THERMOSTAT 0   /* set to 1 for time-dependent thermostat schedule */
#define SIGMA 5.0           /* noise intensity in thermostat */
#define BETA 0.01           /* initial inverse temperature */
#define MU_XI 0.01           /* friction constant in thermostat */
#define KSPRING_BOUNDARY 5.0e9    /* confining harmonic potential outside simulation region */
#define KSPRING_OBSTACLE 5.0e9    /* harmonic potential of obstacles */
#define NBH_DIST_FACTOR 4.5       /* radius in which to count neighbours */
#define GRAVITY 10000.0            /* gravity acting on all particles */
#define INCREASE_GRAVITY 0     /* set to 1 to increase gravity during the simulation */
#define GRAVITY_FACTOR 100.0   /* factor by which to increase gravity */
#define GRAVITY_INITIAL_TIME 500    /* time at start of simulation with constant gravity */
#define GRAVITY_RESTORE_TIME 1000    /* time at end of simulation with gravity restored to initial value */

#define ROTATION 0           /* set to 1 to include rotation of particles */
#define COUPLE_ANGLE_TO_THERMOSTAT 0    /* set to 1 to couple angular degrees of freedom to thermostat */
#define DIMENSION_FACTOR 1.0  /* scaling factor taking into account number of degrees of freedom */  
#define KTORQUE 50.0          /* force constant in angular dynamics */
#define KTORQUE_B 10.0        /* force constant in angular dynamics */
#define KTORQUE_DIFF 150.0    /* force constant in angular dynamics for different particles */
#define DRAW_SPIN 0           /* set to 1 to draw spin vectors of particles */
#define DRAW_SPIN_B 0         /* set to 1 to draw spin vectors of particles */
#define DRAW_CROSS 1          /* set to 1 to draw cross on particles of second type */
#define SPIN_RANGE 7.0       /* range of spin-spin interaction */
#define SPIN_RANGE_B 5.0     /* range of spin-spin interaction for second type of particle */
#define QUADRUPOLE_RATIO 0.6  /* anisotropy in quadrupole potential */ 

#define INCREASE_BETA 0  /* set to 1 to increase BETA during simulation */
#define BETA_FACTOR 20.0   /* factor by which to change BETA during simulation */
#define N_TOSCILLATIONS 1.5   /* number of temperature oscillations in BETA schedule */
#define NO_OSCILLATION 1        /* set to 1 to have exponential BETA change only */
#define FINAL_CONSTANT_PHASE 2000  /* final phase in which temperature is constant */

#define DECREASE_CONTAINER_SIZE 0   /* set to 1 to decrease size of container */
#define SYMMETRIC_DECREASE 0        /* set tp 1 to decrease container symmetrically */
#define COMPRESSION_RATIO 0.3       /* final size of container */
#define RESTORE_CONTAINER_SIZE 1    /* set to 1 to restore container to initial size at end of simulation */
#define RESTORE_TIME 700            /* time before end of sim at which to restore size */

#define MOVE_OBSTACLE 0     /* set to 1 to have a moving obstacle */
#define CENTER_VIEW_ON_OBSTACLE 0   /* set to 1 to center display on moving obstacle */
#define RESAMPLE_Y 0         /* set to 1 to resample y coordinate of moved particles (for shock waves) */
#define NTRIALS 2000         /* number of trials when resampling */
#define OBSTACLE_RADIUS 0.12  /* radius of obstacle for circle boundary conditions */
#define FUNNEL_WIDTH  0.25   /* funnel width for funnel boundary conditions */
#define OBSTACLE_XMIN 0.0    /* initial position of obstacle */
#define OBSTACLE_XMAX 3.0    /* final position of obstacle */
#define RECORD_PRESSURES 0   /* set to 1 to record pressures on obstacle */
#define N_PRESSURES 100      /* number of intervals to record pressure */
#define N_P_AVERAGE 100      /* size of pressure averaging window */
#define N_T_AVERAGE 50       /* size of temperature averaging window */
#define MAX_PRESSURE 3.0e10  /* pressure shown in "hottest" color */
#define PARTIAL_THERMO_COUPLING 0   /* set to 1 to couple only particles to the right of obstacle to thermostat */
#define PARTIAL_THERMO_REGION 1     /* region for partial thermostat coupling (see list in global_ljones.c) */
#define PARTIAL_THERMO_SHIFT 0.5    /* distance from obstacle at the right of which particles are coupled to thermostat */

#define INCREASE_KREPEL 0   /* set to 1 to increase KREPEL during simulation */
#define KREPEL_FACTOR 1000.0   /* factor by which to change KREPEL during simulation */

#define PART_AT_BOTTOM 0     /* set to 1 to include "seed" particles at bottom */
#define MASS_PART_BOTTOM 10000.0 /* mass of particles at bottom */
#define NPART_BOTTOM 100     /* number of particles at the bottom */

#define ADD_PARTICLES 0    /* set to 1 to add particles */
#define ADD_TIME 500       /* time at which to add first particle */
#define ADD_PERIOD 250       /* time interval between adding further particles */
#define FINAL_NOADD_PERIOD 200  /* final period where no particles are added */
#define SAFETY_FACTOR 2.0  /* no particles are added at distance less than MU*SAFETY_FACTOR of other particles */

#define TRACER_PARTICLE 0   /* set to 1 to have a tracer particle */
#define N_TRACER_PARTICLES 3    /* number of tracer particles */
#define TRAJECTORY_LENGTH 8000   /* length of recorded trajectory */
#define TRACER_PARTICLE_MASS 4.0    /* relative mass of tracer particle */
#define TRAJECTORY_WIDTH 3      /* width of tracer particle trajectory */

#define ROTATE_BOUNDARY 1           /* set to 1 to rotate the repelling segments */
#define SMOOTH_ROTATION 1           /* set to 1 to update segments at each time step (rather than at each movie frame) */
#define PERIOD_ROTATE_BOUNDARY 2500 /* period of rotating boundary */
#define ROTATE_INITIAL_TIME 0       /* initial time without rotation */
#define ROTATE_FINAL_TIME 500       /* final time without rotation */
#define PRINT_OMEGA 0               /* set to 1 to print angular speed */
#define PRINT_PARTICLE_SPEEDS 0     /* set to 1 to print average speeds/momenta of particles */

#define MOVE_BOUNDARY 0        /* set to 1 to move repelling segments, due to force from particles */
#define SEGMENTS_MASS 100.0     /* mass of collection of segments */
#define DEACTIVATE_SEGMENT 1    /* set to 1 to deactivate last segment after a certain time */
#define SEGMENT_DEACTIVATION_TIME 1000   /* time at which to deactivate last segment */

#define POSITION_DEPENDENT_TYPE 0   /* set to 1 to make particle type depend on initial position */
#define POSITION_Y_DEPENDENCE 0     /* set to 1 for the separation between particles to be vertical */
#define PRINT_ENTROPY 0     /* set to 1 to compute entropy */

#define REACTION_DIFFUSION 0    /* set to 1 to simulate a chemical reaction (particles may change type) */
#define RD_TYPES 3              /* number of types in reaction-diffusion equation */
#define REACTION_PROB 0.03      /* probability controlling reaction term */ 

#define PRINT_PARTICLE_NUMBER 0     /* set to 1 to print total number of particles */
#define PRINT_LEFT 0        /* set to 1 to print certain parameters at the top left instead of right */

#define EHRENFEST_COPY 0    /* set to 1 to add equal number of larger particles (for Ehrenfest model) */

#define LID_MASS 1000.0     /* mass of lid for BC_RECTANGLE_LID b.c. */
#define LID_WIDTH 0.1       /* width of lid for BC_RECTANGLE_LID b.c. */
#define WALL_MASS 2000.0    /* mass of wall for BC_RECTANGLE_WALL b.c. */
#define WALL_FRICTION 0.0   /* friction on wall for BC_RECTANGLE_WALL b.c. */
#define WALL_WIDTH 0.1      /* width of wall for BC_RECTANGLE_WALL b.c. */
#define WALL_VMAX 100.0     /* max speed of wall */
#define WALL_TIME 500       /* time during which to keep wall */

#define FLOOR_FORCE 1      /* set to 1 to limit force on particle to FMAX */
#define FMAX 1.0e9         /* maximal force */
#define FLOOR_OMEGA 1      /* set to 1 to limit particle momentum to PMAX */
#define PMAX 1000.0        /* maximal force */

#define HASHX 40   /* size of hashgrid in x direction */
#define HASHY 20   /* size of hashgrid in y direction */
#define HASHMAX 100  /* maximal number of particles per hashgrid cell */
#define HASHGRID_PADDING 0.1    /* padding of hashgrid outside simulation window */

#define DRAW_COLOR_SCHEME 0     /* set to 1 to plot the color scheme */
#define COLORBAR_RANGE 8.0    /* scale of color scheme bar */
#define COLORBAR_RANGE_B 12.0    /* scale of color scheme bar for 2nd part */
#define ROTATE_COLOR_SCHEME 0   /* set to 1 to draw color scheme horizontally */

#define NO_WRAP_BC ((BOUNDARY_COND != BC_PERIODIC)&&(BOUNDARY_COND != BC_PERIODIC_CIRCLE)&&(BOUNDARY_COND != BC_PERIODIC_TRIANGLE)&&(BOUNDARY_COND != BC_KLEIN)&&(BOUNDARY_COND != BC_PERIODIC_FUNNEL)&&(BOUNDARY_COND != BC_BOY)&&(BOUNDARY_COND != BC_GENUS_TWO))
#define PERIODIC_BC ((BOUNDARY_COND == BC_PERIODIC)||(BOUNDARY_COND == BC_PERIODIC_CIRCLE)||(BOUNDARY_COND == BC_PERIODIC_FUNNEL)||(BOUNDARY_COND == BC_PERIODIC_TRIANGLE))

double xshift = 0.0;      /* x shift of shown window */
double xspeed = 0.0;      /* x speed of obstacle */
double ylid = 0.9;        /* y coordinate of lid (for BC_RECTANGLE_LID b.c.) */
double vylid = 0.0;       /* y speed coordinate of lid (for BC_RECTANGLE_LID b.c.) */
double xwall = 0.0;       /* x coordinate of wall (for BC_RECTANGLE_WALL b.c.) */
double vxwall = 0.0;      /* x speed of wall (for BC_RECTANGLE_WALL b.c.) */
double angular_speed = 0.0;    /* angular speed of rotating segments */
double xsegments = 0.0;        /* x coordinate of segments (for option MOVE_BOUNDARY) */
double ysegments = 0.0;        /* y coordinate of segments (for option MOVE_BOUNDARY) */
double vxsegments = 0.0;       /* vx coordinate of segments (for option MOVE_BOUNDARY) */
double vysegments = 0.0;       /* vy coordinate of segments (for option MOVE_BOUNDARY) */
int thermostat_on = 1;    /* thermostat switch used when VARY_THERMOSTAT is on */

#define THERMOSTAT_ON ((THERMOSTAT)&&((!VARY_THERMOSTAT)||(thermostat_on)))

#include "global_ljones.c"
#include "sub_lj.c"
#include "sub_hashgrid.c"


/*********************/
/* animation part    */
/*********************/


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
        bexponent = log(BETA_FACTOR)/(double)(NSTEPS - FINAL_CONSTANT_PHASE);
        omega = N_TOSCILLATIONS*DPI/(double)(NSTEPS - FINAL_CONSTANT_PHASE);
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
    double time, acceleration = 40.0;
    double x;
//     static double t1 = 0.5, t2 = 0.75, t3 = 0.875;
    
    if (i < INITIAL_TIME) return(OBSTACLE_XMIN);
    else 
    {
        time = (double)(i-INITIAL_TIME)/(double)(NSTEPS);
        
        x = OBSTACLE_XMIN + 0.5*acceleration*time*time;
        xspeed = acceleration*time;
//         while (x > OBSXMAX) x += OBSXMIN - OBSXMAX;
        return(x);
    }
}

double obstacle_schedule_smooth(int i, int j)
{
    double time, acceleration = 50.0;
    double x;
//     static double t1 = 0.5, t2 = 0.75, t3 = 0.875;
    
    if (i < INITIAL_TIME) return(OBSTACLE_XMIN);
    else 
    {
        time = ((double)(i-INITIAL_TIME) + (double)j/(double)NVID)/(double)(NSTEPS);
        
        x = OBSTACLE_XMIN + 0.5*acceleration*time*time;
        xspeed = acceleration*time;
//         while (x > OBSXMAX) x += OBSXMIN - OBSXMAX;
        return(x);
    }
}

double gravity_schedule(int i, int j)
{
    double time, gravity;
    
    if ((i < INITIAL_TIME + GRAVITY_INITIAL_TIME)||(i > NSTEPS + INITIAL_TIME - GRAVITY_RESTORE_TIME)) return(GRAVITY);
    else 
    {
        time = ((double)(i - INITIAL_TIME - GRAVITY_INITIAL_TIME) 
        + (double)j/(double)NVID)/(double)(NSTEPS - GRAVITY_RESTORE_TIME - GRAVITY_INITIAL_TIME);
        gravity = GRAVITY*(1.0 + time*(GRAVITY_FACTOR - 1.0));
//         printf("i = %i, time = %.3lg, Gravity = %.3lg\n", i, time, gravity); 
        return(gravity);
    }
}

double rotation_angle(double phase)
{
    double omegamax = 15.0;
    
    /* case of rotating hourglass */
    while (phase > DPI) phase -= DPI;
    return(phase - 0.5*sin(2.0*phase));
    
    /* case of centrifuge */
//     while (phase > DPI) phase -= DPI;
//     angular_speed = 0.5*omegamax*(1.0 - cos(phase));
//     return(0.5*omegamax*(phase - sin(phase)));
}

double rotation_schedule(int i)
{
    double phase;
    static int imin = INITIAL_TIME + ROTATE_INITIAL_TIME, imax = INITIAL_TIME + NSTEPS - ROTATE_FINAL_TIME;
    
    if (i < imin) 
    {
        angular_speed = 0.0;
        return(0.0);
    }
    else
    {
        if (i > imax) i = imax;
        phase = (DPI/(double)PERIOD_ROTATE_BOUNDARY)*(double)(i - imin);
        return(rotation_angle(phase));
    }
}

double rotation_schedule_smooth(int i, int j)
{
    double phase;
    static int imin = INITIAL_TIME + ROTATE_INITIAL_TIME, imax = INITIAL_TIME + NSTEPS - ROTATE_FINAL_TIME;
    
    if (i < imin) 
    {
        angular_speed = 0.0;
        return(0.0);
    }
    else
    {
        if (i > imax) i = imax;
        phase = (DPI/(double)PERIOD_ROTATE_BOUNDARY)*((double)(i - imin) + (double)j/(double)NVID);
        return(rotation_angle(phase));
    }
}

int thermostat_schedule(int i)
{
    if (i < INITIAL_TIME) return(1);
    else return(0);
}

double evolve_particles(t_particle particle[NMAXCIRCLES], t_hashgrid hashgrid[HASHX*HASHY], 
                        double qx[NMAXCIRCLES], double qy[NMAXCIRCLES], double qangle[NMAXCIRCLES],
                        double px[NMAXCIRCLES], double py[NMAXCIRCLES], double pangle[NMAXCIRCLES], 
                        double beta, int *nactive, int *nsuccess, int *nmove)
{
    double a, totalenergy = 0.0;  
    static double b = 0.25*SIGMA*SIGMA*DT_PARTICLE/MU_XI, xi = 0.0;
    int j, move;
    
    #pragma omp parallel for private(j,xi,totalenergy,a,move)
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
        
        if ((THERMOSTAT_ON)&&(particle[j].thermostat))
        {
            px[j] *= exp(- 0.5*DT_PARTICLE*xi);
            py[j] *= exp(- 0.5*DT_PARTICLE*xi);
        }
        if ((COUPLE_ANGLE_TO_THERMOSTAT)&&(particle[j].thermostat)) 
            pangle[j] *= exp(- 0.5*DT_PARTICLE*xi);
    }
                
    /* compute kinetic energy */
//     *nactive = 0;
    for (j=0; j<ncircles; j++)
        if ((particle[j].active)&&(particle[j].thermostat))
        {
            totalenergy += particle[j].energy;
//             *nactive++;
        }
    totalenergy *= DIMENSION_FACTOR;    /* normalize energy to take number of degrees of freedom into account */
    if (THERMOSTAT_ON)
    {
        a = DT_PARTICLE*(totalenergy - (double)*nactive/beta)/MU_XI;
        a += SIGMA*sqrt(DT_PARTICLE)*gaussian();
        xi = (xi + a - b*xi)/(1.0 + b);
    }
    
    move = 0;
    for (j=0; j<ncircles; j++) if (particle[j].active) 
    {
        if ((THERMOSTAT_ON)&&(particle[j].thermostat))
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
        
//         particle[j].vx = px[j] + 0.5*DT_PARTICLE*particle[j].fx;
//         particle[j].vy = py[j] + 0.5*DT_PARTICLE*particle[j].fy;
//         particle[j].omega = pangle[j] + 0.5*DT_PARTICLE*particle[j].torque;

        /* TO DO: move this to function wrap_particle */
        if ((BOUNDARY_COND == BC_PERIODIC_CIRCLE)||(BOUNDARY_COND == BC_PERIODIC_FUNNEL)||(BOUNDARY_COND == BC_PERIODIC_TRIANGLE))
        {
//                       if (particle[j].xc < BCXMIN) 
            if (particle[j].xc < xshift + BCXMIN) 
            {
                particle[j].xc += BCXMAX - BCXMIN;
                if (RESAMPLE_Y) 
                {
                    *nmove++;
                    if (resample_particle(j, NTRIALS, particle) == 1) 
                    {
                        px[j] = particle[j].vx;
                        py[j] = particle[j].vy;
                        update_hashgrid(particle, hashgrid, 0);
                        *nsuccess++;
                    }
                }
            }
//                       else if (particle[j].xc > BCXMAX) 
            else if (particle[j].xc > xshift + BCXMAX) 
            {
                particle[j].xc += BCXMIN - BCXMAX;
            }
            if (particle[j].yc > BCYMAX) particle[j].yc += BCYMIN - BCYMAX;
            else if (particle[j].yc < BCYMIN) particle[j].yc += BCYMAX - BCYMIN;
        }
        else if (!NO_WRAP_BC)
        {
            move += wrap_particle(&particle[j], &px[j], &py[j]);
        }
//                 if (move > 0) 
//                 {
//                     compute_relative_positions(particle, hashgrid);
//                     update_hashgrid(particle, hashgrid, 0);   /* REDUNDANT ? */
//                 }
    }
    
    return(totalenergy);
}


void evolve_lid(double fboundary)
{
    double force;
    
    force = fboundary - GRAVITY*LID_MASS;
    if (ylid > BCYMAX + LID_WIDTH) force -= KSPRING_BOUNDARY*(ylid - BCYMAX - LID_WIDTH);
    vylid += force*DT_PARTICLE/LID_MASS;
    ylid += vylid*DT_PARTICLE;
}

void evolve_wall(double fboundary)
{
    double force;
    
    force = fboundary;
    if (xwall > BCYMAX - WALL_WIDTH) force -= KSPRING_BOUNDARY*(xwall - BCYMAX + WALL_WIDTH);
    else if (xwall < BCYMIN + WALL_WIDTH) force += KSPRING_BOUNDARY*(BCYMIN + WALL_WIDTH - xwall);
    
    force -= vxwall*WALL_FRICTION;
    
    vxwall += fboundary*DT_PARTICLE/WALL_MASS;
    
    if (vxwall > WALL_VMAX) vxwall = WALL_VMAX;
    else if (vxwall < -WALL_VMAX) vxwall = -WALL_VMAX;
    
    xwall += vxwall*DT_PARTICLE;
//     printf("fboundary = %.3lg, xwall = %.3lg, vxwall = %.3lg\n", fboundary, xwall, vxwall);
}


void evolve_segments(t_segment segment[NMAXSEGMENTS])
{
    int i, nactive = 0;
    double fx = 0.0, fy = 0.0;
    
    for (i=0; i<nsegments; i++) if (segment[i].active)
    {
        fx += segment[i].fx;
        fy += segment[i].fy;
        nactive++;
    }
    if (nactive > 0)
    {
        fx = fx/(double)nactive;
        fy = fy/(double)nactive;
    }
    if (FLOOR_FORCE)
    {   
        if (fx > FMAX) fx = FMAX;
        else if (fx < -FMAX) fx = -FMAX;
        if (fy > FMAX) fy = FMAX;
        else if (fy < -FMAX) fy = -FMAX;
    }
    vxsegments += fx*DT_PARTICLE/(SEGMENTS_MASS);
    vysegments += fy*DT_PARTICLE/(SEGMENTS_MASS);
    xsegments += vxsegments*DT_PARTICLE;
    ysegments += vysegments*DT_PARTICLE;

    /* to avoid numerical instabilities */
    if (xsegments + 1.0 > BCXMAX) 
    {
        xsegments = BCXMAX - 1.0;
        vxsegments = 0.0;
    }
        
    translate_segments(segment, xsegments, ysegments);
}


void animation()
{
    double time, scale, diss, rgb[3], dissip, gradient[2], x, y, dx, dy, dt, xleft, xright, a, b, 
            length, fx, fy, force[2], totalenergy = 0.0, krepel = KREPEL, pos[2], prop, vx, 
            beta = BETA, xi = 0.0, xmincontainer = BCXMIN, xmaxcontainer = BCXMAX, torque, torque_ij, 
            fboundary = 0.0, pleft = 0.0, pright = 0.0, entropy[2], mean_energy, gravity = GRAVITY;
    double *qx, *qy, *px, *py, *qangle, *pangle, *pressure;
    int i, j, k, n, m, s, ij[2], i0, iplus, iminus, j0, jplus, jminus, p, q, p1, q1, p2, q2, total_neighbours = 0, 
        min_nb, max_nb, close, wrapx = 0, wrapy = 0, nactive = 0, nadd_particle = 0, nmove = 0, nsuccess = 0, 
        tracer_n[N_TRACER_PARTICLES], traj_position = 0, traj_length = 0, move = 0, old, m0, floor, nthermo, wall = 0;
    static int imin, imax;
    static short int first = 1;
    t_particle *particle;
    t_obstacle *obstacle;
    t_segment *segment;
    t_tracer *trajectory;
    t_hashgrid *hashgrid;
    char message[100];


    particle = (t_particle *)malloc(NMAXCIRCLES*sizeof(t_particle));    /* particles */  
    if (ADD_FIXED_OBSTACLES) obstacle = (t_obstacle *)malloc(NMAXOBSTACLES*sizeof(t_obstacle));    /* circular obstacles */  
    if (ADD_FIXED_SEGMENTS) segment = (t_segment *)malloc(NMAXSEGMENTS*sizeof(t_segment));         /* linear obstacles */  

    if (TRACER_PARTICLE) trajectory = (t_tracer *)malloc(TRAJECTORY_LENGTH*N_TRACER_PARTICLES*sizeof(t_tracer));

    hashgrid = (t_hashgrid *)malloc(HASHX*HASHY*sizeof(t_hashgrid));    /* hashgrid */      
            
    qx = (double *)malloc(NMAXCIRCLES*sizeof(double));
    qy = (double *)malloc(NMAXCIRCLES*sizeof(double));
    px = (double *)malloc(NMAXCIRCLES*sizeof(double));
    py = (double *)malloc(NMAXCIRCLES*sizeof(double));
    qangle = (double *)malloc(NMAXCIRCLES*sizeof(double));
    pangle = (double *)malloc(NMAXCIRCLES*sizeof(double));
    pressure = (double *)malloc(N_PRESSURES*sizeof(double));
    
    /* initialise positions and radii of circles */
    init_particle_config(particle);
    init_hashgrid(hashgrid);
    
    xshift = OBSTACLE_XMIN;
    
    if (ADD_FIXED_OBSTACLES) init_obstacle_config(obstacle);
    if (ADD_FIXED_SEGMENTS) init_segment_config(segment);
    
    if (RECORD_PRESSURES) for (i=0; i<N_PRESSURES; i++) pressure[i] = 0.0;

//     printf("1\n");

    nactive = initialize_configuration(particle, hashgrid, obstacle, px, py, pangle, tracer_n);
     
//     xi = 0.0;
    
//     for (i=0; i<ncircles; i++)
//     {
//         printf("Particle %i at (%.3f, %.3f) of energy %.3f\n", i, particle[i].xc, particle[i].yc, particle[i].energy);
//     }
    sleep(1);
        
    update_hashgrid(particle, hashgrid, 1);
    compute_relative_positions(particle, hashgrid);
    
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
        if ((ROTATE_BOUNDARY)&&(!SMOOTH_ROTATION)) rotate_segments(segment, rotation_schedule(i));
        if (VARY_THERMOSTAT) 
        {
            thermostat_on = thermostat_schedule(i);
            printf("Termostat: %i\n", thermostat_on);
        }
        if ((DEACTIVATE_SEGMENT)&&(i > INITIAL_TIME + SEGMENT_DEACTIVATION_TIME))
            segment[nsegments-1].active = 0;
        
        blank();
        
        fboundary = 0.0;
        pleft = 0.0;
        pright = 0.0;
        if (RECORD_PRESSURES) for (j=0; j<N_PRESSURES; j++) pressure[j] = 0.0;
        
        for(n=0; n<NVID; n++)
        {
            if (MOVE_OBSTACLE) 
            {
                xmincontainer = obstacle_schedule_smooth(i, n);
                xshift = xmincontainer;
            }
            if ((ROTATE_BOUNDARY)&&(SMOOTH_ROTATION)) rotate_segments(segment, rotation_schedule_smooth(i,n));
            
            if (INCREASE_GRAVITY) gravity = gravity_schedule(i,n); 
            if ((BOUNDARY_COND == BC_RECTANGLE_WALL)&&(i < INITIAL_TIME + WALL_TIME)) wall = 1;
            else wall = 0;
            
            if (MOVE_BOUNDARY) for (j=0; j<nsegments; j++)
            {
                segment[j].fx = 0.0;
                segment[j].fy = 0.0;
            }
            
            compute_relative_positions(particle, hashgrid);
            update_hashgrid(particle, hashgrid, 0);
            
            /* compute forces on particles */
            for (j=0; j<ncircles; j++) if (particle[j].active)
            {
                particle[j].fx = 0.0;
                particle[j].fy = 0.0;
                particle[j].torque = 0.0;
                
                /* compute force from other particles */
                compute_particle_force(j, krepel, particle, hashgrid);

                /* take care of boundary conditions */
                fboundary += compute_boundary_force(j, particle, obstacle, segment, xmincontainer, xmaxcontainer, &pleft, &pright, pressure, wall);

                /* add gravity */
                if (INCREASE_GRAVITY) particle[j].fy -= gravity;
                else particle[j].fy -= GRAVITY;
                                                
                if (FLOOR_FORCE)
                {
                    if (particle[j].fx > FMAX) particle[j].fx = FMAX;
                    if (particle[j].fx < -FMAX) particle[j].fx = -FMAX;
                    if (particle[j].fy > FMAX) particle[j].fy = FMAX;
                    if (particle[j].fy < -FMAX) particle[j].fy = -FMAX;
                    if (particle[j].torque > FMAX) particle[j].torque = FMAX;
                    if (particle[j].torque < -FMAX) particle[j].torque = -FMAX;
                }
            }
            
            /* timestep of thermostat algorithm */
            totalenergy = evolve_particles(particle, hashgrid, qx, qy, qangle, px, py, pangle, beta, &nactive, &nsuccess, &nmove);
            
            /* evolution of lid coordinate */
            if (BOUNDARY_COND == BC_RECTANGLE_LID) evolve_lid(fboundary);
            if (BOUNDARY_COND == BC_RECTANGLE_WALL) 
            {
                if (i < INITIAL_TIME + WALL_TIME) evolve_wall(fboundary);
                else xwall = 0.0;
            }
            if (MOVE_BOUNDARY) evolve_segments(segment);
        } /* end of for (n=0; n<NVID; n++) */
        
        if (MOVE_BOUNDARY) printf("segments position (%.3lg, %.3lg), speed (%.3lg, %.3lg)\n", xsegments, ysegments, vxsegments, vysegments);
        
//         if ((PARTIAL_THERMO_COUPLING)) 
        if ((PARTIAL_THERMO_COUPLING)&&(i>N_T_AVERAGE)) 
        {
            nthermo = partial_thermostat_coupling(particle, xshift + PARTIAL_THERMO_SHIFT);
            printf("%i particles coupled to thermostat out of %i active\n", nthermo, nactive);
            mean_energy = compute_mean_energy(particle);
        }
        else mean_energy = totalenergy/(double)ncircles;
        
        if (CENTER_PX) center_momentum(px);
        if (CENTER_PY) center_momentum(py);
        if (CENTER_PANGLE) center_momentum(pangle);
        if (FLOOR_OMEGA) floor = floor_momentum(pangle);
            
//         printf("pressure left %.5lg, pressure right %.5lg\n", pleft, pright);
//         for (j=0; j<N_PRESSURES; j++) printf("pressure[%i] = %.5lg\n", j, pressure[j]);
        
        /* reset angular values to [0, 2 Pi) */
        if (ROTATION) for (j=0; j<ncircles; j++) 
        {
            while (particle[j].angle > DPI) particle[j].angle -= DPI;
            while (particle[j].angle < 0.0) particle[j].angle += DPI;
        }
        
            
        /* update tracer particle trajectory */
        if ((TRACER_PARTICLE)&&(i > INITIAL_TIME)) 
        {
            for (j=0; j<N_TRACER_PARTICLES; j++)
            {
                trajectory[j*TRAJECTORY_LENGTH + traj_position].xc = particle[tracer_n[j]].xc;
                trajectory[j*TRAJECTORY_LENGTH + traj_position].yc = particle[tracer_n[j]].yc;
            }
            traj_position++;
            if (traj_position >= TRAJECTORY_LENGTH) traj_position = 0;
            traj_length++;
            if (traj_length >= TRAJECTORY_LENGTH) traj_length = TRAJECTORY_LENGTH - 1;
//             for (j=0; j<traj_length; j++)
//                 printf("Trajectory[%i] = (%.3lg, %.3lg)\n", j, trajectory[j].xc, trajectory[j].yc);
        }
                
        printf("Mean kinetic energy: %.3f\n", totalenergy/(double)ncircles); 
        printf("Boundary force: %.3f\n", fboundary/(double)(ncircles*NVID)); 
        if (RESAMPLE_Y) printf("%i succesful moves out of %i trials\n", nsuccess, nmove);
        if (INCREASE_GRAVITY) printf("Gravity: %.3f\n", gravity);
        
        total_neighbours = 0;
        min_nb = 100;
        max_nb = 0;
        for (j=0; j<ncircles; j++) if (particle[j].active)
        {
            total_neighbours += particle[j].neighb;
            if (particle[j].neighb > max_nb) max_nb = particle[j].neighb; 
            if (particle[j].neighb < min_nb) min_nb = particle[j].neighb; 
        }
//         printf("Mean number of neighbours: %.3f\n", (double)total_neighbours/(double)ncircles); 
//         printf("Min number of neighbours: %i\n", min_nb); 
//         printf("Max number of neighbours: %i\n", max_nb); 
        
        if (TRACER_PARTICLE) draw_trajectory(trajectory, traj_position, traj_length);
        draw_particles(particle, PLOT);
        draw_container(xmincontainer, xmaxcontainer, obstacle, segment, wall);

        /* add a particle */
        if ((ADD_PARTICLES)&&(i > ADD_TIME)&&((i - INITIAL_TIME - ADD_TIME + 1)%ADD_PERIOD == 0)&&(i < NSTEPS - FINAL_NOADD_PERIOD))
            nadd_particle = add_particles(particle, px, py, nadd_particle);
        
        /* case of reaction-diffusion equation */
        if (REACTION_DIFFUSION) update_types(particle); 
 
        update_hashgrid(particle, hashgrid, 1);
        
        print_parameters(beta, mean_energy, krepel, xmaxcontainer - xmincontainer, 
                         fboundary/(double)(ncircles*NVID), PRINT_LEFT, pressure, gravity);
        if ((BOUNDARY_COND == BC_EHRENFEST)||(BOUNDARY_COND == BC_RECTANGLE_WALL)) 
            print_ehrenfest_parameters(particle, pleft, pright);
        else if (PRINT_PARTICLE_NUMBER) print_particle_number(ncircles);
        
        if ((i > INITIAL_TIME + WALL_TIME)&&(PRINT_ENTROPY)) 
        {
            compute_entropy(particle, entropy);
            printf("Entropy 1 = %.5lg, Entropy 2 = %.5lg\n", entropy[0], entropy[1]);
            print_entropy(entropy);
        }
        
        if (PRINT_OMEGA) print_omega(angular_speed); 
        else if (PRINT_PARTICLE_SPEEDS) print_particles_speeds(particle);

	glutSwapBuffers();


	if (MOVIE)
        {
            if (i >= INITIAL_TIME) 
            {
                save_frame_lj();
                if ((TIME_LAPSE)&&((i - INITIAL_TIME)%TIME_LAPSE_FACTOR == 0)&&(!DOUBLE_MOVIE))
                {
                    save_frame_lj_counter(NSTEPS + END_FRAMES + (i - INITIAL_TIME)/TIME_LAPSE_FACTOR);
                }
            }
            else printf("Initial phase time %i of %i\n", i, INITIAL_TIME);
            
            
            if ((i >= INITIAL_TIME)&&(DOUBLE_MOVIE))
            {
                if (TRACER_PARTICLE) draw_trajectory(trajectory, traj_position, traj_length);
                draw_particles(particle, PLOT_B);
                draw_container(xmincontainer, xmaxcontainer, obstacle, segment, wall);
                print_parameters(beta, mean_energy, krepel, xmaxcontainer - xmincontainer, 
                         fboundary/(double)(ncircles*NVID), PRINT_LEFT, pressure, gravity);
                if (BOUNDARY_COND == BC_EHRENFEST) print_ehrenfest_parameters(particle, pleft, pright);
                else if (PRINT_PARTICLE_NUMBER) print_particle_number(ncircles);
                if (PRINT_OMEGA) print_omega(angular_speed); 
                else if (PRINT_PARTICLE_SPEEDS) print_particles_speeds(particle);
                glutSwapBuffers();
                save_frame_lj_counter(NSTEPS + MID_FRAMES + 1 + counter);
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
            if (TRACER_PARTICLE) draw_trajectory(trajectory, traj_position, traj_length);
            draw_particles(particle, PLOT); 
            draw_container(xmincontainer, xmaxcontainer, obstacle, segment, wall);
            print_parameters(beta, mean_energy, krepel, xmaxcontainer - xmincontainer, 
                         fboundary/(double)(ncircles*NVID), PRINT_LEFT, pressure, gravity);
            if (BOUNDARY_COND == BC_EHRENFEST) print_ehrenfest_parameters(particle, pleft, pright);
            else if (PRINT_PARTICLE_NUMBER) print_particle_number(ncircles);
            if (PRINT_OMEGA) print_omega(angular_speed); 
            else if (PRINT_PARTICLE_SPEEDS) print_particles_speeds(particle);
            glutSwapBuffers();
        }
        for (i=0; i<MID_FRAMES; i++) save_frame_lj();
        if (DOUBLE_MOVIE) 
        {
            if (TRACER_PARTICLE) draw_trajectory(trajectory, traj_position, traj_length);
            draw_particles(particle, PLOT_B);
            draw_container(xmincontainer, xmaxcontainer, obstacle, segment, wall);
            print_parameters(beta, mean_energy, krepel, xmaxcontainer - xmincontainer, 
                         fboundary/(double)(ncircles*NVID), PRINT_LEFT, pressure, gravity);
            if (BOUNDARY_COND == BC_EHRENFEST) print_ehrenfest_parameters(particle, pleft, pright);
            else if (PRINT_PARTICLE_NUMBER) print_particle_number(ncircles);
            if (PRINT_OMEGA) print_omega(angular_speed); 
            else if (PRINT_PARTICLE_SPEEDS) print_particles_speeds(particle);
            glutSwapBuffers();
        }
        for (i=0; i<END_FRAMES; i++) save_frame_lj_counter(NSTEPS + MID_FRAMES + 1 + counter + i);
        if ((TIME_LAPSE)&&(!DOUBLE_MOVIE)) 
            for (i=0; i<END_FRAMES; i++) save_frame_lj_counter(NSTEPS + END_FRAMES + NSTEPS/TIME_LAPSE_FACTOR + i);

        s = system("mv lj*.tif tif_ljones/");
    }
    
    printf("%i active particles\n", nactive);  
    
    free(particle);
    if (ADD_FIXED_OBSTACLES) free(obstacle);
    if (ADD_FIXED_SEGMENTS) free(segment);
    if (TRACER_PARTICLE) free(trajectory);
    free(hashgrid);
    free(qx);
    free(qy);
    free(px);
    free(py);
    free(qangle);
    free(pangle);
    free(pressure);
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
    glutCreateWindow("Particles with Lennard-Jones interaction in a planar domain");

    init();

    glutDisplayFunc(display);

    glutMainLoop();

    return 0;
}

