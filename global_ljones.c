/* Global variables and parameters for lennardjones */

/* Basic math */

#define PI 	3.141592654
#define DPI 	6.283185307
#define PID 	1.570796327

/* shape of domain */

#define D_CIRCLES 20    /* several circles */
#define D_CIRCLES_IN_RECT 201   /* several circles in a rectangle */

#define NMAXCIRCLES 20000       /* total number of circles/polygons (must be at least NCX*NCY for square grid) */
#define MAXNEIGH 20         /* max number of neighbours kept in memory */
#define NMAXOBSTACLES 100   /* max number of obstacles */

#define C_SQUARE 0          /* square grid of circles */
#define C_HEX 1             /* hexagonal/triangular grid of circles */
#define C_RAND_DISPLACED 2  /* randomly displaced square grid */
#define C_RAND_PERCOL 3     /* random percolation arrangement */
#define C_RAND_POISSON 4    /* random Poisson point process */
#define C_CLOAK 5           /* invisibility cloak */
#define C_CLOAK_A 6         /* first optimized invisibility cloak */
#define C_LASER 7           /* laser fight in a room of mirrors */
#define C_POISSON_DISC 8    /* Poisson disc sampling */

#define C_GOLDEN_MEAN 10    /* pattern based on vertical shifts by golden mean */
#define C_GOLDEN_SPIRAL 11  /* spiral pattern based on golden mean */
#define C_SQUARE_HEX 12     /* alternating between square and hexagonal/triangular */

#define C_ONE 97            /* one single circle, as for Sinai */
#define C_TWO 98            /* two concentric circles of different type */
#define C_NOTHING 99        /* no circle at all, for comparisons */

/* pattern of additional obstacles */
#define O_CORNERS 0         /* obstacles in the corners (for Boy b.c.) */
#define O_GALTON_BOARD 1    /* Galton board pattern */

/* particle interaction */

#define I_COULOMB 0         /* Coulomb force */
#define I_LENNARD_JONES 1   /* Lennard-Jones force */
#define I_LJ_DIRECTIONAL 2  /* Lennard-Jones with direction dependence of square symmetry */
#define I_LJ_PENTA 3        /* Lennard-Jones with pentagonal symmetry */
#define I_GOLDENRATIO 4     /* Lennard-Jones type with equilibria related by golden ratio */
#define I_LJ_DIPOLE 5       /* Lennard-Jones with a dipolar angle dependence */
#define I_LJ_QUADRUPOLE 6   /* Lennard-Jones with a quadropolar angle dependence */
#define I_LJ_WATER 7        /* model for water molecule */

/* Boundary conditions */

#define BC_SCREEN 0         /* harmonic boundary conditions outside screen area */
#define BC_RECTANGLE 1      /* harmonic boundary conditions on a resizeable rectangle */
#define BC_CIRCLE 2         /* harmonic boundary conditions outside a moving circle */
#define BC_PERIODIC 3       /* periodic boundary conditions */
#define BC_PERIODIC_CIRCLE 4  /* periodic boundary conditions and harmonic b.c. outside moving circle */
#define BC_EHRENFEST 5      /* Ehrenfest urn-type configuration */
#define BC_PERIODIC_FUNNEL 6    /* funnel with periodic boundary conditions */
#define BC_RECTANGLE_LID 7  /* rectangular container with moving lid */
#define BC_PERIODIC_TRIANGLE 8  /* periodic boundary conditions and harmonic b.c. outside moving triangle */
#define BC_KLEIN 11         /* Klein bottle (periodic with twisted vertical parts) */
#define BC_SCREEN_BINS 12   /* harmonic boundary conditions outside screen area plus "bins" (for Galton board) */
#define BC_BOY 13           /* Boy surface/projective plane (periodic with twisted horizontal and vertical parts) */

/* Plot types */

#define P_KINETIC 0       /* colors represent kinetic energy of particles */
#define P_NEIGHBOURS 1    /* colors represent number of neighbours */
#define P_HEALTH 2        /* colors represent health (for SIR model) */
#define P_BONDS 3         /* draw lattice based on neighbours */
#define P_ANGLE 4         /* colors represent angle/spin of particle */
#define P_TYPE 5          /* colors represent type of particle */
#define P_DIRECTION 6     /* colors represent direction of velocity */
#define P_ANGULAR_SPEED 7 /* colors represent angular speed */

/* Color schemes */

#define C_LUM 0          /* color scheme modifies luminosity (with slow drift of hue) */
#define C_HUE 1          /* color scheme modifies hue */
#define C_PHASE 2        /* color scheme shows phase */
#define C_ONEDIM 3       /* use preset 1d color scheme (for Turbo, Viridis, Magma, Inferno, Plasma) */
#define C_ONEDIM_LINEAR 4   /* use preset 1d color scheme with linear scale */

/* Color palettes */

#define COL_JET 0       /* JET color palette */
#define COL_HSLUV 1     /* HSLUV color palette (perceptually uniform) */
#define COL_GRAY 2      /* grayscale */

#define COL_TURBO 10     /* TURBO color palette (by Anton Mikhailov) */
#define COL_VIRIDIS 11   /* Viridis color palette */
#define COL_MAGMA 12     /* Magma color palette */
#define COL_INFERNO 13   /* Inferno color palette */
#define COL_PLASMA 14    /* Plasma color palette */
#define COL_CIVIDIS 15   /* Cividis color palette */
#define COL_PARULA 16    /* Parula color palette */
#define COL_TWILIGHT 17  /* Twilight color palette */
#define COL_TWILIGHT_SHIFTED 18  /* Shifted twilight color palette */

#define COL_TURBO_CYCLIC 101    /* TURBO color palette (by Anton Mikhailov) corrected to be cyclic, beta */


typedef struct
{
    double xc, yc, radius;      /* center and radius of circle */
    double angle;               /* angle of particle's "spin" */
    short int active;           /* circle is active */
    double energy;              /* dissipated energy */
    double vx;                  /* x velocity of particle */
    double vy;                  /* y velocity of particle */
    double omega;               /* angular velocity of particle */
    double mass_inv;            /* inverse of particle mass */
    double inertia_moment_inv;  /* inverse of moment of inertia */
    double fx;                  /* x component of force on particle */
    double fy;                  /* y component of force on particle */
    double torque;              /* torque on particle */
    short int thermostat;       /* whether particle is coupled to thermostat */
    int hashcell;               /* hash cell in which particle is located */
    int neighb;                 /* number of neighbours within given distance */
    int hash_nneighb;           /* number of neighbours in hashgrid */
    int hashneighbour[9*HASHMAX];   /* particle numbers of neighbours in hashgrid */
    double deltax[9*HASHMAX];   /* relative position of neighbours */
    double deltay[9*HASHMAX];   /* relative position of neighbours */
    short int type;             /* type of particle, for mixture simulations */
    short int interaction;      /* type of interaction */
    double eq_dist;             /* equilibrium distance */
    double spin_range;          /* range of spin-spin interaction */
    double spin_freq;           /* angular frequency of spin-spin interaction */
} t_particle;

typedef struct
{
    int number;                 /* total number of particles in cell */
    int particles[HASHMAX];     /* numbers of particles in cell */
    int nneighb;                /* number of neighbouring cells */
    int neighbour[9];           /* numbers of neighbouring cells */
} t_hashgrid;

typedef struct
{
    double xc, yc, radius;      /* center and radius of circle */
    short int active;           /* circle is active */
    double energy;              /* dissipated energy */
    double vx;                  /* x velocity of particle */
    double vy;                  /* y velocity of particle */
    double mass_inv;            /* inverse of particle mass */
    double fx;                  /* x component of force on particle */
    double fy;                  /* y component of force on particle */
    int hashx;                  /* hash grid positions of particles */
    int hashy;                  /* hash grid positions of particles */
    int neighb;                 /* number of neighbours */
    int health;                 /* 0 = healthy, 1 = infected, 2 = recovered */
    double infected_time;       /* time since infected */
    int protected;              /* 0 = not protected, 1 = protected */
} t_person;

typedef struct
{
    double xc, yc, radius;      /* center and radius of circle */
    short int active;           /* circle is active */
} t_obstacle;

typedef struct
{
    double xc, yc;              /* center of circle */
} t_tracer;


int ncircles, nobstacles, counter = 0;

