/* global variables and definitions used by sub_wave_3d.c */

/* Choice of simulated reaction-diffusion equation in rde.c */

#define E_HEAT 0            /* heat equation */
#define E_ALLEN_CAHN 1      /* Allen-Cahn equation */
#define E_CAHN_HILLIARD 2   /* Cahn-Hilliard equation */
#define E_FHN 3             /* FitzHugh-Nagumo equation */
#define E_RPS 4             /* rock-paper-scissors equation */
#define E_RPSLZ 41          /* rock-paper-scissors-lizard-Spock equation */
#define E_SCHRODINGER 5     /* Schrodinger equation */
#define E_EULER_INCOMP 6    /* incompressible Euler equation */
#define E_EULER_COMP 7      /* compressible Euler equation */
#define E_SHALLOW_WATER 8   /* shallow water equation */

/* Choice of potential */

#define POT_HARMONIC 1      /* harmonic oscillator */
#define POT_COULOMB 2       /* Coulomb (1/r) potential */
#define POT_PERIODIC 3      /* periodic potential */
#define POT_DOUBLE_COULOMB 4      /* sum of Coulomb potentials located at focal points of ellipse */
#define POT_FERMIONS 5      /* two interacting 1D fermions */
#define POT_FERMIONS_PERIODIC 6      /* two interacting 1D fermions on the circle */
#define POT_MAZE 7          /* higher potential on walls of a maze */
#define POT_IOR 10          /* index of refraction, for z coordinate of wave equation */

/* Choice of vector potential */

#define VPOT_CONSTANT_FIELD 100  /* constant magnetic field */
#define VPOT_AHARONOV_BOHM 101   /* single flux line for Aharonov-Bohm effect */

/* Choice of force field in compressible Euler equation */

#define GF_VERTICAL 0       /* gravity acting vertically */
#define GF_CIRCLE 1         /* repelling circle */
#define GF_ELLIPSE 2        /* repelling ellipse */
#define GF_AIRFOIL 3        /* curved repelling ellipse */
#define GF_WING 4           /* wing shape */
#define GF_COMPUTE_FROM_BC 5    /* compute force field as gradient of bc_field2 */
#define GF_EARTH 6          /* field depends on altitude on continents */
#define GF_MARS 7           /* field depends on altitude on Mars */
#define GF_VENUS 8          /* field depends on altitude on Venus */

/* Choice of water depth for shallow water equation */

#define SH_CIRCLE 1     /* circular shallower obstacle */
#define SH_CIRCLES 2    /* shallow obstacle specified by CIRCLE_PATTERN */
#define SH_COAST 3      /* depth varying with x-coordinate */
#define SH_COAST_MONOTONE 4      /* depth decreasing with x-coordinate */
#define SH_SPHERE_CUBE 5    /* cube embedded in sphere */
#define SH_SPHERE_OCTAHEDRON 6   /* octahedron embedded in sphere */
#define SH_SPHERE_DODECAHEDRON 7 /* dodecahedron embedded in sphere */
#define SH_SPHERE_ICOSAHEDRON 8  /* icosahedron embedded in sphere */
#define SH_EARTH 10         /* depth of Earth oceans */

/* Type of rotating viewpoint */

#define VP_HORIZONTAL 0     /* rotate in a horizontal plane (constant latitude) */
#define VP_ORBIT 1          /* rotate in a plane containing the origin */
#define VP_ORBIT2 11        /* rotate in a plane specified by max latitude */
#define VP_POLAR 2          /* polar orbit */

/* Type of digital elevation model */

#define DEM_EARTH 0         /* DEM of Earth */
#define DEM_MARS 1          /* DEM of Mars */
#define DEM_MOON 2          /* DEM of the Moon */
#define DEM_VENUS 3         /* DEM of Venus */
#define DEM_MERCURY 4       /* DEM of Mercury */

/* macros to avoid unnecessary computations in 3D plots */

#define COMPUTE_THETA ((cplot == Z_POLAR)||(cplot == Z_NORM_GRADIENT)||(cplot == Z_ANGLE_GRADIENT)||(cplot == Z_NORM_GRADIENT_INTENSITY)||(cplot == Z_VORTICITY)||(cplot == Z_VORTICITY_ABS))
#define COMPUTE_THETAZ ((zplot == Z_POLAR)||(zplot == Z_NORM_GRADIENT)||(zplot == Z_ANGLE_GRADIENT)||(zplot == Z_NORM_GRADIENT_INTENSITY)||(zplot == Z_VORTICITY)||(zplot == Z_VORTICITY_ABS))

#define COMPUTE_ENERGY ((zplot == P_3D_ENERGY)||(cplot == P_3D_ENERGY)||(zplot == P_3D_LOG_ENERGY)||(cplot == P_3D_LOG_ENERGY)||(zplot == P_3D_TOTAL_ENERGY)||(cplot == P_3D_TOTAL_ENERGY)||(zplot == P_3D_LOG_TOTAL_ENERGY)||(cplot == P_3D_LOG_TOTAL_ENERGY)||(zplot == P_3D_MEAN_ENERGY)||(cplot == P_3D_MEAN_ENERGY)||(zplot == P_3D_LOG_MEAN_ENERGY)||(cplot == P_3D_LOG_MEAN_ENERGY)||(ZPLOT == P_3D_FLUX_INTENSITY)||(CPLOT == P_3D_FLUX_INTENSITY)||(ZPLOT_B == P_3D_FLUX_INTENSITY)||(CPLOT_B == P_3D_FLUX_INTENSITY)||(ZPLOT == P_3D_FLUX_DIRECTION)||(CPLOT == P_3D_FLUX_DIRECTION)||(ZPLOT_B == P_3D_FLUX_DIRECTION)||(CPLOT_B == P_3D_FLUX_DIRECTION))

#define COMPUTE_LOG_TOTAL_ENERGY ((ZPLOT == P_3D_LOG_TOTAL_ENERGY)||(CPLOT == P_3D_LOG_TOTAL_ENERGY)||(ZPLOT_B == P_3D_LOG_TOTAL_ENERGY)||(CPLOT_B == P_3D_LOG_TOTAL_ENERGY))

#define COMPUTE_LOG_MEAN_ENERGY ((ZPLOT == P_3D_LOG_MEAN_ENERGY)||(CPLOT == P_3D_LOG_MEAN_ENERGY)||(ZPLOT_B == P_3D_LOG_MEAN_ENERGY)||(CPLOT_B == P_3D_LOG_MEAN_ENERGY))

#define COMPUTE_LOG_ENERGY ((ZPLOT == P_3D_LOG_TOTAL_ENERGY)||(CPLOT == P_3D_LOG_TOTAL_ENERGY)||(ZPLOT_B == P_3D_LOG_TOTAL_ENERGY)||(CPLOT_B == P_3D_LOG_TOTAL_ENERGY)||(ZPLOT == P_3D_LOG_MEAN_ENERGY)||(CPLOT == P_3D_LOG_MEAN_ENERGY)||(ZPLOT_B == P_3D_LOG_MEAN_ENERGY)||(CPLOT_B == P_3D_LOG_MEAN_ENERGY))

#define COMPUTE_MEAN_ENERGY ((ZPLOT == P_3D_MEAN_ENERGY)||(CPLOT == P_3D_MEAN_ENERGY)||(ZPLOT_B == P_3D_MEAN_ENERGY)||(CPLOT_B == P_3D_MEAN_ENERGY)||(ZPLOT == P_3D_LOG_MEAN_ENERGY)||(CPLOT == P_3D_LOG_MEAN_ENERGY)||(ZPLOT_B == P_3D_LOG_MEAN_ENERGY)||(CPLOT_B == P_3D_LOG_MEAN_ENERGY))

#define COMPUTE_ENERGY_FLUX ((ZPLOT == P_3D_FLUX_INTENSITY)||(CPLOT == P_3D_FLUX_INTENSITY)||(ZPLOT_B == P_3D_FLUX_INTENSITY)||(CPLOT_B == P_3D_FLUX_INTENSITY)||(ZPLOT == P_3D_FLUX_DIRECTION)||(CPLOT == P_3D_FLUX_DIRECTION)||(ZPLOT_B == P_3D_FLUX_DIRECTION)||(CPLOT_B == P_3D_FLUX_DIRECTION))

#define COMPUTE_TOTAL_ENERGY ((ZPLOT == P_3D_TOTAL_ENERGY)||(CPLOT == P_3D_TOTAL_ENERGY)||(ZPLOT == P_3D_LOG_TOTAL_ENERGY)||(CPLOT == P_3D_LOG_TOTAL_ENERGY)||(ZPLOT == P_3D_MEAN_ENERGY)||(CPLOT == P_3D_MEAN_ENERGY)||(ZPLOT == P_3D_LOG_MEAN_ENERGY)||(CPLOT == P_3D_LOG_MEAN_ENERGY)||(ZPLOT_B == P_3D_TOTAL_ENERGY)||(CPLOT_B == P_3D_TOTAL_ENERGY)||(ZPLOT_B == P_3D_LOG_TOTAL_ENERGY)||(CPLOT_B == P_3D_LOG_TOTAL_ENERGY)||(ZPLOT_B == P_3D_MEAN_ENERGY)||(CPLOT_B == P_3D_MEAN_ENERGY)||(ZPLOT_B == P_3D_LOG_MEAN_ENERGY)||(CPLOT_B == P_3D_LOG_MEAN_ENERGY))

#define PLANET ((B_DOMAIN == D_SPHERE_EARTH)||(B_DOMAIN == D_SPHERE_MARS)||(B_DOMAIN == D_SPHERE_MOON)||(B_DOMAIN == D_SPHERE_VENUS)||(B_DOMAIN == D_SPHERE_MERCURY))
#define OTHER_PLANET ((B_DOMAIN == D_SPHERE_MARS)||(B_DOMAIN == D_SPHERE_MOON)||(B_DOMAIN == D_SPHERE_VENUS)||(B_DOMAIN == D_SPHERE_MERCURY))

#define RDE_PLANET (((ADAPT_STATE_TO_BC)&&((OBSTACLE_GEOMETRY == D_SPHERE_EARTH)||(OBSTACLE_GEOMETRY == D_SPHERE_MARS)||(OBSTACLE_GEOMETRY == D_SPHERE_VENUS))))
// #define RDE_PLANET (((ADAPT_STATE_TO_BC)&&((OBSTACLE_GEOMETRY == D_SPHERE_EARTH)||(OBSTACLE_GEOMETRY == D_SPHERE_MARS)||(OBSTACLE_GEOMETRY == D_SPHERE_VENUS)))||((RDE_EQUATION == E_SHALLOW_WATER)&&(SWATER_DEPTH == SH_EARTH)))

#define NMAXCIRC_SPHERE 100     /* max number of circles on sphere */
#define NMAX_TRACER_PTS 20       /* max number of tracer points recorded per cell */

int global_time = 0;
double max_depth = 1.0;
int moon_position;

/* structure used for color and height representations */
/* possible extra fields: zfield, cfield, interpolated coordinates */

typedef struct
{
    double energy;              /* wave energy */
    double phase;               /* wave phase */
    double log_energy;          /* log of wave energy */
    double total_energy;        /* total energy since beginning of simulation */
    double log_total_energy;    /* log of total energy since beginning of simulation */
    double mean_energy;         /* energy averaged since beginning of simulation */
    double log_mean_energy;     /* log of energy averaged since beginning of simulation */
    double cos_angle;           /* cos of angle between normal vector and direction of light */
    double flux_intensity;      /* intensity of energy flux */
    double flux_direction;      /* direction of energy flux */
    double flux_int_table[FLUX_WINDOW];   /* table of energy flux intensities (for averaging) */
    short int flux_counter;     /* counter for averaging of energy flux */
    double rgb[3];              /* RGB color code */
    double *potential;          /* pointer to "potential" to add to z-coordinate */
    double *p_zfield[2];        /* pointers to z field (second pointer for option DOUBLE_MOVIE) */
    double *p_cfield[4];        /* pointers to color field (second pointer for option DOUBLE_MOVIE) */
                                /* third and fourth pointer for color luminosity (for energy flux) */
} t_wave;


typedef struct
{
    double theta;               /* angle for Rock-Paper-Scissors equation */
    double nablax;              /* gradient of first field */
    double nablay;              /* gradient of first field */
    double field_norm;          /* norm of field or gradient */
    double field_arg;           /* argument of field or gradient */
    double curl;                /* curl of field */
    double cos_angle;           /* cos of angle between normal vector and direction of light */
    double log_vorticity;       /* logarithm of vorticity (for Euler equation) */
    double Lpressure;           /* Laplacian of pressure (for Euler equation) */
    double height;              /* wave height (for shallow wave equation) */
    double dxu, dyu, dxv, dyv;  /* gradient of velocity field (for compressible Euler equation) */
    double rgb[3];              /* RGB color code */
    double *p_zfield[2];        /* pointers to z field (second pointer for option DOUBLE_MOVIE) */
    double *p_cfield[2];        /* pointers to color field (second pointer for option DOUBLE_MOVIE) */
    double depth;               /* water depth */
    double cos_depth_angle;     /* cos of angle of depth profile */
    double gradx, grady;        /* gradient of water depth */
    short int tracer;           /* has value 1 if cell contains a tracer */
    short int n_tracer_pts;     /* number of recorded tracer points per cell */
    double tracerx[NMAX_TRACER_PTS], tracery[NMAX_TRACER_PTS];    /* coordinates of tracer */
    int prev_cell;              /* cell where tracer was previously */
} t_rde;


typedef struct
{
    double depth;               /* water depth */
    double gradx, grady;        /* gradient of water depth */
} t_swater_depth;


typedef struct
{
    double phi, theta;          /* phi, theta angles */
    double cphi, sphi;          /* cos and sin of phi */
    double ctheta, stheta, cottheta;   /* cos, sin and cotangent of theta */
    double reg_cottheta;        /* regularized cotangent of theta */
    double x, y, z;             /* x, y, z coordinates of point on sphere */
    double radius;              /* radius with wave height */
    double radius_dem;          /* radius with digital elevation model */
    double r, g, b;             /* RGB values for image */
    short int indomain;         /* has value 1 if lattice point is in domain */
    short int draw_wave;        /* has value 1 if wave instead of DEM is drawn */
    short int evolve_wave;      /* has value 1 where there is wave evolution */
    double x2d, y2d;            /* x and y coordinates for 2D representation */
    double altitude;            /* altitude in case of Earth with digital elevation model */
    double cos_angle;           /* cosine of light angle */
    double cos_angle_sphere;    /* cosing of light angle for perfect sphere */
    double force;             /* external forcing */
} t_wave_sphere;


typedef struct 
{
    double phi, theta;          /* longitude, latitude */
    double radius;              /* radius */
    double x, y, z;             /* x, y, z coordinates of point on sphere */
} t_circles_sphere;


t_circles_sphere circ_sphere[NMAXCIRC_SPHERE];      /* circular scatterers on sphere */
