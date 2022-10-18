/* global variables and definitions used by sub_wave_3d.c */

/* plot types used by wave_3d */

#define P_3D_AMPLITUDE  101     /* height/color depends on amplitude - DEPRECATED, instead use set SHADE_3D to 0 */
#define P_3D_ANGLE 102          /* height/color depends on angle with fixed direction - TO DO */
#define P_3D_AMP_ANGLE 103      /* height/color depends on amplitude, luminosity depends on angle */
#define P_3D_ENERGY 104         /* height/color depends on energy, luminosity depends on angle */
#define P_3D_LOG_ENERGY 105     /* height/color depends on logarithm of energy, luminosity depends on angle */
#define P_3D_TOTAL_ENERGY 106   /* height/color depends on total energy over time, luminosity depends on angle */
#define P_3D_LOG_TOTAL_ENERGY 107 /* height/color depends on log on total energy over time, luminosity depends on angle */
#define P_3D_MEAN_ENERGY 108      /* height/color depends on energy averaged over time, luminosity depends on angle */
#define P_3D_LOG_MEAN_ENERGY 109  /* height/color depends on log on energy averaged over time, luminosity depends on angle */

#define P_3D_PHASE 111          /* phase of wave */

/* Choice of simulated reaction-diffusion equation in rde.c */

#define E_HEAT 0            /* heat equation */
#define E_ALLEN_CAHN 1      /* Allen-Cahn equation */
#define E_CAHN_HILLIARD 2   /* Cahn-Hilliard equation */
#define E_FHN 3             /* FitzHugh-Nagumo equation */
#define E_RPS 4             /* rock-paper-scissors equation */
#define E_RPSLZ 41          /* rock-paper-scissors-lizard-Spock equation */
#define E_SCHRODINGER 5     /* Schrodinger equation */
#define E_EULER_INCOMP 6    /* incompressigle Euler equation */

/* Choice of potential */

#define POT_HARMONIC 1      /* harmonic oscillator */
#define POT_COULOMB 2       /* Coulomb (1/r) potential */
#define POT_PERIODIC 3      /* periodic potential */
#define POT_DOUBLE_COULOMB 4      /* sum of Coulomb potentials located at focal points of ellipse */
#define POT_FERMIONS 5      /* two interacting 1D fermions */
#define POT_FERMIONS_PERIODIC 6      /* two interacting 1D fermions on the circle */
#define POT_MAZE 7          /* higher potential on walls of a maze */

/* Choice of vector potential */

#define VPOT_CONSTANT_FIELD 100  /* constant magnetic field */
#define VPOT_AHARONOV_BOHM 101   /* single flux line for Aharonov-Bohm effect */

/* plot types used by rde */

#define Z_AMPLITUDE 0   /* amplitude of first field */
#define Z_RGB 20        /* RGB plot */
#define Z_POLAR 21      /* polar angle associated with RBG plot */
#define Z_NORM_GRADIENT 22      /* gradient of polar angle */
#define Z_ANGLE_GRADIENT 221    /* direction of polar angle */
#define Z_NORM_GRADIENTX 23     /* norm of gradient of u */
#define Z_ANGLE_GRADIENTX 231   /* direction of gradient of u */
#define Z_NORM_GRADIENT_INTENSITY 24  /* gradient and intensity of polar angle */
#define Z_VORTICITY 25  /* curl of polar angle */
#define Z_VORTICITY_ABS 251  /* absolute value of curl of polar angle */

/* for Schrodinger equation */
#define Z_MODULE 30       /* module squared of first two fields */
#define Z_ARGUMENT 31     /* argument of first two fields, with luminosity depending on module */
#define Z_REALPART 32     /* first field, with luminosity depending on module */

/* for RPSLZ equation */
#define Z_MAXTYPE_RPSLZ 40      /* color of type with maximal density */
#define Z_THETA_RPSLZ 41        /* polar angle */
#define Z_NORM_GRADIENT_RPSLZ 42    /* gradient of polar angle */

/* for Euler equation */
#define Z_EULER_VORTICITY 50    /* vorticity of velocity */
#define Z_EULER_LOG_VORTICITY 51    /* log of vorticity of velocity */
#define Z_EULER_VORTICITY_ASYM 52   /* vorticity of velocity */

/* macros to avoid unnecessary computations in 3D plots */

#define COMPUTE_THETA ((cplot == Z_POLAR)||(cplot == Z_NORM_GRADIENT)||(cplot == Z_ANGLE_GRADIENT)||(cplot == Z_NORM_GRADIENT_INTENSITY)||(cplot == Z_VORTICITY)||(cplot == Z_VORTICITY_ABS))
#define COMPUTE_THETAZ ((zplot == Z_POLAR)||(zplot == Z_NORM_GRADIENT)||(zplot == Z_ANGLE_GRADIENT)||(zplot == Z_NORM_GRADIENT_INTENSITY)||(zplot == Z_VORTICITY)||(zplot == Z_VORTICITY_ABS))

#define COMPUTE_ENERGY ((zplot == P_3D_ENERGY)||(cplot == P_3D_ENERGY)||(zplot == P_3D_LOG_ENERGY)||(cplot == P_3D_LOG_ENERGY)||(zplot == P_3D_TOTAL_ENERGY)||(cplot == P_3D_TOTAL_ENERGY)||(zplot == P_3D_LOG_TOTAL_ENERGY)||(cplot == P_3D_LOG_TOTAL_ENERGY)||(zplot == P_3D_MEAN_ENERGY)||(cplot == P_3D_MEAN_ENERGY)||(zplot == P_3D_LOG_MEAN_ENERGY)||(cplot == P_3D_LOG_MEAN_ENERGY))

#define COMPUTE_LOG_TOTAL_ENERGY ((ZPLOT == P_3D_LOG_TOTAL_ENERGY)||(CPLOT == P_3D_LOG_TOTAL_ENERGY)||(ZPLOT_B == P_3D_LOG_TOTAL_ENERGY)||(CPLOT_B == P_3D_LOG_TOTAL_ENERGY))

#define COMPUTE_LOG_MEAN_ENERGY ((ZPLOT == P_3D_LOG_MEAN_ENERGY)||(CPLOT == P_3D_LOG_MEAN_ENERGY)||(ZPLOT_B == P_3D_LOG_MEAN_ENERGY)||(CPLOT_B == P_3D_LOG_MEAN_ENERGY))

#define COMPUTE_LOG_ENERGY ((ZPLOT == P_3D_LOG_TOTAL_ENERGY)||(CPLOT == P_3D_LOG_TOTAL_ENERGY)||(ZPLOT_B == P_3D_LOG_TOTAL_ENERGY)||(CPLOT_B == P_3D_LOG_TOTAL_ENERGY)||(ZPLOT == P_3D_LOG_MEAN_ENERGY)||(CPLOT == P_3D_LOG_MEAN_ENERGY)||(ZPLOT_B == P_3D_LOG_MEAN_ENERGY)||(CPLOT_B == P_3D_LOG_MEAN_ENERGY))

#define COMPUTE_MEAN_ENERGY ((ZPLOT == P_3D_MEAN_ENERGY)||(CPLOT == P_3D_MEAN_ENERGY)||(ZPLOT_B == P_3D_MEAN_ENERGY)||(CPLOT_B == P_3D_MEAN_ENERGY)||(ZPLOT == P_3D_LOG_MEAN_ENERGY)||(CPLOT == P_3D_LOG_MEAN_ENERGY)||(ZPLOT_B == P_3D_LOG_MEAN_ENERGY)||(CPLOT_B == P_3D_LOG_MEAN_ENERGY))

#define COMPUTE_TOTAL_ENERGY ((ZPLOT == P_3D_TOTAL_ENERGY)||(CPLOT == P_3D_TOTAL_ENERGY)||(ZPLOT == P_3D_LOG_TOTAL_ENERGY)||(CPLOT == P_3D_LOG_TOTAL_ENERGY)||(ZPLOT == P_3D_MEAN_ENERGY)||(CPLOT == P_3D_MEAN_ENERGY)||(ZPLOT == P_3D_LOG_MEAN_ENERGY)||(CPLOT == P_3D_LOG_MEAN_ENERGY)||(ZPLOT_B == P_3D_TOTAL_ENERGY)||(CPLOT_B == P_3D_TOTAL_ENERGY)||(ZPLOT_B == P_3D_LOG_TOTAL_ENERGY)||(CPLOT_B == P_3D_LOG_TOTAL_ENERGY)||(ZPLOT_B == P_3D_MEAN_ENERGY)||(CPLOT_B == P_3D_MEAN_ENERGY)||(ZPLOT_B == P_3D_LOG_MEAN_ENERGY)||(CPLOT_B == P_3D_LOG_MEAN_ENERGY))


int global_time = 0;

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
    double rgb[3];              /* RGB color code */
    double *p_zfield[2];        /* pointers to z field (second pointer for option DOUBLE_MOVIE) */
    double *p_cfield[2];        /* pointers to color field (second pointer for option DOUBLE_MOVIE) */
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
    double rgb[3];              /* RGB color code */
    double *p_zfield[2];        /* pointers to z field (second pointer for option DOUBLE_MOVIE) */
    double *p_cfield[2];        /* pointers to color field (second pointer for option DOUBLE_MOVIE) */
} t_rde;



