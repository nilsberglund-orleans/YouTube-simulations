/* global variables and definition used by sub_wave_3d.c */

/* plot types used by wave_3d */

#define P_3D_AMPLITUDE  101     /* color depends on amplitude */
#define P_3D_ANGLE 102          /* color depends on angle with fixed direction */
#define P_3D_AMP_ANGLE 103      /* color depends on amplitude, luminosity depends on angle */
#define P_3D_ENERGY 104         /* color depends on energy, luminosity depends on angle */
#define P_3D_LOG_ENERGY 105     /* color depends on logarithm of energy, luminosity depends on angle */

#define P_3D_PHASE 111          /* phase of wave */

/* plot types used by rde */

#define Z_AMPLITUDE 0   /* amplitude of first field */
#define Z_NORM_GRADIENT 22   /* gradient of polar angle */
#define Z_NORM_GRADIENTX 23  /* direction of gradient of u */
#define Z_NORM_GRADIENT_INTENSITY 24  /* gradient and intensity of polar angle */
#define Z_VORTICITY 25  /* curl of polar angle */


#define COMPUTE_THETA ((plot == P_POLAR)||(plot == P_GRADIENT)||(plot == P_GRADIENTX)||(plot == P_GRADIENT_INTENSITY)||(plot == P_VORTICITY))
#define COMPUTE_THETAZ ((zplot == P_POLAR)||(zplot == P_GRADIENT)||(zplot == P_GRADIENTX)||(zplot == P_GRADIENT_INTENSITY)||(zplot == P_VORTICITY))


/* structure used for color and height representations */
/* possible extra fields: zfield, cfield, interpolated coordinates */

typedef struct
{
    double energy;         /* wave energy */
    double phase;          /* wave phase */
    double log_energy;     /* log of wave energy */
    double total_energy;   /* total energy since beginning of simulation */
    double cos_angle;      /* cos of angle between normal vector and direction of light */
    double rgb[3];         /* RGB color code */
    double *p_zfield;      /* pointer to z field */
    double *p_cfield;      /* pointer to color field */
} t_wave;
