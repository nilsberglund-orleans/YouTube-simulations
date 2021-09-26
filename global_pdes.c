/* Global variables and parameters for wave_billiard, heat and schrodinger */

/* Basic math */

#define PI 	3.141592654
#define DPI 	6.283185307
#define PID 	1.570796327

/* shape of domain */

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

/* The following 3 types are superseded by D_CIRCLES and can be removed later */
#define D_DISK_GRID 12  /* grid of disks */
#define D_DISK_HEX 13   /* haxagonl grid of disks */
#define D_DISK_PERCOL 14    /* random grid of percolation type */

#define D_MENGER 15     /* Menger-Sierpinski carpet */ 
#define D_JULIA_INT 16  /* interior of Julia set */ 
#define D_MENGER_ROTATED 17  /* rotated Menger-Sierpinski carpet */
#define D_PARABOLA 18   /* parabolic domain */
#define D_TWO_PARABOLAS 19   /* two facing parabolic antennas */

#define D_CIRCLES 20    /* several circles */

#define D_FOUR_PARABOLAS 31     /* four parabolas with axes in NSEW directions */
#define D_POLY_PARABOLAS 32     /* polygon with parabolic sides */
#define D_PENROSE 33            /* Penrose illumination problem */
#define D_HYPERBOLA 34          /* one branch of hyperbola */

#define NMAXCIRCLES 1000        /* total number of circles (must be at least NCX*NCY for square grid) */

#define C_SQUARE 0          /* square grid of circles */
#define C_HEX 1             /* hexagonal/triangular grid of circles */
#define C_RAND_DISPLACED 2  /* randomly displaced square grid */
#define C_RAND_PERCOL 3     /* random percolation arrangement */
#define C_RAND_POISSON 4    /* random Poisson point process */
#define C_CLOAK 5           /* invisibility cloak */
#define C_CLOAK_A 6         /* first optimized invisibility cloak */

#define C_POISSON_DISC 8    /* Poisson disc sampling */

#define C_GOLDEN_MEAN 10    /* pattern based on vertical shifts by golden mean */
#define C_GOLDEN_SPIRAL 11  /* spiral pattern based on golden mean */
#define C_SQUARE_HEX 12     /* alternating between square and hexagonal/triangular */

#define C_ONE 97            /* one single circle, as for Sinai */
#define C_TWO 98            /* two concentric circles of different type */
#define C_NOTHING 99        /* no circle at all, for comparisons */


/* Billiard tables for heat equation */

#define D_ANNULUS_HEATED 21 /* annulus with different temperatures */
#define D_MENGER_HEATED 22  /* Menger gasket with different temperatures */
#define D_MENGER_H_OPEN 23  /* Menger gasket with different temperatures and larger domain */
#define D_MANDELBROT 24     /* Mandelbrot set */
#define D_JULIA 25          /* Julia set */
#define D_MANDELBROT_CIRCLE 26     /* Mandelbrot set with circular conductor */

/* Boundary conditions */

#define BC_DIRICHLET 0   /* Dirichlet boundary conditions */
#define BC_PERIODIC 1    /* periodic boundary conditions */
#define BC_ABSORBING 2   /* absorbing boundary conditions (beta version) */
#define BC_VPER_HABS 3   /* vertically periodic and horizontally absorbing boundary conditions */
#define BC_ABS_REFLECT 4   /* absorbing boundary conditions, except reflecting at y=0, for comparisons */
// #define BC_OSCILL_ABSORB 5  /* oscillating boundary condition on the left, absorbing on other walls */ 

/* For debugging purposes only */
// #define FLOOR 0         /* set to 1 to limit wave amplitude to VMAX */
// #define VMAX 10.0       /* max value of wave amplitude */

/* Plot types */

/* For wave equation */
#define P_AMPLITUDE 0    /* plot amplitude of wave */
#define P_ENERGY 1       /* plot energy of wave */
#define P_MIXED 2        /* plot amplitude in upper half, energy in lower half */
#define P_MEAN_ENERGY 3  /* energy averaged over time */

/* For Schrodinger equation */
#define P_MODULE 10        /* plot module of wave function squared */
#define P_PHASE 11         /* plot phase of wave function */
#define P_REAL 12          /* plot real part */
#define P_IMAGINARY 13     /* plot imaginary part */


/* Color schemes */

#define C_LUM 0          /* color scheme modifies luminosity (with slow drift of hue) */
#define C_HUE 1          /* color scheme modifies hue */
#define C_PHASE 2        /* color scheme shows phase */


double circlex[NMAXCIRCLES], circley[NMAXCIRCLES], circlerad[NMAXCIRCLES];      /* position and radius of circular scatterers */
short int circleactive[NMAXCIRCLES];                                      /* tells which circular scatters are active */
int ncircles = NMAXCIRCLES;            /* actual number of circles, can be decreased e.g. for random patterns */

double julia_x = -0.5, julia_y = 0.5;    /* parameters for Julia sets */
// double julia_x = 0.33267, julia_y = 0.06395;    /* parameters for Julia sets */
// double julia_x = 0.37468, julia_y = 0.21115;    /* parameters for Julia sets */
