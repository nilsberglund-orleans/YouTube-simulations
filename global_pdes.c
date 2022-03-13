/* Global variables and parameters for wave_billiard, heat and schrodinger */

// #include "hsluv.c"


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
#define D_CIRCLES_IN_RECT 201   /* several circles in a rectangle */

#define D_FOUR_PARABOLAS 31     /* four parabolas with axes in NSEW directions */
#define D_POLY_PARABOLAS 32     /* polygon with parabolic sides */
#define D_PENROSE 33            /* Penrose unilluminable room */
#define D_HYPERBOLA 34          /* one branch of hyperbola */
#define D_TOKARSKY 35           /* Tokarsky unilluminable room */
#define D_TOKA_PRIME 36         /* Tokarsky room made of 86 triangles */
#define D_ISOSPECTRAL 37        /* isospectral billiards */
#define D_HOMOPHONIC 38         /* homophonic billiards */

#define D_POLYGONS 40           /* several polygons */
#define D_VONKOCH 41            /* von Koch snowflake fractal */
#define D_STAR 42               /* star shape */
#define D_FRESNEL 43            /* Fresnel lens */
#define D_NOISEPANEL 44         /* zigzag noise insulating panel */
#define D_DOUBLE_FRESNEL 45     /* two facing Fresnel lenses */
#define D_QRD 46                /* quadratic resonance diffuser */
#define D_CIRCLE_SEGMENT 47     /* lens-shaped circular segment */

#define NMAXCIRCLES 10000       /* total number of circles/polygons (must be at least NCX*NCY for square grid) */
#define NMAXPOLY 50000          /* maximal number of vertices of polygonal lines (for von Koch et al) */
// #define NMAXCIRCLES 10000        /* total number of circles/polygons (must be at least NCX*NCY for square grid) */

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
#define C_HEX_NONUNIF 13    /* triangular grid with non-constant column distance */

#define C_RINGS 20          /* obstacles arranged in concentric rings */
#define C_RINGS_T 201       /* obstacles arranged in concentric rings, triangular lattice */
#define C_RINGS_SPIRAL 202  /* obstacles arranged on a "subflower" spiral, similar to C_GOLDEN_SPIRAL */

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
#define D_VONKOCH_HEATED 27 /* von Koch snowflake in larger circle */

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
#define P_LOG_ENERGY 4  /* log of energy averaged over time */
#define P_LOG_MEAN_ENERGY 5  /* log of energy averaged over time */

/* For Schrodinger equation */
#define P_MODULE 10        /* plot module of wave function squared */
#define P_PHASE 11         /* plot phase of wave function */
#define P_REAL 12          /* plot real part */
#define P_IMAGINARY 13     /* plot imaginary part */


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
    short int active, top;      /* circle is active, circle is in top half */
} t_circle;

typedef struct
{
    double xc, yc, radius, angle;      /* center, radius and angle of polygon */
    int nsides;                         /* number of sides of polygon */
    short int active, top;              /* polygon is active, polygon is in top half */
} t_polygon;

typedef struct
{
    double x, y;           /* (x,y) coordinates of vertex */
    double posi, posj;     /* (i,j) coordinates of vertex */
} t_vertex;


// double circlex[NMAXCIRCLES], circley[NMAXCIRCLES], circlerad[NMAXCIRCLES];      /* position and radius of circular scatterers */
// short int circleactive[NMAXCIRCLES];                                      /* tells which circular scatters are active */
int ncircles = NMAXCIRCLES;         /* actual number of circles, can be decreased e.g. for random patterns */
int npolyline = NMAXPOLY;           /* actual length of polyline */

t_circle circles[NMAXCIRCLES];      /* circular scatterers */
t_polygon polygons[NMAXCIRCLES];    /* polygonal scatterers */
t_vertex polyline[NMAXPOLY];        /* vertices of polygonal line */

double julia_x = -0.5, julia_y = 0.5;    /* parameters for Julia sets */
// double julia_x = 0.33267, julia_y = 0.06395;    /* parameters for Julia sets */
// double julia_x = 0.37468, julia_y = 0.21115;    /* parameters for Julia sets */
