// double circlex[NMAXCIRCLES], circley[NMAXCIRCLES], circlerad[NMAXCIRCLES];      /* position and radius of circular scatters */
// short int circleactive[NMAXCIRCLES];     /* tells which circular scatters are active */
// short int newcircle[NMAXCIRCLES];        /* takes value 1 when circle has just been hit */ 
// int circlecolor[NMAXCIRCLES];          /* color of circular scatterer */
int ncircles = NMAXCIRCLES;            /* actual number of circles, can be decreased e.g. for random patterns */
int nsides = NMAXPOLY;                 /* actual number of sides of polygonal line */
int narcs = NMAXCIRCLES;              /* actual number of arcs */

typedef struct
{
    double xc, yc, radius;                  /* center and radius of circle */
    short int active, new, double_circle;   /* circle is active, has just been hit, has a partner (for torus) */
    int color, partner;                     /* circle color, number of partner */
} t_circle;

typedef struct
{
    double x1, x2, y1, y2, length, angle;
    int color;
} t_segment;

typedef struct
{
    double xc, yc, radius, angle1, dangle;
    int color;
} t_arc;

typedef struct
{
    short int left, right;
} t_exit;

t_circle circles[NMAXCIRCLES];
t_segment polyline[NMAXPOLY];
t_arc arcs[NMAXCIRCLES];

double x_shooter = -0.2, y_shooter = -0.6, x_target = 0.4, y_target = 0.7;    
/* shooter and target positions for "laser in room of mirrors" simulations, with default values for square domain */

/* some basic math */

#define PI 	3.141592654
#define DPI 	6.283185307
#define PID 	1.570796327

/* Choice of the billiard table */

#define D_RECTANGLE 0   /* rectangular domain */
#define D_ELLIPSE 1     /* elliptical domain */
#define D_STADIUM 2     /* stadium-shaped domain */
#define D_SINAI 3       /* Sinai billiard */
#define D_DIAMOND 4     /* diamond-shaped billiard */
#define D_TRIANGLE 5    /* triangular billiard */
#define D_ANNULUS 7     /* annulus */
#define D_POLYGON 8     /* polygon */
#define D_REULEAUX 9    /* Reuleaux and star shapes */
#define D_FLOWER 10     /* Bunimovich flower */
#define D_ALT_REU 11    /* alternating between star and Reuleaux */
#define D_ANGLE 12      /* angular sector */
#define D_LSHAPE 13     /* L-shaped billiard for conical singularity */
#define D_GENUSN 14     /* polygon with identified opposite sides */
#define D_PARABOLAS 15  /* polygon with parabolic sides */
#define D_PENROSE 16    /* Penrose solution to illumination problem */

#define D_CIRCLES 20     /* several circles */
#define D_CIRCLES_IN_RECT 21     /* several circles inside a rectangle */
#define D_CIRCLES_IN_GENUSN 22   /* several circles in polygon with identified opposite sides */
#define D_CIRCLES_IN_TORUS 23    /* several circles in a rectangle with periodic boundary conditions */

#define C_FOUR_CIRCLES 0  /* four circles almost touching each other */
#define C_SQUARE 1        /* square grid of circles */
#define C_HEX 2           /* hexagonal/triangular grid of circles */
#define C_TRI 21          /* equilateral triangular grid of circles */
#define C_GOLDEN_MEAN 3   /* golden mean grid */
#define C_GOLDEN_SPIRAL 4   /* golden spiral (sunflower) grid */
#define C_RAND_DISPLACED 5  /* randomly displaced square grid */
#define C_RAND_POISSON 6    /* random Poisson point process */
#define C_POISSON_DISC 7    /* Poisson disc sampling */

#define C_LASER 11          /* laser fight in a room of mirrors */
#define C_LASER_GENUSN 12   /* laser fight in a translation surface */

#define D_POLYLINE 30       /* polygonal line */
#define D_POLYLINE_ARCS 31  /* polygonal line and circular arcs */

#define P_RECTANGLE 0     /* rectangle (for test purposes) */
#define P_TOKARSKY 1      /* Tokarsky unilluminable room */
#define P_POLYRING 2      /* polygonal ring */
#define P_SIERPINSKI 3    /* sierpinski carpet */
#define P_VONKOCH 4       /* von Koch curve */
#define P_POLYGON 5       /* regular polygon, alternative for D_POLYGON */
#define P_TOKA_PRIME 6    /* Tokarsky room made of 86 triangles */
#define P_TREE 7          /* pine tree */
#define P_TOKA_NONSELF 8  /* Tokarsky non-self-unilluminable room */
#define P_MAZE 10         /* maze */
#define P_MAZE_DIAG 11    /* maze with 45 degrees angles */
#define P_MAZE_RANDOM 12  /* maze with randomized wall positions */
#define P_MAZE_CIRCULAR 13  /* circular maze */
#define P_MAZE_CIRC_SCATTERER 14  /* circular maze with scatterers */
#define P_MAZE_HEX 15     /* hexagonal maze */
#define P_MAZE_OCT 16     /* maze with octagonal and square cells */

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


