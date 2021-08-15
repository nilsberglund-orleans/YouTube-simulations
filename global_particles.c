double circlex[NMAXCIRCLES], circley[NMAXCIRCLES], circlerad[NMAXCIRCLES];      /* position and radius of circular scatters */
short int circleactive[NMAXCIRCLES];                                      /* tells which circular scatters are active */
int ncircles = NMAXCIRCLES;            /* actual number of circles, can be decreased e.g. for random patterns */

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
#define D_GENUSN 14     /* polygon with identifies opposite sides */

#define D_CIRCLES 20     /* several circles */
#define D_CIRCLES_IN_RECT 21     /* several circles inside a rectangle */
#define D_CIRCLES_IN_GENUSN 22   /* several circles in polygon with identified opposite sides */

#define C_FOUR_CIRCLES 0  /* four circles almost touching each other */
#define C_SQUARE 1        /* square grid of circles */
#define C_HEX 2           /* hexagonal/triangular grid of circles */
#define C_GOLDEN_MEAN 3   /* golden mean grid */

#define C_LASER 11          /* laser fight in a room of mirrors */
#define C_LASER_GENUSN 12   /* laser fight in a translation surface */