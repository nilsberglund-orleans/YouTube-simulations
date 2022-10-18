/* Global variables and parameters for wave_billiard, heat and schrodinger */

/* Basic math */

#define PI 	3.141592654
#define DPI 	6.283185307
#define PID 	1.570796327

/* Lattice types */

#define BC_SQUARE_DIRICHLET 0   /* square lattice, Dirichlet boundary conditions */
#define BC_SQUARE_PERIODIC 1    /* square lattice, periodic boundary conditions */
#define BC_SQUARE_BOND_DIRICHLET 2  /* square lattice, bond percolation, Dirichlet b.c. */
#define BC_HEX_SITE_DIRICHLET 10    /* hexagonal lattice, site percolation, Dirichlet b.c. */
#define BC_HEX_BOND_DIRICHLET 12    /* hexagonal lattice, bond percolation, Dirichlet b.c. */
#define BC_TRIANGLE_SITE_DIRICHLET 20   /* triangular lattice, site percolation, Dirichlet b.c. */
#define BC_POISSON_DISC 50      /* Poisson disc (blue noise) lattice */

#define BC_CUBIC_DIRICHLET 100  /* cubic lattice */

/* Plot types */

#define PLOT_SQUARES 0      /* plot squares */
#define PLOT_SQUARE_BONDS 1 /* plot edges of square lattice */
#define PLOT_HEX 2          /* plot hexagons */
#define PLOT_TRIANGLE 3     /* plot triangles */
#define PLOT_HEX_BONDS 4    /* plot edges of hexagonal lattice */
#define PLOT_POISSON_DISC 5 /* plot Poisson disc process */

#define PLOT_CUBES 10       /* plot cubes */

/* 3D representation */

#define REP_AXO_3D 0        /* linear projection (axonometry) */
#define REP_PROJ_3D 1       /* projection on plane orthogonal to observer line of sight */


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


#define NSEG 50         /* number of segments of drawn circles */
#define MAX_NEIGHB 10   /* maximal number of neighbours */
#define NPOISSON_CANDIDATES 30  /* number of candidates in Poisson disc process */

typedef struct
{
    double proba;           /* probability of being open */
    short int nneighb;      /* number of neighbours */
    int nghb[MAX_NEIGHB];   /* indices of neighbours */
    short int open;         /* vertex is open, flooded */
    short int flooded;      /* vertex is open, flooded */
    short int active;       /* used in cluster algorithm */
    int cluster;            /* number of cluster, for search of all clusters */
    int previous_cluster;   /* number of cluster of previous p, to limit number of color changes */
    int tested;             /* 1 if the site has been tested for belonging to a cluster */
    double x, y;            /* coordinates of center, for Poisson disc process */
} t_perco;



