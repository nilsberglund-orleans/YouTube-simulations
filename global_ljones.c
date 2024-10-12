/* Global variables and parameters for lennardjones */

/* Basic math */

#define PI 	3.141592654
#define DPI 	6.283185307
#define PID 	1.570796327

/* shape of domain */

#define D_CIRCLES 20    /* several circles */
#define D_CIRCLES_IN_RECT 201   /* several circles in a rectangle */

#define NMAXCIRCLES 200000       /* total number of circles/polygons (must be at least NCX*NCY for square grid) */
#define MAXNEIGH 20         /* max number of neighbours kept in memory */
#define NMAXOBSTACLES 1000  /* max number of obstacles */
#define NMAXSEGMENTS 1000   /* max number of repelling segments */
#define NMAXGROUPS 50       /* max number of groups of segments */
#define NMAXCOLLISIONS 200000   /* max number of collisions */
#define NMAXPARTNERS 30     /* max number of partners in molecule */
#define NMAXPARTNERMOLECULES 30 /* max number of partners of a molecule */
#define NMAXPARTINCLUSTER 500   /* max number of particles in cluster */
#define NMAXBELTS 10        /* max number of conveyor belts */
#define NMAXSHOVELS 50      /* max number of shovels */

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

#define C_POOL_TABLE 20     /* pool table initial position */

#define C_ONE 97            /* one single circle, as for Sinai */
#define C_TWO 98            /* two concentric circles of different type */
#define C_NOTHING 99        /* no circle at all, for comparisons */

/* pattern of additional obstacles */
#define O_CORNERS 0         /* obstacles in the corners (for Boy b.c.) */
#define O_GALTON_BOARD 1    /* Galton board pattern */
#define O_GENUS_TWO 2       /* obstacles in corners of L-shape domeain (for genus 2 b.c.) */
#define O_POOL_TABLE 3      /* obstacles around pockets of pool table */
#define O_HLINE_HOLE_SPOKES 181    /* tips of spokes for S_HLINE_HOLE_SPOKES segment pattern */
#define O_CIRCLE 4          /* one circle at the origin */
#define O_FOUR_CIRCLES 5    /* four circles  */
#define O_HEX 6             /* hexagonal lattice */
#define O_SIDES 7           /* grid along the sides of the simulation rectangle */
#define O_SIDES_B 71        /* finer grid along the sides of the simulation rectangle */
#define O_SIEVE 8           /* obstacles form a sieve */
#define O_SIEVE_B 81        /* obstacles form a sieve, v2 */
#define O_SIEVE_LONG 82     /* obstacles form a long sieve */
#define O_SIEVE_LONG_B 83   /* obstacles form a long sieve, version with varying spacing */

/* pattern of additional repelling segments */
#define S_RECTANGLE 0       /* segments forming a rectangle */
#define S_CUP 1             /* segments forming a cup (for increasing gravity) */
#define S_HOURGLASS 2       /* segments forming an hour glass */
#define S_PENTA 3           /* segments forming a pentagon with 3 angles of 120° and 2 right angles */
#define S_CENTRIFUGE 4      /* segments forming "centrifuge" (polygon with radial segments) */
#define S_POLY_ELLIPSE 5    /* segments forming a polygonal approximation of an ellipse */
#define S_POOL_TABLE 6      /* pool table with pockets */
#define S_CENTRIFUGE_RND 7  /* segments forming centrifuge with more rounded bins */
#define S_CENTRIFUGE_LEAKY 8  /* segments forming centrifuge with rounded bins and holes */
#define S_CIRCLE_EXT 9      /* segments forming a repelling cicle */
#define S_ROCKET_NOZZLE 10  /* segments forming a rocket with bell-shaped nozzle */
#define S_ROCKET_NOZZLE_ROTATED 101 /* rotated version of rocket with bell-shaped nozzle */
#define S_TWO_ROCKETS 102     /* two different rockets, with nozzles specified by NOZZLE_SHAPE and NOZZLE_SHAPE_B */
#define S_TWO_CIRCLES_EXT 11  /* segments forming two repelling cicle */
#define S_DAM 12              /* segments forming a dam that can break */
#define S_DAM_WITH_HOLE 13    /* segments forming a dam in which a hole can open */
#define S_DAM_WITH_HOLE_AND_RAMP 14    /* segments forming a dam in which a hole can open */
#define S_MAZE 15           /* segments forming a maze */
#define S_MAZE_DIAG 151     /* segments forming a maze with diagonally opposed exits */
#define S_EXT_RECTANGLE 16  /* particles outside a rectangle */
#define S_DAM_BRICKS 17     /* dam made of several bricks */
#define S_HLINE_HOLE 18    /* horizontal line with a hole in the bottom */
#define S_HLINE_HOLE_SPOKES 181    /* horizontal line with a hole in the bottom and extra spokes */
#define S_HLINE_HOLE_SLOPED 182 /* slanted lines with a hole */
#define S_EXT_CIRCLE_RECT 19    /* particles outside a circle and a rectangle */
#define S_BIN_OPENING 20        /* bin containing particles opening at deactivation time */
#define S_BIN_LARGE 201         /* larger bin */
#define S_POLYGON_EXT 21        /* exterior of a regular polygon */
#define S_WEDGE_EXT 22          /* exterior of a wedge */ 
#define S_MIXER 23              /* exterior of a blender made of rectangles */
#define S_AIRFOIL 24            /* exterior of an air foil */
#define S_COANDA 25             /* wall for Coanda effect */
#define S_COANDA_SHORT 26       /* shorter wall for Coanda effect */
#define S_CYLINDER 27           /* walls at top and bottom, for cylindrical b.c. */
#define S_TREE 28               /* Christmas tree(s) */
#define S_CONE 29               /* cone */
#define S_CONVEYOR_BELT 30      /* conveyor belt */
#define S_TWO_CONVEYOR_BELTS 31 /* two angled conveyor belts */
#define S_PERIODIC_CONVEYORS 32 /* one wrapping belt, and one short horizontal belt */
#define S_TEST_CONVEYORS 321    /* test */
#define S_CONVEYOR_SHOVELS 33   /* conveyor belt with shovels */
#define S_CONVEYOR_MIXED 34     /* multiple conveyor belts with and without shovels */
#define S_CONVEYOR_SIEVE 35     /* conveyor belts for polygon sieve */
#define S_CONVEYOR_SIEVE_B 351  /* conveyor belts for polygon sieve, v2 with backward top conveyor */
#define S_CONVEYOR_SIEVE_LONG 352  /* conveyor belts for long polygon sieve */
#define S_MASS_SPECTROMETER 36  /* bins for mass spectrometer */
#define S_WIND_FORCE 361        /* bins for sorting by wind force */

/* particle interaction */

#define I_COULOMB 0         /* Coulomb force */
#define I_LENNARD_JONES 1   /* Lennard-Jones force */
#define I_LJ_DIRECTIONAL 2  /* Lennard-Jones with direction dependence of square symmetry */
#define I_LJ_PENTA 3        /* Lennard-Jones with pentagonal symmetry */
#define I_GOLDENRATIO 4     /* Lennard-Jones type with equilibria related by golden ratio */
#define I_LJ_DIPOLE 5       /* Lennard-Jones with a dipolar angle dependence */
#define I_LJ_QUADRUPOLE 6   /* Lennard-Jones with a quadropolar angle dependence */
#define I_LJ_WATER 7        /* model for water molecule */
#define I_VICSEK 8          /* Vicsek-type interaction */
#define I_VICSEK_REPULSIVE 9  /* Vicsek-type interaction with harmonic repulsion */
#define I_VICSEK_SPEED 10   /* Vicsek-type interaction with speed adjustment */
#define I_VICSEK_SHARK 11   /* Vicsek-type interaction with speed adjustment, and one shark */
#define I_COULOMB_LJ 12     /* Coulomb force regularised by Lennard-Jones repulsion */
#define I_COULOMB_PENTA 13  /* Lennard-Jones force with or without pentagonal symmetry depending on charge */
#define I_COULOMB_IMAGINARY 14  /* Coulomb interaction with "imaginary charge" */
#define I_DNA_CHARGED 15    /* Coulomb-type interaction between end points of DNA nucleotides */
#define I_DNA_CHARGED_B 151    /* stronger Coulomb-type interaction between end points of DNA nucleotides */
#define I_SEGMENT 16           /* harmonic interaction between segments */
#define I_SEGMENT_CHARGED 161  /* harmonic interaction between segments and Coulomb interaction between ends*/
#define I_POLYGON 17           /* harmonic interaction between regular polygons */
#define I_POLYGON_ALIGN 171    /* harmonic interaction between polygons with an aligning torque */

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
#define BC_RECTANGLE_WALL 9 /* rectangular container with vertical movable wall */
#define BC_KLEIN 11         /* Klein bottle (periodic with twisted vertical parts) */
#define BC_SCREEN_BINS 12   /* harmonic boundary conditions outside screen area plus "bins" (for Galton board) */
#define BC_BOY 13           /* Boy surface/projective plane (periodic with twisted horizontal and vertical parts) */
#define BC_GENUS_TWO 14     /* surface of genus 2, obtained by identifying opposite sides of an L shape */
#define BC_ABSORBING 20     /* "no-return" boundary conditions outside BC area */
#define BC_REFLECT_ABS 21   /* reflecting on lower boundary, and "no-return" boundary conditions outside BC area */
#define BC_REFLECT_ABS_BOTTOM 22    /* absorbing on lower boundary, and reflecting elsewhere */

/* Regions for partial thermostat couplings */

#define TH_VERTICAL 0       /* only particles at the right of x = PARTIAL_THERMO_SHIFT are coupled */
#define TH_INSEGMENT 1      /* only particles in region defined by segments are coupled */
#define TH_INBOX 2          /* only particles in a given box are coupled */
#define TH_LAYER 3          /* only particles above PARTIAL_THERMO_HEIGHT are coupled */
#define TH_LAYER_TYPE2 4    /* only particles above highest type 2 particle are coupled */
#define TH_RING 5           /* only particles outside disc of radius PARTIAL_THERMO_WIDTH are coupled */
#define TH_RING_EXPAND 6    /* only particles outside disc of radius changing from PARTIAL_THERMO_RIN to PARTIAL_THERMO_RFIN are coupled */
#define TH_INIT 7           /* only particles in region defined by INITXMIN, etc are coupled */
#define TH_THERMO 8         /* only particles in region defined by THERMOXMIN, etc are coupled */
#define TH_CONE 9           /* cone defined by S_CONE */

/* temperature schedules */

#define TS_EXPONENTIAL 0     /* temperature follows an exponential in time */
#define TS_CYCLING 1         /* temperature cycling */
#define TS_PERIODIC 2        /* periodic time dependence */
#define TS_LINEAR 3          /* linear time dependence */
#define TS_COSINE 4          /* periodic time dependence, cosine */
#define TS_EXPCOS 5          /* periodic time dependence, exponential of cosine */
#define TS_ASYM_EXPCOS 6     /* periodic time dependence, asymmetric exponential of cosine */
#define TS_ATAN 7            /* atan approaching asymptotic value */
#define TS_TANH 8            /* tanh approaching asymptotic value */

/* Gravity schedules */

#define G_INCREASE_RELEASE 1    /* slow increase and instant release */
#define G_INCREASE_DECREASE 2   /* slow increase and decrease */

/* Rocket shapes */

#define RCK_DISC 0      /* disc-shaped rocket */
#define RCK_RECT 1      /* rectangular rocket */
#define RCK_RECT_HAT 2  /* rectangular rocket with a hat */
#define RCK_RECT_BAR 3  /* rectangular rocket with a hat and a separating bar */

/* Nozzle shapes */

#define NZ_STRAIGHT 0   /* straight nozzle */
#define NZ_BELL 1       /* bell-shaped nozzle */
#define NZ_GLAS 2       /* glas-shaped nozzle */
#define NZ_CONE 3       /* cone-shaped nozzle */
#define NZ_TRUMPET 4    /* trumpet-shaped nozzle */
#define NZ_BROAD 5      /* broad straight nozzle */
#define NZ_DELAVAL 6    /* a type of de Laval nozzle */
#define NZ_NONE 99      /* no nozzle */

/* Types of chemical reactions */

#define CHEM_RPS 0      /* rock-paper-scissors reaction */
#define CHEM_AAB 1      /* reaction A + A -> B */
#define CHEM_ABC 2      /* reaction A + B -> C */
#define CHEM_A2BC 3      /* reaction 2A + B -> C */
#define CHEM_CATALYSIS 4    /* reaction 2A + C -> B + C */
#define CHEM_BAA 5      /* reaction B -> A + A (dissociation) */
#define CHEM_AABAA 6    /* reaction A + A <-> B (reversible) */
#define CHEM_POLYMER 7  /* reaction A + B -> C, A + C -> D, etc */
#define CHEM_POLYMER_DISS 8  /* polimerisation with dissociation */
#define CHEM_POLYMER_STEP 9  /* step growth polimerisation with dissociation */
#define CHEM_AUTOCATALYSIS 10   /* autocatalytic reaction 2A + B -> 2B */
#define CHEM_CATALYTIC_A2D 11   /* catalytic reaction A + B -> C, A + C -> B + D */
#define CHEM_ABCAB 12       /* reaction A + B <-> C (reversible) */
#define CHEM_ABCDABC 13     /* reactions A + B <-> C, A + C <-> D */
#define CHEM_BZ 14          /* simplified Belousov-Zhabotinski reaction with 6 types (Oregonator) */
#define CHEM_BRUSSELATOR 15 /* Brusselator oscillating reaction */
#define CHEM_ABDACBE 16     /* A + B -> D, A + C -> B + E */
#define CHEM_H2O_H_OH 20    /* H2O <-> H+ + OH- */
#define CHEM_2H2O_H3O_OH 21 /* 2 H2O <-> H3O+ + OH- */
#define CHEM_AGGREGATION 22 /* agregation of molecules coming close */
#define CHEM_AGGREGATION_CHARGE 23 /* agregation of charged molecules coming close */
#define CHEM_AGGREGATION_NNEIGH 24 /* agregation of molecules with limitation on neighbours */
#define CHEM_DNA 25         /* aggregation of DNA molecules */
#define CHEM_DNA_ALT 251    /* aggregation of DNA molecules with constraints on connections */
#define CHEM_DNA_DOUBLE 252 /* aggregation of DNA molecules with different ends */
#define CHEM_DNA_DSPLIT 253 /* aggregation/splitting of DNA molecules with different ends */
#define CHEM_DNA_BASE_SPLIT 254 /* aggregation/splitting of DNA molecules when base pairs don't match */
#define CHEM_DNA_ENZYME 255 /* aggregation/splitting of DNA molecules in presence of enzymes */
#define CHEM_DNA_ENZYME_REPAIR 256 /* aggregation/splitting of DNA molecules in presence of enzymes and additional repairing of bad connections */
#define CHEM_POLYGON_AGGREGATION 26 /* aggregation of polygons */
#define CHEM_POLYGON_CLUSTER 261    /* clustering of polygons into new clusters */
#define CHEM_POLYGON_ONECLUSTER 262 /* clustering of polygons, with only one cluster allowed */

/* Initial conditions for chemical reactions */

#define IC_NOTHING 99      /* do not change particle types */
#define IC_UNIFORM 0       /* all particles have type 1 */
#define IC_UNIFORM2 20     /* all particles have type 2 */
#define IC_RANDOM_UNIF 1   /* particle type chosen uniformly at random */
#define IC_RANDOM_TWO 2    /* particle type chosen randomly between 1 and 2, with TYPE_PROPORTION */
#define IC_CIRCLE 3        /* type 1 in a disc */
#define IC_CATALYSIS 4     /* mix of 1 and 2 in left half, only 1 in right half */
#define IC_LAYERS 5        /* layer of 2 below 1 */
#define IC_BZ 6            /* initial state for BZ reaction */
#define IC_SIGNX 7         /* type 1 or 2 depending on sign of x */
#define IC_TWOROCKETS 8    /* type 1 or 2 depending on rocket position */
#define IC_TWOROCKETS_TWOFUELS 9    /* type 1 and 2 or 1 and 3 depending on rocket */
#define IC_DNA_POLYMERASE 10    /* initial condition for DNA polymerase */
#define IC_DNA_POLYMERASE_REC 11    /* initial condition for DNA polymerase with recombination */

/* Initial conditions for option TWO_TYPES */

#define TTC_RANDOM 0        /* assign types randomly */
#define TTC_CHESSBOARD 1    /* assign types according to chessboard, works with hex initial config */
#define TTC_COANDA 2        /* type 1 in a band of width LAMBDA */

/* Initial speed distribution */

#define VI_RANDOM 0         /* random (Gaussian) initial speed distribution */
#define VI_COANDA 1         /* nonzero speed in a band of width LAMBDA */

/* Plot types */

#define P_KINETIC 0       /* colors represent kinetic energy of particles */
#define P_NEIGHBOURS 1    /* colors represent number of neighbours */
#define P_HEALTH 2        /* colors represent health (for SIR model) */
#define P_BONDS 3         /* draw lattice based on neighbours */
#define P_ANGLE 4         /* colors represent angle/spin of particle */
#define P_TYPE 5          /* colors represent type of particle */
#define P_DIRECTION 6     /* colors represent direction of velocity */
#define P_ANGULAR_SPEED 7 /* colors represent angular speed */
#define P_DIRECT_ENERGY 8 /* hue represents direction, luminosity represents energy */
#define P_DIFF_NEIGHB 9   /* colors represent number of neighbours of different type */
#define P_THERMOSTAT 10   /* colors show which particles are coupled to the thermostat */
#define P_INITIAL_POS 11  /* colors depend on initial position of particle */
#define P_NUMBER 12       /* colors depend on particle number */
#define P_EMEAN 13        /* averaged kinetic energy (with exponential damping) */
#define P_LOG_EMEAN 131   /* log of averaged kinetic energy (with exponential damping) */
#define P_EMEAN_DENSITY 132 /* averaged kinetic energy divided by the cluster size */
#define P_DIRECT_EMEAN 14 /* averaged version of P_DIRECT_ENERGY */
#define P_NOPARTICLE 15   /* particles are not drawn (only the links between them) */
#define P_NPARTNERS 16    /* number of partners */
#define P_CHARGE 17       /* colors represent charge */
#define P_MOL_ANGLE 18    /* orientation of molecule defined by partners */
#define P_CLUSTER 19      /* colors depend on connected component */
#define P_CLUSTER_SIZE 20 /* colors depend on size of connected component */
#define P_CLUSTER_SELECTED 21   /* colors show which clusters are slected for growth */
#define P_COLLISION 22    /* colors depend on number of collision/reaction */
#define P_RADIUS 23       /* colors depend on particle radius */

/* Rotation schedules */

#define ROT_SPEEDUP_SLOWDOWN 0  /* rotation speeds up and then slows down to zero */
#define ROT_BACK_FORTH 1        /* rotation goes in one direction and then back */

/* Initial position dependence types */

#define IP_X 0      /* color depends on x coordinate of initial position */
#define IP_Y 1      /* color depends on y coordinate of initial position */

/* Space dependence of magnetic field */

#define BF_CONST 0  /* constant magnetic field */
#define BF_SQUARE 1 /* magnetic field concentrated in square */

/* Interaction types for polyatomic molecules */

#define POLY_STAR 0     /* star-shaped graph (central molecule attracts outer ones) */
#define POLY_ALL 1      /* all-to-all coupling */
#define POLY_STAR_CHARGED 11     /* star-shaped graph with charged molecules */
#define POLY_POLYGON 12 /* polygonal shape */
#define POLY_WATER 2    /* star-shaped with a 120° separation between anions */
#define POLY_SOAP 3     /* polymers with all-to-all coupling and polar end */
#define POLY_SOAP_B 4   /* polymers with pairwise coupling and polar end */
#define POLY_SOAP_N 41  /* polymers with pairwise coupling and neutral polar end */
#define POLY_SOAP_NMIX 42   /* polymers mixing neutral polar and neutral end */
#define POLY_PLUSMINUS 5    /* polymers with ends of opposite charge */
#define POLY_HYDRA 6        /* star-shaped with longer arms */
#define POLY_HYDRA_RIGID 61 /* star-shaped with longer arms and rigid first ring */
#define POLY_DNA 7          /* simplified model for DNA */
#define POLY_DNA_ALT 71     /* simplified model for DNA with different short ends */
#define POLY_DNA_DOUBLE 72  /* simplified model for DNA with double ends for rigidity */
#define POLY_DNA_FLEX 73    /* simplified model for DNA with less backbone rigidity (beta) */
#define POLY_KITE 8         /* kite for kites and darts quasicrystal */
#define POLY_DART 81        /* dart for kites and darts quasicrystal */
#define POLY_SEG_POLYGON 9  /* polygon of segments */

/* Background color schemes */

#define BG_NONE 0       /* no background color */
#define BG_DENSITY 1    /* background color depends on number of particles */
#define BG_CHARGE 2     /* background color depends on charge density */
#define BG_EKIN 3       /* background color depends on kinetic energy */
#define BG_FORCE 4      /* background color depends on total force */

/* Particle add regions */

#define ADD_RECTANGLE 0     /* rectangular region, defined by ADDXMIN, etc */
#define ADD_RING 1          /* ring_shaped region, defined by ADDRMIN, ADDRMAX */

/* Color schemes */

#define C_LUM 0          /* color scheme modifies luminosity (with slow drift of hue) */
#define C_HUE 1          /* color scheme modifies hue */
#define C_PHASE 2        /* color scheme shows phase */
#define C_ONEDIM 3       /* use preset 1d color scheme (for Turbo, Viridis, Magma, Inferno, Plasma, Twilight) */
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

#define VICSEK_INT ((INTERACTION == I_VICSEK)||(INTERACTION == I_VICSEK_REPULSIVE)||(INTERACTION == I_VICSEK_SPEED)||(INTERACTION == I_VICSEK_SHARK))

typedef struct
{
    double xc, yc, radius;      /* center and radius of circle */
    double angle;               /* angle of particle's "spin" */
    short int active;           /* circle is active */
    double energy;              /* dissipated energy */
    double emean;               /* mean energy */
    double vx;                  /* x velocity of particle */
    double vy;                  /* y velocity of particle */
    double omega;               /* angular velocity of particle */
    double mass_inv;            /* inverse of particle mass */
    double charge;              /* electric charge */
    double inertia_moment_inv;  /* inverse of moment of inertia */
    double fx;                  /* x component of force on particle */
    double fy;                  /* y component of force on particle */
    double torque;              /* torque on particle */
    double damping;             /* factor in front of damping coefficient */
    double dirmean;             /* time averaged direction */
    int close_to_boundary;      /* has value 1 if particle is close to a boundary */
    short int thermostat;       /* whether particle is coupled to thermostat */
    int hashcell;               /* hash cell in which particle is located */
    int neighb;                 /* number of neighbours within given distance */
    int diff_neighb;            /* number of neighbours of different type */
    int hash_nneighb;           /* number of neighbours in hashgrid */
    int hashneighbour[9*HASHMAX];   /* particle numbers of neighbours in hashgrid */
    double deltax[9*HASHMAX];   /* relative position of neighbours */
    double deltay[9*HASHMAX];   /* relative position of neighbours */
    short int type;             /* type of particle, for mixture simulations */
    short int interaction;      /* type of interaction */
    double eq_dist;             /* equilibrium distance */
    double spin_range;          /* range of spin-spin interaction */
    double spin_freq;           /* angular frequency of spin-spin interaction */
    double color_hue;           /* color hue of particle, for P_INITIAL_POS plot type */
    int color_rgb[3];           /* RGB colors code of particle, for use in ljones_movie.c */
    int partner[NMAXPARTNERS];  /* partner particles for option PAIR_PARTICLES */
    short int npartners;        /* number of partner particles */
    double partner_eqd[NMAXPARTNERS];   /* equilibrium distances between partners */
    double partner_eqa[NMAXPARTNERS];   /* equilibrium angle between partners */
    int p0, p1;                 /* numbers of two first partners (for P_MOL_ANGLE color scheme) */
//     short int mol_angle;        /* for color scheme P_MOL_ANGLE */
    int cluster;                /* number of cluster */
    int cluster_color;          /* color of cluster */
    int cluster_size;           /* size of cluster */
    int molecule;               /* number of molecule */
    short int tested, cactive;  /* for cluster search */
    short int coulomb;          /* has value 1 if DNA-Coulomb interaction is attractive */
    short int added;            /* has value 1 if particle has been added */
    short int reactive;         /* has value 1 if particle can react */
    short int paired;           /* has value 1 if belongs to base-paired molecule */
    short int flip;             /* keeps track of which particles in a cluster are flipped by PI */
    int partner_molecule;       /* number of partner molecule */
    int collision;              /* number of collision */
} t_particle;

typedef struct
{
    short int active;           /* has value 1 if cluster is active */
    short int thermostat;       /* has value 1 if cluster is coupled to thermostat */
    short int selected;         /* has value 1 if cluster is selected to be able to grow */
    double xg, yg;              /* center of gravity */
    double vx, vy;              /* velocity of center of gravity */
    double angle;               /* orientation of cluster */
    double omega;               /* angular velocity of cluster */
    double mass, mass_inv;      /* mass of cluster and its inverse */
    double inertia_moment, inertia_moment_inv;   /* moment of inertia */
    double fx, fy, torque;      /* force and torque */
    double energy, emean;       /* energy and averaged energy */
    double dirmean;             /* time-averaged direction */
    int nparticles;             /* number of particles in cluster */
    int particle[NMAXPARTINCLUSTER];    /* list of particles in cluster */
//     int angle_ref;              /* reference particle for orientation */
} t_cluster;

typedef struct
{
    int number;                 /* total number of particles in cell */
    int particles[HASHMAX];     /* numbers of particles in cell */
    int nneighb;                /* number of neighbouring cells */
    int neighbour[9];           /* numbers of neighbouring cells */
    double x1, y1, x2, y2;      /* coordinates of hashcell corners */
    double hue1, hue2;          /* color hues */
    double charge;              /* charge of fixed obstacles */
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
    double xc0, yc0;            /* center of oscillation for option RATTLE_OBSTACLES */
    short int active;           /* circle is active */
    double charge;              /* charge of obstacle, for EM simulations */
    double omega0, omega;       /* speed of rotation */
    double angle;               /* angle of obstacle */
    short int oscillate;        /* has value 1 if the obstacles oscillates over time */
    int period;                 /* oscillation period */
    double amplitude, phase;    /* amplitude and phase of oscillation */
} t_obstacle;

typedef struct
{
    double x1, x2, y1, y2;      /* extremities of segment */
    double xc, yc;              /* mid-point of segment */
    double nx, ny;              /* normal vector */
    double nangle;              /* angle of normal vector */
    double c;                   /* constant term in cartesian eq nx*x + ny*y = c */
    double length;              /* length of segment */
    short int concave;          /* corner is concave, to add extra repelling force */
    short int cycle;            /* set to 1 if (x2, y2) is equal to (x1, y1) of next segment */
    short int group;            /* group to which segment belongs (for several obstacles) */
    double angle1, angle2;      /* angles in which concave corners repel */
    short int active;           /* segment is active */
    double x01, x02, y01, y02;  /* initial values of extremities, in case of rotation/translation */
    double nx0, ny0;            /* initial normal vector */
    double angle01, angle02;    /* initial values of angles in which concave corners repel */
    double fx, fy;              /* x and y-components of force on segment */
    double torque;              /* torque on segment with respect to its center */
    double pressure;            /* pressure acting on segement */
    double avrg_pressure;       /* time-averaged pressure */
    short int inactivate;       /* set to 1 for segment to become inactive at time SEGMENT_DEACTIVATION_TIME */
    short int conveyor;         /* set to 1 for segment to exert lateral force */
    double conveyor_speed;      /* speed of conveyor belt */
    short int align_torque;     /* set to 1 for segment to exert aligning torque */
} t_segment;

typedef struct
{
    double xc, yc;              /* center of circle */
} t_tracer;

typedef struct
{
    double xc, yc;              /* center of mass of obstacle */
    double angle;               /* orientation of obstacle */
    double vx, vy;              /* velocity of center of mass */
    double omega;               /* angular velocity */ 
    double mass;                /* mass of obstacle */
    double moment_inertia;      /* moment of inertia */
} t_group_segments;

typedef struct
{
    double xc, yc;              /* coordinates of centers of mass */
    double vx, vy;              /* velocities */
    double omega;               /* angular velocity */
} t_group_data;


typedef struct
{
    double x, y;                /* location of collision */
    int time;                   /* time since collision */
    int color;                  /* color hue in case of different collisions */
} t_collision;

typedef struct
{
    int nparticles;             /* number of particles */
    int particle[2*NPARTNERS+1];  /* list of particles */
    int npartners;              /* number of partner molecules */
    int partner[NMAXPARTNERMOLECULES];  /* list of partner molecules */
    int connection_type[NMAXPARTNERMOLECULES];  /* types of particles in connection */
    short int added;            /* has value 1 if molecule has been added */
} t_molecule;

typedef struct 
{
    double x1, y1, x2, y2;      /* positions of extremities */
    double width;               /* width of belt */
    double speed;               /* speed of conveyor belt */
    double position;            /* position of belt (needed for display of rotating parts) */
    double length;              /* distance between (x1,x2) and (y1,y2) */
    double angle;               /* angle of (x1,x2) - (y1,y2) */
    double tx, ty;              /* coordinates of tangent vector */
    int nshovels;               /* number of shovels */
    double shovel_pos[NMAXSHOVELS];    /* position od each shovel */
    int shovel_segment[NMAXSHOVELS];   /* first segment of each shovel */
} t_belt;

typedef struct
{
    int nactive;                /* number of active particles */
    double beta;                /* inverse temperature */
    double mean_energy;         /* mean energy */
    double krepel;              /* force constant */
    double xmincontainer, xmaxcontainer;      /* container size */
    double fboundary;           /* boundary force */
    double pressure;            /* pressure */
    double gravity;             /* gravity */
    double radius;              /* particle radius */
    double angle;               /* orientation of obstacle */
    double omega;               /* angular speed of obstacle */
    double bdry_fx, bdry_fy;    /* components of boundary force */
    double efield, bfield;      /* electric and magnetic field */
    double prop;                /* proportion of types */
    double thermo_radius;       /* radius of thermostat region */
} t_lj_parameters;



int frame_time = 0, ncircles, nobstacles, nsegments, ngroups = 1, counter = 0, nmolecules = 0, nbelts = 0, n_tracers = 0;
FILE *lj_log;

