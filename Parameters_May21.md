### Parameter values for YouTube simulations ###

Created by **Nils Berglund** and optimized by **Marco Mancini**

C code for videos on YouTube Channel https://www.youtube.com/c/NilsBerglund

Below are parameter values used for different simulations, as well as initial conditions used in 
function animation. Some simulations use variants of the published code. The list is going to be 
updated gradually. Some constants have been moved to global files in later versions.


### 31 May 21 - Video #100: 100 000 particles in a concave hectogon billiard ###

**Program:** `particle_billiard.c`

**Initial condition in function `animation()`:** `init_drop_config(0.0, 0.5, 0.0, DPI, configs);`

```
#define MOVIE 1         /* set to 1 to generate movie */

#define WINWIDTH 	1280  /* window width */
#define WINHEIGHT 	720   /* window height */

#define XMIN -2.0
#define XMAX 2.0	/* x interval */
#define YMIN -1.125
#define YMAX 1.125	/* y interval for 9/16 aspect ratio */

/* Choice of the billiard table */

#define B_DOMAIN 9      /* choice of domain shape */

#define D_RECTANGLE 0   /* rectangular domain */
#define D_ELLIPSE 1     /* elliptical domain */
#define D_STADIUM 2     /* stadium-shaped domain */
#define D_SINAI 3       /* Sinai billiard */
#define D_DIAMOND 4     /* diamond-shaped billiard */
#define D_TRIANGLE 5    /* triangular billiard */
#define D_ANNULUS 7     /* annulus */
#define D_POLYGON 8     /* polygon */
#define D_REULEAUX 9    /* Reuleaux and star shapes */

#define LAMBDA 0.04588546 
#define MU 0.1          /* second parameter controlling shape of billiard */
#define FOCI 1          /* set to 1 to draw focal points of ellipse */
#define NPOLY 100             /* number of sides of polygon */
#define APOLY 1.0           /* angle by which to turn polygon, in units of Pi/2 */ 

#define RESAMPLE 0      /* set to 1 if particles should be added when dispersion too large */
#define DEBUG 0         /* draw trajectories, for debugging purposes */

#define NPART 100000         /* number of particles */
#define NPARTMAX 100000	/* maximal number of particles after resampling */
#define LMAX 0.01       /* minimal segment length triggering resampling */ 
#define DMIN 0.02       /* minimal distance to boundary for triggering resampling */ 
#define CYCLE 1         /* set to 1 for closed curve (start in all directions) */

#define NSTEPS 4000     /* number of frames of movie */
#define TIME 125       /* time between movie frames, for fluidity of real-time simulation */ 
#define DPHI 0.00001    /* integration step */
#define NVID 150         /* number of iterations between images displayed on screen */

#define NCOLORS -7      /* number of colors */
#define COLORSHIFT 220    /* hue of initial color */ 
#define NSEG 100        /* number of segments of boundary */
#define LENGTH 0.003     /* length of velocity vectors */
#define BILLIARD_WIDTH 4    /* width of billiard */
#define PARTICLE_WIDTH 2       /* width of particles */
#define FRONT_WIDTH 3       /* width of wave front */

#define BLACK 1             /* set to 1 for black background */
#define COLOR_OUTSIDE 0     /* set to 1 for colored outside */ 
#define OUTER_COLOR 270.0   /* color outside billiard */
#define PAINT_INT 1         /* set to 1 to paint interior in other color (for polygon) */

/* Decreasing TIME accelerates the animation and the movie                               */
/* For constant speed of movie, TIME*DPHI should be kept constant                        */
/* However, increasing DPHI too much deterioriates quality of simulation                 */
/* NVID tells how often a picture is drawn in the animation, increase it for faster anim */
/* For a good quality movie, take for instance TIME = 400, DPHI = 0.00005, NVID = 100    */

#define PAUSE 1000       /* number of frames after which to pause */
#define PSLEEP 1         /* sleep time during pause */
#define SLEEP1  1        /* initial sleeping time */
#define SLEEP2  1000      /* final sleeping time */

#define PI 	3.141592654
#define DPI 	6.283185307
#define PID 	1.570796327

```

### 30 May 21 - The quantum Ehrenfest model ###

**Program:** `schrodinger.c`

**Initial condition in function `animation()`:** `init_coherent_state(-1.2, 0.0, 8.0, 8.0, 0.08, phi, psi, xy_in);`

```
#define MOVIE 1         /* set to 1 to generate movie */

/* General geometrical parameters */

#define WINWIDTH 	1280  /* window width */
#define WINHEIGHT 	720   /* window height */

#define NX 640          /* number of grid points on x axis */
#define NY 360          /* number of grid points on y axis */

/* setting NX to WINWIDTH and NY to WINHEIGHT increases resolution */
/* but will multiply run time by 4                                 */

#define XMIN -2.0
#define XMAX 2.0	/* x interval */
#define YMIN -1.125
#define YMAX 1.125	/* y interval for 9/16 aspect ratio */

/* Choice of the billiard table */

#define B_DOMAIN 11     /* choice of domain shape */

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

#define LAMBDA 0.8	    /* parameter controlling the dimensions of domain */
#define MU 0.05	            /* parameter controlling the dimensions of domain */
#define NPOLY 6             /* number of sides of polygon */
#define APOLY 1.0           /* angle by which to turn polygon, in units of Pi/2 */ 
#define FOCI 1              /* set to 1 to draw focal points of ellipse */

/* You can add more billiard tables by adapting the functions */
/* xy_in_billiard and draw_billiard below */

/* Physical patameters of wave equation */

#define DT 0.00000002
#define HBAR 1.0

/* For debugging purposes only */
#define FLOOR 0         /* set to 1 to limit wave amplitude to VMAX */
#define VMAX 10.0       /* max value of wave amplitude */

/* Parameters for length and speed of simulation */

#define NSTEPS 5230     //7500      /* number of frames of movie */
#define NVID 1000          /* number of iterations between images displayed on screen */
#define NSEG 100         /* number of segments of boundary */

#define PAUSE 1000         /* number of frames after which to pause */
#define PSLEEP 1         /* sleep time during pause */
#define SLEEP1  1        /* initial sleeping time */
#define SLEEP2  1   /* final sleeping time */

/* Plot type */

#define PLOT 0

#define P_MODULE 0        /* plot module of wave function squared */
#define P_PHASE 1         /* plot phase of wave function */
#define P_REAL 2          /* plot real part */
#define P_IMAGINARY 3     /* plot imaginary part */

/* Boundary conditions */

#define B_COND 1

#define BC_DIRICHLET 0   /* Dirichlet boundary conditions */
#define BC_PERIODIC 1    /* periodic boundary conditions */
#define BC_ABSORBING 2   /* absorbing boundary conditions */

/* Color schemes */

#define BLACK 1          /* background */

#define COLOR_SCHEME 1   /* choice of color scheme */

#define C_LUM 0          /* color scheme modifies luminosity (with slow drift of hue) */
#define C_HUE 1          /* color scheme modifies hue */
#define C_PHASE 2        /* color scheme shows phase */

#define SCALE 1          /* set to 1 to adjust color scheme to variance of field */
#define SLOPE 1.0        /* sensitivity of color on wave amplitude */
#define ATTENUATION 0.0  /* exponential attenuation coefficient of contrast with time */

#define COLORHUE 260     /* initial hue of water color for scheme C_LUM */
#define COLORDRIFT 0.0   /* how much the color hue drifts during the whole simulation */
#define LUMMEAN 0.5      /* amplitude of luminosity variation for scheme C_LUM */
#define LUMAMP 0.3       /* amplitude of luminosity variation for scheme C_LUM */
#define HUEMEAN 150.0    /* mean value of hue for color scheme C_HUE */
#define HUEAMP -150.0      /* amplitude of variation of hue for color scheme C_HUE */

/* Basic math */

#define PI 	3.141592654
#define DPI 	6.283185307
#define PID 	1.570796327

```

### 29 May 21 - A quantum Sinai billiard: probability density ###

**Program:** `schrodinger.c`

**Initial condition in function `animation()`:** `init_coherent_state(-1.2, 0.0, 20.0, 0.0, 0.2, phi, psi, xy_in);`

```
#define MOVIE 1         /* set to 1 to generate movie */

/* General geometrical parameters */

#define WINWIDTH 	1280  /* window width */
#define WINHEIGHT 	720   /* window height */

#define NX 640          /* number of grid points on x axis */
#define NY 360          /* number of grid points on y axis */

/* setting NX to WINWIDTH and NY to WINHEIGHT increases resolution */
/* but will multiply run time by 4                                 */

#define XMIN -2.0
#define XMAX 2.0	/* x interval */
#define YMIN -1.125
#define YMAX 1.125	/* y interval for 9/16 aspect ratio */

/* Choice of the billiard table */

#define B_DOMAIN 3      /* choice of domain shape */

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

#define LAMBDA 0.2	    /* parameter controlling the dimensions of domain */
#define MU 0.05	            /* parameter controlling the dimensions of domain */
#define NPOLY 6             /* number of sides of polygon */
#define APOLY 1.0           /* angle by which to turn polygon, in units of Pi/2 */ 
#define FOCI 1              /* set to 1 to draw focal points of ellipse */

/* You can add more billiard tables by adapting the functions */
/* xy_in_billiard and draw_billiard below */

/* Physical patameters of wave equation */

#define DT 0.00000002
#define HBAR 1.0
/* The Courant number is given by c*DT/DX, where DT is the time step and DX the lattice spacing */
/* The physical damping coefficient is given by GAMMA/(DT)^2 */
/* Increasing COURANT speeds up the simulation, but decreases accuracy */
/* For similar wave forms, COURANT^2*GAMMA should be kept constant */

/* For debugging purposes only */
#define FLOOR 0         /* set to 1 to limit wave amplitude to VMAX */
#define VMAX 10.0       /* max value of wave amplitude */

/* Parameters for length and speed of simulation */

#define NSTEPS 3500     //7500      /* number of frames of movie */
#define NVID 750          /* number of iterations between images displayed on screen */
#define NSEG 100         /* number of segments of boundary */

#define PAUSE 1000         /* number of frames after which to pause */
#define PSLEEP 1         /* sleep time during pause */
#define SLEEP1  1        /* initial sleeping time */
#define SLEEP2  1   /* final sleeping time */

/* Plot type */

#define PLOT 1

#define P_MODULE 0        /* plot module of wave function squared */
#define P_PHASE 1         /* plot phase of wave function */
#define P_REAL 2          /* plot real part */
#define P_IMAGINARY 3     /* plot imaginary part */

/* Boundary conditions */

#define B_COND 1

#define BC_DIRICHLET 0   /* Dirichlet boundary conditions */
#define BC_PERIODIC 1    /* periodic boundary conditions */
#define BC_ABSORBING 2   /* absorbing boundary conditions */

/* Color schemes */

#define BLACK 1          /* background */

#define COLOR_SCHEME 1   /* choice of color scheme */

#define C_LUM 0          /* color scheme modifies luminosity (with slow drift of hue) */
#define C_HUE 1          /* color scheme modifies hue */
#define C_PHASE 2        /* color scheme shows phase */

#define SCALE 1          /* set to 1 to adjust color scheme to variance of field */
#define SLOPE 1.2        /* sensitivity of color on wave amplitude */
#define ATTENUATION 0.0  /* exponential attenuation coefficient of contrast with time */

#define COLORHUE 260     /* initial hue of water color for scheme C_LUM */
#define COLORDRIFT 0.0   /* how much the color hue drifts during the whole simulation */
#define LUMMEAN 0.5      /* amplitude of luminosity variation for scheme C_LUM */
#define LUMAMP 0.3       /* amplitude of luminosity variation for scheme C_LUM */
#define HUEMEAN 130.0    /* mean value of hue for color scheme C_HUE */
#define HUEAMP -150.0      /* amplitude of variation of hue for color scheme C_HUE */

/* Basic math */

#define PI 	3.141592654
#define DPI 	6.283185307
#define PID 	1.570796327

```

### 28 May 21 - The Vitruvian pond: wave front in a 5-pointed star ###

**Program:** `drop_billiard.c`

**Initial condition in function `animation()`:** `init_drop_config(0.0, 0.1, 0.0, DPI, configs);`

```
#define MOVIE 1         /* set to 1 to generate movie */

#define WINWIDTH 	1280  /* window width */
#define WINHEIGHT 	720   /* window height */

#define XMIN -1.8
#define XMAX 1.8	/* x interval */
#define YMIN -0.91
#define YMAX 1.115	/* y interval for 9/16 aspect ratio */

/* Choice of the billiard table */

#define B_DOMAIN 9      /* choice of domain shape */

#define D_RECTANGLE 0   /* rectangular domain */
#define D_ELLIPSE 1     /* elliptical domain */
#define D_STADIUM 2     /* stadium-shaped domain */
#define D_SINAI 3       /* Sinai billiard */
#define D_DIAMOND 4     /* diamond-shaped billiard */
#define D_TRIANGLE 5    /* triangular billiard */
#define D_ANNULUS 7     /* annulus */
#define D_POLYGON 8     /* polygon */
#define D_REULEAUX 9    /* Reuleaux and star shapes */

#define LAMBDA 1.124950941	/* sin(36°)/sin(31.5°) for 5-star shape with 45° angles */
#define MU 0.1          /* second parameter controlling shape of billiard */
#define FOCI 1          /* set to 1 to draw focal points of ellipse */
#define NPOLY 5             /* number of sides of polygon */
#define APOLY -1.0           /* angle by which to turn polygon, in units of Pi/2 */ 

#define RESAMPLE 0      /* set to 1 if particles should be added when dispersion too large */

#define NPART 300000	/* number of particles */
#define NPARTMAX 500000	/* maximal number of particles after resampling */
#define LMAX 0.01       /* minimal segment length triggering resampling */ 
#define LPERIODIC 1.0   /* lines longer than this are not drawn (useful for Sinai billiard) */
#define LCUT 2000.0        /* controls the max size of segments not considered as being cut */
#define DMIN 0.02       /* minimal distance to boundary for triggering resampling */ 
#define CYCLE 0         /* set to 1 for closed curve (start in all directions) */
#define ORDER_COLORS 1  /* set to 1 if colors should be drawn in order */ 

#define NSTEPS 4000         /* number of frames of movie */
#define TIME 15             /* time between movie frames, for fluidity of real-time simulation */ 
#define DPHI 0.0001         /* integration step */
#define NVID 10             /* number of iterations between images displayed on screen */
#define NCOLORS 10          /* number of colors */
#define COLORSHIFT 200      /* hue of initial color */ 
#define NSEG 100            /* number of segments of boundary */
#define BILLIARD_WIDTH 6    /* width of billiard */
#define FRONT_WIDTH 3       /* width of wave front */

#define BLACK 1             /* set to 1 for black background */
#define COLOR_OUTSIDE 0     /* set to 1 for colored outside */ 
#define OUTER_COLOR 300.0   /* color outside billiard */
#define PAINT_INT 1         /* set to 1 to paint interior in other color (for polygon) */

/* Decreasing TIME accelerates the animation and the movie               */
/* For constant speed of movie, TIME*DPHI should be kept constant        */
/* However, increasing DPHI too much deterioriates quality of simulation */
/* For a good quality movie, take for instance TIME = 50, DPHI = 0.0002  */

#define PAUSE 1000          /* number of frames after which to pause */
#define PSLEEP 1         /* sleep time during pause */
#define SLEEP1  1        /* initial sleeping time */
#define SLEEP2  100      /* final sleeping time */

#define PI 	3.141592654
#define DPI 	6.283185307
#define PID 	1.570796327

```

### 27 May 21 - A quantum Sinai billiard, phase evolution ###

**Program:** `schrodinger.c`

**Initial condition in function `animation()`:** `init_coherent_state(-1.2, 0.0, 20.0, 0.0, 0.2, phi, psi, xy_in);`

```
#define MOVIE 1         /* set to 1 to generate movie */

/* General geometrical parameters */

#define WINWIDTH 	1280  /* window width */
#define WINHEIGHT 	720   /* window height */

#define NX 640          /* number of grid points on x axis */
#define NY 360          /* number of grid points on y axis */

/* setting NX to WINWIDTH and NY to WINHEIGHT increases resolution */
/* but will multiply run time by 4                                 */

#define XMIN -2.0
#define XMAX 2.0	/* x interval */
#define YMIN -1.125
#define YMAX 1.125	/* y interval for 9/16 aspect ratio */

/* Choice of the billiard table */

#define B_DOMAIN 3      /* choice of domain shape */

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

#define LAMBDA 0.2	    /* parameter controlling the dimensions of domain */
#define MU 0.05	            /* parameter controlling the dimensions of domain */
#define NPOLY 6             /* number of sides of polygon */
#define APOLY 1.0           /* angle by which to turn polygon, in units of Pi/2 */ 
#define FOCI 1              /* set to 1 to draw focal points of ellipse */

/* You can add more billiard tables by adapting the functions */
/* xy_in_billiard and draw_billiard below */

/* Physical patameters of wave equation */

#define DT 0.00000002
#define HBAR 1.0
/* The Courant number is given by c*DT/DX, where DT is the time step and DX the lattice spacing */
/* The physical damping coefficient is given by GAMMA/(DT)^2 */
/* Increasing COURANT speeds up the simulation, but decreases accuracy */
/* For similar wave forms, COURANT^2*GAMMA should be kept constant */

/* For debugging purposes only */
#define FLOOR 0         /* set to 1 to limit wave amplitude to VMAX */
#define VMAX 10.0       /* max value of wave amplitude */

/* Parameters for length and speed of simulation */

#define NSTEPS 3500      /* number of frames of movie */
#define NVID 750         /* number of iterations between images displayed on screen */
#define NSEG 100         /* number of segments of boundary */

#define PAUSE 1000         /* number of frames after which to pause */
#define PSLEEP 1         /* sleep time during pause */
#define SLEEP1  1        /* initial sleeping time */
#define SLEEP2  1   /* final sleeping time */

/* Plot type */

#define PLOT 1

#define P_MODULE 0        /* plot module of wave function squared */
#define P_PHASE 1         /* plot phase of wave function */
#define P_REAL 2          /* plot real part */
#define P_IMAGINARY 3     /* plot imaginary part */

/* Boundary conditions */

#define B_COND 1

#define BC_DIRICHLET 0   /* Dirichlet boundary conditions */
#define BC_PERIODIC 1    /* periodic boundary conditions */
#define BC_ABSORBING 2   /* absorbing boundary conditions */

/* Color schemes */

#define BLACK 1          /* background */

#define COLOR_SCHEME 1   /* choice of color scheme */

#define C_LUM 0          /* color scheme modifies luminosity (with slow drift of hue) */
#define C_HUE 1          /* color scheme modifies hue */
#define C_PHASE 2        /* color scheme shows phase */

#define SCALE 1          /* set to 1 to adjust color scheme to variance of field */
#define SLOPE 1.2        /* sensitivity of color on wave amplitude */
#define ATTENUATION 0.0  /* exponential attenuation coefficient of contrast with time */

#define COLORHUE 260     /* initial hue of water color for scheme C_LUM */
#define COLORDRIFT 0.0   /* how much the color hue drifts during the whole simulation */
#define LUMMEAN 0.5      /* amplitude of luminosity variation for scheme C_LUM */
#define LUMAMP 0.3       /* amplitude of luminosity variation for scheme C_LUM */
#define HUEMEAN 130.0    /* mean value of hue for color scheme C_HUE */
#define HUEAMP -150.0      /* amplitude of variation of hue for color scheme C_HUE */

/* Basic math */

#define PI 	3.141592654
#define DPI 	6.283185307
#define PID 	1.570796327

```

### 26 May 21 - 50 pence for your thoughts - Billiard in a 7-sided Reuleaux t̶r̶i̶a̶n̶g̶l̶e̶ polygon ###

**Program:** `particle_billiard.c`

**Initial condition in function `animation()`:** `alpha = 3.0*PI/7.0;  init_sym_drop_config(-0.99, 0.0, -alpha, alpha, configs);`

```
#define MOVIE 0         /* set to 1 to generate movie */

#define WINWIDTH 	1280  /* window width */
#define WINHEIGHT 	720   /* window height */

#define XMIN -2.0
#define XMAX 2.0	/* x interval */
#define YMIN -1.125
#define YMAX 1.125	/* y interval for 9/16 aspect ratio */

/* Choice of the billiard table */

#define B_DOMAIN 9      /* choice of domain shape */

#define D_RECTANGLE 0   /* rectangular domain */
#define D_ELLIPSE 1     /* elliptical domain */
#define D_STADIUM 2     /* stadium-shaped domain */
#define D_SINAI 3       /* Sinai billiard */
#define D_DIAMOND 4     /* diamond-shaped billiard */
#define D_TRIANGLE 5    /* triangular billiard */
#define D_ANNULUS 7     /* annulus */
#define D_POLYGON 8     /* polygon */
#define D_REULEAUX 9    /* Reuleaux and star shapes */

#define LAMBDA -1.949855824	/* 7-sided Reuleaux triangle */
#define MU 0.1          /* second parameter controlling shape of billiard */
#define FOCI 1          /* set to 1 to draw focal points of ellipse */
#define NPOLY 7             /* number of sides of polygon */
#define APOLY 0.0           /* angle by which to turn polygon, in units of Pi/2 */ 

#define RESAMPLE 0      /* set to 1 if particles should be added when dispersion too large */
#define DEBUG 0         /* draw trajectories, for debugging purposes */

/* Simulation parameters */

#define NPART 10000     /* number of particles */
#define NPARTMAX 100000	/* maximal number of particles after resampling */
#define LMAX 0.01       /* minimal segment length triggering resampling */ 
#define DMIN 0.02       /* minimal distance to boundary for triggering resampling */ 
#define CYCLE 1         /* set to 1 for closed curve (start in all directions) */

#define NSTEPS 9000     /* number of frames of movie */
#define TIME 125       /* time between movie frames, for fluidity of real-time simulation */ 
#define DPHI 0.00001    /* integration step */
#define NVID 150         /* number of iterations between images displayed on screen */

/* Decreasing TIME accelerates the animation and the movie                               */
/* For constant speed of movie, TIME*DPHI should be kept constant                        */
/* However, increasing DPHI too much deterioriates quality of simulation                 */
/* NVID tells how often a picture is drawn in the animation, increase it for faster anim */
/* For a good quality movie, take for instance TIME = 400, DPHI = 0.00005, NVID = 100    */

/* Colors and other graphical parameters */

#define NCOLORS -7      /* number of colors */
#define COLORSHIFT 220    /* hue of initial color */ 
#define NSEG 100        /* number of segments of boundary */
#define LENGTH 0.01     /* length of velocity vectors */
#define BILLIARD_WIDTH 4    /* width of billiard */
#define PARTICLE_WIDTH 2       /* width of particles */
#define FRONT_WIDTH 3       /* width of wave front */

#define BLACK 1             /* set to 1 for black background */
#define COLOR_OUTSIDE 0     /* set to 1 for colored outside */ 
#define OUTER_COLOR 270.0   /* color outside billiard */
#define PAINT_INT 1         /* set to 1 to paint interior in other color (for polygon) */


#define PAUSE 1000       /* number of frames after which to pause */
#define PSLEEP 1         /* sleep time during pause */
#define SLEEP1  1        /* initial sleeping time */
#define SLEEP2  1000      /* final sleeping time */

#define PI 	3.141592654
#define DPI 	6.283185307
#define PID 	1.570796327

```

### 25 May 21 - Young's double-slit experiment for a quantum particle, phase representation ###

**Program:** `schrodinger.c`

**Initial condition in function `animation()`:** `init_coherent_state(-1.2, 0.0, 20.0, 0.0, 0.25, phi, psi, xy_in);`

```
#define MOVIE 1         /* set to 1 to generate movie */

/* General geometrical parameters */

#define WINWIDTH 	1280  /* window width */
#define WINHEIGHT 	720   /* window height */

#define NX 640          /* number of grid points on x axis */
#define NY 360          /* number of grid points on y axis */

/* setting NX to WINWIDTH and NY to WINHEIGHT increases resolution */
/* but will multiply run time by 4                                 */

#define XMIN -2.0
#define XMAX 2.0	/* x interval */
#define YMIN -1.125
#define YMAX 1.125	/* y interval for 9/16 aspect ratio */

/* Choice of the billiard table */

#define B_DOMAIN 9      /* choice of domain shape */

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

#define LAMBDA 0.3	    /* parameter controlling the dimensions of domain */
#define MU 0.05	            /* parameter controlling the dimensions of domain */
#define NPOLY 6             /* number of sides of polygon */
#define APOLY 1.0           /* angle by which to turn polygon, in units of Pi/2 */ 
#define FOCI 1              /* set to 1 to draw focal points of ellipse */

/* You can add more billiard tables by adapting the functions */
/* xy_in_billiard and draw_billiard below */

/* Physical patameters of wave equation */

#define DT 0.00000002
#define HBAR 1.0
/* The Courant number is given by c*DT/DX, where DT is the time step and DX the lattice spacing */
/* The physical damping coefficient is given by GAMMA/(DT)^2 */
/* Increasing COURANT speeds up the simulation, but decreases accuracy */
/* For similar wave forms, COURANT^2*GAMMA should be kept constant */

/* For debugging purposes only */
#define FLOOR 0         /* set to 1 to limit wave amplitude to VMAX */
#define VMAX 10.0       /* max value of wave amplitude */

/* Parameters for length and speed of simulation */

#define NSTEPS 4000     //7500      /* number of frames of movie */
#define NVID 750          /* number of iterations between images displayed on screen */
#define NSEG 100         /* number of segments of boundary */

#define PAUSE 1000         /* number of frames after which to pause */
#define PSLEEP 1         /* sleep time during pause */
#define SLEEP1  1        /* initial sleeping time */
#define SLEEP2  1   /* final sleeping time */

/* Plot type */

#define PLOT 1

#define P_MODULE 0        /* plot module of wave function squared */
#define P_PHASE 1         /* plot phase of wave function */
#define P_REAL 2          /* plot real part */
#define P_IMAGINARY 3     /* plot imaginary part */

/* Boundary conditions */

#define B_COND 1

#define BC_DIRICHLET 0   /* Dirichlet boundary conditions */
#define BC_PERIODIC 1    /* periodic boundary conditions */
#define BC_ABSORBING 2   /* absorbing boundary conditions */

/* Color schemes */

#define BLACK 1          /* background */

#define COLOR_SCHEME 1   /* choice of color scheme */

#define C_LUM 0          /* color scheme modifies luminosity (with slow drift of hue) */
#define C_HUE 1          /* color scheme modifies hue */
#define C_PHASE 2        /* color scheme shows phase */

#define SCALE 1          /* set to 1 to adjust color scheme to variance of field */
#define SLOPE 1.5        /* sensitivity of color on wave amplitude */
#define ATTENUATION 0.0  /* exponential attenuation coefficient of contrast with time */

#define COLORHUE 260     /* initial hue of water color for scheme C_LUM */
#define COLORDRIFT 0.0   /* how much the color hue drifts during the whole simulation */
#define LUMMEAN 0.5      /* amplitude of luminosity variation for scheme C_LUM */
#define LUMAMP 0.3       /* amplitude of luminosity variation for scheme C_LUM */
#define HUEMEAN 130.0    /* mean value of hue for color scheme C_HUE */
#define HUEAMP -130.0      /* amplitude of variation of hue for color scheme C_HUE */

/* Basic math */

#define PI 	3.141592654
#define DPI 	6.283185307
#define PID 	1.570796327

```

### 24 May 21 - The redemption of the pentagon ###

**Program:** `drop_billiard.c`

**Initial condition in function `animation()`:** `init_drop_config(0.0, 0.0, PID, DPI+PID, configs);`

```
#define MOVIE 1         /* set to 1 to generate movie */

#define WINWIDTH 	1280  /* window width */
#define WINHEIGHT 	720   /* window height */

#define XMIN -1.8
#define XMAX 1.8	/* x interval */
#define YMIN -0.91
#define YMAX 1.115	/* y interval for 9/16 aspect ratio */

/* Choice of the billiard table */

#define B_DOMAIN 9      /* choice of domain shape */

#define D_RECTANGLE 0   /* rectangular domain */
#define D_ELLIPSE 1     /* elliptical domain */
#define D_STADIUM 2     /* stadium-shaped domain */
#define D_SINAI 3       /* Sinai billiard */
#define D_DIAMOND 4     /* diamond-shaped billiard */
#define D_TRIANGLE 5    /* triangular billiard */
#define D_ANNULUS 7     /* annulus */
#define D_POLYGON 8     /* polygon */
#define D_REULEAUX 9    /* Reuleaux and star shapes */

#define LAMBDA 3.75738973	/* sin(36°)/sin(9°) for 5-star shape with 90° angles */
#define MU 0.1          /* second parameter controlling shape of billiard */
#define FOCI 1          /* set to 1 to draw focal points of ellipse */
#define NPOLY 5             /* number of sides of polygon */
#define APOLY -1.0           /* angle by which to turn polygon, in units of Pi/2 */ 

#define RESAMPLE 0      /* set to 1 if particles should be added when dispersion too large */

#define NPART 100000	/* number of particles */
#define NPARTMAX 100000	/* maximal number of particles after resampling */
#define LMAX 0.01       /* minimal segment length triggering resampling */ 
#define LPERIODIC 1.0   /* lines longer than this are not drawn (useful for Sinai billiard) */
#define LCUT 1000.0        /* controls the max size of segments not considered as being cut */
#define DMIN 0.02       /* minimal distance to boundary for triggering resampling */ 
#define CYCLE 0         /* set to 1 for closed curve (start in all directions) */

#define NSTEPS 4500         /* number of frames of movie */
#define TIME 20             /* time between movie frames, for fluidity of real-time simulation */ 
#define DPHI 0.0001         /* integration step */
#define NVID 10             /* number of iterations between images displayed on screen */
#define NCOLORS 10          /* number of colors */
#define COLORSHIFT 270      /* hue of initial color */ 
#define NSEG 100            /* number of segments of boundary */
#define BILLIARD_WIDTH 4    /* width of billiard */
#define FRONT_WIDTH 3       /* width of wave front */

#define BLACK 1             /* set to 1 for black background */
#define COLOR_OUTSIDE 0     /* set to 1 for colored outside */ 
#define OUTER_COLOR 300.0   /* color outside billiard */
#define PAINT_INT 1         /* set to 1 to paint interior in other color (for polygon) */

/* Decreasing TIME accelerates the animation and the movie               */
/* For constant speed of movie, TIME*DPHI should be kept constant        */
/* However, increasing DPHI too much deterioriates quality of simulation */
/* For a good quality movie, take for instance TIME = 50, DPHI = 0.0002  */

#define PAUSE 1000          /* number of frames after which to pause */
#define PSLEEP 1         /* sleep time during pause */
#define SLEEP1  1        /* initial sleeping time */
#define SLEEP2  100      /* final sleeping time */

#define PI 	3.141592654
#define DPI 	6.283185307
#define PID 	1.570796327

```

### 23 May 21 - Young's double-slit experiment for a quantum particle ###

**Program:** `schrodinger.c`

**Initial condition in function `animation()`:** `init_coherent_state(-1.2, 0.0, 20.0, 0.0, 0.25, phi, psi, xy_in);`

```
#define MOVIE 1         /* set to 1 to generate movie */

/* General geometrical parameters */

#define WINWIDTH 	1280  /* window width */
#define WINHEIGHT 	720   /* window height */

#define NX 640          /* number of grid points on x axis */
#define NY 360          /* number of grid points on y axis */

/* setting NX to WINWIDTH and NY to WINHEIGHT increases resolution */
/* but will multiply run time by 4                                 */

#define XMIN -2.0
#define XMAX 2.0	/* x interval */
#define YMIN -1.125
#define YMAX 1.125	/* y interval for 9/16 aspect ratio */

/* Choice of the billiard table */

#define B_DOMAIN 9      /* choice of domain shape */

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

#define LAMBDA 0.3	    /* parameter controlling the dimensions of domain */
#define MU 0.05	            /* parameter controlling the dimensions of domain */
#define NPOLY 6             /* number of sides of polygon */
#define APOLY 1.0           /* angle by which to turn polygon, in units of Pi/2 */ 
#define FOCI 1              /* set to 1 to draw focal points of ellipse */

/* You can add more billiard tables by adapting the functions */
/* xy_in_billiard and draw_billiard below */

/* Physical patameters of wave equation */

#define DT 0.00000002
#define HBAR 1.0
/* The Courant number is given by c*DT/DX, where DT is the time step and DX the lattice spacing */
/* The physical damping coefficient is given by GAMMA/(DT)^2 */
/* Increasing COURANT speeds up the simulation, but decreases accuracy */
/* For similar wave forms, COURANT^2*GAMMA should be kept constant */

/* For debugging purposes only */
#define FLOOR 0         /* set to 1 to limit wave amplitude to VMAX */
#define VMAX 10.0       /* max value of wave amplitude */

/* Parameters for length and speed of simulation */

#define NSTEPS 4000      /* number of frames of movie */
#define NVID 750         /* number of iterations between images displayed on screen */
#define NSEG 100         /* number of segments of boundary */

#define PAUSE 1000         /* number of frames after which to pause */
#define PSLEEP 1         /* sleep time during pause */
#define SLEEP1  1        /* initial sleeping time */
#define SLEEP2  1   /* final sleeping time */

/* Plot type */

#define PLOT 0

#define P_MODULE 0        /* plot module of wave function squared */
#define P_PHASE 1         /* plot phase of wave function */
#define P_REAL 2          /* plot real part */
#define P_IMAGINARY 3     /* plot imaginary part */

/* Boundary conditions */

#define B_COND 1

#define BC_DIRICHLET 0   /* Dirichlet boundary conditions */
#define BC_PERIODIC 1    /* periodic boundary conditions */
#define BC_ABSORBING 2   /* absorbing boundary conditions */

/* Color schemes */

#define BLACK 1          /* background */

#define COLOR_SCHEME 1   /* choice of color scheme */

#define C_LUM 0          /* color scheme modifies luminosity (with slow drift of hue) */
#define C_HUE 1          /* color scheme modifies hue */
#define C_PHASE 2        /* color scheme shows phase */

#define SCALE 1          /* set to 1 to adjust color scheme to variance of field */
#define SLOPE 1.5        /* sensitivity of color on wave amplitude */
#define ATTENUATION 0.0  /* exponential attenuation coefficient of contrast with time */

#define COLORHUE 260     /* initial hue of water color for scheme C_LUM */
#define COLORDRIFT 0.0   /* how much the color hue drifts during the whole simulation */
#define LUMMEAN 0.5      /* amplitude of luminosity variation for scheme C_LUM */
#define LUMAMP 0.3       /* amplitude of luminosity variation for scheme C_LUM */
#define HUEMEAN 130.0    /* mean value of hue for color scheme C_HUE */
#define HUEAMP -130.0      /* amplitude of variation of hue for color scheme C_HUE */

/* Basic math */

#define PI 	3.141592654
#define DPI 	6.283185307
#define PID 	1.570796327

```

### 22 May 21 - Billiard in a Reuleaux triangle ###

**Program:** `particle_billiard.c`

**Initial condition in function `animation()`:** `init_drop_config(0.0, 0.0, 0.0, DPI, configs);`

```
#define MOVIE 1         /* set to 1 to generate movie */

#define WINWIDTH 	1280  /* window width */
#define WINHEIGHT 	720   /* window height */

#define XMIN -1.8
#define XMAX 1.8	/* x interval */
#define YMIN -0.91
#define YMAX 1.115	/* y interval for 9/16 aspect ratio */

/* Choice of the billiard table */

#define B_DOMAIN 9      /* choice of domain shape */

#define D_RECTANGLE 0   /* rectangular domain */
#define D_ELLIPSE 1     /* elliptical domain */
#define D_STADIUM 2     /* stadium-shaped domain */
#define D_SINAI 3       /* Sinai billiard */
#define D_DIAMOND 4     /* diamond-shaped billiard */
#define D_TRIANGLE 5    /* triangular billiard */
#define D_ANNULUS 7     /* annulus */
#define D_POLYGON 8     /* polygon */
#define D_REULEAUX 9    /* Reuleaux and star shapes */

#define LAMBDA -1.73205080756888	/* -sqrt(3) for triangle tiling plane */
#define MU 0.1          /* second parameter controlling shape of billiard */
#define FOCI 1          /* set to 1 to draw focal points of ellipse */
#define NPOLY 3             /* number of sides of polygon */
#define APOLY 0.3333333333333333333           /* angle by which to turn polygon, in units of Pi/2 */ 

#define RESAMPLE 0      /* set to 1 if particles should be added when dispersion too large */
#define DEBUG 0         /* draw trajectories, for debugging purposes */

#define NPART 10000         /* number of particles */
#define NPARTMAX 50000	/* maximal number of particles after resampling */
#define LMAX 0.01       /* minimal segment length triggering resampling */ 
#define DMIN 0.02       /* minimal distance to boundary for triggering resampling */ 
#define CYCLE 1         /* set to 1 for closed curve (start in all directions) */

#define NSTEPS 9000     /* number of frames of movie */
#define TIME 600       /* time between movie frames, for fluidity of real-time simulation */ 
#define DPHI 0.00001    /* integration step */
#define NVID 75         /* number of iterations between images displayed on screen */

#define NCOLORS 14      /* number of colors */
#define COLORSHIFT 0    /* hue of initial color */ 
#define NSEG 100        /* number of segments of boundary */
#define LENGTH 0.007     /* length of velocity vectors */

#define BLACK 1             /* set to 1 for black background */
#define COLOR_OUTSIDE 0     /* set to 1 for colored outside */ 
#define OUTER_COLOR 270.0   /* color outside billiard */
#define PAINT_INT 1         /* set to 1 to paint interior in other color (for polygon) */

/* Decreasing TIME accelerates the animation and the movie                               */
/* For constant speed of movie, TIME*DPHI should be kept constant                        */
/* However, increasing DPHI too much deterioriates quality of simulation                 */
/* NVID tells how often a picture is drawn in the animation, increase it for faster anim */
/* For a good quality movie, take for instance TIME = 400, DPHI = 0.00005, NVID = 100    */

#define PAUSE 1000       /* number of frames after which to pause */
#define PSLEEP 1         /* sleep time during pause */
#define SLEEP1  1        /* initial sleeping time */
#define SLEEP2  1000      /* final sleeping time */

#define PI 	3.141592654
#define DPI 	6.283185307
#define PID 	1.570796327

```

### 21 May 21 - Phasers locked on target: Phase evolution in the quantum hexagonal billiard ###

**Program:** `schrodinger.c`

**Initial condition in function `animation()`:** `init_coherent_state(0.0, 0.0, 0.0, 10.0, 0.1, phi, psi, xy_in);`

```
#define MOVIE 1         /* set to 1 to generate movie */

/* General geometrical parameters */

#define WINWIDTH 	1280  /* window width */
#define WINHEIGHT 	720   /* window height */

#define NX 640          /* number of grid points on x axis */
#define NY 360          /* number of grid points on y axis */

/* setting NX to WINWIDTH and NY to WINHEIGHT increases resolution */
/* but will multiply run time by 4                                 */

#define XMIN -2.0
#define XMAX 2.0	/* x interval */
#define YMIN -1.125
#define YMAX 1.125	/* y interval for 9/16 aspect ratio */

/* Choice of the billiard table */

#define B_DOMAIN 8      /* choice of domain shape */

#define D_RECTANGLE 0   /* rectangular domain */
#define D_ELLIPSE 1     /* elliptical domain */
#define D_STADIUM 2     /* stadium-shaped domain */
#define D_SINAI 3       /* Sinai billiard */
#define D_DIAMOND 4     /* diamond-shaped billiard */
#define D_TRIANGLE 5    /* triangular billiard */
#define D_FLAT 6        /* flat interface */
#define D_ANNULUS 7     /* annulus */
#define D_POLYGON 8     /* polygon */

#define LAMBDA 1.5	    /* parameter controlling the dimensions of domain */
#define NPOLY 6             /* number of sides of polygon */
#define APOLY 1.0           /* angle by which to turn polygon, in units of Pi/2 */ 
#define FOCI 1              /* set to 1 to draw focal points of ellipse */

/* You can add more billiard tables by adapting the functions */
/* xy_in_billiard and draw_billiard below */

/* Physical patameters of wave equation */

#define DT 0.00000002
#define HBAR 1.0
/* The Courant number is given by c*DT/DX, where DT is the time step and DX the lattice spacing */
/* The physical damping coefficient is given by GAMMA/(DT)^2 */
/* Increasing COURANT speeds up the simulation, but decreases accuracy */
/* For similar wave forms, COURANT^2*GAMMA should be kept constant */

/* For debugging purposes only */
#define FLOOR 0         /* set to 1 to limit wave amplitude to VMAX */
#define VMAX 10.0       /* max value of wave amplitude */

/* Parameters for length and speed of simulation */

#define NSTEPS 6000     //7500      /* number of frames of movie */
#define NVID 750          /* number of iterations between images displayed on screen */
#define NSEG 100         /* number of segments of boundary */

#define PAUSE 1000         /* number of frames after which to pause */
#define PSLEEP 1         /* sleep time during pause */
#define SLEEP1  1        /* initial sleeping time */
#define SLEEP2  1   /* final sleeping time */

/* Plot type */

#define PLOT 1

#define P_MODULE 0        /* plot module of wave function squared */
#define P_PHASE 1         /* plot phase of wave function */
#define P_REAL 2          /* plot real part */
#define P_IMAGINARY 3     /* plot imaginary part */

/* Color schemes */

#define BLACK 1          /* background */

#define COLOR_SCHEME 1   /* choice of color scheme */

#define C_LUM 0          /* color scheme modifies luminosity (with slow drift of hue) */
#define C_HUE 1          /* color scheme modifies hue */
#define C_PHASE 2        /* color scheme shows phase */

#define SCALE 1          /* set to 1 to adjust color scheme to variance of field */
#define SLOPE 1.0        /* sensitivity of color on wave amplitude */
#define ATTENUATION 0.0  /* exponential attenuation coefficient of contrast with time */

#define COLORHUE 260     /* initial hue of water color for scheme C_LUM */
#define COLORDRIFT 0.0   /* how much the color hue drifts during the whole simulation */
#define LUMMEAN 0.5      /* amplitude of luminosity variation for scheme C_LUM */
#define LUMAMP 0.3       /* amplitude of luminosity variation for scheme C_LUM */
#define HUEMEAN 200.0    /* mean value of hue for color scheme C_HUE */
#define HUEAMP -120.0      /* amplitude of variation of hue for color scheme C_HUE */

/* Basic math */

#define PI 	3.141592654
#define DPI 	6.283185307
#define PID 	1.570796327

```

### 20 May 21 - How to distinguish an icosikaipentagon from a circle ###

**Program:** `drop_billiard.c`

**Initial condition in function `animation()`:** `init_drop_config(0.0, 0.0, 0.0, DPI, configs);`

```
#define MOVIE 1         /* set to 1 to generate movie */

#define WINWIDTH 	1280  /* window width */
#define WINHEIGHT 	720   /* window height */

#define XMIN -2.0
#define XMAX 2.0	/* x interval */
#define YMIN -1.125
#define YMAX 1.125	/* y interval for 9/16 aspect ratio */

/* Choice of the billiard table */

#define B_DOMAIN 8      /* choice of domain shape */

#define D_RECTANGLE 0   /* rectangular domain */
#define D_ELLIPSE 1     /* elliptical domain */
#define D_STADIUM 2     /* stadium-shaped domain */
#define D_SINAI 3       /* Sinai billiard */
#define D_DIAMOND 4     /* diamond-shaped billiard */
#define D_TRIANGLE 5    /* triangular billiard */
#define D_ANNULUS 7     /* annulus */
#define D_POLYGON 8     /* polygon */
#define D_REULEAUX 9    /* Reuleaux and star shapes */

#define LAMBDA 1.0	/* parameter controlling shape of billiard */
#define MU 0.1          /* second parameter controlling shape of billiard */
#define FOCI 1          /* set to 1 to draw focal points of ellipse */
#define NPOLY 25             /* number of sides of polygon */
#define APOLY 1.0           /* angle by which to turn polygon, in units of Pi/2 */ 

#define RESAMPLE 0      /* set to 1 if particles should be added when dispersion too large */

#define NPART 40000	/* number of particles */
#define NPARTMAX 100000	/* maximal number of particles after resampling */
#define LMAX 0.01       /* minimal segment length triggering resampling */ 
#define LPERIODIC 1.0   /* lines longer than this are not drawn (useful for Sinai billiard) */
#define LCUT 2.0        /* controls the max size of segments not considered as being cut */
#define DMIN 0.02       /* minimal distance to boundary for triggering resampling */ 
#define CYCLE 0         /* set to 1 for closed curve (start in all directions) */

#define NSTEPS 8000     /* number of frames of movie */
#define TIME 50         /* time between movie frames, for fluidity of real-time simulation */ 
#define DPHI 0.0001     /* integration step */
#define NVID 25         /* number of iterations between images displayed on screen */
#define NCOLORS -6      /* number of colors */
#define COLORSHIFT 240  /* hue of initial color */ 
#define NSEG 100        /* number of segments of boundary */

#define BLACK 1             /* set to 1 for black background */
#define COLOR_OUTSIDE 1     /* set to 1 for colored outside */ 
#define OUTER_COLOR 300.0   /* color outside billiard */
#define PAINT_INT 1         /* set to 1 to paint interior in other color (for polygon) */

/* Decreasing TIME accelerates the animation and the movie               */
/* For constant speed of movie, TIME*DPHI should be kept constant        */
/* However, increasing DPHI too much deterioriates quality of simulation */
/* For a good quality movie, take for instance TIME = 50, DPHI = 0.0002  */

#define PAUSE 1000          /* number of frames after which to pause */
#define PSLEEP 1         /* sleep time during pause */
#define SLEEP1  1        /* initial sleeping time */
#define SLEEP2  100      /* final sleeping time */

#define PI 	3.141592654
#define DPI 	6.283185307
#define PID 	1.570796327

```

### 19 May 21 - The hexagonal quantum billiard: probability density ###

**Program:** `schrodinger.c`

**Initial condition in function `animation()`:** `init_coherent_state(0.0, 0.0, 10.0, 0.0, 0.1, phi, psi, xy_in);`

```
#define MOVIE 0         /* set to 1 to generate movie */

/* General geometrical parameters */

#define WINWIDTH 	1280  /* window width */
#define WINHEIGHT 	720   /* window height */

#define NX 640          /* number of grid points on x axis */
#define NY 360          /* number of grid points on y axis */

/* setting NX to WINWIDTH and NY to WINHEIGHT increases resolution */
/* but will multiply run time by 4                                 */

#define XMIN -2.0
#define XMAX 2.0	/* x interval */
#define YMIN -1.125
#define YMAX 1.125	/* y interval for 9/16 aspect ratio */

/* Choice of the billiard table */

#define B_DOMAIN 8      /* choice of domain shape */

#define D_RECTANGLE 0   /* rectangular domain */
#define D_ELLIPSE 1     /* elliptical domain */
#define D_STADIUM 2     /* stadium-shaped domain */
#define D_SINAI 3       /* Sinai billiard */
#define D_DIAMOND 4     /* diamond-shaped billiard */
#define D_TRIANGLE 5    /* triangular billiard */
#define D_FLAT 6        /* flat interface */
#define D_ANNULUS 7     /* annulus */
#define D_POLYGON 8     /* polygon */

#define LAMBDA 1.5	    /* parameter controlling the dimensions of domain */
#define NPOLY 6             /* number of sides of polygon */
#define APOLY 1.0           /* angle by which to turn polygon, in units of Pi/2 */ 
#define FOCI 1              /* set to 1 to draw focal points of ellipse */

/* You can add more billiard tables by adapting the functions */
/* xy_in_billiard and draw_billiard below */

/* Physical patameters of wave equation */

#define DT 0.00000005
#define HBAR 1.0
/* The Courant number is given by c*DT/DX, where DT is the time step and DX the lattice spacing */
/* The physical damping coefficient is given by GAMMA/(DT)^2 */
/* Increasing COURANT speeds up the simulation, but decreases accuracy */
/* For similar wave forms, COURANT^2*GAMMA should be kept constant */

/* For debugging purposes only */
#define FLOOR 0         /* set to 1 to limit wave amplitude to VMAX */
#define VMAX 10.0       /* max value of wave amplitude */

/* Parameters for length and speed of simulation */

#define NSTEPS 6000     //7500      /* number of frames of movie */
#define NVID 750          /* number of iterations between images displayed on screen */
#define NSEG 100         /* number of segments of boundary */

#define PAUSE 1000         /* number of frames after which to pause */
#define PSLEEP 1         /* sleep time during pause */
#define SLEEP1  1        /* initial sleeping time */
#define SLEEP2  1   /* final sleeping time */

/* Plot type */

#define PLOT 0

#define P_MODULE 0        /* plot module of wave function squared */
#define P_PHASE 1         /* plot phase of wave function */
#define P_REAL 2          /* plot real part */
#define P_IMAGINARY 3     /* plot imaginary part */

/* Color schemes */

#define BLACK 1          /* background */

#define COLOR_SCHEME 1   /* choice of color scheme */

#define C_LUM 0          /* color scheme modifies luminosity (with slow drift of hue) */
#define C_HUE 1          /* color scheme modifies hue */
#define C_PHASE 2        /* color scheme shows phase */

#define SCALE 1          /* set to 1 to adjust color scheme to variance of field */
#define SLOPE 1.0        /* sensitivity of color on wave amplitude */
#define ATTENUATION 0.0  /* exponential attenuation coefficient of contrast with time */

#define COLORHUE 260     /* initial hue of water color for scheme C_LUM */
#define COLORDRIFT 0.0   /* how much the color hue drifts during the whole simulation */
#define LUMMEAN 0.5      /* amplitude of luminosity variation for scheme C_LUM */
#define LUMAMP 0.3       /* amplitude of luminosity variation for scheme C_LUM */
#define HUEMEAN 200.0    /* mean value of hue for color scheme C_HUE */
#define HUEAMP -120.0      /* amplitude of variation of hue for color scheme C_HUE */

/* Basic math */

#define PI 	3.141592654
#define DPI 	6.283185307
#define PID 	1.570796327

```

### 18 May 21 - Diamonds are a nerd's best friend - A finite horizon Sinai billiard ###

**Program:** `particle_billiard.c`

**Initial condition in function `animation()`:** `init_line_config(0.0, 0.1, 0.0, 0.3, 0.0, configs);`

```
#define MOVIE 1         /* set to 1 to generate movie */

#define WINWIDTH 	1280  /* window width */
#define WINHEIGHT 	720   /* window height */

#define XMIN -2.0
#define XMAX 2.0	/* x interval */
#define YMIN -1.125
#define YMAX 1.125	/* y interval for 9/16 aspect ratio */

/* Choice of the billiard table */

#define B_DOMAIN 9      /* choice of domain shape */

#define D_RECTANGLE 0   /* rectangular domain */
#define D_ELLIPSE 1     /* elliptical domain */
#define D_STADIUM 2     /* stadium-shaped domain */
#define D_SINAI 3       /* Sinai billiard */
#define D_DIAMOND 4     /* diamond-shaped billiard */
#define D_TRIANGLE 5    /* triangular billiard */
#define D_ANNULUS 7     /* annulus */
#define D_POLYGON 8     /* polygon */
#define D_REULEAUX 9    /* Reuleaux and star shapes */

#define LAMBDA 10.0	/* parameter controlling shape of billiard */
#define MU 0.1          /* second parameter controlling shape of billiard */
#define FOCI 1          /* set to 1 to draw focal points of ellipse */
#define NPOLY 4             /* number of sides of polygon */
#define APOLY 0.5           /* angle by which to turn polygon, in units of Pi/2 */ 

#define RESAMPLE 0      /* set to 1 if particles should be added when dispersion too large */
#define DEBUG 0         /* draw trajectories, for debugging purposes */

#define NPART 5000      /* number of particles */
#define NPARTMAX 50000	/* maximal number of particles after resampling */
#define LMAX 0.01       /* minimal segment length triggering resampling */ 
#define DMIN 0.02       /* minimal distance to boundary for triggering resampling */ 
#define CYCLE 1         /* set to 1 for closed curve (start in all directions) */

#define NSTEPS 10000     /* number of frames of movie */
#define TIME 500       /* time between movie frames, for fluidity of real-time simulation */ 
#define DPHI 0.00001    /* integration step */
#define NVID 75         /* number of iterations between images displayed on screen */

#define NCOLORS 14      /* number of colors */
#define COLORSHIFT 0    /* hue of initial color */ 
#define NSEG 100        /* number of segments of boundary */
#define LENGTH 0.007     /* length of velocity vectors */

#define BLACK 1             /* set to 1 for black background */
#define COLOR_OUTSIDE 0     /* set to 1 for colored outside */ 
#define OUTER_COLOR 270.0   /* color outside billiard */
#define PAINT_INT 1         /* set to 1 to paint interior in other color (for polygon) */

/* Decreasing TIME accelerates the animation and the movie                               */
/* For constant speed of movie, TIME*DPHI should be kept constant                        */
/* However, increasing DPHI too much deterioriates quality of simulation                 */
/* NVID tells how often a picture is drawn in the animation, increase it for faster anim */
/* For a good quality movie, take for instance TIME = 400, DPHI = 0.00005, NVID = 100    */

#define PAUSE 1000       /* number of frames after which to pause */
#define PSLEEP 1         /* sleep time during pause */
#define SLEEP1  1        /* initial sleeping time */
#define SLEEP2  1000      /* final sleeping time */

#define PI 	3.141592654
#define DPI 	6.283185307
#define PID 	1.570796327

```

### 17 May 21 - The hexagonal quantum billiard: real part of the wave function ###

**Program:** `schrodinger.c`

**Initial condition in function `animation()`:** `init_coherent_state(0.0, 0.0, 0.0, 10.0, 0.1, phi, psi, xy_in);`

```
#define MOVIE 1         /* set to 1 to generate movie */

/* General geometrical parameters */

#define WINWIDTH 	1280  /* window width */
#define WINHEIGHT 	720   /* window height */

#define NX 640          /* number of grid points on x axis */
#define NY 360          /* number of grid points on y axis */

/* setting NX to WINWIDTH and NY to WINHEIGHT increases resolution */
/* but will multiply run time by 4                                 */

#define XMIN -2.0
#define XMAX 2.0	/* x interval */
#define YMIN -1.125
#define YMAX 1.125	/* y interval for 9/16 aspect ratio */

/* Choice of the billiard table */

#define B_DOMAIN 8      /* choice of domain shape */

#define D_RECTANGLE 0   /* rectangular domain */
#define D_ELLIPSE 1     /* elliptical domain */
#define D_STADIUM 2     /* stadium-shaped domain */
#define D_SINAI 3       /* Sinai billiard */
#define D_DIAMOND 4     /* diamond-shaped billiard */
#define D_TRIANGLE 5    /* triangular billiard */
#define D_FLAT 6        /* flat interface */
#define D_ANNULUS 7     /* annulus */
#define D_POLYGON 8     /* polygon */

#define LAMBDA 1.5	    /* parameter controlling the dimensions of domain */
#define NPOLY 6             /* number of sides of polygon */
#define APOLY 1.0           /* angle by which to turn polygon, in units of Pi/2 */ 
#define FOCI 1              /* set to 1 to draw focal points of ellipse */

/* You can add more billiard tables by adapting the functions */
/* xy_in_billiard and draw_billiard below */

/* Physical patameters of wave equation */

#define DT 0.00000002
#define HBAR 1.0
/* The Courant number is given by c*DT/DX, where DT is the time step and DX the lattice spacing */
/* The physical damping coefficient is given by GAMMA/(DT)^2 */
/* Increasing COURANT speeds up the simulation, but decreases accuracy */
/* For similar wave forms, COURANT^2*GAMMA should be kept constant */

/* For debugging purposes only */
#define FLOOR 0         /* set to 1 to limit wave amplitude to VMAX */
#define VMAX 10.0       /* max value of wave amplitude */

/* Parameters for length and speed of simulation */

#define NSTEPS 6000     //7500      /* number of frames of movie */
#define NVID 750          /* number of iterations between images displayed on screen */
#define NSEG 100         /* number of segments of boundary */

#define PAUSE 1000         /* number of frames after which to pause */
#define PSLEEP 1         /* sleep time during pause */
#define SLEEP1  1        /* initial sleeping time */
#define SLEEP2  1   /* final sleeping time */

/* Plot type */

#define PLOT 2

#define P_MODULE 0        /* plot module of wave function squared */
#define P_PHASE 1         /* plot phase of wave function */
#define P_REAL 2          /* plot real part */
#define P_IMAGINARY 3     /* plot imaginary part */

/* Color schemes */

#define BLACK 1          /* background */

#define COLOR_SCHEME 1   /* choice of color scheme */

#define C_LUM 0          /* color scheme modifies luminosity (with slow drift of hue) */
#define C_HUE 1          /* color scheme modifies hue */
#define C_PHASE 2        /* color scheme shows phase */

#define SCALE 1          /* set to 1 to adjust color scheme to variance of field */
#define SLOPE 1.0        /* sensitivity of color on wave amplitude */
#define ATTENUATION 0.0  /* exponential attenuation coefficient of contrast with time */

#define COLORHUE 260     /* initial hue of water color for scheme C_LUM */
#define COLORDRIFT 0.0   /* how much the color hue drifts during the whole simulation */
#define LUMMEAN 0.5      /* amplitude of luminosity variation for scheme C_LUM */
#define LUMAMP 0.3       /* amplitude of luminosity variation for scheme C_LUM */
#define HUEMEAN 200.0    /* mean value of hue for color scheme C_HUE */
#define HUEAMP -120.0      /* amplitude of variation of hue for color scheme C_HUE */

/* Basic math */

#define PI 	3.141592654
#define DPI 	6.283185307
#define PID 	1.570796327

```

### 16 May 21 - Hexing the circle: are hexagons really bestagons? ###

**Program:** `drop_billiard.c`

**Initial condition in function `animation()`:** `init_drop_config(0.0, 0.0, 0.0, DPI, configs);`

```
#define MOVIE 1         /* set to 1 to generate movie */

#define WINWIDTH 	1280  /* window width */
#define WINHEIGHT 	720   /* window height */

#define XMIN -2.0
#define XMAX 2.0	/* x interval */
#define YMIN -1.125
#define YMAX 1.125	/* y interval for 9/16 aspect ratio */

/* Choice of the billiard table */

#define B_DOMAIN 8      /* choice of domain shape */

#define D_RECTANGLE 0   /* rectangular domain */
#define D_ELLIPSE 1     /* elliptical domain */
#define D_STADIUM 2     /* stadium-shaped domain */
#define D_SINAI 3       /* Sinai billiard */
#define D_DIAMOND 4     /* diamond-shaped billiard */
#define D_TRIANGLE 5    /* triangular billiard */
#define D_ANNULUS 7     /* annulus */
#define D_POLYGON 8     /* polygon */

#define LAMBDA 1.0	/* parameter controlling shape of billiard */
#define MU 0.1          /* second parameter controlling shape of billiard */
#define FOCI 1          /* set to 1 to draw focal points of ellipse */
#define NPOLY 6             /* number of sides of polygon */
#define APOLY 1.0           /* angle by which to turn polygon, in units of Pi/2 */ 

#define RESAMPLE 0      /* set to 1 if particles should be added when dispersion too large */

#define NPART 40000	/* number of particles */
#define NPARTMAX 100000	/* maximal number of particles after resampling */
#define LMAX 0.01       /* minimal segment length triggering resampling */ 
#define LPERIODIC 1.0   /* lines longer than this are not drawn (useful for Sinai billiard) */
#define LCUT 2.0        /* controls the max size of segments not considered as being cut */
#define DMIN 0.02       /* minimal distance to boundary for triggering resampling */ 
#define CYCLE 0         /* set to 1 for closed curve (start in all directions) */

#define NSTEPS 5175     /* number of frames of movie */
#define TIME 50         /* time between movie frames, for fluidity of real-time simulation */ 
#define DPHI 0.0001     /* integration step */
#define NVID 25         /* number of iterations between images displayed on screen */
#define NCOLORS -12      /* number of colors */
#define COLORSHIFT 240  /* hue of initial color */ 
#define NSEG 100        /* number of segments of boundary */

#define BLACK 1             /* set to 1 for black background */
#define COLOR_OUTSIDE 1     /* set to 1 for colored outside */ 
#define OUTER_COLOR 270.0   /* color outside billiard */
#define PAINT_INT 1         /* set to 1 to paint interior in other color (for polygon) */

/* Decreasing TIME accelerates the animation and the movie               */
/* For constant speed of movie, TIME*DPHI should be kept constant        */
/* However, increasing DPHI too much deterioriates quality of simulation */
/* For a good quality movie, take for instance TIME = 50, DPHI = 0.0002  */

#define PAUSE 1000          /* number of frames after which to pause */
#define PSLEEP 1         /* sleep time during pause */
#define SLEEP1  1        /* initial sleeping time */
#define SLEEP2  100      /* final sleeping time */

#define PI 	3.141592654
#define DPI 	6.283185307
#define PID 	1.570796327

```

### 15 May 21 - Scarred for life, no compensation - Schrödinger's equation in a stadium ###

**Program:** `schrodinger.c`

**Initial condition in function `animation()`:** `init_coherent_state(0.0, 0.0, 0.0, 5.0, 0.03, phi, psi, xy_in);`

```
#define MOVIE 1         /* set to 1 to generate movie */

/* General geometrical parameters */

#define WINWIDTH 	1280  /* window width */
#define WINHEIGHT 	720   /* window height */

#define NX 640          /* number of grid points on x axis */
#define NY 360          /* number of grid points on y axis */

/* setting NX to WINWIDTH and NY to WINHEIGHT increases resolution */
/* but will multiply run time by 4                                 */

#define XMIN -2.0
#define XMAX 2.0	/* x interval */
#define YMIN -1.125
#define YMAX 1.125	/* y interval for 9/16 aspect ratio */

/* Choice of the billiard table */

#define B_DOMAIN 2      /* choice of domain shape */

#define D_RECTANGLE 0   /* rectangular domain */
#define D_ELLIPSE 1     /* elliptical domain */
#define D_STADIUM 2     /* stadium-shaped domain */
#define D_SINAI 3       /* Sinai billiard */
#define D_DIAMOND 4     /* diamond-shaped billiard */
#define D_TRIANGLE 5    /* triangular billiard */
#define D_FLAT 6        /* flat interface */
#define D_ANNULUS 7     /* annulus */
#define D_POLYGON 8     /* polygon */

#define LAMBDA 1.5	    /* parameter controlling the dimensions of domain */
#define NPOLY 6             /* number of sides of polygon */
#define APOLY 1.0           /* angle by which to turn polygon, in units of Pi/2 */ 
#define FOCI 1              /* set to 1 to draw focal points of ellipse */

/* You can add more billiard tables by adapting the functions */
/* xy_in_billiard and draw_billiard below */

/* Physical patameters of wave equation */

#define DT 0.00000002
#define HBAR 1.0
/* The Courant number is given by c*DT/DX, where DT is the time step and DX the lattice spacing */
/* The physical damping coefficient is given by GAMMA/(DT)^2 */
/* Increasing COURANT speeds up the simulation, but decreases accuracy */
/* For similar wave forms, COURANT^2*GAMMA should be kept constant */

/* For debugging purposes only */
#define FLOOR 0         /* set to 1 to limit wave amplitude to VMAX */
#define VMAX 10.0       /* max value of wave amplitude */

/* Parameters for length and speed of simulation */

#define NSTEPS 6000     //7500      /* number of frames of movie */
#define NVID 500          /* number of iterations between images displayed on screen */
#define NSEG 100         /* number of segments of boundary */

#define PAUSE 1000         /* number of frames after which to pause */
#define PSLEEP 1         /* sleep time during pause */
#define SLEEP1  1        /* initial sleeping time */
#define SLEEP2  1   /* final sleeping time */

/* Plot type */

#define PLOT 0

#define P_MODULE 0        /* plot module of wave function squared */
#define P_PHASE 1         /* plot phase of wave function */
#define P_REAL 2          /* plot real part */
#define P_IMAGINARY 3     /* plot imaginary part */

/* Color schemes */

#define BLACK 1          /* background */

#define COLOR_SCHEME 1   /* choice of color scheme */

#define C_LUM 0          /* color scheme modifies luminosity (with slow drift of hue) */
#define C_HUE 1          /* color scheme modifies hue */
#define C_PHASE 2        /* color scheme shows phase */

#define SCALE 1          /* set to 1 to adjust color scheme to variance of field */
#define SLOPE 1.0        /* sensitivity of color on wave amplitude */
#define ATTENUATION 0.0  /* exponential attenuation coefficient of contrast with time */

#define COLORHUE 260     /* initial hue of water color for scheme C_LUM */
#define COLORDRIFT 0.0   /* how much the color hue drifts during the whole simulation */
#define LUMMEAN 0.5      /* amplitude of luminosity variation for scheme C_LUM */
#define LUMAMP 0.3       /* amplitude of luminosity variation for scheme C_LUM */
#define HUEMEAN 160.0    /* mean value of hue for color scheme C_HUE */
#define HUEAMP -120.0      /* amplitude of variation of hue for color scheme C_HUE */

/* Basic math */

#define PI 	3.141592654
#define DPI 	6.283185307
#define PID 	1.570796327

```

### 14 May 21 - What's going on with the pentagon? ###

**Program:** `drop_billiard.c`

**Initial condition in function `animation()`:** `init_drop_config(0.0, 0.0, 0.0, DPI, configs);`

```
#define MOVIE 1         /* set to 1 to generate movie */

#define WINWIDTH 	1280  /* window width */
#define WINHEIGHT 	720   /* window height */

#define XMIN -2.0
#define XMAX 2.0	/* x interval */
#define YMIN -1.125
#define YMAX 1.125	/* y interval for 9/16 aspect ratio */

/* Choice of the billiard table */

#define B_DOMAIN 8      /* choice of domain shape */

#define D_RECTANGLE 0   /* rectangular domain */
#define D_ELLIPSE 1     /* elliptical domain */
#define D_STADIUM 2     /* stadium-shaped domain */
#define D_SINAI 3       /* Sinai billiard */
#define D_DIAMOND 4     /* diamond-shaped billiard */
#define D_TRIANGLE 5    /* triangular billiard */
#define D_ANNULUS 7     /* annulus */
#define D_POLYGON 8     /* polygon */

#define LAMBDA 1.0	/* parameter controlling shape of billiard */
#define MU 0.1          /* second parameter controlling shape of billiard */
#define FOCI 1          /* set to 1 to draw focal points of ellipse */
#define NPOLY 5             /* number of sides of polygon */
#define APOLY 1.0           /* angle by which to turn polygon, in units of Pi/2 */ 

#define RESAMPLE 0      /* set to 1 if particles should be added when dispersion too large */

#define NPART 10000	/* number of particles */
#define NPARTMAX 20000	/* maximal number of particles after resampling */
#define LMAX 0.01       /* minimal segment length triggering resampling */ 
#define LPERIODIC 1.0   /* lines longer than this are not drawn (useful for Sinai billiard) */
#define LCUT 2.0        /* controls the max size of segments not considered as being cut */
#define DMIN 0.02       /* minimal distance to boundary for triggering resampling */ 
#define CYCLE 0         /* set to 1 for closed curve (start in all directions) */

#define NSTEPS 5000     /* number of frames of movie */
#define TIME 50         /* time between movie frames, for fluidity of real-time simulation */ 
#define DPHI 0.0001     /* integration step */
#define NVID 25         /* number of iterations between images displayed on screen */
#define NCOLORS 8       /* number of colors */
#define COLORSHIFT 180  /* hue of initial color */ 
#define NSEG 100        /* number of segments of boundary */

#define BLACK 1         /* set to 1 for black background */
#define PAINT_INT 1     /* set to 1 to paint interior in other color (for polygon) */

/* Decreasing TIME accelerates the animation and the movie               */
/* For constant speed of movie, TIME*DPHI should be kept constant        */
/* However, increasing DPHI too much deterioriates quality of simulation */
/* For a good quality movie, take for instance TIME = 50, DPHI = 0.0002  */

#define PAUSE 1000          /* number of frames after which to pause */
#define PSLEEP 1         /* sleep time during pause */
#define SLEEP1  1        /* initial sleeping time */
#define SLEEP2  100      /* final sleeping time */

#define PI 	3.141592654
#define DPI 	6.283185307
#define PID 	1.570796327

```

### 13 May 21 - Quantum of focus: Schrödinger's equation in an elliptical domain ###

**Program:** `schrodinger.c`

**Initial condition in function `animation()`:** `init_coherent_state(-0.5, 0.0, 1.0, 1.0, 0.03, phi, psi, xy_in);`

```
#define MOVIE 0         /* set to 1 to generate movie */

/* General geometrical parameters */

#define WINWIDTH 	1280  /* window width */
#define WINHEIGHT 	720   /* window height */

#define NX 640          /* number of grid points on x axis */
#define NY 360          /* number of grid points on y axis */

/* setting NX to WINWIDTH and NY to WINHEIGHT increases resolution */
/* but will multiply run time by 4                                 */

#define XMIN -2.0
#define XMAX 2.0	/* x interval */
#define YMIN -1.125
#define YMAX 1.125	/* y interval for 9/16 aspect ratio */

/* Choice of the billiard table */

#define B_DOMAIN 1      /* choice of domain shape */

#define D_RECTANGLE 0   /* rectangular domain */
#define D_ELLIPSE 1     /* elliptical domain */
#define D_STADIUM 2     /* stadium-shaped domain */
#define D_SINAI 3       /* Sinai billiard */
#define D_DIAMOND 4     /* diamond-shaped billiard */
#define D_TRIANGLE 5    /* triangular billiard */
#define D_FLAT 6        /* flat interface */
#define D_ANNULUS 7     /* annulus */
#define D_POLYGON 8     /* polygon */

#define LAMBDA 1.5	    /* parameter controlling the dimensions of domain */
#define NPOLY 6             /* number of sides of polygon */
#define APOLY 1.0           /* angle by which to turn polygon, in units of Pi/2 */ 
#define FOCI 1              /* set to 1 to draw focal points of ellipse */

/* You can add more billiard tables by adapting the functions */
/* xy_in_billiard and draw_billiard below */

/* Physical patameters of wave equation */

#define DT 0.00000002
#define HBAR 1.0
/* The Courant number is given by c*DT/DX, where DT is the time step and DX the lattice spacing */
/* The physical damping coefficient is given by GAMMA/(DT)^2 */
/* Increasing COURANT speeds up the simulation, but decreases accuracy */
/* For similar wave forms, COURANT^2*GAMMA should be kept constant */

/* For debugging purposes only */
#define FLOOR 0         /* set to 1 to limit wave amplitude to VMAX */
#define VMAX 10.0       /* max value of wave amplitude */

/* Parameters for length and speed of simulation */

#define NSTEPS 4000     //7500      /* number of frames of movie */
#define NVID 500          /* number of iterations between images displayed on screen */
#define NSEG 100         /* number of segments of boundary */

#define PAUSE 1000         /* number of frames after which to pause */
#define PSLEEP 1         /* sleep time during pause */
#define SLEEP1  1        /* initial sleeping time */
#define SLEEP2  1   /* final sleeping time */

/* Plot type */

#define PLOT 0

#define P_MODULE 0        /* plot module of wave function squared */
#define P_PHASE 1         /* plot phase of wave function */
#define P_REAL 2          /* plot real part */
#define P_IMAGINARY 3     /* plot imaginary part */

/* Color schemes */

#define BLACK 1          /* background */

#define COLOR_SCHEME 1   /* choice of color scheme */

#define C_LUM 0          /* color scheme modifies luminosity (with slow drift of hue) */
#define C_HUE 1          /* color scheme modifies hue */
#define C_PHASE 2        /* color scheme shows phase */

#define SCALE 1          /* set to 1 to adjust color scheme to variance of field */
#define SLOPE 1.0        /* sensitivity of color on wave amplitude */
#define ATTENUATION 0.0  /* exponential attenuation coefficient of contrast with time */

#define COLORHUE 260     /* initial hue of water color for scheme C_LUM */
#define COLORDRIFT 0.0   /* how much the color hue drifts during the whole simulation */
#define LUMMEAN 0.5      /* amplitude of luminosity variation for scheme C_LUM */
#define LUMAMP 0.3       /* amplitude of luminosity variation for scheme C_LUM */
#define HUEMEAN 120.0    /* mean value of hue for color scheme C_HUE */
#define HUEAMP -120.0      /* amplitude of variation of hue for color scheme C_HUE */

/* Basic math */

#define PI 	3.141592654
#define DPI 	6.283185307
#define PID 	1.570796327

```

### 12 May 21 - The eccentric annular billiard, or how to mix a Pan Galactic Gargle Blaster ###

**Program:** `particle_billiard.c`

**Initial condition in function `animation()`:** `init_line_config(-0.7, 0.2, -0.7, 0.7, 0.0, configs);`

```
#define MOVIE 0         /* set to 1 to generate movie */

#define WINWIDTH 	1280  /* window width */
#define WINHEIGHT 	720   /* window height */

#define XMIN -2.0
#define XMAX 2.0	/* x interval */
#define YMIN -1.125
#define YMAX 1.125	/* y interval for 9/16 aspect ratio */

/* Choice of the billiard table */

#define B_DOMAIN 7      /* choice of domain shape */

#define D_RECTANGLE 0   /* rectangular domain */
#define D_ELLIPSE 1     /* elliptical domain */
#define D_STADIUM 2     /* stadium-shaped domain */
#define D_SINAI 3       /* Sinai billiard */
#define D_DIAMOND 4     /* diamond-shaped billiard */
#define D_TRIANGLE 5    /* triangular billiard */
#define D_ANNULUS 7     /* annulus */
#define D_POLYGON 8     /* polygon */

#define LAMBDA 0.5	/* parameter controlling shape of billiard */
// #define LAMBDA 1.73205080756888	/* sqrt(3) for triangle tiling plane */
#define MU 0.1          /* second parameter controlling shape of billiard */
#define FOCI 1          /* set to 1 to draw focal points of ellipse */
#define NPOLY 8             /* number of sides of polygon */
#define APOLY 1.0           /* angle by which to turn polygon, in units of Pi/2 */ 

#define RESAMPLE 0      /* set to 1 if particles should be added when dispersion too large */
#define DEBUG 0         /* draw trajectories, for debugging purposes */

#define NPART 50      /* number of particles */
// #define NPART 50      /* number of particles */
#define NPARTMAX 50000	/* maximal number of particles after resampling */
#define LMAX 0.01       /* minimal segment length triggering resampling */ 
#define DMIN 0.02       /* minimal distance to boundary for triggering resampling */ 
#define CYCLE 1         /* set to 1 for closed curve (start in all directions) */

#define NSTEPS 7500     /* number of frames of movie */
#define TIME 750       /* time between movie frames, for fluidity of real-time simulation */ 
// #define TIME 1500       /* time between movie frames, for fluidity of real-time simulation */ 
#define DPHI 0.00001    /* integration step */
#define NVID 150        /* number of iterations between images displayed on screen */
// #define NVID 500        /* number of iterations between images displayed on screen */

#define NCOLORS 5      /* number of colors */
// #define NCOLORS 14      /* number of colors */
#define COLORSHIFT 0    /* hue of initial color */ 
#define NSEG 100        /* number of segments of boundary */
#define LENGTH 0.05     /* length of velocity vectors */

#define BLACK 1         /* set to 1 for black background */

/* Decreasing TIME accelerates the animation and the movie                               */
/* For constant speed of movie, TIME*DPHI should be kept constant                        */
/* However, increasing DPHI too much deterioriates quality of simulation                 */
/* NVID tells how often a picture is drawn in the animation, increase it for faster anim */
/* For a good quality movie, take for instance TIME = 400, DPHI = 0.00005, NVID = 100    */

#define PAUSE 1000       /* number of frames after which to pause */
#define PSLEEP 1         /* sleep time during pause */
#define SLEEP1  1        /* initial sleeping time */
#define SLEEP2  1000      /* final sleeping time */

#define PI 	3.141592654
#define DPI 	6.283185307
#define PID 	1.570796327

```

### 11 May 21 - Diffraction of a wave on a grating ###

**Program:** `wave_billiard.c`

**Initial condition in function `animation()`:** `init_wave(-1.0, 0.0, phi, psi, xy_in);`

```
#define MOVIE 1         /* set to 1 to generate movie */

/* General geometrical parameters */

#define WINWIDTH 	1280  /* window width */
#define WINHEIGHT 	720   /* window height */

#define NX 640          /* number of grid points on x axis */
#define NY 360          /* number of grid points on y axis */

/* setting NX to WINWIDTH and NY to WINHEIGHT increases resolution */
/* but will multiply run time by 4                                 */

#define XMIN -2.0
#define XMAX 2.0	/* x interval */
#define YMIN -1.125
#define YMAX 1.125	/* y interval for 9/16 aspect ratio */

/* Choice of the billiard table */

#define B_DOMAIN 10      /* choice of domain shape */

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

#define LAMBDA 0.08	    /* parameter controlling the dimensions of domain */
#define MU 0.02	            /* parameter controlling the dimensions of domain */
#define NPOLY 8             /* number of sides of polygon */
#define APOLY 1.0           /* angle by which to turn polygon, in units of Pi/2 */ 
#define FOCI 1              /* set to 1 to draw focal points of ellipse */

/* You can add more billiard tables by adapting the functions */
/* xy_in_billiard and draw_billiard below */

/* Physical patameters of wave equation */

#define OMEGA 0.7           /* frequency of periodic excitation */
#define COURANT 0.01       /* Courant number */
#define GAMMA 0.0      /* damping factor in wave equation */
#define KAPPA 0.0       /* "elasticity" term enforcing oscillations */
/* The Courant number is given by c*DT/DX, where DT is the time step and DX the lattice spacing */
/* The physical damping coefficient is given by GAMMA/(DT)^2 */
/* Increasing COURANT speeds up the simulation, but decreases accuracy */
/* For similar wave forms, COURANT^2*GAMMA should be kept constant */

/* For debugging purposes only */
#define FLOOR 0         /* set to 1 to limit wave amplitude to VMAX */
#define VMAX 10.0       /* max value of wave amplitude */

/* Parameters for length and speed of simulation */

#define NSTEPS 4500     //7500      /* number of frames of movie */
#define NVID 50          /* number of iterations between images displayed on screen */
#define NSEG 100         /* number of segments of boundary */

#define PAUSE 100         /* number of frames after which to pause */
#define PSLEEP 1         /* sleep time during pause */
#define SLEEP1  1        /* initial sleeping time */
#define SLEEP2  1   /* final sleeping time */

/* Color schemes */

#define COLOR_SCHEME 0   /* choice of color scheme */

#define C_LUM 0          /* color scheme modifies luminosity (with slow drift of hue) */
#define C_HUE 1          /* color scheme modifies hue */

#define SCALE 1          /* set to 1 to adjust color scheme to variance of field */
#define SLOPE 1.0        /* sensitivity of color on wave amplitude */
#define ATTENUATION 0.0  /* exponential attenuation coefficient of contrast with time */

#define COLORHUE 270     /* initial hue of water color for scheme C_LUM */
#define COLORDRIFT 0.0   /* how much the color hue drifts during the whole simulation */
#define LUMMEAN 0.5      /* amplitude of luminosity variation for scheme C_LUM */
#define LUMAMP 0.3       /* amplitude of luminosity variation for scheme C_LUM */
#define HUEMEAN 220.0    /* mean value of hue for color scheme C_HUE */
#define HUEAMP 60.0      /* amplitude of variation of hue for color scheme C_HUE */

/* Basic math */

#define PI 	3.141592654
#define DPI 	6.283185307
#define PID 	1.570796327

```

### 10 May 21 - A perfect wave front in the Bermuda triangle billiard ###

**Program:** `drop_billiard.c`

**Initial condition in function `animation()`:** `init_drop_config(0.0, 0.0, 0.0, DPI, configs);`

```
#define MOVIE 0         /* set to 1 to generate movie */

#define WINWIDTH 	1280  /* window width */
#define WINHEIGHT 	720   /* window height */

#define XMIN -2.0
#define XMAX 2.0	/* x interval */
#define YMIN -1.125
#define YMAX 1.125	/* y interval for 9/16 aspect ratio */

/* Choice of the billiard table */

#define B_DOMAIN 8      /* choice of domain shape */

#define D_RECTANGLE 0   /* rectangular domain */
#define D_ELLIPSE 1     /* elliptical domain */
#define D_STADIUM 2     /* stadium-shaped domain */
#define D_SINAI 3       /* Sinai billiard */
#define D_DIAMOND 4     /* diamond-shaped billiard */
#define D_TRIANGLE 5    /* triangular billiard */
#define D_ANNULUS 7     /* annulus */
#define D_POLYGON 8     /* polygon */

#define LAMBDA 1.0	/* parameter controlling shape of billiard */
#define MU 0.1          /* second parameter controlling shape of billiard */
#define FOCI 1          /* set to 1 to draw focal points of ellipse */
#define NPOLY 3             /* number of sides of polygon */
#define APOLY 1.0           /* angle by which to turn polygon, in units of Pi/2 */ 

#define RESAMPLE 0      /* set to 1 if particles should be added when dispersion too large */

#define NPART 10000	/* number of particles */
#define NPARTMAX 20000	/* maximal number of particles after resampling */
#define LMAX 0.01       /* minimal segment length triggering resampling */ 
#define LPERIODIC 0.01   /* lines longer than this are not drawn (useful for Sinai billiard) */
#define DMIN 0.02       /* minimal distance to boundary for triggering resampling */ 
#define MARGIN 0.02      /* distance above which points of curve are not drawn */
#define CYCLE 0         /* set to 1 for closed curve (start in all directions) */

#define NSTEPS 5000     /* number of frames of movie */
#define TIME 150        /* time between movie frames, for fluidity of real-time simulation */ 
#define DPHI 0.0001     /* integration step */
#define NVID 25        /* number of iterations between images displayed on screen */
#define NCOLORS 10      /* number of colors */
#define COLORSHIFT 60   /* hue of initial color */ 
#define NSEG 100        /* number of segments of boundary */

#define BLACK 1         /* set to 1 for black background */

/* Decreasing TIME accelerates the animation and the movie               */
/* For constant speed of movie, TIME*DPHI should be kept constant        */
/* However, increasing DPHI too much deterioriates quality of simulation */
/* For a good quality movie, take for instance TIME = 50, DPHI = 0.0002  */

#define PAUSE 1000          /* number of frames after which to pause */
#define PSLEEP 1         /* sleep time during pause */
#define SLEEP1  1        /* initial sleeping time */
#define SLEEP2  100      /* final sleeping time */

#define PI 	3.141592654
#define DPI 	6.283185307
#define PID 	1.570796327

```

### 9 May 21 - Waves around a conical singularity ###

**Program:** Special code `wave_conical.c`, not yet incorporated 

Main changes: 

```
#define D_CONICAL 999     /* L-shape with conical singularity */

#define LAMBDA 0.3	    /* parameter controlling the dimensions of domain */
#define NPOLY 6             /* number of sides of polygon */
#define APOLY 1.0           /* angle by which to turn polygon, in units of Pi/2 */ 
#define FOCI 1              /* set to 1 to draw focal points of ellipse */

/* You can add more billiard tables by adapting the functions */
/* xy_in_billiard and draw_billiard below */

/* Physical patameters of wave equation */

#define COURANT 0.01       /* Courant number */
#define GAMMA 0.0      /* damping factor in wave equation */
#define KAPPA 5.0e-7       /* "elasticity" term enforcing oscillations */
/* The Courant number is given by c*DT/DX, where DT is the time step and DX the lattice spacing */
/* The physical damping coefficient is given by GAMMA/(DT)^2 */
/* Increasing COURANT speeds up the simulation, but decreases accuracy */
/* For similar wave forms, COURANT^2*GAMMA should be kept constant */

/* For debugging purposes only */
#define FLOOR 0         /* set to 1 to limit wave amplitude to VMAX */
#define VMAX 10.0       /* max value of wave amplitude */

/* Parameters for length and speed of simulation */

#define NSTEPS 7500     //7500      /* number of frames of movie */
#define NVID 50          /* number of iterations between images displayed on screen */
#define NSEG 100         /* number of segments of boundary */

#define PAUSE 1000         /* number of frames after which to pause */
#define PSLEEP 1         /* sleep time during pause */
#define SLEEP1  1        /* initial sleeping time */
#define SLEEP2  1   /* final sleeping time */

/* Color schemes */

#define COLOR_SCHEME 1   /* choice of color scheme */

#define C_LUM 0          /* color scheme modifies luminosity (with slow drift of hue) */
#define C_HUE 1          /* color scheme modifies hue */

#define SCALE 1          /* set to 1 to adjust color scheme to variance of field */
#define SLOPE 1.0        /* sensitivity of color on wave amplitude */
#define ATTENUATION 0.0  /* exponential attenuation coefficient of contrast with time */

#define COLORHUE 260     /* initial hue of water color for scheme C_LUM */
#define COLORDRIFT 0.0   /* how much the color hue drifts during the whole simulation */
#define LUMMEAN 0.5      /* amplitude of luminosity variation for scheme C_LUM */
#define LUMAMP 0.3       /* amplitude of luminosity variation for scheme C_LUM */
#define HUEMEAN 120.0    /* mean value of hue for color scheme C_HUE */
#define HUEAMP 120.0      /* amplitude of variation of hue for color scheme C_HUE */

/* Basic math */

#define PI 	3.141592654
#define DPI 	6.283185307
#define PID 	1.570796327

int xy_in_billiard(x, y)
/* returns 1 if (x,y) represents a point in the billiard */
double x, y;
{
    double l2, r2, omega, c, angle;
    int k, condition;

    switch (B_DOMAIN) {
         case D_CONICAL:
        {
            if ((vabs(x) < 1.0)&&(y > -1.0)&&(y < 0.0)) return(1);
            else if ((vabs(y) < 1.0)&&(x > -1.0)&&(x < 0.0)) return(1);
            else return(0);
        }
}

void draw_segment(x1, y1, x2, y2, h, s, l)
/* draw line segment (x1,y1)-(x2,y2) in color (h,s,l) */
double x1, y1, x2, y2, h, s, l;
{
    double rgb[3], pos[2];
    
    glBegin(GL_LINE_STRIP);
    hsl_to_rgb(h, s, l, rgb);
    glColor3f(rgb[0], rgb[1], rgb[2]);
    xy_to_pos(x1, y1, pos);
    glVertex2d(pos[0], pos[1]);
    xy_to_pos(x2, y2, pos);
    glVertex2d(pos[0], pos[1]);
    glEnd();
}

void draw_billiard()      /* draws the billiard boundary */
{
    double x0, x, y, phi, r = 0.01, pos[2], pos1[2], alpha, dphi, omega, rgb[3];
    int i;

    glColor3f(0.0, 0.0, 0.0);
    glLineWidth(5);

    glEnable(GL_LINE_SMOOTH);

    switch (B_DOMAIN) {
        case (D_CONICAL):
        {
            draw_segment(-1.0, -1.0, 0.0, -1.0, 0.0, 1.0, 0.5);
            draw_segment(-1.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.5);
            draw_segment(-1.0, -1.0, -1.0, 0.0, 220.0, 1.0, 0.5);
            draw_segment(1.0, -1.0, 1.0, 0.0, 220.0, 1.0, 0.5);
            draw_segment(0.0, -1.0, 1.0, -1.0, 60.0, 1.0, 0.5);
            draw_segment(0.0, 0.0, 1.0, 0.0, 60.0, 1.0, 0.5);
            draw_segment(-1.0, 0.0, -1.0, 1.0, 180.0, 1.0, 0.5);
            draw_segment(0.0, 0.0, 0.0, 1.0, 180.0, 1.0, 0.5);
            break;
        }
    }
}

void evolve_wave(phi, psi, xy_in)
/* time step of field evolution */
/* phi is value of field at time t, psi at time t-1 */
    double *phi[NX], *psi[NX]; short int *xy_in[NX];
{
    int i, j, iplus, iminus, jplus, jminus, i1, i2, i3, j1, j2, j3, ij[2];
    double delta, x, y;

    /* boundaries */
    xy_to_ij(-1.0, -1.0, ij);
    i1 = ij[0];     j1 = ij[1];
    xy_to_ij(0.0, 0.0, ij);
    i2 = ij[0];     j2 = ij[1];
    xy_to_ij(1.0, 1.0, ij);
    i3 = ij[0];     j3 = ij[1];
    
    #pragma omp parallel for private(i,j,iplus,iminus,delta,x,y)
    for (i=0; i<NX; i++){
        for (j=0; j<NY; j++){
            if (xy_in[i][j]){
                /* boundary conditions */
                iplus = i+1;
                if ((iplus == i2)&&(j > j2)) iplus = i1+1;
                else if ((iplus == i3)&&(j <= j2)) iplus = i1+1;
                
                iminus = i-1;
                if ((iminus == i1)&&(j > j2)) iminus = i2-1;
                else if ((iminus == i1)&&(j <= j2)) iminus = i3-1;
                
                jplus = j+1;
                if ((jplus == j2)&&(i > i2)) jplus = j1+1;
                else if ((jplus == j3)&&(i <= i2)) jplus = j1+1;
                
                jminus = j-1;
                if ((jminus == j1)&&(i > i2)) jminus = j2-1;
                else if ((jminus == j1)&&(i <= i2)) jminus = j3-1;

                /* discretized Laplacian */
                delta = phi[iplus][j] + phi[iminus][j] + phi[i][jplus] + phi[i][jminus] - 4.0*phi[i][j];

                x = phi[i][j];
		y = psi[i][j];

                /* evolve phi */
                phi[i][j] = -y + 2*x + courant2*delta - KAPPA*x - GAMMA*(x-y);

                psi[i][j] = x;

                if (FLOOR)
                {
                    if (phi[i][j] > VMAX) phi[i][j] = VMAX;
                    if (phi[i][j] < -VMAX) phi[i][j] = -VMAX;
                    if (psi[i][j] > VMAX) psi[i][j] = VMAX;
                    if (psi[i][j] < -VMAX) psi[i][j] = -VMAX;
                }
            }
        }
    }
}
```

### 8 May 21 - If life gives you a lemon billiard, use it to make lemonade! ###

**Program:** `particle_billiard.c`

**Initial condition in function `animation()`:** `init_line_config(0.0, -0.7, 0.0, 0.7, 0.0, configs);`

```
#define MOVIE 1         /* set to 1 to generate movie */

#define WINWIDTH 	1280  /* window width */
#define WINHEIGHT 	720   /* window height */

#define XMIN -2.0
#define XMAX 2.0	/* x interval */
#define YMIN -1.125
#define YMAX 1.125	/* y interval for 9/16 aspect ratio */

/* Choice of the billiard table */

#define B_DOMAIN 2      /* choice of domain shape */

#define D_RECTANGLE 0   /* rectangular domain */
#define D_ELLIPSE 1     /* elliptical domain */
#define D_STADIUM 2     /* stadium-shaped domain */
#define D_SINAI 3       /* Sinai billiard */
#define D_DIAMOND 4     /* diamond-shaped billiard */
#define D_TRIANGLE 5    /* triangular billiard */
#define D_ANNULUS 7     /* annulus */
#define D_POLYGON 8     /* polygon */

#define LAMBDA -0.8	/* parameter controlling shape of billiard */
#define MU 0.1          /* second parameter controlling shape of billiard */
#define FOCI 1          /* set to 1 to draw focal points of ellipse */
#define NPOLY 8             /* number of sides of polygon */
#define APOLY 1.0           /* angle by which to turn polygon, in units of Pi/2 */ 

#define RESAMPLE 0      /* set to 1 if particles should be added when dispersion too large */
#define DEBUG 0         /* draw trajectories, for debugging purposes */

#define NPART 10000     /* number of particles */
#define NPARTMAX 50000	/* maximal number of particles after resampling */
#define LMAX 0.01       /* minimal segment length triggering resampling */ 
#define DMIN 0.02       /* minimal distance to boundary for triggering resampling */ 
#define CYCLE 1         /* set to 1 for closed curve (start in all directions) */

#define NSTEPS 5000     /* number of frames of movie */
#define TIME 750       /* time between movie frames, for fluidity of real-time simulation */ 
#define DPHI 0.00001    /* integration step */
#define NVID 750        /* number of iterations between images displayed on screen */

#define NCOLORS 10      /* number of colors */
#define COLORSHIFT 0    /* hue of initial color */ 
#define NSEG 100        /* number of segments of boundary */
#define LENGTH 0.05     /* length of velocity vectors */

#define BLACK 1         /* set to 1 for black background */

/* Decreasing TIME accelerates the animation and the movie                               */
/* For constant speed of movie, TIME*DPHI should be kept constant                        */
/* However, increasing DPHI too much deterioriates quality of simulation                 */
/* NVID tells how often a picture is drawn in the animation, increase it for faster anim */
/* For a good quality movie, take for instance TIME = 400, DPHI = 0.00005, NVID = 100    */

#define PAUSE 1000       /* number of frames after which to pause */
#define PSLEEP 1         /* sleep time during pause */
#define SLEEP1  1        /* initial sleeping time */
#define SLEEP2  1000      /* final sleeping time */

#define PI 	3.141592654
#define DPI 	6.283185307
#define PID 	1.570796327

```

### 7 May 21 - Diffraction - Young's double-slit experiment with waves ###

**Program:** `wave_billiard.c`

**Initial condition in function `animation()`:** `init_wave(-1.0, 0.0, phi, psi, xy_in);`

```
#define MOVIE 1         /* set to 1 to generate movie */

/* General geometrical parameters */

#define WINWIDTH 	1280  /* window width */
#define WINHEIGHT 	720   /* window height */

#define NX 640          /* number of grid points on x axis */
#define NY 360          /* number of grid points on y axis */

/* setting NX to WINWIDTH and NY to WINHEIGHT increases resolution */
/* but will multiply run time by 4                                 */

#define XMIN -2.0
#define XMAX 2.0	/* x interval */
#define YMIN -1.125
#define YMAX 1.125	/* y interval for 9/16 aspect ratio */

/* Choice of the billiard table */

#define B_DOMAIN 9      /* choice of domain shape */

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

#define LAMBDA 0.5	    /* parameter controlling the dimensions of domain */
#define MU 0.04	            /* parameter controlling the dimensions of domain */
#define NPOLY 8             /* number of sides of polygon */
#define APOLY 1.0           /* angle by which to turn polygon, in units of Pi/2 */ 
#define FOCI 1              /* set to 1 to draw focal points of ellipse */

/* You can add more billiard tables by adapting the functions */
/* xy_in_billiard and draw_billiard below */

/* Physical patameters of wave equation */

#define OMEGA 0.7           /* frequency of periodic excitation */
#define COURANT 0.01       /* Courant number */
#define GAMMA 0.0      /* damping factor in wave equation */
#define KAPPA 0.0       /* "elasticity" term enforcing oscillations */
/* The Courant number is given by c*DT/DX, where DT is the time step and DX the lattice spacing */
/* The physical damping coefficient is given by GAMMA/(DT)^2 */
/* Increasing COURANT speeds up the simulation, but decreases accuracy */
/* For similar wave forms, COURANT^2*GAMMA should be kept constant */

/* For debugging purposes only */
#define FLOOR 0         /* set to 1 to limit wave amplitude to VMAX */
#define VMAX 10.0       /* max value of wave amplitude */

/* Parameters for length and speed of simulation */

#define NSTEPS 4000     //7500      /* number of frames of movie */
#define NVID 50          /* number of iterations between images displayed on screen */
#define NSEG 100         /* number of segments of boundary */

#define PAUSE 100         /* number of frames after which to pause */
#define PSLEEP 1         /* sleep time during pause */
#define SLEEP1  1        /* initial sleeping time */
#define SLEEP2  1   /* final sleeping time */

/* Color schemes */

#define COLOR_SCHEME 1   /* choice of color scheme */

#define C_LUM 0          /* color scheme modifies luminosity (with slow drift of hue) */
#define C_HUE 1          /* color scheme modifies hue */

#define SCALE 1          /* set to 1 to adjust color scheme to variance of field */
#define SLOPE 1.0        /* sensitivity of color on wave amplitude */
#define ATTENUATION 0.0  /* exponential attenuation coefficient of contrast with time */

#define COLORHUE 260     /* initial hue of water color for scheme C_LUM */
#define COLORDRIFT 0.0   /* how much the color hue drifts during the whole simulation */
#define LUMMEAN 0.5      /* amplitude of luminosity variation for scheme C_LUM */
#define LUMAMP 0.3       /* amplitude of luminosity variation for scheme C_LUM */
#define HUEMEAN 220.0    /* mean value of hue for color scheme C_HUE */
#define HUEAMP 60.0      /* amplitude of variation of hue for color scheme C_HUE */

/* Basic math */

#define PI 	3.141592654
#define DPI 	6.283185307
#define PID 	1.570796327

```

### 6 May 21 - Wave fronts in a square: Dr Jekyll shouting from Hyde Park Corner ###

**Program:** `drop_billiard.c`

**Initial condition in function `animation()`:** `init_drop_config(1.0, 1.0, PI, 3.0*DPI, configs);`

```
#define MOVIE 0         /* set to 1 to generate movie */

#define WINWIDTH 	1280  /* window width */
#define WINHEIGHT 	720   /* window height */

#define XMIN -2.0
#define XMAX 2.0	/* x interval */
#define YMIN -1.125
#define YMAX 1.125	/* y interval for 9/16 aspect ratio */

/* Choice of the billiard table */

#define B_DOMAIN 0      /* choice of domain shape */

#define D_RECTANGLE 0   /* rectangular domain */
#define D_ELLIPSE 1     /* elliptical domain */
#define D_STADIUM 2     /* stadium-shaped domain */
#define D_SINAI 3       /* Sinai billiard */
#define D_DIAMOND 4     /* diamond-shaped billiard */
#define D_TRIANGLE 5    /* triangular billiard */
#define D_ANNULUS 7     /* annulus */
#define D_POLYGON 8     /* polygon */

#define LAMBDA 1.0	/* parameter controlling shape of billiard */
#define MU 0.1          /* second parameter controlling shape of billiard */
#define FOCI 1          /* set to 1 to draw focal points of ellipse */
#define NPOLY 8             /* number of sides of polygon */
#define APOLY 1.0           /* angle by which to turn polygon, in units of Pi/2 */ 

#define RESAMPLE 0      /* set to 1 if particles should be added when dispersion too large */

#define NPART 20000	/* number of particles */
#define NPARTMAX 20000	/* maximal number of particles after resampling */
#define LMAX 0.01       /* minimal segment length triggering resampling */ 
#define LPERIODIC 2.0   /* lines longer than this are not drawn (useful for Sinai billiard) */
#define DMIN 0.02       /* minimal distance to boundary for triggering resampling */ 
#define MARGIN 1.0      /* distance above which points of curve are not drawn */
#define CYCLE 0         /* set to 1 for closed curve (start in all directions) */

#define NSTEPS 10000    /* number of frames of movie */
#define TIME 100        /* time between movie frames, for fluidity of real-time simulation */ 
#define DPHI 0.0001     /* integration step */
#define NVID 100        /* number of iterations between images displayed on screen */
#define NCOLORS 20      /* number of colors */
#define COLORSHIFT 90   /* hue of initial color */ 
#define NSEG 100        /* number of segments of boundary */

#define BLACK 1         /* set to 1 for black background */

/* Decreasing TIME accelerates the animation and the movie               */
/* For constant speed of movie, TIME*DPHI should be kept constant        */
/* However, increasing DPHI too much deterioriates quality of simulation */
/* For a good quality movie, take for instance TIME = 50, DPHI = 0.0002  */

#define PAUSE 1000          /* number of frames after which to pause */
#define PSLEEP 1         /* sleep time during pause */
#define SLEEP1  1        /* initial sleeping time */
#define SLEEP2  100      /* final sleeping time */

#define PI 	3.141592654
#define DPI 	6.283185307
#define PID 	1.570796327

```

### 5 May 21 - Winter is coming ###

**Program:** `wave_billiard.c`

**Initial condition in function `animation()`:** `init_wave(0.0, 0.0, phi, psi, xy_in);`

```
#define MOVIE 1         /* set to 1 to generate movie */

/* General geometrical parameters */

#define WINWIDTH 	1280  /* window width */
#define WINHEIGHT 	720   /* window height */

#define NX 640          /* number of grid points on x axis */
#define NY 360          /* number of grid points on y axis */

/* setting NX to WINWIDTH and NY to WINHEIGHT increases resolution */
/* but will multiply run time by 4                                 */

#define XMIN -2.0
#define XMAX 2.0	/* x interval */
#define YMIN -1.125
#define YMAX 1.125	/* y interval for 9/16 aspect ratio */

/* Choice of the billiard table */

#define B_DOMAIN 8      /* choice of domain shape */

#define D_RECTANGLE 0   /* rectangular domain */
#define D_ELLIPSE 1     /* elliptical domain */
#define D_STADIUM 2     /* stadium-shaped domain */
#define D_SINAI 3       /* Sinai billiard */
#define D_DIAMOND 4     /* diamond-shaped billiard */
#define D_TRIANGLE 5    /* triangular billiard */
#define D_FLAT 6        /* flat interface */
#define D_ANNULUS 7     /* annulus */
#define D_POLYGON 8     /* polygon */

#define LAMBDA 1.0	    /* parameter controlling the dimensions of domain */
#define NPOLY 6             /* number of sides of polygon */
#define FOCI 1              /* set to 1 to draw focal points of ellipse */

/* You can add more billiard tables by adapting the functions */
/* xy_in_billiard and draw_billiard below */

/* Physical patameters of wave equation */

#define COURANT 0.01       /* Courant number */
#define GAMMA 0.0      /* damping factor in wave equation */
#define KAPPA 5.0e-7       /* "elasticity" term enforcing oscillations */
/* The Courant number is given by c*DT/DX, where DT is the time step and DX the lattice spacing */
/* The physical damping coefficient is given by GAMMA/(DT)^2 */
/* Increasing COURANT speeds up the simulation, but decreases accuracy */
/* For similar wave forms, COURANT^2*GAMMA should be kept constant */

/* For debugging purposes only */
#define FLOOR 0         /* set to 1 to limit wave amplitude to VMAX */
#define VMAX 10.0       /* max value of wave amplitude */

/* Parameters for length and speed of simulation */

#define NSTEPS 4575     //7500      /* number of frames of movie */
#define NVID 50          /* number of iterations between images displayed on screen */
#define NSEG 100         /* number of segments of boundary */

#define PAUSE 1000         /* number of frames after which to pause */
#define PSLEEP 1         /* sleep time during pause */
#define SLEEP1  1        /* initial sleeping time */
#define SLEEP2  1   /* final sleeping time */

/* Color schemes */

#define COLOR_SCHEME 1   /* choice of color scheme */

#define C_LUM 0          /* color scheme modifies luminosity (with slow drift of hue) */
#define C_HUE 1          /* color scheme modifies hue */

#define SCALE 1          /* set to 1 to adjust color scheme to variance of field */
#define SLOPE 1.0        /* sensitivity of color on wave amplitude */
#define ATTENUATION 0.0  /* exponential attenuation coefficient of contrast with time */

#define COLORHUE 260     /* initial hue of water color for scheme C_LUM */
#define COLORDRIFT 0.0   /* how much the color hue drifts during the whole simulation */
#define LUMMEAN 0.5      /* amplitude of luminosity variation for scheme C_LUM */
#define LUMAMP 0.3       /* amplitude of luminosity variation for scheme C_LUM */
#define HUEMEAN 235.0    /* mean value of hue for color scheme C_HUE */
#define HUEAMP 60.0      /* amplitude of variation of hue for color scheme C_HUE */

/* Basic math */

#define PI 	3.141592654
#define DPI 	6.283185307
#define PID 	1.570796327

```

### 4 May 21 - Why the Bunimovich stadium is chaotic ###

**Program:** `particle_billiard.c` (adapted to August 1 version)

**Initial condition in function `animation()`:** `init_line_config(0.0, -0.05, 0.0, 0.05, 0.0, configs);`

```
#define MOVIE 0         /* set to 1 to generate movie */

#define WINWIDTH 	1280  /* window width */
#define WINHEIGHT 	720   /* window height */

#define XMIN -2.0
#define XMAX 2.0	/* x interval */
#define YMIN -1.125
#define YMAX 1.125	/* y interval for 9/16 aspect ratio */

#define SCALING_FACTOR 1.0       /* scaling factor of drawing, needed for flower billiards, otherwise set to 1.0 */

/* Choice of the billiard table, see global_particles.c */

#define B_DOMAIN 2      /* choice of domain shape */

#define CIRCLE_PATTERN 2    /* pattern of circles */

#define NMAXCIRCLES 1000        /* total number of circles (must be at least NCX*NCY for square grid) */
#define NCX 15            /* number of circles in x direction */
#define NCY 20            /* number of circles in y direction */

#define LAMBDA 0.5	/* parameter controlling shape of domain */
#define MU 0.035          /* second parameter controlling shape of billiard */
#define FOCI 1          /* set to 1 to draw focal points of ellipse */
#define NPOLY 9             /* number of sides of polygon */
#define APOLY -1.0           /* angle by which to turn polygon, in units of Pi/2 */ 
#define DRAW_BILLIARD 1     /* set to 1 to draw billiard */
#define DRAW_CONSTRUCTION_LINES 0   /* set to 1 to draw additional construction lines for billiard */
#define PERIODIC_BC 0       /* set to 1 to enforce periodic boundary conditions when drawing particles */

#define RESAMPLE 0      /* set to 1 if particles should be added when dispersion too large */
#define DEBUG 0         /* draw trajectories, for debugging purposes */

/* Simulation parameters */

#define NPART 5000       /* number of particles */
#define NPARTMAX 100000	/* maximal number of particles after resampling */
#define LMAX 0.01       /* minimal segment length triggering resampling */ 
#define DMIN 0.02       /* minimal distance to boundary for triggering resampling */ 
#define CYCLE 1         /* set to 1 for closed curve (start in all directions) */

#define NSTEPS 3000     /* number of frames of movie */
#define TIME 750       /* time between movie frames, for fluidity of real-time simulation */ 
#define DPHI 0.00002    /* integration step */
#define NVID 150         /* number of iterations between images displayed on screen */

/* Decreasing TIME accelerates the animation and the movie                               */
/* For constant speed of movie, TIME*DPHI should be kept constant                        */
/* However, increasing DPHI too much deterioriates quality of simulation                 */
/* NVID tells how often a picture is drawn in the animation, increase it for faster anim */
/* For a good quality movie, take for instance TIME = 400, DPHI = 0.00005, NVID = 100    */

/* Colors and other graphical parameters */

#define NCOLORS 10      /* number of colors */
#define COLORSHIFT 0     /* hue of initial color */ 
#define RAINBOW_COLOR 0  /* set to 1 to use different colors for all particles */
#define FLOWER_COLOR 0   /* set to 1 to adapt initial colors to flower billiard (tracks vs core) */
#define NSEG 100         /* number of segments of boundary */
#define LENGTH 0.05      /* length of velocity vectors */
#define BILLIARD_WIDTH 3    /* width of billiard */
#define PARTICLE_WIDTH 2    /* width of particles */
#define FRONT_WIDTH 3       /* width of wave front */

#define BLACK 1             /* set to 1 for black background */
#define COLOR_OUTSIDE 0     /* set to 1 for colored outside */ 
#define OUTER_COLOR 270.0   /* color outside billiard */
#define PAINT_INT 0         /* set to 1 to paint interior in other color (for polygon/Reuleaux) */


#define PAUSE 1000       /* number of frames after which to pause */
#define PSLEEP 1         /* sleep time during pause */
#define SLEEP1  1        /* initial sleeping time */
#define SLEEP2  1000      /* final sleeping time */

```

### 3 May 21 - The Revelation of Sonmi 451 ###

**Program:** `wave_billiard.c`

**Initial condition in function `animation()`:** `init_wave(0.0, 0.6, phi, psi, xy_in);`

```
#define MOVIE 0         /* set to 1 to generate movie */

/* General geometrical parameters */

#define WINWIDTH 	1280  /* window width */
#define WINHEIGHT 	720   /* window height */

#define NX 640          /* number of grid points on x axis */
#define NY 360          /* number of grid points on y axis */

/* setting NX to WINWIDTH and NY to WINHEIGHT increases resolution */
/* but will multiply run time by 4                                 */

#define XMIN -2.0
#define XMAX 2.0	/* x interval */
#define YMIN -1.125
#define YMAX 1.125	/* y interval for 9/16 aspect ratio */

/* Choice of the billiard table */

#define B_DOMAIN 7      /* choice of domain shape */

#define D_RECTANGLE 0   /* rectangular domain */
#define D_ELLIPSE 1     /* elliptical domain */
#define D_STADIUM 2     /* stadium-shaped domain */
#define D_SINAI 3       /* Sinai billiard */
#define D_DIAMOND 4     /* diamond-shaped billiard */
#define D_TRIANGLE 5    /* triangular billiard */
#define D_FLAT 6        /* flat interface */
#define D_ANNULUS 7     /* annulus */
#define D_POLYGON 8     /* polygon */

#define LAMBDA 0.5	    /* parameter controlling the dimensions of domain */
#define NPOLY 6             /* number of sides of polygon */
#define FOCI 1              /* set to 1 to draw focal points of ellipse */

/* You can add more billiard tables by adapting the functions */
/* xy_in_billiard and draw_billiard below */

/* Physical patameters of wave equation */

#define COURANT 0.01       /* Courant number */
#define GAMMA 0.0      /* damping factor in wave equation */
#define KAPPA 1.0e-8       /* "elasticity" term enforcing oscillations */
/* The Courant number is given by c*DT/DX, where DT is the time step and DX the lattice spacing */
/* The physical damping coefficient is given by GAMMA/(DT)^2 */
/* Increasing COURANT speeds up the simulation, but decreases accuracy */
/* For similar wave forms, COURANT^2*GAMMA should be kept constant */

/* For debugging purposes only */
#define FLOOR 0         /* set to 1 to limit wave amplitude to VMAX */
#define VMAX 10.0       /* max value of wave amplitude */

/* Parameters for length and speed of simulation */

#define NSTEPS 7900     //7500      /* number of frames of movie */
#define NVID 75          /* number of iterations between images displayed on screen */
#define NSEG 100         /* number of segments of boundary */

#define PAUSE 1000         /* number of frames after which to pause */
#define PSLEEP 1         /* sleep time during pause */
#define SLEEP1  1        /* initial sleeping time */
#define SLEEP2  1   /* final sleeping time */

/* Color schemes */

#define COLOR_SCHEME 1   /* choice of color scheme */

#define C_LUM 0          /* color scheme modifies luminosity (with slow drift of hue) */
#define C_HUE 1          /* color scheme modifies hue */

#define SCALE 1          /* set to 1 to adjust color scheme to variance of field */
#define SLOPE 1.0        /* sensitivity of color on wave amplitude */
#define ATTENUATION 0.0  /* exponential attenuation coefficient of contrast with time */

#define COLORHUE 260     /* initial hue of water color for scheme C_LUM */
#define COLORDRIFT 0.0   /* how much the color hue drifts during the whole simulation */
#define LUMMEAN 0.5      /* amplitude of luminosity variation for scheme C_LUM */
#define LUMAMP 0.3       /* amplitude of luminosity variation for scheme C_LUM */
#define HUEMEAN 310.0    /* mean value of hue for color scheme C_HUE */
#define HUEAMP 110.0      /* amplitude of variation of hue for color scheme C_HUE */

/* Basic math */

#define PI 	3.141592654
#define DPI 	6.283185307
#define PID 	1.570796327

```

### 2 May 21 - Triangulating the circle ###

**Program:** `drop_billiard.c` (old version)

**Initial condition in function `animation()`:** `init_drop_config(0.0, -0.5, 0.0, DPI, configs);`

```
#define WINWIDTH 	1280  /* window width */
#define WINHEIGHT 	720   /* window height */

/* Choice of the billiard table */

#define B_DOMAIN 5      /* choice of domain shape */

#define D_RECTANGLE 0   /* rectangular domain */
#define D_ELLIPSE 1     /* elliptical domain */
#define D_STADIUM 2     /* stadium-shaped domain */
#define D_SINAI 3       /* Sinai billiard */
#define D_DIAMOND 4     /* diamond-shaped billiard */
#define D_TRIANGLE 5    /* triangular billiard */


#define LAMBDA 1.73205080756888	/* parameter controlling shape of billiard */
#define FOCI 1          /* set to 1 to draw focal points of ellipse */

#define MOVIE 1         /* set to 1 to generate movie */
#define RESAMPLE 0      /* set to 1 if particles should be added when dispersion too large */

#define XMIN -2.0
#define XMAX 2.0	/* x interval */
#define YMIN -1.125
#define YMAX 1.125	/* y interval for 9/16 aspect ratio */

#define NPART 10000	/* number of particles */
#define NPARTMAX 20000	/* maximal number of particles after resampling */
#define LMAX 0.01       /* minimal segment length triggering resampling */ 
#define DMIN 0.02       /* minimal distance to boundary for triggering resampling */ 
#define MARGIN 1.0     /* distance above which points of curve are not drawn */
#define CYCLE 0         /* set to 1 for closed curve (start in all directions) */

#define NSTEPS 5000     /* number of frames of movie */
#define TIME 75         /* time between movie frames, for fluidity of real-time simulation */ 
#define DPHI 0.0001     /* integration step */
#define NVID 50         /* number of iterations between images displayed on screen */
                        /* TIME shoud be a multiple of NVID */
#define NCOLORS 10      /* number of colors */
#define COLORSHIFT 200  /* hue of initial color */ 
#define NSEG 100        /* number of segments of boundary */

#define BLACK 1         /* set to 1 for black background */


/* Decreasing TIME accelerates the animation and the movie               */
/* For constant speed of movie, TIME*DPHI should be kept constant        */
/* However, increasing DPHI too much deterioriates quality of simulation */
/* For a good quality movie, take for instance TIME = 50, DPHI = 0.0002  */

#define PAUSE 1000          /* number of frames after which to pause */
#define PSLEEP 1         /* sleep time during pause */
#define SLEEP1  1        /* initial sleeping time */
#define SLEEP2  100      /* final sleeping time */

#define PI 	3.141592654
#define DPI 	6.283185307
#define PID 	1.570796327

```

### 1 May 21 - Snell's law of refraction ###

**Program:** `wave_billiard.c` (old version, should be remade by using option `TWOSPEED 1`)

**Initial condition in function `animation()`:** `init_wave(-1.0, 0.5, phi, psi, xy_in);`

```
#define MOVIE 0         /* set to 1 to generate movie */

/* General geometrical parameters */

#define WINWIDTH 	1280  /* window width */
#define WINHEIGHT 	720   /* window height */

#define NX 640          /* number of grid points on x axis */
#define NY 360          /* number of grid points on y axis */

/* setting NX to WINWIDTH and NY to WINHEIGHT increases resolution */
/* but will multiply run time by 4                                 */

#define XMIN -2.0
#define XMAX 2.0	/* x interval */
#define YMIN -1.125
#define YMAX 1.125	/* y interval for 9/16 aspect ratio */

/* Choice of the billiard table */

#define B_DOMAIN 6      /* choice of domain shape */

#define D_RECTANGLE 0   /* rectangular domain */
#define D_ELLIPSE 1     /* elliptical domain */
#define D_STADIUM 2     /* stadium-shaped domain */
#define D_SINAI 3       /* Sinai billiard */
#define D_DIAMOND 4     /* diamond-shaped billiard */
#define D_TRIANGLE 5    /* triangular billiard */
#define D_FLAT 6        /* flat interface */

#define LAMBDA 1.0	    /* parameter controlling the dimensions of domain */
#define FOCI 1              /* set to 1 to draw focal points of ellipse */

/* You can add more billiard tables by adapting the functions */
/* xy_in_billiard and draw_billiard below */

/* Physical patameters of wave equation */

#define COURANT 0.01       /* Courant number */
#define COURANT1 0.005       /* Courant number in second domain */
#define GAMMA 0.0      /* damping factor in wave equation */
#define KAPPA 5.0e-9       /* "elasticity" term enforcing oscillations */
/* The Courant number is given by c*DT/DX, where DT is the time step and DX the lattice spacing */
/* The physical damping coefficient is given by GAMMA/(DT)^2 */
/* Increasing COURANT speeds up the simulation, but decreases accuracy */
/* For similar wave forms, COURANT^2*GAMMA should be kept constant */

/* For debugging purposes only */
#define FLOOR 0         /* set to 1 to limit wave amplitude to VMAX */
#define VMAX 10.0       /* max value of wave amplitude */

/* Parameters for length and speed of simulation */

#define NSTEPS 8000     //7500      /* number of frames of movie */
#define NVID 75          /* number of iterations between images displayed on screen */
#define NSEG 100         /* number of segments of boundary */

#define PAUSE 100         /* number of frames after which to pause */
#define PSLEEP 1         /* sleep time during pause */
#define SLEEP1  1        /* initial sleeping time */
#define SLEEP2  1   /* final sleeping time */

/* Color schemes */

#define COLOR_SCHEME 0   /* choice of color scheme */

#define C_LUM 0          /* color scheme modifies luminosity (with slow drift of hue) */
#define C_HUE 1          /* color scheme modifies hue */

#define SCALE 1          /* set to 1 to adjust color scheme to variance of field */
#define SLOPE 1.0        /* sensitivity of color on wave amplitude */
#define ATTENUATION 0.0  /* exponential attenuation coefficient of contrast with time */

#define COLORHUE 260     /* initial hue of water color for scheme C_LUM */
#define COLORDRIFT 0.0   /* how much the color hue drifts during the whole simulation */
#define LUMMEAN 0.5      /* amplitude of luminosity variation for scheme C_LUM */
#define LUMAMP 0.3       /* amplitude of luminosity variation for scheme C_LUM */
#define HUEMEAN 330.0    /* mean value of hue for color scheme C_HUE */
#define HUEAMP 90.0      /* amplitude of variation of hue for color scheme C_HUE */

/* Basic math */

#define PI 	3.141592654
#define DPI 	6.283185307
#define PID 	1.570796327

```

