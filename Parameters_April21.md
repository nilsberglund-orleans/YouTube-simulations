### Parameter values for YouTube simulations ###

Created by **Nils Berglund** and optimized by **Marco Mancini**

C code for videos on YouTube Channel https://www.youtube.com/c/NilsBerglund

Below are parameter values used for different simulations, as well as initial conditions used in 
function animation. Some simulations use variants of the published code. The list is going to be 
updated gradually. Some constants have been moved to global files in later versions.
Some parameters in early versions may be off, since not all versions were backed up. 


### 30 April 21 - Why the Sinai billiard is chaotic ###

**Program:** `vid_part_billiard.c` (old version of `part_billiard.c`)

**Initial condition in function `animation()`:** `init_drop_config(-1.5, 0.0, 0.0, 0.2, configs);`

```
#define WINWIDTH 	1280  /* window width */
#define WINHEIGHT 	720   /* window height */

/* Choice of the billiard table */

#define B_DOMAIN 3      /* choice of domain shape */

#define D_RECTANGLE 0   /* rectangular domain */
#define D_ELLIPSE 1     /* elliptical domain */
#define D_STADIUM 2     /* stadium-shaped domain */
#define D_SINAI 3       /* Sinai billiard */
#define D_DIAMOND 4     /* diamond-shaped billiard */
#define D_TRIANGLE 5    /* triangular billiard */


#define LAMBDA 0.5	/* parameter controlling shape of billiard */
#define FOCI 1          /* set to 1 to draw focal points of ellipse */

#define MOVIE 0         /* set to 1 to generate movie */
#define RESAMPLE 0      /* set to 1 if particles should be added when dispersion too large */

#define XMIN -2.0
#define XMAX 2.0	/* x interval */
#define YMIN -1.125
#define YMAX 1.125	/* y interval for 9/16 aspect ratio */

#define NPART 10000     /* number of particles */
#define NPARTMAX 50000	/* maximal number of particles after resampling */
#define LMAX 0.01       /* minimal segment length triggering resampling */ 
#define DMIN 0.02       /* minimal distance to boundary for triggering resampling */ 
#define CYCLE 1         /* set to 1 for closed curve (start in all directions) */

#define NSTEPS 5000     /* number of frames of movie */
#define TIME 2000       /* time between movie frames, for fluidity of real-time simulation */ 
#define DPHI 0.00001    /* integration step */
#define NVID 500       /* number of iterations between images displayed on screen */
#define NCOLORS 20      /* number of colors */
#define COLORSHIFT 220  /* hue of initial color */ 
#define NSEG 100        /* number of segments of boundary */
#define LENGTH 0.05     /* length of velocity vectors */

#define BLACK 1         /* set to 1 for black background */

/* Decreasing TIME accelerates the animation and the movie                               */
/* For constant speed of movie, TIME*DPHI should be kept constant                        */
/* However, increasing DPHI too much deterioriates quality of simulation                 */
/* NVID tells how often a picture is drawn in the animation, increase it for faster anim */
/* For a good quality movie, take for instance TIME = 400, DPHI = 0.00005, NVID = 100    */

#define PAUSE 100        /* number of frames after which to pause */
#define PSLEEP 1         /* sleep time during pause */
#define SLEEP1  1        /* initial sleeping time */
#define SLEEP2  100      /* final sleeping time */

#define PI 	3.141592654
#define DPI 	6.283185307
#define PID 	1.570796327

```

### 29 April 21 - The dark side of the pool ###

**Program:** `wave_billiard.c` (parameters reconstructed for 1 August version)

**Initial condition in function `animation()`:** `init_wave(LAMBDA, 0.0, phi, psi, xy_in);`

```
#define MOVIE 1         /* set to 1 to generate movie */

/* General geometrical parameters */

#define WINWIDTH 	1280  /* window width */
#define WINHEIGHT 	720   /* window height */

// #define NX 1280          /* number of grid points on x axis */
// #define NY 720          /* number of grid points on y axis */
#define NX 640          /* number of grid points on x axis */
#define NY 360          /* number of grid points on y axis */

/* setting NX to WINWIDTH and NY to WINHEIGHT increases resolution */
/* but will multiply run time by 4                                 */

#define XMIN -2.0
#define XMAX 2.0	/* x interval */
#define YMIN -1.125
#define YMAX 1.125	/* y interval for 9/16 aspect ratio */

#define JULIA_SCALE 1.0 /* scaling for Julia sets */

/* Choice of the billiard table */

#define B_DOMAIN 1      /* choice of domain shape, see list in global_pdes.c */

#define CIRCLE_PATTERN 10    /* pattern of circles, see list in global_pdes.c */

#define P_PERCOL 0.25       /* probability of having a circle in C_RAND_PERCOL arrangement */
#define NPOISSON 300        /* number of points for Poisson C_RAND_POISSON arrangement */

#define LAMBDA 1.5	    /* parameter controlling the dimensions of domain */
#define MU 0.03	    /* parameter controlling the dimensions of domain */
#define NPOLY 3             /* number of sides of polygon */
#define APOLY 1.0           /* angle by which to turn polygon, in units of Pi/2 */ 
#define MDEPTH 4            /* depth of computation of Menger gasket */
#define MRATIO 3            /* ratio defining Menger gasket */
#define MANDELLEVEL 1000      /* iteration level for Mandelbrot set */
#define MANDELLIMIT 10.0     /* limit value for approximation of Mandelbrot set */
#define FOCI 1              /* set to 1 to draw focal points of ellipse */
#define NGRIDX 15            /* number of grid point for grid of disks */
#define NGRIDY 20           /* number of grid point for grid of disks */

/* You can add more billiard tables by adapting the functions */
/* xy_in_billiard and draw_billiard below */

/* Physical parameters of wave equation */

#define TWOSPEEDS 0         /* set to 1 to replace hardcore boundary by medium with different speed */

#define OMEGA 0.9           /* frequency of periodic excitation */
#define COURANT 0.01       /* Courant number */
#define COURANTB 0.0075      /* Courant number in medium B */
#define GAMMA 0.0      /* damping factor in wave equation */
#define GAMMA_SIDES 1.0e-4      /* damping factor on boundary */
#define GAMMA_TOPBOT 1.0e-6      /* damping factor on boundary */
#define KAPPA 0.0       /* "elasticity" term enforcing oscillations */
#define KAPPA_SIDES 5.0e-4       /* "elasticity" term on absorbing boundary */
#define KAPPA_TOPBOT 0.0       /* "elasticity" term on absorbing boundary */
/* The Courant number is given by c*DT/DX, where DT is the time step and DX the lattice spacing */
/* The physical damping coefficient is given by GAMMA/(DT)^2 */
/* Increasing COURANT speeds up the simulation, but decreases accuracy */
/* For similar wave forms, COURANT^2*GAMMA should be kept constant */

/* Boundary conditions, see list in global_pdes.c  */

#define B_COND 3

/* Parameters for length and speed of simulation */

#define NSTEPS 8000      /* number of frames of movie */
#define NVID 40          /* number of iterations between images displayed on screen */
#define NSEG 100         /* number of segments of boundary */
#define INITIAL_TIME 0    /* time after which to start saving frames */
#define BOUNDARY_WIDTH 3    /* width of billiard boundary */

#define PAUSE 1000         /* number of frames after which to pause */
#define PSLEEP 1         /* sleep time during pause */
#define SLEEP1  1        /* initial sleeping time */
#define SLEEP2  1   /* final sleeping time */

/* Plot type, see list in global_pdes.c  */

#define PLOT 0

/* Color schemes */

#define BLACK 0          /* background */

#define COLOR_SCHEME 0   /* choice of color scheme, see list in global_pdes.c  */

#define SCALE 0          /* set to 1 to adjust color scheme to variance of field */
#define SLOPE 2.0        /* sensitivity of color on wave amplitude */
#define ATTENUATION 0.0  /* exponential attenuation coefficient of contrast with time */
#define E_SCALE 750.0     /* scaling factor for energy representation */

#define COLORHUE 260     /* initial hue of water color for scheme C_LUM */
#define COLORDRIFT 0.0   /* how much the color hue drifts during the whole simulation */
#define LUMMEAN 0.5      /* amplitude of luminosity variation for scheme C_LUM */
#define LUMAMP 0.3       /* amplitude of luminosity variation for scheme C_LUM */
#define HUEMEAN 220.0    /* mean value of hue for color scheme C_HUE */
#define HUEAMP -220.0      /* amplitude of variation of hue for color scheme C_HUE */

```

### 28 April 21 - Squaring the circle ###

**Program:** `vid_drop_billiard.c` (old version of `drop_billiard.c`)

**Initial condition in function `animation()`:** `init_drop_config(0.0, 0.0, 0.0, DPI, configs);`

```
#define WINWIDTH 	1280  /* window width */
#define WINHEIGHT 	720   /* window height */

/* Choice of the billiard table */

#define B_DOMAIN 0      /* choice of domain shape */

#define D_RECTANGLE 0   /* rectangular domain */
#define D_ELLIPSE 1     /* elliptical domain */
#define D_STADIUM 2     /* stadium-shaped domain */
#define D_SINAI 3       /* Sinai billiard */
#define D_DIAMOND 4     /* diamond-shaped billiard */
#define D_TRIANGLE 5    /* triangular billiard */


#define LAMBDA 1.0	/* parameter controlling shape of billiard */
#define FOCI 1          /* set to 1 to draw focal points of ellipse */

#define MOVIE 0         /* set to 1 to generate movie */
#define RESAMPLE 0      /* set to 1 if particles should be added when dispersion too large */

#define XMIN -2.0
#define XMAX 2.0	/* x interval */
#define YMIN -1.125
#define YMAX 1.125	/* y interval for 9/16 aspect ratio */

#define NPART 10000	/* number of particles */
#define NPARTMAX 20000	/* maximal number of particles after resampling */
#define LMAX 0.01       /* minimal segment length triggering resampling */ 
#define DMIN 0.02       /* minimal distance to boundary for triggering resampling */ 
#define CYCLE 0         /* set to 1 for closed curve (start in all directions) */

#define NSTEPS 7500     /* number of frames of movie */
#define TIME 50         /* time between movie frames, for fluidity of real-time simulation */ 
#define DPHI 0.0001     /* integration step */
#define NVID 20         /* number of iterations between images displayed on screen */
                        
#define NCOLORS 10      /* number of colors */
#define COLORSHIFT 200  /* hue of initial color */ 
#define NSEG 100        /* number of segments of boundary */

/* Decreasing TIME accelerates the animation and the movie               */
/* For constant speed of movie, TIME*DPHI should be kept constant        */
/* However, increasing DPHI too much deterioriates quality of simulation */
/* For a good quality movie, take for instance TIME = 50, DPHI = 0.0002  */

#define PAUSE 500          /* number of frames after which to pause */
#define PSLEEP 1         /* sleep time during pause */
#define SLEEP1  1        /* initial sleeping time */
#define SLEEP2  100      /* final sleeping time */

#define PI 	3.141592654
#define DPI 	6.283185307
#define PID 	1.570796327

```

### 27 April 21 - Sunset at Research Triangle Park ###

**Program:** `vid_wave_billiard.c` (old version of `wave_billiard.c`)

**Initial condition in function `animation()`:** `init_wave(0.0, -0.5, phi, psi);`

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

#define B_DOMAIN 5      /* choice of domain shape */

#define D_RECTANGLE 0   /* rectangular domain */
#define D_ELLIPSE 1     /* elliptical domain */
#define D_STADIUM 2     /* stadium-shaped domain */
#define D_SINAI 3       /* Sinai billiard */
#define D_DIAMOND 4     /* diamond-shaped billiard */
#define D_TRIANGLE 5    /* triangular billiard */

#define LAMBDA 1.732050808	    /* parameter controlling the dimensions of domain */
#define FOCI 1              /* set to 1 to draw focal points of ellipse */

/* You can add more billiard tables by adapting the functions */
/* xy_in_billiard and draw_billiard below */

/* Physical parameters of wave equation */

#define COURANT 0.01       /* Courant number */
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

#define NSTEPS 7500      /* number of frames of movie */
#define NVID 75          /* number of iterations between images displayed on screen */
#define NSEG 100         /* number of segments of boundary */

#define PAUSE 10         /* number of frames after which to pause */
#define PSLEEP 1         /* sleep time during pause */
#define SLEEP1  1        /* initial sleeping time */
#define SLEEP2  100000   /* final sleeping time */

/* Color schemes */

#define COLOR_SCHEME 1   /* choice of color scheme */

#define C_LUM 0          /* color scheme modifies luminosity (with slow drift of hue) */
#define C_HUE 1          /* color scheme modifies hue */

#define SCALE 1          /* set to 1 to adjust color scheme to variance of field */
#define SLOPE 1.0        /* sensitivity of color on wave amplitude */
#define ATTENUATION 0.0  /* exponential attenuation coefficient of contrast with time */

#define COLORHUE 220     /* initial hue of water color for scheme C_LUM */ 
#define COLORDRIFT -40.0 /* how much the color hue drifts during the whole simulation */
#define LUMMEAN 0.5      /* amplitude of luminosity variation for scheme C_LUM */
#define LUMAMP 0.2       /* amplitude of luminosity variation for scheme C_LUM */
#define HUEMEAN 330.0    /* mean value of hue for color scheme C_HUE */
#define HUEAMP 90.0      /* amplitude of variation of hue for color scheme C_HUE */

/* Basic math */

#define PI 	3.141592654
#define DPI 	6.283185307
#define PID 	1.570796327

```

### 26 April 21 - Is the billiard in a stadium reversible? ###

**Program:** Variant of `vid_part_billiard.c`

**Initial condition in function `animation()`:** `init_drop_config(0.5, -1.0, 1.8, 2.0, configs);`

```
#define WINWIDTH 	1280  /* window width */
#define WINHEIGHT 	720   /* window height */

/* Choice of the billiard table */

#define B_DOMAIN 2      /* choice of domain shape */

#define D_RECTANGLE 0   /* rectangular domain */
#define D_ELLIPSE 1     /* elliptical domain */
#define D_STADIUM 2     /* stadium-shaped domain */
#define D_SINAI 3       /* Sinai billiard */
#define D_DIAMOND 4     /* diamond-shaped billiard */

#define LAMBDA 1.0	/* aspect ratio of ellipse */
#define FOCI 1          /* set to 1 to draw focal points of ellipse */

#define MOVIE 1         /* set to 1 to generate movie */
#define RESAMPLE 0      /* set to 1 if particles should be added when dispersion too large */

#define XMIN -2.0
#define XMAX 2.0	/* x interval */
#define YMIN -1.125
#define YMAX 1.125	/* y interval for 9/16 aspect ratio */

#define NPART 200       /* number of particles */
#define NPARTMAX 50000	/* maximal number of particles after resampling */
#define LMAX 0.01       /* minimal segment length triggering resampling */ 
#define DMIN 0.02       /* minimal distance to boundary for triggering resampling */ 
#define CYCLE 1         /* set to 1 for closed curve (start in all directions) */

#define NSTEPS 3000     /* number of frames of movie */
#define TIME 1000       /* time between movie frames, for fluidity of real-time simulation */ 
#define DPHI 0.00001    /* integration step */
#define NVID 100        /* number of iterations between images displayed on screen */
#define NCOLORS 20      /* number of colors */
#define COLORSHIFT 220  /* hue of initial color */ 
#define NSEG 100        /* number of segments of boundary */
#define LENGTH 0.05     /* length of velocity vectors */

/* Decreasing TIME accelerates the animation and the movie                               */
/* For constant speed of movie, TIME*DPHI should be kept constant                        */
/* However, increasing DPHI too much deterioriates quality of simulation                 */
/* NVID tells how often a picture is drawn in the animation, increase it for faster anim */
/* For a good quality movie, take for instance TIME = 400, DPHI = 0.00005, NVID = 100    */

#define PAUSE 10          /* number of frames after which to pause */
#define PSLEEP 1         /* sleep time during pause */
#define SLEEP1  1        /* initial sleeping time */
#define SLEEP2  100      /* final sleeping time */

#define PI 	3.141592654
#define DPI 	6.283185307
#define PID 	1.570796327

```

### 25 April 21 - Like tears in rain ###

**Program:** `vid_wave_billiard.c` (old version of `wave_billiard.c`)

**Initial condition in function `animation()`:** `init_wave(0.4, 0.4, phi, psi);`

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

#define B_DOMAIN 1      /* choice of domain shape */

#define D_RECTANGLE 0   /* rectangular domain */
#define D_ELLIPSE 1     /* elliptical domain */
#define D_STADIUM 2     /* stadium-shaped domain */
#define D_SINAI 3       /* Sinai billiard */
#define D_DIAMOND 4     /* diamond-shaped billiard */

#define LAMBDA 1.0	    /* parameter controlling the dimensions of domain */
#define FOCI 1              /* set to 1 to draw focal points of ellipse */

/* You can add more billiard tables by adapting the functions */
/* xy_in_billiard and draw_billiard below */

/* Physical patameters of wave equation */

#define COURANT 0.01       /* Courant number */
#define GAMMA 5.0e-10      /* damping factor in wave equation */
#define KAPPA 5.0e-8       /* "elasticity" term enforcing oscillations */
/* The Courant number is given by c*DT/DX, where DT is the time step and DX the lattice spacing */
/* The physical damping coefficient is given by GAMMA/(DT)^2 */
/* Increasing COURANT speeds up the simulation, but decreases accuracy */
/* For similar wave forms, COURANT^2*GAMMA should be kept constant */

/* For debugging purposes only */
#define FLOOR 0         /* set to 1 to limit wave amplitude to VMAX */
#define VMAX 10.0       /* max value of wave amplitude */

/* Parameters for length and speed of simulation */

#define NSTEPS 6000      /* number of frames of movie */
#define NVID 50          /* number of iterations between images displayed on screen */
#define NSEG 100         /* number of segments of boundary */

#define PAUSE 10         /* number of frames after which to pause */
#define PSLEEP 1         /* sleep time during pause */
#define SLEEP1  1        /* initial sleeping time */
#define SLEEP2  100000   /* final sleeping time */

/* Color schemes */

#define COLOR_SCHEME 1   /* choice of color scheme */

#define C_LUM 0          /* color scheme modifies luminosity (with slow drift of hue) */
#define C_HUE 1          /* color scheme modifies hue */

#define SCALE 1          /* set to 1 to adjust color scheme to variance of field */
#define SLOPE 1.0        /* sensitivity of color on wave amplitude */
#define ATTENUATION 0.0  /* exponential attenuation coefficient of contrast with time */

#define COLORHUE 220     /* initial hue of water color for scheme C_LUM */ 
#define COLORDRIFT -40.0 /* how much the color hue drifts during the whole simulation */
#define LUMMEAN 0.5      /* amplitude of luminosity variation for scheme C_LUM */
#define LUMAMP 0.2       /* amplitude of luminosity variation for scheme C_LUM */
#define HUEMEAN 150.0    /* mean value of hue for color scheme C_HUE */
#define HUEAMP -120.0      /* amplitude of variation of hue for color scheme C_HUE */

/* Basic math */

#define PI 	3.141592654
#define DPI 	6.283185307
#define PID 	1.570796327

```

### 24 April 21 - Contrails in the elliptical billiard ###

**Program:** Variant of `vid_part_billiard.c` (contrail part not yet published)

**Initial condition in function `animation()`:** `init_boundary_config(-0.35, -0.15, PID, PID, configs);`

```
#define WINWIDTH 	1280  /* window width */
#define WINHEIGHT 	720   /* window height */

#define LAMBDA 1.5	/* aspect ratio of ellipse */
#define FOCI 1          /* set to 1 to draw focal points of ellipse */

#define MOVIE 0         /* set to 1 to generate movie */
#define RESAMPLE 0      /* set to 1 if particles should be added when dispersion too large */

#define XMIN -2.0
#define XMAX 2.0	/* x interval */
#define YMIN -1.125
#define YMAX 1.125	/* y interval for 9/16 aspect ratio */

#define NPART 501       /* number of particles */
#define NPARTMAX 50000	/* maximal number of particles after resampling */
#define LMAX 0.01       /* minimal segment length triggering resampling */ 
#define DMIN 0.02       /* minimal distance to boundary for triggering resampling */ 
#define CYCLE 1         /* set to 1 for closed curve (start in all directions) */

#define NSTEPS 4500     /* number of frames of movie */
#define TIME 2000       /* time between movie frames, for fluidity of real-time simulation */ 
#define DPHI 0.00003    /* integration step */
#define NVID 500        /* number of iterations between images displayed on screen */
#define NCOLORS 20      /* number of colors */
#define COLORSHIFT 220  /* hue of initial color */ 
#define NSEG 100        /* number of segments of boundary */
#define LENGTH 0.05     /* length of velocity vectors */

/* Decreasing TIME accelerates the animation and the movie                               */
/* For constant speed of movie, TIME*DPHI should be kept constant                        */
/* However, increasing DPHI too much deterioriates quality of simulation                 */
/* NVID tells how often a picture is drawn in the animation, increase it for faster anim */
/* For a good quality movie, take for instance TIME = 400, DPHI = 0.00005, NVID = 100    */

#define PAUSE 10         /* number of frames after which to pause */
#define PSLEEP 1         /* sleep time during pause */
#define SLEEP1  1        /* initial sleeping time */
#define SLEEP2  100      /* final sleeping time */

#define PI 	3.141592654
#define DPI 	6.283185307
#define PID 	1.570796327

```

### 23 April 21 - Low sea in the square with diamonds ###

**Program:** `vid_wave_billiard.c` (old version of `wave_billiard.c`)

**Initial condition in function `animation()`:** `init_wave(0.0, 0.0, phi, psi);`

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

#define B_DOMAIN 0      /* choice of domain shape */

#define D_RECTANGLE 0   /* rectangular domain */
#define D_ELLIPSE 1     /* elliptical domain */
#define D_STADIUM 2     /* stadium-shaped domain */
#define D_SINAI 3       /* Sinai billiard */
#define D_DIAMOND 4     /* diamond-shaped billiard */

#define LAMBDA 1.0	    /* parameter controlling the dimensions of domain */
#define FOCI 1              /* set to 1 to draw focal points of ellipse */

/* You can add more billiard tables by adapting the functions */
/* xy_in_billiard and draw_billiard below */

/* Physical patameters of wave equation */

#define COURANT 0.01       /* Courant number */
#define GAMMA 5.0e-06      /* damping factor in wave equation */
/* The Courant number is given by c*DT/DX, where DT is the time step and DX the lattice spacing */
/* The physical damping coefficient is given by GAMMA/(DT)^2 */
/* Increasing COURANT speeds up the simulation, but decreases accuracy */
/* For similar wave forms, COURANT^2*GAMMA should be kept constant */

/* For debugging purposes only */
#define FLOOR 0         /* set to 1 to limit wave amplitude to VMAX */
#define VMAX 10.0       /* max value of wave amplitude */

/* Parameters for length and speed of simulation */

#define NSTEPS 6000      /* number of frames of movie */
#define NVID 35          /* number of iterations between images displayed on screen */
#define NSEG 100         /* number of segments of boundary */

#define PAUSE 10         /* number of frames after which to pause */
#define PSLEEP 1         /* sleep time during pause */
#define SLEEP1  1        /* initial sleeping time */
#define SLEEP2  100000   /* final sleeping time */

/* Color schemes */

#define COLOR_SCHEME 1   /* choice of color scheme */

#define C_LUM 0          /* color scheme modifies luminosity (with slow drift of hue) */
#define C_HUE 1          /* color scheme modifies hue */

#define SCALE 1          /* set to 1 to adjust color scheme to variance of field */
#define SLOPE 1.0        /* sensitivity of color on wave amplitude */
#define ATTENUATION 0.0  /* exponential attenuation coefficient of contrast with time */

#define COLORHUE 220     /* initial hue of water color for scheme C_LUM */ 
#define COLORDRIFT -40.0 /* how much the color hue drifts during the whole simulation */
#define LUMMEAN 0.5      /* amplitude of luminosity variation for scheme C_LUM */
#define LUMAMP 0.2       /* amplitude of luminosity variation for scheme C_LUM */
#define HUEMEAN 260.0    /* mean value of hue for color scheme C_HUE */
#define HUEAMP 80.0      /* amplitude of variation of hue for color scheme C_HUE */

/* Basic math */

#define PI 	3.141592654
#define DPI 	6.283185307
#define PID 	1.570796327

```

### 22 April 21 - Just another drop in the pool ###

**Program:** `vid_wave_billiard.c` (old version of `wave_billiard.c`)

**Initial condition in function `animation()`:** `init_wave(sqrt(LAMBDA*LAMBDA - 1.0), 0.0, phi, psi); add_drop_to_wave(-1.0,-sqrt(LAMBDA*LAMBDA - 1.0), 0.0, phi, psi);`

```
#define WINWIDTH 	1280  /* window width */
#define WINHEIGHT 	720   /* window height */

#define NX 640          /* number of grid points on x axis */
#define NY 360          /* number of grid points on y axis */

#define B_DOMAIN 1      /* choice of domain shape */

#define D_RECTANGLE 0   /* rectangular domain */
#define D_ELLIPSE 1     /* elliptical domain */
#define D_STADIUM 2     /* stadium-shaped domain */
#define D_SINAI 3       /* Sinai billiard */
#define D_DIAMOND 4     /* diamond-shaped billiard */

#define LAMBDA 1.5	    /* parameter controlling the dimensions of domain */
#define FOCI 1              /* set to 1 to draw focal points of ellipse */
#define COURANT 0.01       /* Courant number */
#define GAMMA 5.0e-05         /* damping factor in wave equation */
// The Courant number is given by c*DT/DX, where DT is the time step and DX the lattice spacing
// The physical damping coefficient is given by GAMMA/(DT)^2

#define MOVIE 1         /* set to 1 to generate movie */

#define XMIN -2.0
#define XMAX 2.0	/* x interval */
#define YMIN -1.125
#define YMAX 1.125	/* y interval for 9/16 aspect ratio */

#define VMAX 10.0       /* max value of wave amplitude */
#define FLOOR 0         /* set to 1 to limit wave amplitude to VMAX */

#define NSTEPS 4500      /* number of frames of movie */
#define NVID 35          /* number of iterations between images displayed on screen */
#define NCOLORS 20       /* number of colors */
#define COLORHUE 220     /* initial hue of water color */ 
#define COLORDRIFT -40.0 /* how much the color hue drifts during the whole simulation */
#define NSEG 100         /* number of segments of boundary */
#define SLOPE 1.0        /* for color scheme: sensitivity of color on wave amplitude */
#define ATTENUATION 0.0  /* exponential attenuation coefficient of contrat with time */
#define LUMMEAN 0.5      /* amplitude of luminosity variation */
#define LUMAMP 0.2       /* amplitude of luminosity variation */
#define SCALE 1          /* set to 1 to adjust color scheme to variance of field */

/* Decreasing TIME accelerates the animation and the movie                               */
/* For constant speed of movie, TIME*DPHI should be kept constant                        */
/* However, increasing DPHI too much deterioriates quality of simulation                 */
/* NVID tells how often a picture is drawn in the animation, increase it for faster anim */
/* For a good quality movie, take for instance TIME = 400, DPHI = 0.00005, NVID = 100    */

#define PAUSE 10         /* number of frames after which to pause */
#define PSLEEP 1         /* sleep time during pause */
#define SLEEP1  1        /* initial sleeping time */
#define SLEEP2  100000   /* final sleeping time */

#define PI 	3.141592654
#define DPI 	6.283185307
#define PID 	1.570796327

```

### 21 April 21 - Rain on my crazy diamond ###

**Program:** `vid_wave_billiard.c` (old version of `wave_billiard.c`)

**Initial condition in function `animation()`:** `init_wave(0.0, 0.0, phi, psi);`

```
#define WINWIDTH 	1280  /* window width */
#define WINHEIGHT 	720   /* window height */

#define NX 640          /* number of grid points on x axis */
#define NY 360          /* number of grid points on y axis */

#define B_DOMAIN 4      /* choice of domain shape */

#define D_RECTANGLE 0   /* rectangular domain */
#define D_ELLIPSE 1     /* elliptical domain */
#define D_STADIUM 2     /* stadium-shaped domain */
#define D_SINAI 3       /* Sinai billiard */
#define D_DIAMOND 4     /* diamond-shaped billiard */

#define LAMBDA 2.0	    /* parameter controlling the dimensions of domain */
#define FOCI 1              /* set to 1 to draw focal points of ellipse */
#define COURANT 0.01       /* Courant number */
#define GAMMA 5.0e-05         /* damping factor in wave equation */
// The Courant number is given by c*DT/DX, where DT is the time step and DX the lattice spacing
// The physical damping coefficient is given by GAMMA/(DT)^2


#define MOVIE 1         /* set to 1 to generate movie */

#define XMIN -2.0
#define XMAX 2.0	/* x interval */
#define YMIN -1.125
#define YMAX 1.125	/* y interval for 9/16 aspect ratio */

#define VMAX 10.0       /* max value of wave amplitude */
#define FLOOR 0         /* set to 1 to limit wave amplitude to VMAX */

#define NSTEPS 7500      /* number of frames of movie */
#define NVID 25          /* number of iterations between images displayed on screen */
#define NCOLORS 20       /* number of colors */
#define COLORHUE 220     /* initial hue of water color */ 
#define COLORDRIFT -40.0 /* how much the color hue drifts during the whole simulation */
#define NSEG 100         /* number of segments of boundary */
#define SLOPE 1.0        /* for color scheme: sensitivity of color on wave amplitude */
#define ATTENUATION 0.0  /* exponential attenuation coefficient of contrat with time */
#define LUMMEAN 0.5      /* amplitude of luminosity variation */
#define LUMAMP 0.2       /* amplitude of luminosity variation */
#define SCALE 1          /* set to 1 to adjust color scheme to variance of field */

/* Decreasing TIME accelerates the animation and the movie                               */
/* For constant speed of movie, TIME*DPHI should be kept constant                        */
/* However, increasing DPHI too much deterioriates quality of simulation                 */
/* NVID tells how often a picture is drawn in the animation, increase it for faster anim */
/* For a good quality movie, take for instance TIME = 400, DPHI = 0.00005, NVID = 100    */

#define PAUSE 10         /* number of frames after which to pause */
#define PSLEEP 1         /* sleep time during pause */
#define SLEEP1  1        /* initial sleeping time */
#define SLEEP2  100000   /* final sleeping time */

#define PI 	3.141592654
#define DPI 	6.283185307
#define PID 	1.570796327

```

### 20 April 21 - A saucerful of droplets ###

**Program:** Variant of `vid_wave_billiard.c` (old version of `vid_wave_billiard.c`)

**Initial condition in function `animation()`:** Raindrops are generated by 
```
while (!in_bill)
    {
        x = (double)rand() / RAND_MAX;
        y = (double)rand() / RAND_MAX;
        x = LAMBDA*(2.0*x - 1.0);
        y = 2.0*y - 1.0;
        in_bill = xy_in_billiard(x, y);
    }
    init_wave(x, y, phi, psi);
```
    

```
#define WINWIDTH 	1280  /* window width */
#define WINHEIGHT 	720   /* window height */

#define NX 640          /* number of grid points on x axis */
#define NY 360          /* number of grid points on y axis */

#define B_DOMAIN 1      /* choice of domain shape */

#define D_RECTANGLE 0   /* rectangular domain */
#define D_ELLIPSE 1     /* elliptical domain */
#define D_STADIUM 2     /* stadium-shaped domain */
#define D_SINAI 3       /* Sinai billiard */

#define LAMBDA 1.3	    /* parameter controlling the dimensions of domain */
#define FOCI 0              /* set to 1 to draw focal points of ellipse */
#define COURANT 0.003125    /* Courant number */
#define GAMMA 1e-06         /* damping factor in wave equation */
// The Courant number is given by c*DT/DX, where DT is the time step and DX the lattice spacing
// The physical damping coefficient is given by GAMMA/(DT)^2


#define MOVIE 1         /* set to 1 to generate movie */

#define XMIN -2.0
#define XMAX 2.0	/* x interval */
#define YMIN -1.125
#define YMAX 1.125	/* y interval for 9/16 aspect ratio */

#define VMAX 10.0       /* max value of wave amplitude */
#define FLOOR 0         /* set to 1 to limit wave amplitude to VMAX */

#define NSTEPS 3000      /* number of frames of movie */
#define NVID 60          /* number of iterations between images displayed on screen */
#define NCOLORS 20       /* number of colors */
#define COLORHUE 200     /* initial hue of water color */ 
#define COLORDRIFT -20.0 /* how much the color hue drifts during the whole simulation */
#define NSEG 100         /* number of segments of boundary */
#define SLOPE 0.75       /* for color scheme: sensitivity of color on wave amplitude */
#define ATTENUATION 0.0  /* exponential attenuation coefficient of contrat with time */
#define LUMMEAN 0.5      /* amplitude of luminosity variation */
#define LUMAMP 0.2       /* amplitude of luminosity variation */
#define SCALE 1          /* set to 1 to adjust color scheme to variance of field */

/* Decreasing TIME accelerates the animation and the movie                               */
/* For constant speed of movie, TIME*DPHI should be kept constant                        */
/* However, increasing DPHI too much deterioriates quality of simulation                 */
/* NVID tells how often a picture is drawn in the animation, increase it for faster anim */
/* For a good quality movie, take for instance TIME = 400, DPHI = 0.00005, NVID = 100    */

#define PAUSE 10         /* number of frames after which to pause */
#define PSLEEP 1         /* sleep time during pause */
#define SLEEP1  1        /* initial sleeping time */
#define SLEEP2  100000   /* final sleeping time */

#define PI 	3.141592654
#define DPI 	6.283185307
#define PID 	1.570796327

```

### 19 April 21 - Raindrop in the Sinai billiard ###

**Program:** `vid_wave_billiard.c` (old version of `wave_billiard.c`)

**Initial condition in function `animation()`:** `init_wave(0.6, 0.9, phi, psi);`

```
#define WINWIDTH 	1280  /* window width */
#define WINHEIGHT 	720   /* window height */

// #define NX 1280         /* number of grid points on x axis */
// #define NY 720          /* number of grid points on y axis */
#define NX 640          /* number of grid points on x axis */
#define NY 360          /* number of grid points on y axis */

#define B_DOMAIN 3      /* choice of domain shape */

#define D_RECTANGLE 0   /* rectangular domain */
#define D_ELLIPSE 1     /* elliptical domain */
#define D_STADIUM 2     /* stadium-shaped domain */
#define D_SINAI 3       /* Sinai billiard */

#define LAMBDA 0.7	/* parameter controlling the dimensions of domain */
#define FOCI 1          /* set to 1 to draw focal points of ellipse */
#define K 20.0          /* coupling constant */
#define GAMMA 0.00005   /* damping factor in wave equation */
#define DT 0.05         /* time step */ 

#define MOVIE 0         /* set to 1 to generate movie */

#define XMIN -2.0
#define XMAX 2.0	/* x interval */
#define YMIN -1.125
#define YMAX 1.125	/* y interval for 9/16 aspect ratio */

#define VMAX 10.0       /* max value of wave amplitude */
#define FLOOR 0         /* set to 1 to limit wave amplitude to VMAX */

#define NSTEPS 7500      /* number of frames of movie */
#define NVID 60          /* number of iterations between images displayed on screen */
#define NCOLORS 20       /* number of colors */
#define COLORHUE 180     /* hue of water color */ 
#define COLORDRIFT 180.0 /* how much the color hue drifts during the whole simulation */
#define NSEG 100         /* number of segments of boundary */
#define SLOPE 0.75       /* for color scheme: sensitivity of color on wave amplitude */
#define ATTENUATION 0.0  /* exponential attenuation coefficient of contrat with time */
#define LUMAMP 0.2       /* amplitude of luminosity variation */
#define SCALE 1          /* set to 1 to adjust color scheme to variance of field */

/* Decreasing TIME accelerates the animation and the movie                               */
/* For constant speed of movie, TIME*DPHI should be kept constant                        */
/* However, increasing DPHI too much deterioriates quality of simulation                 */
/* NVID tells how often a picture is drawn in the animation, increase it for faster anim */
/* For a good quality movie, take for instance TIME = 400, DPHI = 0.00005, NVID = 100    */

#define PAUSE 10         /* number of frames after which to pause */
#define PSLEEP 1         /* sleep time during pause */
#define SLEEP1  1        /* initial sleeping time */
#define SLEEP2  100000   /* final sleeping time */

#define PI 	3.141592654
#define DPI 	6.283185307
#define PID 	1.570796327

```

### 18 April 21 - Dropping a stone in a flooded stadium ###

**Program:** `xxx.c`

**Initial condition in function `animation()`:** `xxx`

```
#define WINWIDTH 	1280  /* window width */
#define WINHEIGHT 	720   /* window height */

#define NX 640          /* number of grid points on x axis */
#define NY 360          /* number of grid points on y axis */

#define B_DOMAIN 1        /* choice of domain shape */

#define D_ELLIPSE 0     /* elliptical domain */
#define D_STADIUM 1     /* stadium-shaped domain */

#define LAMBDA 1.0	/* parameter controlling the dimensions of domain */
#define FOCI 1          /* set to 1 to draw focal points of ellipse */
#define K 20.0          /* coupling constant */
#define GAMMA 0.00005   /* damping factor in wave equation */
#define DT 0.05         /* time step */ 

#define MOVIE 0         /* set to 1 to generate movie */

#define XMIN -2.0
#define XMAX 2.0	/* x interval */
#define YMIN -1.125
#define YMAX 1.125	/* y interval for 9/16 aspect ratio */

#define VMAX 10.0       /* max value of wave amplitude */
#define FLOOR 0         /* set to 1 to limit wave amplitude to VMAX */

#define NSTEPS 7500     /* number of frames of movie */
#define NVID 60         /* number of iterations between images displayed on screen */
#define NCOLORS 20      /* number of colors */
#define COLORHUE 200    /* hue of water color */ 
#define NSEG 100        /* number of segments of boundary */
#define SLOPE 1.25      /* for color scheme: sensitivity of color on wave amplitude */
#define ATTENUATION 0.0 /* exponential attenuation coefficient of contrat with time */
#define LUMAMP 0.25     /* amplitude of luminosity variation */

/* Decreasing TIME accelerates the animation and the movie                               */
/* For constant speed of movie, TIME*DPHI should be kept constant                        */
/* However, increasing DPHI too much deterioriates quality of simulation                 */
/* NVID tells how often a picture is drawn in the animation, increase it for faster anim */
/* For a good quality movie, take for instance TIME = 400, DPHI = 0.00005, NVID = 100    */

#define PAUSE 10         /* number of frames after which to pause */
#define PSLEEP 1         /* sleep time during pause */
#define SLEEP1  1        /* initial sleeping time */
#define SLEEP2  100000   /* final sleeping time */

#define PI 	3.141592654
#define DPI 	6.283185307
#define PID 	1.570796327

```

### 17 April 21 - Waves in an elliptical swimming pool 3: Starting from a focus ###

**Program:** `vid_wave_ellipse.c` (old version of `wave_billiard.c`)

**Initial condition in function `animation()`:** 
`init_drop_config(sqrt(LAMBDA*LAMBDA-1.0) - 0.1,0.0, 0.0, DPI, configs);`

```
#define WINWIDTH 	1280  /* window width */
#define WINHEIGHT 	720   /* window height */

#define LAMBDA 1.5	/* aspect ratio of ellipse */
#define FOCI 1          /* set to 1 to draw focal points of ellipse */

#define MOVIE 1         /* set to 1 to generate movie */
#define RESAMPLE 1      /* set to 1 if particles should be added when dispersion too large */

#define XMIN -2.0
#define XMAX 2.0	/* x interval */
#define YMIN -1.125
#define YMAX 1.125	/* y interval for 9/16 aspect ratio */

#define NPART 2000	/* number of particles */
#define NPARTMAX 200000	/* maximal number of particles after resampling */
#define LMAX 0.01       /* minimal segment length triggering resampling */ 
#define DMIN 0.02       /* minimal distance to boundary for triggering resampling */ 
#define CYCLE 1         /* set to 1 for closed curve (start in all directions) */

#define NSTEPS 130      /* number of frames of movie */
#define TIME 25         /* time between movie frames, for fluidity of real-time simulation */ 
#define DPHI 0.0004     /* integration step */
#define NVID 10         /* number of iterations between images displayed on screen */
#define NCOLORS 10      /* number of colors */
#define COLORSHIFT 0    /* hue of initial color */ 
#define NSEG 100        /* number of segments of boundary */

/* Decreasing TIME accelerates the animation and the movie               */
/* For constant speed of movie, TIME*DPHI should be kept constant        */
/* However, increasing DPHI too much deterioriates quality of simulation */
/* For a good quality movie, take for instance TIME = 50, DPHI = 0.0002  */

#define PAUSE 1          /* number of frames after which to pause */
#define PSLEEP 1         /* sleep time during pause */
#define SLEEP1  1        /* initial sleeping time */
#define SLEEP2  100      /* final sleeping time */

#define PI 	3.141592654
#define DPI 	6.283185307
#define PID 	1.570796327

```

### 17 April 21 - Waves in an elliptical swimming pool 2 ###

**Program:** `vid_wave_ellipse.c` (old version of `wave_billiard.c`)

**Initial condition in function `animation()`:** `init_wave(0.8, 0.6, phi, psi);`

```
#include <math.h>
#include <string.h>
#include <GL/glut.h>
#include <GL/glu.h>
#include <unistd.h>
#include <sys/types.h>
#include <tiffio.h>     /* Sam Leffler's libtiff library. */

#define WINWIDTH 	1280  /* window width */
#define WINHEIGHT 	720   /* window height */

#define NX 640          /* number of grid points on x axis */
#define NY 360          /* number of grid points on y axis */

#define B_DOMAIN 1        /* choice of domain shape */

#define D_ELLIPSE 0     /* elliptical domain */
#define D_STADIUM 1     /* stadium-shaped domain */

#define LAMBDA 1.0	/* parameter controlling the dimensions of domain */
#define FOCI 1          /* set to 1 to draw focal points of ellipse */
#define K 20.0          /* coupling constant */
#define GAMMA 0.00005   /* damping factor in wave equation */
#define DT 0.05         /* time step */ 

#define MOVIE 0         /* set to 1 to generate movie */

#define XMIN -2.0
#define XMAX 2.0	/* x interval */
#define YMIN -1.125
#define YMAX 1.125	/* y interval for 9/16 aspect ratio */

#define VMAX 10.0       /* max value of wave amplitude */
#define FLOOR 0         /* set to 1 to limit wave amplitude to VMAX */

#define NSTEPS 4500     /* number of frames of movie */
#define NVID 60         /* number of iterations between images displayed on screen */
#define NCOLORS 20      /* number of colors */
#define COLORHUE 200    /* hue of water color */ 
#define NSEG 100        /* number of segments of boundary */
#define SLOPE 1.25      /* for color scheme: sensitivity of color on wave amplitude */
#define ATTENUATION 0.0 /* exponential attenuation coefficient of contrat with time */
#define LUMMEAN 0.5      /* amplitude of luminosity variation */
#define LUMAMP 0.2       /* amplitude of luminosity variation */
#define SCALE 1          /* set to 1 to adjust color scheme to variance of field */

/* Decreasing TIME accelerates the animation and the movie                               */
/* For constant speed of movie, TIME*DPHI should be kept constant                        */
/* However, increasing DPHI too much deterioriates quality of simulation                 */
/* NVID tells how often a picture is drawn in the animation, increase it for faster anim */
/* For a good quality movie, take for instance TIME = 400, DPHI = 0.00005, NVID = 100    */

#define PAUSE 10         /* number of frames after which to pause */
#define PSLEEP 1         /* sleep time during pause */
#define SLEEP1  1        /* initial sleeping time */
#define SLEEP2  100000   /* final sleeping time */

#define PI 	3.141592654
#define DPI 	6.283185307
#define PID 	1.570796327

```

### 16 April 21 - Waves in an elliptical swimming pool ###

**Program:** `vid_wave_billiard.c` (old version of `wave_billiard.c`)

**Initial condition in function `animation()`:** `init_wave(0.8, 0.6, phi, psi);`

```
#define WINWIDTH 	1280  /* window width */
#define WINHEIGHT 	720   /* window height */

#define NX 640          /* number of grid points on x axis */
#define NY 360          /* number of grid points on y axis */

#define LAMBDA 1.5	/* aspect ratio of ellipse */
#define FOCI 1          /* set to 1 to draw focal points of ellipse */
#define K 20.0          /* coupling constant */
#define GAMMA 0.00005       /* damping */
#define DT 0.05         /* time step */ 
#define PHIMAX 1.0      /* max value of field */

#define MOVIE 1         /* set to 1 to generate movie */

#define XMIN -2.0
#define XMAX 2.0	/* x interval */
#define YMIN -1.125
#define YMAX 1.125	/* y interval for 9/16 aspect ratio */

#define VMAX 2.0        /* max value for color scheme */

#define NSTEPS 4500     /* number of frames of movie */
#define NVID 50         /* number of iterations between images displayed on screen */
#define NCOLORS 20      /* number of colors */
#define COLORSHIFT 220  /* hue of initial color */ 
#define NSEG 100        /* number of segments of boundary */

/* Decreasing TIME accelerates the animation and the movie                               */
/* For constant speed of movie, TIME*DPHI should be kept constant                        */
/* However, increasing DPHI too much deterioriates quality of simulation                 */
/* NVID tells how often a picture is drawn in the animation, increase it for faster anim */
/* For a good quality movie, take for instance TIME = 400, DPHI = 0.00005, NVID = 100    */

#define PAUSE 10         /* number of frames after which to pause */
#define PSLEEP 1         /* sleep time during pause */
#define SLEEP1  1        /* initial sleeping time */
#define SLEEP2  100      /* final sleeping time */

#define PI 	3.141592654
#define DPI 	6.283185307
#define PID 	1.570796327

```

### 15 April 21 - Starting near a focus in the elliptical billiard ###

**Program:** `vid_particle_ellipse.c` (Old version of `particle_billiard.c`)

**Initial condition in function `animation()`:** `init_drop_config(sqrt(LAMBDA*LAMBDA-1.0)-0.01,0.0, 0.0, DPI);`

```
#define WINWIDTH 	1280  /* window width */
#define WINHEIGHT 	720   /* window height */

#define LAMBDA 1.2	/* aspect ratio of ellipse */
#define FOCI 1          /* set to 1 to draw focal points of ellipse */

#define MOVIE 1         /* set to 1 to generate movie */
#define RESAMPLE 0      /* set to 1 if particles should be added when dispersion too large */

#define XMIN -2.0
#define XMAX 2.0	/* x interval */
#define YMIN -1.125
#define YMAX 1.125	/* y interval for 9/16 aspect ratio */

#define NPART 500       /* number of particles */
#define NPARTMAX 50000	/* maximal number of particles after resampling */
#define LMAX 0.01       /* minimal segment length triggering resampling */ 
#define DMIN 0.02       /* minimal distance to boundary for triggering resampling */ 
#define CYCLE 1         /* set to 1 for closed curve (start in all directions) */

#define NSTEPS 4500     /* number of frames of movie */
#define TIME 2000       /* time between movie frames, for fluidity of real-time simulation */ 
#define DPHI 0.00001    /* integration step */
#define NVID 500        /* number of iterations between images displayed on screen */
#define NCOLORS 20      /* number of colors */
#define COLORSHIFT 220  /* hue of initial color */ 
#define NSEG 100        /* number of segments of boundary */
#define LENGTH 0.05     /* length of velocity vectors */

/* Decreasing TIME accelerates the animation and the movie                               */
/* For constant speed of movie, TIME*DPHI should be kept constant                        */
/* However, increasing DPHI too much deterioriates quality of simulation                 */
/* NVID tells how often a picture is drawn in the animation, increase it for faster anim */
/* For a good quality movie, take for instance TIME = 400, DPHI = 0.00005, NVID = 100    */

#define PAUSE 10         /* number of frames after which to pause */
#define PSLEEP 1         /* sleep time during pause */
#define SLEEP1  1        /* initial sleeping time */
#define SLEEP2  100      /* final sleeping time */

#define PI 	3.141592654
#define DPI 	6.283185307
#define PID 	1.570796327

```

### 14 April 21 - Billiard in an ellipse ###

**Program:** `vid_particle_ellipse.c` (Old version of `particle_billiard.c`)

**Initial condition in function `animation()`:** ?

```
#define WINWIDTH 	1280  /* window width */
#define WINHEIGHT 	720   /* window height */

#define LAMBDA 1.2	/* aspect ratio of ellipse */
#define FOCI 1          /* set to 1 to draw focal points of ellipse */

#define MOVIE 1         /* set to 1 to generate movie */
#define RESAMPLE 0      /* set to 1 if particles should be added when dispersion too large */

#define XMIN -2.0
#define XMAX 2.0	/* x interval */
#define YMIN -1.125
#define YMAX 1.125	/* y interval for 9/16 aspect ratio */

#define NPART 500       /* number of particles */
#define NPARTMAX 50000	/* maximal number of particles after resampling */
#define LMAX 0.01       /* minimal segment length triggering resampling */ 
#define DMIN 0.02       /* minimal distance to boundary for triggering resampling */ 
#define CYCLE 1         /* set to 1 for closed curve (start in all directions) */

#define NSTEPS 4500     /* number of frames of movie */
#define TIME 2000       /* time between movie frames, for fluidity of real-time simulation */ 
#define DPHI 0.00001    /* integration step */
#define NVID 500        /* number of iterations between images displayed on screen */
#define NCOLORS 20      /* number of colors */
#define COLORSHIFT 220  /* hue of initial color */ 
#define NSEG 100        /* number of segments of boundary */
#define LENGTH 0.05     /* length of velocity vectors */

/* Decreasing TIME accelerates the animation and the movie                               */
/* For constant speed of movie, TIME*DPHI should be kept constant                        */
/* However, increasing DPHI too much deterioriates quality of simulation                 */
/* NVID tells how often a picture is drawn in the animation, increase it for faster anim */
/* For a good quality movie, take for instance TIME = 400, DPHI = 0.00005, NVID = 100    */

#define PAUSE 10         /* number of frames after which to pause */
#define PSLEEP 1         /* sleep time during pause */
#define SLEEP1  1        /* initial sleeping time */
#define SLEEP2  100      /* final sleeping time */

#define PI 	3.141592654
#define DPI 	6.283185307
#define PID 	1.570796327

```

### 14 April 21 - Billiard in an ellipse ###

**Program:** `vid_particle_ellipse.c` (Old version of `particle_billiard.c`)

**Initial condition in function `animation()`:** ?

```
#define WINWIDTH 	1280  /* window width */
#define WINHEIGHT 	720   /* window height */

#define LAMBDA 1.2	/* aspect ratio of ellipse */
#define FOCI 1          /* set to 1 to draw focal points of ellipse */

#define MOVIE 1         /* set to 1 to generate movie */
#define RESAMPLE 0      /* set to 1 if particles should be added when dispersion too large */

#define XMIN -2.0
#define XMAX 2.0	/* x interval */
#define YMIN -1.125
#define YMAX 1.125	/* y interval for 9/16 aspect ratio */

#define NPART 500       /* number of particles */
#define NPARTMAX 50000	/* maximal number of particles after resampling */
#define LMAX 0.01       /* minimal segment length triggering resampling */ 
#define DMIN 0.02       /* minimal distance to boundary for triggering resampling */ 
#define CYCLE 1         /* set to 1 for closed curve (start in all directions) */

#define NSTEPS 4500     /* number of frames of movie */
#define TIME 2000       /* time between movie frames, for fluidity of real-time simulation */ 
#define DPHI 0.00001    /* integration step */
#define NVID 500        /* number of iterations between images displayed on screen */
#define NCOLORS 20      /* number of colors */
#define COLORSHIFT 220  /* hue of initial color */ 
#define NSEG 100        /* number of segments of boundary */
#define LENGTH 0.05     /* length of velocity vectors */

/* Decreasing TIME accelerates the animation and the movie                               */
/* For constant speed of movie, TIME*DPHI should be kept constant                        */
/* However, increasing DPHI too much deterioriates quality of simulation                 */
/* NVID tells how often a picture is drawn in the animation, increase it for faster anim */
/* For a good quality movie, take for instance TIME = 400, DPHI = 0.00005, NVID = 100    */

#define PAUSE 10         /* number of frames after which to pause */
#define PSLEEP 1         /* sleep time during pause */
#define SLEEP1  1        /* initial sleeping time */
#define SLEEP2  100      /* final sleeping time */

#define PI 	3.141592654
#define DPI 	6.283185307
#define PID 	1.570796327

```

### 14 April 21 - Billiard in a circle (very long and HD-ish) ###

**Program:** `xxx.c`

**Initial condition in function `animation()`:** `xxx`

```

```

### 11 April 21 - Drop in an elliptic pond ###

**Program:** `vid_drop_ellipse.c` (Old version of `drop_billiard.c`)

**Initial condition in function `animation()`:** `init_boundary_config(0.0, 0.0, 0.0, PI, configs);`

```
#define WINWIDTH 	1280  /* window width */
#define WINHEIGHT 	720   /* window height */

#define LAMBDA 1.5	/* aspect ratio of ellipse */
#define FOCI 1          /* set to 1 to draw focal points of ellipse */

#define MOVIE 0         /* set to 1 to generate movie */
#define RESAMPLE 1      /* set to 1 if particles should be added when dispersion too large */

#define XMIN -2.0
#define XMAX 2.0	/* x interval */
#define YMIN -1.125
#define YMAX 1.125	/* y interval for 9/16 aspect ratio */

#define NPART 2000	/* number of particles */
#define NPARTMAX 50000	/* maximal number of particles after resampling */
#define LMAX 0.01       /* minimal segment length triggering resampling */ 
#define DMIN 0.02       /* minimal distance to boundary for triggering resampling */ 
#define CYCLE 1         /* set to 1 for closed curve (start in all directions) */

#define NSTEPS 2000     /* number of frames of movie */
#define TIME 25         /* time between movie frames, for fluidity of real-time simulation */ 
#define DPHI 0.0004     /* integration step */
#define NVID 10         /* number of iterations between images displayed on screen */
#define NCOLORS 10      /* number of colors */
#define COLORSHIFT 0    /* hue of initial color */ 
#define NSEG 100        /* number of segments of boundary */

/* Decreasing TIME accelerates the animation and the movie               */
/* For constant speed of movie, TIME*DPHI should be kept constant        */
/* However, increasing DPHI too much deterioriates quality of simulation */
/* For a good quality movie, take for instance TIME = 50, DPHI = 0.0002  */

#define PAUSE 10         /* number of frames after which to pause */
#define PSLEEP 1         /* sleep time during pause */
#define SLEEP1  1        /* initial sleeping time */
#define SLEEP2  100      /* final sleeping time */

#define PI 	3.141592654
#define DPI 	6.283185307
#define PID 	1.570796327

```

### 11 April 21 - Drop in an elliptic pond ###

**Program:** `vid_drop_ellipse.c` (Old version of `drop_billiard.c`)

**Initial condition in function `animation()`:** `init_drop_config(0.8, 0.6, 0.0, DPI, configs);`

```
#define WINWIDTH 	1280  /* window width */
#define WINHEIGHT 	720   /* window height */

#define LAMBDA 1.5	/* aspect ratio of ellipse */

#define XMIN -2.0
#define XMAX 2.0	/* x interval */
#define YMIN -1.125
#define YMAX 1.125	/* y interval for 9/16 aspect ratio */

#define NPART 2000	/* number of particles */
#define NPARTMAX 50000	/* maximal number of particles after resampling */
#define LMAX 0.01       /* minimal segment length triggering resampling */ 
#define DMIN 0.02       /* minimal distance to boundary for triggering resampling */ 
#define CYCLE 1         /* set to 1 for closed curve (start in all directions) */

#define NSTEPS 2000     /* number of frames of movie */
#define TIME 25         /* time between movie frames, for fluidity of real-time simulation */ 
#define DPHI 0.0004     /* integration step */
#define NCOLORS 10      /* number of colors */
#define COLORSHIFT 0    /* hue of initial color */ 
#define NSEG 100        /* number of segments of boundary */

/* Decreasing TIME accelerates the animation and the movie               */
/* For constant speed of movie, TIME*DPHI should be kept constant        */
/* However, increasing DPHI too much deterioriates quality of simulation */


#define SLEEP1  1        /* initial sleeping time */
#define SLEEP2  100      /* final sleeping time */

#define PI 	3.141592654
#define DPI 	6.283185307
#define PID 	1.570796327

```

### 11 April 21 - Drop in an elliptic pond, starting from a focus ###

**Program:** `vid_drop_ellipse.c` (Old version of `drop_billiard.c`)

**Initial condition in function `animation()`:** `init_drop_config(sqrt(LAMBDA*LAMBDA-1.0),0.0, 0.0, DPI, configs);`

```
#define WINWIDTH 	1280  /* window width */
#define WINHEIGHT 	720   /* window height */

#define LAMBDA 1.5	/* aspect ratio of ellipse */

#define XMIN -2.0
#define XMAX 2.0	/* x interval */
#define YMIN -1.125
#define YMAX 1.125	/* y interval for 9/16 aspect ratio */

#define NPART 2000	/* number of particles */
#define NPARTMAX 50000	/* maximal number of particles after resampling */
#define LMAX 0.01       /* minimal segment length triggering resampling */ 
#define DMIN 0.02       /* minimal distance to boundary for triggering resampling */ 
#define CYCLE 1         /* set to 1 for closed curve (start in all directions) */

#define NSTEPS 2000     /* number of frames of movie */
#define TIME 25         /* time between movie frames, for fluidity of real-time simulation */ 
#define DPHI 0.0004     /* integration step */
#define NCOLORS 10      /* number of colors */
#define COLORSHIFT 0    /* hue of initial color */ 
#define NSEG 100        /* number of segments of boundary */

/* Decreasing TIME accelerates the animation and the movie               */
/* For constant speed of movie, TIME*DPHI should be kept constant        */
/* However, increasing DPHI too much deterioriates quality of simulation */


#define SLEEP1  1        /* initial sleeping time */
#define SLEEP2  100      /* final sleeping time */

#define PI 	3.141592654
#define DPI 	6.283185307
#define PID 	1.570796327

```

### 11 April 21 - Drop in an elliptic pond ###

**Program:** `vid_drop_ellipse.c` (Old version of `drop_billiard.c`)

**Initial condition in function `animation()`:** `init_drop_config(0.0, 0.0, 0.0, DPI, configs);`

```
#define WINWIDTH 	1280  /* window width */
#define WINHEIGHT 	720   /* window height */

#define LAMBDA 1.5	/* aspect ratio of ellipse */

#define XMIN -2.0
#define XMAX 2.0	/* x interval */
#define YMIN -1.125
#define YMAX 1.125	/* y interval for 9/16 aspect ratio */

#define NPART 2000	/* number of particles */
#define NPARTMAX 50000	/* maximal number of particles after resampling */
#define LMAX 0.01       /* minimal segment length triggering resampling */ 
#define DMIN 0.02       /* minimal distance to boundary for triggering resampling */ 
#define CYCLE 1         /* set to 1 for closed curve (start in all directions) */

#define NSTEPS 2000     /* number of frames of movie */
#define TIME 25         /* time between movie frames, for fluidity of real-time simulation */ 
#define DPHI 0.0004     /* integration step */
#define NCOLORS 10      /* number of colors */
#define COLORSHIFT 0    /* hue of initial color */ 
#define NSEG 100        /* number of segments of boundary */

/* Decreasing TIME accelerates the animation and the movie               */
/* For constant speed of movie, TIME*DPHI should be kept constant        */
/* However, increasing DPHI too much deterioriates quality of simulation */


#define SLEEP1  1        /* initial sleeping time */
#define SLEEP2  100      /* final sleeping time */

#define PI 	3.141592654
#define DPI 	6.283185307
#define PID 	1.570796327

```


### 9 April 21 - Drop in a circular water bowl (higher resolution) ###

**Program:** `vid_drop_circle.c` (ancestor of `drop_billiard.c`)

```
#define NPART 5000	/* number of particles */

#define LAMBDA 0.0	/* aspect ratio of ellipse */
#define MU 2.0	/* Larmor radius for billiard in magnetic field */

#define XMIN -2.0
#define XMAX 2.0		/* x interval */
#define YMIN -1.125
#define YMAX 1.125		/* y interval */

#define NSTEPS 3200
#define TIME 400
#define DPHI 0.00003
#define LPART 0.001
#define NCOLORS 15
#define COLORSHIFT 220

#define SLEEP1  1       	   /* initial sleeping time */
#define SLEEP2  100000		   /* final sleeping time */

```

### 9 April 21 - Droplet in a stadium-shaped billiard (higher res version) ###

**Program:** `vid_drop_stadium.c` (ancestor of `drop_billiard.c`)

**Initial condition in function `animation()`:** `init_drop_config(0.0,0.0);`

```
#define NPART 20000	/* number of particles */

#define LAMBDA 0.75	/* dimensions of stadium */

#define XMIN -2.0
#define XMAX 2.0		/* x interval */
#define YMIN -1.125
#define YMAX 1.125		/* y interval */

#define NSTEPS 1600
#define TIME 400
#define DPHI 0.00003
#define LPART 0.001
#define NCOLORS 10

#define SLEEP1  1       	   /* initial sleeping time */
#define SLEEP2  100000		   /* final sleeping time */

```

### 7 April 21 - Droplet in a stadium-shaped billiard ###

**Program:** `xxx.c`

**Initial condition in function `animation()`:** `xxx`

```

```

### 6 April 21 - Drop in a circular water bowl ###

**Program:** `xxx.c`

**Initial condition in function `animation()`:** `xxx`

```

```

### 6 April 21 - Drop in a circular water bowl ###

**Program:** `xxx.c`

**Initial condition in function `animation()`:** `xxx`

```

```

### 4 April 21 - Billiard in a circle (long version) ###

**Program:** `vid_circle_new.c` (Longer version of old 2012 code)

```
#define NPART 500		/* number of particles */

#define LAMBDA 0.0	/* aspect ratio of ellipse */
#define MU 2.0	/* Larmor radius for billiard in magnetic field */

#define XMIN -2.0
#define XMAX 2.0		/* x interval */
#define YMIN -1.125
#define YMAX 1.125		/* y interval */

#define NSTEPS 3200
#define TIME 400
#define DPHI 0.0001
#define LPART 0.05

#define SLEEP1  1       	   /* initial sleeping time */
#define SLEEP2  100000		   /* final sleeping time */

```


