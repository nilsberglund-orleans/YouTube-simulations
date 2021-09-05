### Parameter values for YouTube simulations ###

Created by **Nils Berglund** and optimized by **Marco Mancini**

C code for videos on YouTube Channel https://www.youtube.com/c/NilsBerglund

Below are parameter values used for different simulations, as well as initial conditions used in 
function animation. Some simulations use variants of the published code. The list is going to be 
updated gradually. 


### 5 September 21 -  ###

**Program:** `particle_pinball.c` 

**Initial condition in function `animation()`:** `init_drop_config(0.0, 0.0, -0.5*PID, 0.5*PID, configs); `

```
#define MOVIE 1         /* set to 1 to generate movie */

#define WINWIDTH 	1280  /* window width */
#define WINHEIGHT 	720   /* window height */

#define XMIN -2.0
#define XMAX 2.0	/* x interval */
#define YMIN -1.125
#define YMAX 1.125	/* y interval for 9/16 aspect ratio */

#define SCALING_FACTOR 1.0       /* scaling factor of drawing, needed for flower billiards, otherwise set to 1.0 */

/* Choice of the billiard table, see global_particles.c */

#define B_DOMAIN 21      /* choice of domain shape */

#define CIRCLE_PATTERN 21    /* pattern of circles */

#define ABSORBING_CIRCLES 0 /* set to 1 for circular scatterers to be absorbing */

#define NMAXCIRCLES 2000        /* total number of circles (must be at least NCX*NCY for square grid) */
#define NCX 50            /* number of circles in x direction */
#define NCY 26            /* number of circles in y direction */
#define NPOISSON 300        /* number of points for Poisson C_RAND_POISSON arrangement */

#define LAMBDA 1.68	/* parameter controlling shape of domain */
#define MU 0.036         /* second parameter controlling shape of billiard */
#define FOCI 1          /* set to 1 to draw focal points of ellipse */
#define NPOLY 6             /* number of sides of polygon */
#define APOLY 0.5           /* angle by which to turn polygon, in units of Pi/2 */ 
#define DRAW_BILLIARD 1     /* set to 1 to draw billiard */
#define DRAW_CONSTRUCTION_LINES 0   /* set to 1 to draw additional construction lines for billiard */
#define PERIODIC_BC 0       /* set to 1 to enforce periodic boundary conditions when drawing particles */

#define RESAMPLE 0      /* set to 1 if particles should be added when dispersion too large */
#define DEBUG 0         /* draw trajectories, for debugging purposes */

/* Simulation parameters */

#define NPART 1         /* number of particles */
#define NPARTMAX 100000	/* maximal number of particles after resampling */
#define LMAX 0.01       /* minimal segment length triggering resampling */ 
#define DMIN 0.02       /* minimal distance to boundary for triggering resampling */ 
#define CYCLE 1         /* set to 1 for closed curve (start in all directions) */
#define SHOWTRAILS 1    /* set to 1 to keep trails of the particles */
#define TEST_ACTIVE 0   /* set to 1 to test whether particle is in billiard */

#define NSTEPS 13000     /* number of frames of movie */
#define TIME 2500         /* time between movie frames, for fluidity of real-time simulation */ 
#define DPHI 0.00004     /* integration step */
#define NVID 150         /* number of iterations between images displayed on screen */

/* Decreasing TIME accelerates the animation and the movie                               */
/* For constant speed of movie, TIME*DPHI should be kept constant                        */
/* However, increasing DPHI too much deterioriates quality of simulation                 */
/* NVID tells how often a picture is drawn in the animation, increase it for faster anim */
/* For a good quality movie, take for instance TIME = 400, DPHI = 0.00005, NVID = 100    */

/* Colors and other graphical parameters */

#define NCOLORS 128      /* number of colors */
#define COLORSHIFT 0     /* hue of initial color */ 
#define RAINBOW_COLOR 1  /* set to 1 to use different colors for all particles */
#define SINGLE_COLOR 1   /* set to 1 to make all particles a single color */
#define FLOWER_COLOR 0   /* set to 1 to adapt initial colors to flower billiard (tracks vs core) */
#define NSEG 100         /* number of segments of boundary */
#define LENGTH 0.045      /* length of velocity vectors */
#define BILLIARD_WIDTH 2    /* width of billiard */
#define PARTICLE_WIDTH 2    /* width of particles */
#define FRONT_WIDTH 3       /* width of wave front */
#define COLOR_TRAJECTORY 90    /* hue for single color */

#define BLACK 1             /* set to 1 for black background */
#define COLOR_OUTSIDE 0     /* set to 1 for colored outside */ 
#define OUTER_COLOR 270.0   /* color outside billiard */
#define PAINT_INT 0         /* set to 1 to paint interior in other color (for polygon/Reuleaux) */
#define ERASE_OUTSIDE 1     /* set to 1 to erase outside of rectangular billiard (beta) */

#define PAUSE 1000       /* number of frames after which to pause */
#define PSLEEP 1         /* sleep time during pause */
#define SLEEP1  1        /* initial sleeping time */
#define SLEEP2  1000      /* final sleeping time */
#define END_FRAMES 100   /* number of still frames at end of movie */

```

### 4 September 21 - Energy evolution in a square grid wave protection ###

**Program:** `wave_energy.c` (variant of `wave_comparison.c`) 

**Initial condition in function `animation()`:** `int_planar_wave_comp(XMIN + 0.015, 0.0, phi, psi, xy_in);`

```
#define MOVIE 1         /* set to 1 to generate movie */

#define WINWIDTH 	1280  /* window width */
#define WINHEIGHT 	720   /* window height */

#define NX 1280          /* number of grid points on x axis */
#define NY 720          /* number of grid points on y axis */

#define XMIN -2.0
#define XMAX 2.0	/* x interval */
#define YMIN -1.125
#define YMAX 1.125	/* y interval for 9/16 aspect ratio */

#define JULIA_SCALE 1.0 /* scaling for Julia sets */

/* Choice of the billiard table */

#define B_DOMAIN 20      /* choice of domain shape, see list in global_pdes.c */
#define B_DOMAIN_B 20    /* choice of domain shape, see list in global_pdes.c */

#define CIRCLE_PATTERN 0       /* pattern of circles, see list in global_pdes.c */
#define CIRCLE_PATTERN_B 0     /* pattern of circles, see list in global_pdes.c */

#define P_PERCOL 0.25       /* probability of having a circle in C_RAND_PERCOL arrangement */
#define NPOISSON 300        /* number of points for Poisson C_RAND_POISSON arrangement */

#define LAMBDA 0.75	    /* parameter controlling the dimensions of domain */
#define MU 0.03 	    /* parameter controlling the dimensions of domain */
#define MUB 0.03	    /* parameter controlling the dimensions of domain */
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
#define OSCILLATE_LEFT 0    /* set to 1 to add oscilating boundary condition on the left */
#define OSCILLATE_TOPBOT 0  /* set to 1 to enforce a planar wave on top and bottom boundary */

#define OMEGA 0.0          /* frequency of periodic excitation */
#define AMPLITUDE 0.025      /* amplitude of periodic excitation */ 
#define COURANT 0.02       /* Courant number */
#define COURANTB 0.004    /* Courant number in medium B */
#define GAMMA 0.0      /* damping factor in wave equation */
#define GAMMAB 1.0e-8           /* damping factor in wave equation */
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

#define NSTEPS 2800      /* number of frames of movie */
#define NVID 25          /* number of iterations between images displayed on screen */
#define NSEG 100         /* number of segments of boundary */
#define INITIAL_TIME 200    /* time after which to start saving frames */
#define COMPUTE_ENERGIES 1  /* set to 1 to compute and print energies */
#define BOUNDARY_WIDTH 2    /* width of billiard boundary */

#define PAUSE 1000         /* number of frames after which to pause */
#define PSLEEP 1         /* sleep time during pause */
#define SLEEP1  1        /* initial sleeping time */
#define SLEEP2  1   /* final sleeping time */
#define END_FRAMES 100   /* number of still frames at end of movie */

/* Plot type, see list in global_pdes.c  */

#define PLOT 1

/* Color schemes */

#define BLACK 1          /* background */

#define COLOR_SCHEME 1   /* choice of color scheme, see list in global_pdes.c  */

#define SCALE 0          /* set to 1 to adjust color scheme to variance of field */
#define SLOPE 50.0        /* sensitivity of color on wave amplitude */
#define ATTENUATION 0.0  /* exponential attenuation coefficient of contrast with time */
#define E_SCALE 500.0     /* scaling factor for energy representation */

#define COLORHUE 260     /* initial hue of water color for scheme C_LUM */
#define COLORDRIFT 0.0   /* how much the color hue drifts during the whole simulation */
#define LUMMEAN 0.5      /* amplitude of luminosity variation for scheme C_LUM */
#define LUMAMP 0.3       /* amplitude of luminosity variation for scheme C_LUM */
#define HUEMEAN 220.0    /* mean value of hue for color scheme C_HUE */
#define HUEAMP -220.0      /* amplitude of variation of hue for color scheme C_HUE */

/* For debugging purposes only */
#define FLOOR 0         /* set to 1 to limit wave amplitude to VMAX */
#define VMAX 5.0       /* max value of wave amplitude */

```

### 3 September 21 - Illustrating the butterfly effect with a laser, long version without pauses ###

**Program:** `particle_pinball_angles.c` (variant of `particle_billiard.c`) 

**Initial condition in function `animation()`:** `angle0 = 0.2 + 2.4815e-6, x0 = -1.1, y0 = -0.5625;`

```
#define MOVIE 1         /* set to 1 to generate movie */

#define WINWIDTH 	1280  /* window width */
#define WINHEIGHT 	720   /* window height */

#define XMIN -2.0
#define XMAX 2.0	/* x interval */
#define YMIN -1.125
#define YMAX 1.125	/* y interval for 9/16 aspect ratio */

#define SCALING_FACTOR 1.0       /* scaling factor of drawing, needed for flower billiards, otherwise set to 1.0 */

/* Choice of the billiard table, see global_particles.c */

#define B_DOMAIN 20      /* choice of domain shape */

#define CIRCLE_PATTERN 3    /* pattern of circles */

#define ABSORBING_CIRCLES 0 /* set to 1 for circular scatterers to be absorbing */

#define NMAXCIRCLES 1000        /* total number of circles (must be at least NCX*NCY for square grid) */
#define NCX 3            /* number of circles in x direction */
#define NCY 8            /* number of circles in y direction */
#define NPOISSON 300        /* number of points for Poisson C_RAND_POISSON arrangement */

#define LAMBDA 0.5	/* parameter controlling shape of domain */
#define MU 0.11          /* second parameter controlling shape of billiard */
#define FOCI 1          /* set to 1 to draw focal points of ellipse */
#define NPOLY 6             /* number of sides of polygon */
#define APOLY 0.0           /* angle by which to turn polygon, in units of Pi/2 */ 
#define DRAW_BILLIARD 1     /* set to 1 to draw billiard */
#define DRAW_CONSTRUCTION_LINES 0   /* set to 1 to draw additional construction lines for billiard */
#define PERIODIC_BC 0       /* set to 1 to enforce periodic boundary conditions when drawing particles */

#define RESAMPLE 0      /* set to 1 if particles should be added when dispersion too large */
#define DEBUG 0         /* draw trajectories, for debugging purposes */

/* Simulation parameters */

#define NPART 1        /* number of particles */
#define NPARTMAX 100000	/* maximal number of particles after resampling */
#define LMAX 0.01       /* minimal segment length triggering resampling */ 
#define DMIN 0.02       /* minimal distance to boundary for triggering resampling */ 
#define CYCLE 1         /* set to 1 for closed curve (start in all directions) */
#define SHOWTRAILS 1    /* set to 1 to keep trails of the particles */

#define NSTEPS 6000     /* number of frames of movie */
#define TIME 4200        /* time between movie frames, for fluidity of real-time simulation */ 
#define DPHI 0.00001    /* integration step */
#define NVID 150         /* number of iterations between images displayed on screen */

/* Decreasing TIME accelerates the animation and the movie                               */
/* For constant speed of movie, TIME*DPHI should be kept constant                        */
/* However, increasing DPHI too much deterioriates quality of simulation                 */
/* NVID tells how often a picture is drawn in the animation, increase it for faster anim */
/* For a good quality movie, take for instance TIME = 400, DPHI = 0.00005, NVID = 100    */

/* Colors and other graphical parameters */

#define NCOLORS 128      /* number of colors */
#define COLORSHIFT 0     /* hue of initial color */ 
#define RAINBOW_COLOR 0  /* set to 1 to use different colors for all particles */
#define FLOWER_COLOR 0   /* set to 1 to adapt initial colors to flower billiard (tracks vs core) */
#define NSEG 100         /* number of segments of boundary */
#define LENGTH 0.02      /* length of velocity vectors */
#define BILLIARD_WIDTH 2    /* width of billiard */
#define PARTICLE_WIDTH 3    /* width of particles */
#define FRONT_WIDTH 3       /* width of wave front */

#define BLACK 1             /* set to 1 for black background */
#define COLOR_OUTSIDE 0     /* set to 1 for colored outside */ 
#define OUTER_COLOR 270.0   /* color outside billiard */
#define PAINT_INT 0         /* set to 1 to paint interior in other color (for polygon/Reuleaux) */
#define ERASE_OUTSIDE 0     /* set to 1 to erase outside of rectangular domain */

#define PAUSE 1000       /* number of frames after which to pause */
#define PSLEEP 1         /* sleep time during pause */
#define SLEEP1  1        /* initial sleeping time */
#define SLEEP2  1000      /* final sleeping time */
#define END_FRAMES 50   /* number of still frames at end of movie */

#define LUM_ATTENUATION 0.005  /* luminosity attenuation */

```

### 2 September 21 - An elliptical resonator ###

**Program:** `wave_billiard.c` 

**Initial condition in function `animation()`:** ``

```
#define MOVIE 1         /* set to 1 to generate movie */
#define DOUBLE_MOVIE 0  /* set to 1 to produce movies for wave height and energy simultaneously */

/* General geometrical parameters */

#define WINWIDTH 	720  /* window width */
#define WINHEIGHT 	720   /* window height */

#define NX 720          /* number of grid points on x axis */
#define NY 720          /* number of grid points on y axis */

#define XMIN -1.6
#define XMAX 1.6	/* x interval */
#define YMIN -1.6
#define YMAX 1.6	/* y interval for 9/16 aspect ratio */

#define JULIA_SCALE 1.0 /* scaling for Julia sets */

/* Choice of the billiard table */

#define B_DOMAIN 1      /* choice of domain shape, see list in global_pdes.c */

#define CIRCLE_PATTERN 8   /* pattern of circles, see list in global_pdes.c */

#define P_PERCOL 0.25       /* probability of having a circle in C_RAND_PERCOL arrangement */
#define NPOISSON 300        /* number of points for Poisson C_RAND_POISSON arrangement */

#define LAMBDA 1.5	    /* parameter controlling the dimensions of domain */
#define MU 1.25             /* parameter controlling the dimensions of domain */
#define NPOLY 3             /* number of sides of polygon */
#define APOLY 1.0           /* angle by which to turn polygon, in units of Pi/2 */ 
#define MDEPTH 4            /* depth of computation of Menger gasket */
#define MRATIO 3            /* ratio defining Menger gasket */
#define MANDELLEVEL 1000    /* iteration level for Mandelbrot set */
#define MANDELLIMIT 10.0    /* limit value for approximation of Mandelbrot set */
#define FOCI 1              /* set to 1 to draw focal points of ellipse */
#define NGRIDX 16           /* number of grid point for grid of disks */
#define NGRIDY 20           /* number of grid point for grid of disks */

/* You can add more billiard tables by adapting the functions */
/* xy_in_billiard and draw_billiard below */

/* Physical parameters of wave equation */

#define TWOSPEEDS 0         /* set to 1 to replace hardcore boundary by medium with different speed */
#define OSCILLATE_LEFT 0     /* set to 1 to add oscilating boundary condition on the left */
#define OSCILLATE_TOPBOT 0   /* set to 1 to enforce a planar wave on top and bottom boundary */

#define OMEGA 0.002        /* frequency of periodic excitation */
#define AMPLITUDE 1.0      /* amplitude of periodic excitation */ 
#define COURANT 0.02       /* Courant number */
#define COURANTB 0.01      /* Courant number in medium B */
#define GAMMA 0.0          /* damping factor in wave equation */
#define GAMMAB 1.0e-6           /* damping factor in wave equation */
#define GAMMA_SIDES 1.0e-4      /* damping factor on boundary */
#define GAMMA_TOPBOT 1.0e-7     /* damping factor on boundary */
#define KAPPA 0.0           /* "elasticity" term enforcing oscillations */
#define KAPPA_SIDES 5.0e-4  /* "elasticity" term on absorbing boundary */
#define KAPPA_TOPBOT 0.0    /* "elasticity" term on absorbing boundary */
/* The Courant number is given by c*DT/DX, where DT is the time step and DX the lattice spacing */
/* The physical damping coefficient is given by GAMMA/(DT)^2 */
/* Increasing COURANT speeds up the simulation, but decreases accuracy */
/* For similar wave forms, COURANT^2*GAMMA should be kept constant */

/* Boundary conditions, see list in global_pdes.c  */

#define B_COND 2

/* Parameters for length and speed of simulation */

#define NSTEPS 1400      /* number of frames of movie */
#define NVID 45          /* number of iterations between images displayed on screen */
#define NSEG 100         /* number of segments of boundary */
#define INITIAL_TIME 0      /* time after which to start saving frames */
#define BOUNDARY_WIDTH 2    /* width of billiard boundary */

#define PAUSE 1000       /* number of frames after which to pause */
#define PSLEEP 1         /* sleep time during pause */
#define SLEEP1  1        /* initial sleeping time */
#define SLEEP2  1        /* final sleeping time */
#define MID_FRAMES 20    /* number of still frames between parts of two-part movie */
#define END_FRAMES 100   /* number of still frames at end of movie */

/* Plot type, see list in global_pdes.c  */

#define PLOT 1

#define PLOT_B 0        /* plot type for second movie */

/* Color schemes */

#define BLACK 1          /* background */

#define COLOR_SCHEME 1   /* choice of color scheme, see list in global_pdes.c  */

#define SCALE 0          /* set to 1 to adjust color scheme to variance of field */
#define SLOPE 1.0        /* sensitivity of color on wave amplitude */
#define ATTENUATION 0.0  /* exponential attenuation coefficient of contrast with time */
#define E_SCALE 500.0     /* scaling factor for energy representation */

#define COLORHUE 260     /* initial hue of water color for scheme C_LUM */
#define COLORDRIFT 0.0   /* how much the color hue drifts during the whole simulation */
#define LUMMEAN 0.5      /* amplitude of luminosity variation for scheme C_LUM */
#define LUMAMP 0.3       /* amplitude of luminosity variation for scheme C_LUM */
#define HUEMEAN 220.0    /* mean value of hue for color scheme C_HUE */
#define HUEAMP -220.0      /* amplitude of variation of hue for color scheme C_HUE */

/* For debugging purposes only */
#define FLOOR 0         /* set to 1 to limit wave amplitude to VMAX */
#define VMAX 10.0       /* max value of wave amplitude */

```

### 1 September 21 - A parabolic resonating cavity ###

**Program:** `wave_billiard.c` 

**Initial condition in function `animation()`:** `init_wave(-LAMBDA, 0.0, phi, psi, xy_in);`

**In `init_wave()`:** `if ((xy_in[i][j])||(TWOSPEEDS)) phi[i][j] = 0.2*exp(-dist2/0.001)*cos(-sqrt(dist2)/0.01)`

```
#define MOVIE 1         /* set to 1 to generate movie */
#define DOUBLE_MOVIE 0  /* set to 1 to produce movies for wave height and energy simultaneously */

/* General geometrical parameters */

#define WINWIDTH 	1280  /* window width */
#define WINHEIGHT 	720   /* window height */

#define NX 1280          /* number of grid points on x axis */
#define NY 720          /* number of grid points on y axis */

#define XMIN -2.0
#define XMAX 2.0	/* x interval */
#define YMIN -1.125
#define YMAX 1.125	/* y interval for 9/16 aspect ratio */

#define JULIA_SCALE 1.0 /* scaling for Julia sets */

/* Choice of the billiard table */

#define B_DOMAIN 19      /* choice of domain shape, see list in global_pdes.c */

#define CIRCLE_PATTERN 8   /* pattern of circles, see list in global_pdes.c */

#define P_PERCOL 0.25       /* probability of having a circle in C_RAND_PERCOL arrangement */
#define NPOISSON 300        /* number of points for Poisson C_RAND_POISSON arrangement */

#define LAMBDA 0.0	    /* parameter controlling the dimensions of domain */
#define MU 1.25             /* parameter controlling the dimensions of domain */
#define NPOLY 3             /* number of sides of polygon */
#define APOLY 1.0           /* angle by which to turn polygon, in units of Pi/2 */ 
#define MDEPTH 4            /* depth of computation of Menger gasket */
#define MRATIO 3            /* ratio defining Menger gasket */
#define MANDELLEVEL 1000    /* iteration level for Mandelbrot set */
#define MANDELLIMIT 10.0    /* limit value for approximation of Mandelbrot set */
#define FOCI 1              /* set to 1 to draw focal points of ellipse */
#define NGRIDX 16           /* number of grid point for grid of disks */
#define NGRIDY 20           /* number of grid point for grid of disks */

/* You can add more billiard tables by adapting the functions */
/* xy_in_billiard and draw_billiard below */

/* Physical parameters of wave equation */

#define TWOSPEEDS 0         /* set to 1 to replace hardcore boundary by medium with different speed */
#define OSCILLATE_LEFT 0     /* set to 1 to add oscilating boundary condition on the left */
#define OSCILLATE_TOPBOT 0   /* set to 1 to enforce a planar wave on top and bottom boundary */

#define OMEGA 0.002        /* frequency of periodic excitation */
#define AMPLITUDE 1.0      /* amplitude of periodic excitation */ 
#define COURANT 0.02       /* Courant number */
#define COURANTB 0.01      /* Courant number in medium B */
#define GAMMA 0.0          /* damping factor in wave equation */
#define GAMMAB 1.0e-6           /* damping factor in wave equation */
#define GAMMA_SIDES 1.0e-4      /* damping factor on boundary */
#define GAMMA_TOPBOT 1.0e-7     /* damping factor on boundary */
#define KAPPA 0.0           /* "elasticity" term enforcing oscillations */
#define KAPPA_SIDES 5.0e-4  /* "elasticity" term on absorbing boundary */
#define KAPPA_TOPBOT 0.0    /* "elasticity" term on absorbing boundary */
/* The Courant number is given by c*DT/DX, where DT is the time step and DX the lattice spacing */
/* The physical damping coefficient is given by GAMMA/(DT)^2 */
/* Increasing COURANT speeds up the simulation, but decreases accuracy */
/* For similar wave forms, COURANT^2*GAMMA should be kept constant */

/* Boundary conditions, see list in global_pdes.c  */

#define B_COND 2

/* Parameters for length and speed of simulation */

#define NSTEPS 5000      /* number of frames of movie */
#define NVID 25          /* number of iterations between images displayed on screen */
#define NSEG 100         /* number of segments of boundary */
#define INITIAL_TIME 0      /* time after which to start saving frames */
#define BOUNDARY_WIDTH 2    /* width of billiard boundary */

#define PAUSE 1000       /* number of frames after which to pause */
#define PSLEEP 1         /* sleep time during pause */
#define SLEEP1  1        /* initial sleeping time */
#define SLEEP2  1        /* final sleeping time */
#define MID_FRAMES 20    /* number of still frames between parts of two-part movie */
#define END_FRAMES 100   /* number of still frames at end of movie */

/* Plot type, see list in global_pdes.c  */

#define PLOT 1

#define PLOT_B 0        /* plot type for second movie */

/* Color schemes */

#define BLACK 1          /* background */

#define COLOR_SCHEME 1   /* choice of color scheme, see list in global_pdes.c  */

#define SCALE 0          /* set to 1 to adjust color scheme to variance of field */
#define SLOPE 1.0        /* sensitivity of color on wave amplitude */
#define ATTENUATION 0.0  /* exponential attenuation coefficient of contrast with time */
#define E_SCALE 500.0     /* scaling factor for energy representation */

#define COLORHUE 260     /* initial hue of water color for scheme C_LUM */
#define COLORDRIFT 0.0   /* how much the color hue drifts during the whole simulation */
#define LUMMEAN 0.5      /* amplitude of luminosity variation for scheme C_LUM */
#define LUMAMP 0.3       /* amplitude of luminosity variation for scheme C_LUM */
#define HUEMEAN 220.0    /* mean value of hue for color scheme C_HUE */
#define HUEAMP -220.0      /* amplitude of variation of hue for color scheme C_HUE */

/* For debugging purposes only */
#define FLOOR 0         /* set to 1 to limit wave amplitude to VMAX */
#define VMAX 10.0       /* max value of wave amplitude */

```

