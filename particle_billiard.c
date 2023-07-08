/*********************************************************************************/
/*                                                                               */
/*  Animation of particles in billiard                                           */
/*                                                                               */
/*  N. Berglund, december 2012, april 2021                                       */
/*  UPDATE 14 April 21 : graphics files go to subfolder,                         */
/*  Switch MOVIE to decide whether to create a movie                             */
/*  UPDATE 3 May 21 : new domains                                                */
/*                                                                               */
/*  Feel free to reuse, but if doing so it would be nice to drop a               */
/*  line to nils.berglund@univ-orleans.fr - Thanks!                              */
/*                                                                               */
/*  compile with                                                                 */
/*  gcc -o particle_billiard particle_billiard.c                                 */
/*  -O3 -L/usr/X11R6/lib -ltiff -lm -lGL -lGLU -lX11 -lXmu -lglut                */
/*                                                                               */
/*  OMP acceleration may be more effective after executing, e.g.,                */
/*  export OMP_NUM_THREADS=2 in the shell before running the program             */
/*                                                                               */
/*  To make a video, set MOVIE to 1 and create subfolder tif_part                */
/*  It may be possible to increase parameter PAUSE                               */
/*                                                                               */
/*  create movie using                                                           */
/*  ffmpeg -i part.%05d.tif -vcodec libx264 part.mp4                             */
/*                                                                               */
/*********************************************************************************/
  
#include <math.h>
#include <string.h>
#include <GL/glut.h>
#include <GL/glu.h>
#include <unistd.h>
#include <sys/types.h>
#include <tiffio.h>     /* Sam Leffler's libtiff library. */
#include <omp.h>
#include <time.h>

#define MOVIE 1         /* set to 1 to generate movie */
#define SAVE_MEMORY 1       /* set to 1 to save memory when writing tiff images */
#define INVERT_COUNTER 0    /* set to 1 to save frames in inverse order */

// #define WINWIDTH 	1280  /* window width */
#define WINWIDTH 	720  /* window width */
#define WINHEIGHT 	720   /* window height */

// #define XMIN -1.5
// #define XMAX 2.5	/* x interval */
#define XMIN -1.125
#define XMAX 1.125	/* x interval */
#define YMIN -1.125
#define YMAX 1.125	/* y interval for 9/16 aspect ratio */

#define SCALING_FACTOR 1.0       /* scaling factor of drawing, needed for flower billiards, otherwise set to 1.0 */

/* Choice of the billiard table, see global_particles.c */

#define B_DOMAIN 30     /* choice of domain shape */

#define CIRCLE_PATTERN 1    /* pattern of circles */
#define POLYLINE_PATTERN 10  /* pattern of polyline */

#define ABSORBING_CIRCLES 0 /* set to 1 for circular scatterers to be absorbing */

#define NMAXCIRCLES 100000     /* total number of circles (must be at least NCX*NCY for square grid) */
#define NMAXPOLY 100000        /* total number of sides of polygonal line */   
#define NCX 30            /* number of circles in x direction */
#define NCY 20            /* number of circles in y direction */
#define NPOISSON 500        /* number of points for Poisson C_RAND_POISSON arrangement */
#define NGOLDENSPIRAL 2000  /* max number of points for C_GOLDEN_SPIRAL arrandement */
#define SDEPTH 1            /* Sierpinski gastket depth */

#define LAMBDA 1.5	/* parameter controlling shape of domain */
#define MU 0.005          /* second parameter controlling shape of billiard */
#define FOCI 1          /* set to 1 to draw focal points of ellipse */
#define NPOLY 6             /* number of sides of polygon */
#define APOLY 0.0           /* angle by which to turn polygon, in units of Pi/2 */ 
#define PENROSE_RATIO 2.5    /* parameter controlling the shape of small ellipses in Penrose room */

#define DRAW_BILLIARD 1     /* set to 1 to draw billiard */
#define DRAW_CONSTRUCTION_LINES 0   /* set to 1 to draw additional construction lines for billiard */
#define PERIODIC_BC 0       /* set to 1 to enforce periodic boundary conditions when drawing particles */

#define RESAMPLE 0      /* set to 1 if particles should be added when dispersion too large */
#define DEBUG 0         /* draw trajectories, for debugging purposes */

/* Simulation parameters */

// #define NPART 10      /* number of particles */
#define NPART 50000      /* number of particles */
#define NPARTMAX 100000	/* maximal number of particles after resampling */
#define LMAX 0.01       /* minimal segment length triggering resampling */ 
#define DMIN 0.02       /* minimal distance to boundary for triggering resampling */ 
#define CYCLE 1         /* set to 1 for closed curve (start in all directions) */
#define SHOWTRAILS 0    /* set to 1 to keep trails of the particles */
#define HEATMAP 1       /* set to 1 to show heat map of particles */
#define DRAW_FINAL_HEATMAP 1       /* set to 1 to show final heat map of particles */
#define DRAW_HEATMAP_HISTOGRAM 0   /* set to 1 to draw a histogram of particle distribution in heat map */
#define NBIN_FACTOR 6.0             /* constant controlling number of bins in histogram */
#define DRAW_HEATMAP_PARTICLES 1    /* set to 1 to draw particles in heat map */
#define HEATMAP_MAX_PART_BY_CELL 50     /* set to positive value to draw only limited number of particles in cell */
#define PLOT_HEATMAP_AVERAGE 1      /* set to 1 to plot average number of particles in heat map */
#define SHOWZOOM 0      /* set to 1 to show zoom on specific area */
#define PRINT_PARTICLE_NUMBER 0 /* set to 1 to print number of particles */
#define PRINT_LEFT_RIGHT_PARTICLE_NUMBER 0 /* set to 1 to print number of particles on left and right side */
#define PRINT_CIRCLE_PARTICLE_NUMBER 0 /* set to 1 to print number of particles outside circular maze */
#define PRINT_COLLISION_NUMBER 0 /* set to 1 to print number of collisions */
#define TEST_ACTIVE 1   /* set to 1 to test whether particle is in billiard */

#define TEST_INITIAL_COND 0     /* set to 1 to allow only initial conditions that pass a test */

#define NSTEPS 1300      /* number of frames of movie */
#define TIME 3000        /* time between movie frames, for fluidity of real-time simulation */ 
// #define DPHI 0.000002     /* integration step */
#define DPHI 0.00002     /* integration step */
#define NVID 25          /* number of iterations between images displayed on screen */
#define END_FRAMES 50    /* number of still frames at the end of the movie */

/* Decreasing TIME accelerates the animation and the movie                               */
/* For constant speed of movie, TIME*DPHI should be kept constant                        */
/* However, increasing DPHI too much deterioriates quality of simulation                 */
/* NVID tells how often a picture is drawn in the animation, increase it for faster anim */
/* For a good quality movie, take for instance TIME = 400, DPHI = 0.00005, NVID = 100    */

/* Colors and other graphical parameters */

#define COLOR_PALETTE 11     /* Color palette, see list in global_pdes.c  */

#define NCOLORS 500      /* number of colors */
#define COLORSHIFT 0     /* hue of initial color */ 
#define COLOR_HUEMIN 0   /* minimal color hue */
#define COLOR_HUEMAX 300 /* maximal color hue */
#define RAINBOW_COLOR 1  /* set to 1 to use different colors for all particles */
#define FLOWER_COLOR 0   /* set to 1 to adapt initial colors to flower billiard (tracks vs core) */
#define NSEG 100         /* number of segments of boundary */
#define LENGTH 0.025       /* length of velocity vectors */
#define BILLIARD_WIDTH 2    /* width of billiard */
#define PARTICLE_WIDTH 2    /* width of particles */
#define FRONT_WIDTH 3       /* width of wave front */

#define BLACK 1             /* set to 1 for black background */
#define COLOR_OUTSIDE 0     /* set to 1 for colored outside */ 
#define OUTER_COLOR 270.0   /* color outside billiard */
#define PAINT_INT 0         /* set to 1 to paint interior in other color (for polygon/Reuleaux) */
#define PAINT_EXT 1         /* set to 1 to paint exterior */

#define PAUSE 1000       /* number of frames after which to pause */
#define PSLEEP 2         /* sleep time during pause */
#define SLEEP1  1        /* initial sleeping time */
#define SLEEP2  1       /* final sleeping time */

#define NXMAZE 18       /* width of maze */
#define NYMAZE 18      /* height of maze */
#define MAZE_MAX_NGBH 8     /* max number of neighbours of maze cell */
#define RAND_SHIFT 15        /* seed of random number generator */
#define MAZE_XSHIFT 0.0     /* horizontal shift of maze */
#define MAZE_RANDOM_FACTOR 0.1     /* randomization factor for S_MAZE_RANDOM */
#define MAZE_CORNER_RADIUS 0.5     /* radius of tounded corners in maze */
#define CLOSE_MAZE 1        /* set to 1 to close maze exits */

#include "global_particles.c"
#include "sub_maze.c"
#include "sub_part_billiard.c"

int ncollisions = 0;

/*********************/
/* animation part    */
/*********************/

void init_boundary_config(double smin, double smax, double anglemin, double anglemax, double *configs[NPARTMAX])
/* initialize configuration: drop on the boundary, beta version */
/* WORKS FOR ELLIPSE, HAS TO BE ADAPTED TO GENERAL BILLIARD */
{
    int i;
    double ds, da, s, angle, theta, alpha, pos[2];
  
    if (anglemin <= 0.0) anglemin = PI/((double)NPART);
    if (anglemax >= PI) anglemax = PI*(1.0 - 1.0/((double)NPART));
    ds = (smax - smin)/((double)NPART);
    da = (anglemax - anglemin)/((double)NPART);
    for (i=0; i<NPART; i++) 
    {
        s = smin + ds*((double)i);
        angle = anglemin + da*((double)i),
        pos[0] = LAMBDA*cos(s);
        pos[1] = sin(s);
        theta = argument(-LAMBDA*pos[1], pos[0]/LAMBDA);
        alpha = theta + angle; 
        
        vbilliard_xy(configs[i], alpha, pos);
    }
}

void init_drop_config(double x0, double y0, double angle1, double angle2, double *configs[NPARTMAX])   
/* initialize configuration: drop at (x0,y0) */
{
    int i;
    double dalpha, alpha;
    double conf[2], pos[2];
  
    while (angle2 < angle1) angle2 += DPI;
    if (NPART > 1) dalpha = (angle2 - angle1)/((double)(NPART-1));
    else dalpha = 0.0;
    for (i=0; i<NPART; i++) 
    {
        alpha = angle1 + dalpha*((double)i);  
        
//         printf("alpha=%.5lg\n", alpha);
      
        pos[0] = x0;
        pos[1] = y0;
        vbilliard_xy(configs[i], alpha, pos);
    }
}
 
void init_partial_drop_config(double x0, double y0, double angle1, double angle2, int particle1, int particle2, 
                              int col, double *configs[NPARTMAX], int color[NPARTMAX], int newcolor[NPARTMAX])   
/* initialize configuration: drop at (x0,y0) for a range of particles */
{
    int i;
    double dalpha, alpha;
    double conf[2], pos[2];
  
    while (angle2 < angle1) angle2 += DPI;
    if (particle2 - particle1 > 1) dalpha = (angle2 - angle1)/((double)(particle2 - particle1-1));
    else dalpha = 0.0;
    for (i=particle1; i<particle2; i++) 
    {
        alpha = angle1 + dalpha*((double)i);  
        
//         printf("alpha=%.5lg\n", alpha);
      
        pos[0] = x0;
        pos[1] = y0;
        vbilliard_xy(configs[i], alpha, pos);
        color[i] = col;
        newcolor[i] = col;
    }
}
 
void init_sym_drop_config(double x0, double y0, double angle1, double angle2, double *configs[NPARTMAX])   
/* initialize configuration with two symmetric partial drops */
{
    int i;
    double dalpha, alpha, meanangle;
    double conf[2], pos[2];
  
    while (angle2 < angle1) angle2 += DPI;
    meanangle = 0.5*(angle1 + angle2);
    dalpha = (angle2 - angle1)/((double)(NPART-1));
    for (i=0; i<NPART/2; i++) 
    {
        alpha = meanangle + dalpha*((double)i);  
        pos[0] = x0;
        pos[1] = y0;
        vbilliard_xy(configs[i], alpha, pos);
    }
    for (i=0; i<NPART/2; i++)
    {
        alpha = meanangle - dalpha*((double)i);  
        pos[0] = x0;
        pos[1] = y0;
        vbilliard_xy(configs[NPART/2 + i], alpha, pos);
    }

}
 
void init_line_config(double x0, double y0, double x1, double y1, double angle, double *configs[NPARTMAX])   
/* initialize configuration: line (x0,y0)-(x1,y1) in direction alpha */
{
    int i;
    double dx, dy;
    double conf[2], pos[2];
  
    dx = (x1-x0)/((double)(NPART));
    dy = (y1-y0)/((double)(NPART));
//     dx = (x1-x0)/((double)(NPART-1));
//     dy = (y1-y0)/((double)(NPART-1));
    for (i=0; i<NPART; i++) 
    {
        pos[0] = x0 + ((double)i)*dx;
        pos[1] = y0 + ((double)i)*dy;
        vbilliard_xy(configs[i], angle, pos);
    }
}

void draw_zoom(int color[NPARTMAX], double *configs[NPARTMAX], int active[NPARTMAX], 
               double x_target, double y_target, double width, double shiftx, double shifty, double zoomwidth, int shooter)
/* draw zoom around target (for laser in room of mirrors) */
{
    int i;
    double x1, y1, x2, y2, xb, yb, cosphi, sinphi, rgb[3], tradius, phi;
//     double x1, y1, x2, y2, xb, yb, cosphi, sinphi, rgb[3], shiftx = 0.0, shifty = 0.65, tradius, phi, zoomwidth = 0.4;
    
//     shiftx = 1.65;
//     shifty = 0.75;
//     zoomwidth = 0.3;
    
    glEnable(GL_LINE_SMOOTH);
    glColor3f(1.0, 1.0, 1.0);
    
    /* draw zoomed area */
    glLineWidth(BILLIARD_WIDTH/2);
    
    x1 = x_target - width;
    y1 = y_target - width;
    x2 = x_target + width;
    y2 = y_target + width;
    glBegin(GL_LINE_LOOP);    
    glVertex2d(x1, y1);
    glVertex2d(x2, y1);
    glVertex2d(x2, y2);
    glVertex2d(x1, y2);
    glEnd();

    /* draw zoom boundary */
    glLineWidth(BILLIARD_WIDTH*2);
    
    x1 = shiftx - zoomwidth;
    y1 = shifty - zoomwidth;
    x2 = shiftx + zoomwidth;
    y2 = shifty + zoomwidth;
    glBegin(GL_LINE_LOOP);    
    glVertex2d(x1, y1);
    glVertex2d(x2, y1);
    glVertex2d(x2, y2);
    glVertex2d(x1, y2);
    glEnd();
    
    /* draw billiard boundaries in zoom */
    glLineWidth(BILLIARD_WIDTH*2);
    
    if (y_target + width > 1.0)
    {
        yb = shifty + 0.5*(1.0 - y_target)/width;
        glBegin(GL_LINE_STRIP);
        glVertex2d(x1, yb);
        glVertex2d(x2, yb);
        glVertex2d(x2, yb + 0.02);
        glVertex2d(x1, yb + 0.02);
        glEnd();
    }
    /* other boundaries not yet implemented */
    
    /* draw target in zoom */
    glLineWidth(BILLIARD_WIDTH*2);
    
    if (shooter) glColor3f(1.0, 0.0, 0.0);
    else glColor3f(0.0, 0.8, 0.2);
    tradius = zoomwidth*MU/width;
    draw_circle(shiftx, shifty, tradius, NSEG);
    
//     glLineWidth(PARTICLE_WIDTH*2);
    
    for (i=0; i<nparticles; i++)
    {
        cosphi = (configs[i][6] - configs[i][4])/configs[i][3];
        sinphi = (configs[i][7] - configs[i][5])/configs[i][3];
        x1 = (configs[i][4] + configs[i][2]*cosphi - x_target)/width;
        y1 = (configs[i][5] + configs[i][2]*sinphi - y_target)/width;
        x2 = (configs[i][4] + (configs[i][2] + LENGTH)*cosphi - x_target)/width;
        y2 = (configs[i][5] + (configs[i][2] + LENGTH)*sinphi - y_target)/width;
        
        /* adjusting segments that are partly in the domain */
        if ((vabs(x1) < 1.0)&&(vabs(x2) > 1.0)) 
        {
            if (x1 > 0.0) xb = 1.0;
            else xb = -1.0;
            y2 = y1 + (xb - x1)*(y2 - y1)/(x2 - x1);
            x2 = xb;
        }
        else 
        if ((vabs(x1) > 1.0)&&(vabs(x2) < 1.0)) 
        {
            if (x2 > 0.0) xb = 1.0;
            else xb = -1.0;
            y1 = y2 + (xb - x2)*(y1 - y2)/(x1 - x2);
            x1 = xb;
        }
        if ((vabs(y1) < 1.0)&&(vabs(y2) > 1.0)) 
        {
            if (y1 > 0.0) yb = 1.0;
            else yb = -1.0;
            x2 = x1 + (yb - y1)*(x2 - x1)/(y2 - y1);
            y2 = yb;
        }
        else 
        if ((vabs(y1) > 1.0)&&(vabs(y2) < 1.0)) 
        {
            if (y2 > 0.0) yb = 1.0;
            else yb = -1.0;
            x1 = x2 + (yb - y2)*(x1 - x2)/(y1 - y2);
            y1 = yb;
        }
        
//         if ((active[i])&&(vabs(x1) < 1.0)&&(vabs(y1) < 1.0)&&(vabs(x2) < 1.0)&&(vabs(y2) < 1.0))  
        if (((active[i])&&(vabs(x1) < 1.0)&&(vabs(y1) < 1.0))||((vabs(x2) < 1.0)&&(vabs(y2) < 1.0))) 
        {
            rgb_color_scheme_minmax(color[i], rgb);
            glColor3f(rgb[0], rgb[1], rgb[2]);
        
            glBegin(GL_LINE_STRIP);
            glVertex2d(shiftx + zoomwidth*SCALING_FACTOR*x1, shifty + zoomwidth*SCALING_FACTOR*y1);
            glVertex2d(shiftx + zoomwidth*SCALING_FACTOR*x2, shifty + zoomwidth*SCALING_FACTOR*y2);
            glEnd ();
        }
        
    }
}


void draw_config_showtrails(int color[NPARTMAX], double *configs[NPARTMAX], int active[NPARTMAX])
/* draw the particles */
{
    int i;
    double x0, y0, x1, y1, x2, y2, cosphi, sinphi, rgb[3], len;

    glutSwapBuffers(); 
    if (PAINT_INT) paint_billiard_interior();
      
    glLineWidth(PARTICLE_WIDTH);
    
    glEnable(GL_LINE_SMOOTH);

    for (i=0; i<nparticles; i++) 
//         if (active[i])
    {
//         if (configs[i][2]<0.0) 
//         {    
//             vbilliard(configs[i]);
//             if (!RAINBOW_COLOR)
//             {
//                 color[i]++;
//                 if (color[i] >= NCOLORS) color[i] -= NCOLORS;
//             }
//         }
            
        configs[i][2] += DPHI; 
        
        cosphi = (configs[i][6] - configs[i][4])/configs[i][3];
        sinphi = (configs[i][7] - configs[i][5])/configs[i][3];
        len = configs[i][2] + LENGTH;
        if (len > configs[i][3]) len = configs[i][3];
        
        x0 = configs[i][4];
        y0 = configs[i][5]; 
        x1 = configs[i][4] + configs[i][2]*cosphi;
        y1 = configs[i][5] + configs[i][2]*sinphi;
        x2 = configs[i][4] + len*cosphi;
        y2 = configs[i][5] + len*sinphi;
        
        /* test whether particle does not escape billiard */
        if ((TEST_ACTIVE)&&(active[i])) active[i] = xy_in_billiard(x1, y1);
        
        if (active[i])  
        {
            rgb_color_scheme_minmax(color[i], rgb);
            glColor3f(rgb[0], rgb[1], rgb[2]);
        
            glBegin(GL_LINE_STRIP);
            glVertex2d(SCALING_FACTOR*x0, SCALING_FACTOR*y0);
            glVertex2d(SCALING_FACTOR*x2, SCALING_FACTOR*y2);
            glEnd ();
        }

//         if (configs[i][2] > configs[i][3] - DPHI)
//         {
//             glBegin(GL_LINE_STRIP);
//             glVertex2d(SCALING_FACTOR*x0, SCALING_FACTOR*y0);
//             glVertex2d(SCALING_FACTOR*configs[i][6], SCALING_FACTOR*configs[i][7]);
//             glEnd ();
//         }
    }
    if (DRAW_BILLIARD) draw_billiard();    
    
    if (SHOWZOOM) switch (POLYLINE_PATTERN) {
        case (P_TOKA_PRIME):
        {
            draw_zoom(color, configs, active, x_target, y_target, 0.1, 1.65, 0.75, 0.3, 0);
            draw_zoom(color, configs, active, x_shooter, y_shooter, 0.1, -1.65, 0.75, 0.3, 1);
            break;
        }
        case (P_TOKA_NONSELF):
        {
            draw_zoom(color, configs, active, 0.0, 0.0, 0.1, 1.65, 0.75, 0.3, 0);
            break;
        }
    }
//     if (SHOWZOOM) draw_zoom(color, configs, active, 0.95, 0.0, 0.1);
}


void draw_config_heatmap(double *configs[NPARTMAX], int active[NPARTMAX], int heatmap_number[NXMAZE*NYMAZE+1], int heatmap_total[NXMAZE*NYMAZE+1], short int heatmap_visited[NXMAZE*NYMAZE+1], int draw_particles)
/* draw a heat map of particle distribution (for mazes) */
{
    int i, j, n, part_number;
    double x, y, cosphi, sinphi, rgb[3], len, padding = 0.02;
    double *xtable, *ytable;
    short int *drawtable;
    static int time, first = 1;
    static double minprop;
    
    if (first)
    {
        time = 0;
        first = 0;
        minprop = 0.01;
//         if (PLOT_HEATMAP_AVERAGE) minprop = 0.005;
//         else minprop = 0.01;
    }
    time++;
    
    drawtable = malloc(sizeof(short int)*(NPARTMAX));
    xtable = malloc(sizeof(double)*(NPARTMAX));
    ytable = malloc(sizeof(double)*(NPARTMAX));
    
    glutSwapBuffers(); 
    blank();
    if (PAINT_INT) paint_billiard_interior();
    
    for (i=0; i<NXMAZE*NYMAZE+1; i++) heatmap_number[i] = 0;
    
    for (i=0; i<nparticles; i++) 
//         if ((active[i])&&(configs[i][0] < DUMMY_ABSORBING))
    {
//         configs[i][2] += DPHI; 
        
        cosphi = (configs[i][6] - configs[i][4])/configs[i][3];
        sinphi = (configs[i][7] - configs[i][5])/configs[i][3];
        len = configs[i][2];
        if (len > configs[i][3] - padding) len = configs[i][3] - padding;
        if (len < 1.0e-10) len = 1.0e-10;
//         if (len < 0.0) len = 1.0e-10;
        
        x = configs[i][4] + len*cosphi;
        y = configs[i][5] + len*sinphi;
        
        xtable[i] = x;
        ytable[i] = y;
        
        /* test whether particle does not escape billiard */
        if ((TEST_ACTIVE)&&(active[i])) active[i] = xy_in_billiard(x, y);
        
        if (active[i])  
        {
            n = find_maze_cell(x, y);
        
            if ((n > -1)&&(n < NXMAZE*NYMAZE+1)) 
            {
                heatmap_number[n]++;
                heatmap_total[n]++;
                heatmap_visited[n] = 1;
            }
            
            if (HEATMAP_MAX_PART_BY_CELL > 0)
            {
                drawtable[i] = ((n == -2)||((n >= -1)&&(heatmap_number[n] <= HEATMAP_MAX_PART_BY_CELL)));
            }
            else drawtable[i] = 1;
//             else drawtable[i] = (n >= -1);
        }
        
//         printf("Particle %i is in maze cell %i\n", i, n);
    }
    
    for (n=0; n<NXMAZE*NYMAZE+1; n++) 
    {
        if (PLOT_HEATMAP_AVERAGE) 
        {
            if (heatmap_total[n] == 0) part_number = 0;
            else part_number = 1 + (int)((double)heatmap_total[n]/(double)time);
        }
        else part_number = heatmap_number[n];
        
        if ((part_number == 0)&&(heatmap_visited[n])) part_number = -1;
        draw_maze_cell(n, part_number, minprop);
//         if (part_number != 0) printf("%i particles in maze cell %i\n", part_number, n);
    }
    
    glColor3f(1.0, 1.0, 1.0);
    if (draw_particles) 
        for (i=0; i<nparticles; i++) 
            if ((active[i])&&(drawtable[i]))
                draw_circle(SCALING_FACTOR*xtable[i], SCALING_FACTOR*ytable[i], 0.001, 6);
    
    if (DRAW_BILLIARD) draw_billiard(); 
    
    free(xtable);
    free(ytable);
    free(drawtable);
}

double plot_coord(double x, double xmin, double xmax)
{
    return(xmin + x*(xmax - xmin));
}

void draw_chosen_heatmap_histogram(int heatmap_number[NXMAZE*NYMAZE+1], int normalisation)
{
    int i, j, k, n, bin, nbins, binwidth, maxpart, part_number, maxhist = 0, *histo;
    double xmin, xmax, ymin, ymax, xmid, ymid, dx, dy, plotxmin, plotxmax, plotymin, plotymax, c; 
    double x1, x2, y, y1, y2, hue, rgb[3];
    static int first = 1, prevmaxhist, prevmaxpart;
    static double minprop = 0.01;
    char message[100];
    
    if (first)
    {
        xmin = XMAX - 1.25;
        xmax = XMAX - 0.05;
        ymin = YMAX - 1.25;
        ymax = YMAX - 0.05;
                
        xmid = 0.5*(xmin + xmax);
        ymid = 0.5*(ymin + ymax);
        
        dx = 0.5*(xmax - xmin);
        dy = 0.5*(ymax - ymin);
        
        plotxmin = xmin + 0.15;
        plotxmax = xmax - 0.1;
        plotymin = ymin + 0.07;
        plotymax = ymax - 0.1;
        
        prevmaxhist = 1010;
        prevmaxpart = 10;
        
        first = 0;
    }
    
//     nbins = NXMAZE;
//     nbins = (int)(2.0*sqrt((double)NPART/(double)NXMAZE));
    nbins = (int)(NBIN_FACTOR*sqrt((double)NPART/(double)(NXMAZE*NYMAZE)));
    
    histo = (int *)malloc(nbins*sizeof(int));
    
    for (i=0; i<nbins; i++) histo[i] = 0;
    
    maxpart = 0;
    for (j=0; j<NXMAZE*NYMAZE+1; j++) 
        if (heatmap_number[j]/normalisation > maxpart) maxpart = heatmap_number[j]/normalisation;
    
    if (maxpart < 1010) maxpart = 1010;
        
    c = log((double)maxpart)/(double)nbins;
    
    for (j=0; j<NXMAZE*NYMAZE+1; j++) 
    {
        n = heatmap_number[j]/normalisation;
        if (n > 0)
        {
            i = (int)(log((double)n)/c) - 2;
            if (i < 0) i = 0;
            histo[i]++;
        }
    }
    
//     for (i=0; i<nbins; i++) printf("%i cells in bin %i\n", histo[i], i);
    
    for (i=0; i<nbins; i++) if (histo[i] > maxhist) maxhist = histo[i];
    if (maxhist < 10) maxhist = 10;
    
    /* some smoothing in time */
    maxpart = (maxpart + 3*prevmaxpart)/4;
    prevmaxpart = maxpart;
    maxhist = (maxhist + 3*prevmaxhist)/4;
    prevmaxhist = maxhist;
    
    glColor3f(0.0, 0.0, 0.0);
    glLineWidth(1);
    
    x1 = plotxmin;
    y1 = plotymin;
    
    for (i=0; i<nbins; i++)
    {
        x2 = plot_coord((double)i/(double)nbins, plotxmin, plotxmax);
        y = log((double)(histo[i]+1))/log((double)(maxhist+1));
//         y = (double)(histo[i]+1)/(double)(maxhist+1);
        y2 = plot_coord(y, plotymin, plotymax);
        
        part_number = (int)(exp(c*(double)i));
        rgb_color_scheme_density(part_number, rgb, minprop);
//         printf("Bin %i, part nb %i, cell nb %i\n", i, part_number, histo[i]);
        draw_colored_rectangle_rgb(x1, y1, x2, y2, rgb);
        x1 = x2;
    }
    
    glColor3f(1.0, 1.0, 1.0);
    glLineWidth(2);
    
    draw_line(plotxmin, plotymin, plotxmax + 0.05, plotymin);
    draw_line(plotxmin, plotymin, plotxmin, plotymax + 0.1);
            
    /* graduation of x axis */
    for (j=1; j < (int)(log((double)maxpart)/log(10.0)) + 1; j++)
    {
        n = (int)ipow(10.0, j);
        
        i = (int)(log((double)n)/c) - 2;
        x2 = plot_coord((double)i/(double)nbins, plotxmin, plotxmax);
        
        draw_line(x2, plotymin - 0.02, x2, plotymin + 0.02);
        
        if (n <= 1000) sprintf(message, "%i", n); 
        else sprintf(message, "1e%i", j);
        write_text_fixedwidth(x2 - 0.015 - 0.01*(double)j, plotymin - 0.08, message);  
    }
    
    sprintf(message, "Particles");
    write_text_fixedwidth(plotxmax - 0.2, plotymin - 0.15, message); 
    
    /* graduation of y axis */
    for (j=0; j < (int)(log((double)maxhist)/log(10.0)) + 1; j++) for (k=1; k<10; k++)
    {
        if (j==0) n = k;
        else n = k*(int)ipow(10.0, j);
        y = log((double)(n+1))/log((double)(maxhist+1));
        y2 = plot_coord(y, plotymin, plotymax);
        
        if (y < plotymax) draw_line(plotxmin - 0.02, y2, plotxmin + 0.02, y2);
        
        if (((k < 3)||(k == 5))&&(y < plotymax))
        {
            if (n <= 1000) sprintf(message, "%4d", n); 
            else sprintf(message, "1e%i", j);
            write_text_fixedwidth(plotxmin - 0.18, y2 - 0.015, message);  
        }
    }
    
    sprintf(message, "Cells");
    write_text_fixedwidth(plotxmin + 0.05, plotymax + 0.05, message); 
    
    free(histo); 
}

void draw_heatmap_histogram(int heatmap_number[NXMAZE*NYMAZE+1], int heatmap_total[NXMAZE*NYMAZE+1], int time)
{
    if (PLOT_HEATMAP_AVERAGE) draw_chosen_heatmap_histogram(heatmap_total, time+1);
    else draw_chosen_heatmap_histogram(heatmap_number, 1);
}

void draw_config(int color[NPARTMAX], double *configs[NPARTMAX], int active[NPARTMAX])
/* draw the particles */
{
    int i, c;
    double x1, y1, x2, y2, cosphi, sinphi, rgb[3];

    glutSwapBuffers(); 
    if (!SHOWTRAILS) blank();
//     if (!((SHOWTRAILS)||(HEATMAP))) blank();
    if (PAINT_INT) paint_billiard_interior();
      
    glLineWidth(PARTICLE_WIDTH);
    
    glEnable(GL_LINE_SMOOTH);

    for (i=0; i<nparticles; i++) 
//         if (active[i])
    {
//         if (configs[i][2]<0.0) 
//         {    
//             c = vbilliard(configs[i]);
//             if (!RAINBOW_COLOR)
//             {
//                 color[i]++;
//                 if (color[i] >= NCOLORS) color[i] -= NCOLORS;
//             }
//             if ((ABSORBING_CIRCLES)&&(c < 0)) active[i] = 0;
//         }

//         configs[i][2] += DPHI; 
        
        cosphi = (configs[i][6] - configs[i][4])/configs[i][3];
        sinphi = (configs[i][7] - configs[i][5])/configs[i][3];
        x1 = configs[i][4] + configs[i][2]*cosphi;
        y1 = configs[i][5] + configs[i][2]*sinphi;
        x2 = configs[i][4] + (configs[i][2] + LENGTH)*cosphi;
        y2 = configs[i][5] + (configs[i][2] + LENGTH)*sinphi;
        
        /* test whether particle does not escape billiard */
        if ((TEST_ACTIVE)&&(active[i])) active[i] = xy_in_billiard(x1, y1);
        
        if (active[i])  
        {
            rgb_color_scheme_minmax(color[i], rgb);
            glColor3f(rgb[0], rgb[1], rgb[2]);
        
            glBegin(GL_LINE_STRIP);
            glVertex2d(SCALING_FACTOR*x1, SCALING_FACTOR*y1);
            glVertex2d(SCALING_FACTOR*x2, SCALING_FACTOR*y2);
            glEnd ();
        
            /* taking care of boundary conditions - only needed for periodic boundary conditions */
            if (PERIODIC_BC)
            {
                if (SCALING_FACTOR*x2 > XMAX)
                {
                    glBegin(GL_LINE_STRIP);
                    glVertex2d(SCALING_FACTOR*(x1+XMIN-XMAX), SCALING_FACTOR*y1);
                    glVertex2d(SCALING_FACTOR*(x2+XMIN-XMAX), SCALING_FACTOR*y2);
                    glEnd ();
                }
        
                if (SCALING_FACTOR*x2 < XMIN)
                {
                    glBegin(GL_LINE_STRIP);
                    glVertex2d(SCALING_FACTOR*(x1-XMIN+XMAX), SCALING_FACTOR*y1);
                    glVertex2d(SCALING_FACTOR*(x2-XMIN+XMAX), SCALING_FACTOR*y2);
                    glEnd ();
                }

                if (SCALING_FACTOR*y2 > YMAX)
                {
                    glBegin(GL_LINE_STRIP);
                    glVertex2d(SCALING_FACTOR*x1, SCALING_FACTOR*(y1+YMIN-YMAX));
                    glVertex2d(SCALING_FACTOR*x2, SCALING_FACTOR*(y2+YMIN-YMAX));
                    glEnd ();
                }

                if (SCALING_FACTOR*y2 < YMIN)
                {
                    glBegin(GL_LINE_STRIP);
                    glVertex2d(SCALING_FACTOR*x1, SCALING_FACTOR*(y1+YMAX-YMIN));
                    glVertex2d(SCALING_FACTOR*x2, SCALING_FACTOR*(y2+YMAX-YMIN));
                    glEnd ();
                }
            }
        }
        
        /* draw trajectories, for debugging purpose */
        if (DEBUG)
        {
            glLineWidth(1.0);
            glBegin(GL_LINES);
            glVertex2d(SCALING_FACTOR*configs[i][4], SCALING_FACTOR*configs[i][5]);
            glVertex2d(SCALING_FACTOR*configs[i][6], SCALING_FACTOR*configs[i][7]);
            glEnd ();
            glLineWidth(3.0);
        }
    
//         if (configs[i][2] > configs[i][3] - DPHI) configs[i][2] -= configs[i][3];
    }
    if (DRAW_BILLIARD) draw_billiard();    
    
    if (SHOWZOOM) switch (POLYLINE_PATTERN) {
        case (P_TOKA_PRIME):
        {
            draw_zoom(color, configs, active, x_target, y_target, 0.1, 1.65, 0.75, 0.3, 0);
            draw_zoom(color, configs, active, x_shooter, y_shooter, 0.1, -1.65, 0.75, 0.3, 1);
            break;
        }
        case (P_TOKA_NONSELF):
        {
            draw_zoom(color, configs, active, 0.0, 0.0, 0.1, 0.82, 0.82, 0.25, 0);
            break;
        }
    }
//     if (SHOWZOOM) draw_zoom(color, configs, active, 0.95, 0.0, 0.1);
}


void graph_movie(int time, int color[NPARTMAX], double *configs[NPARTMAX], int active[NPARTMAX])
/* compute next movie frame */
{
    int i, j, c;

    for (j=0; j<time; j++)
    {
        #pragma omp parallel for private(i,c)
        for (i=0; i<nparticles; i++) if (active[i])
        {      
            if (configs[i][2]<0.0) 
            {    
//                 printf("reflecting particle %i\n", i);
                c = vbilliard(configs[i]);
                
                if ((ABSORBING_CIRCLES)&&(c < 0)) active[i] = 0;
                if (c < 0) active[i] = 0;
                else ncollisions++;
//                 if (c>=0) color[i]++;
                if ((!RAINBOW_COLOR)&&(c>=0)) color[i]++;
                if (!RAINBOW_COLOR)
                {
                    color[i]++;
                    if (color[i] >= NCOLORS) color[i] -= NCOLORS;
                }                
            }

            configs[i][2] += DPHI; 
        
            if (configs[i][2] > configs[i][3] - DPHI)
            {
                configs[i][2] -= configs[i][3];
            }
        }
    }
    
//     draw_config(color, configs);
}

void print_part_number(double *configs[NPARTMAX], int active[NPARTMAX], double x, double y)
{
    char message[50];
    int i, n_active_particles = 0;
    double rgb[3];
    
    /* count active particles, using the fact that absorbed particles have been given dummy coordinates */
    for (i=0; i<nparticles; i++)
        if (active[i]) n_active_particles++;
//         if (configs[i][0] != -10.0) n_active_particles++;
        
    hsl_to_rgb(0.0, 0.0, 0.0, rgb);
    erase_area(x, y, 0.5, 0.1, rgb);
    
    glColor3f(1.0, 1.0, 1.0);
    sprintf(message, "%i particles", n_active_particles);
    write_text(x, y, message);
    
}

void print_left_right_part_number(double *configs[NPARTMAX], int active[NPARTMAX], double xl, double yl, double xr, double yr, t_exit exits[NPARTMAX])
{
    char message[50];
    int i, nleft = 0, nright = 0;
    double rgb[3], x1, y1, cosphi, sinphi;
    static int first = 1;
    static double xmin, xmax;
    
    if (first)
    {
        compute_maze_boundaries(POLYLINE_PATTERN, &xmin, &xmax);
        first = 0;
    }
    
    /* count active particles, using the fact that absorbed particles have been given dummy coordinates */
    for (i=0; i<nparticles; i++) if ((configs[i][0] == DUMMY_ABSORBING))
    {
        cosphi = (configs[i][6] - configs[i][4])/configs[i][3];
        sinphi = (configs[i][7] - configs[i][5])/configs[i][3];
        x1 = configs[i][4] + configs[i][2]*cosphi;
        y1 = configs[i][5] + configs[i][2]*cosphi;
        if ((x1 < xmin)&&(x1 > XMIN)&&(y1 < YMAX)&&(y1 > YMIN)) exits[i].left = 1;
        else if ((x1 > xmax)&&(x1 < XMAX)&&(y1 < YMAX)&&(y1 > YMIN)) 
        {
            exits[i].right = 1;
            printf("Detected leaving particle %i at (%.2f, %2f)\n\n\n", i, x1, y1);
        }
    }
    
    for (i=0; i<nparticles; i++) 
    {
        if (exits[i].left) nleft++;
        if (exits[i].right) nright++;
//         printf("particle[%i]: left = %i, right = %i\n", i, exits[i].left, exits[i].right);
    }
    
    hsl_to_rgb(0.0, 0.0, 0.0, rgb);
    
    erase_area(xl, yl - 0.03, 0.3, 0.12, rgb);
//     erase_area(xl, yl - 0.03, 0.25, 0.12, rgb);
    erase_area(xr, yr - 0.03, 0.35, 0.12, rgb);
    
    glColor3f(1.0, 1.0, 1.0);
//     if (nleft > 1) sprintf(message, "%i particles", nleft);
//     else sprintf(message, "%i particle", nleft);
    if (nleft > 1) sprintf(message, "%i part.", nleft);
    else sprintf(message, "%i part.", nleft);
    write_text_fixedwidth(xl, yl, message);
    if (nright > 1) sprintf(message, "%i particles", nright);
    else sprintf(message, "%i particle", nright);
    write_text_fixedwidth(xr, yr, message);
}

void print_circle_part_number(double *configs[NPARTMAX], int active[NPARTMAX], double xr, double yr)
{
    char message[50];
    int i, npart = 0;
    double rgb[3], x1, y1, cosphi, sinphi;
    
    /* count active particles, using the fact that absorbed particles have been given dummy coordinates */
    for (i=0; i<nparticles; i++) if (configs[i][0] >= DUMMY_ABSORBING) npart++;
        
    hsl_to_rgb(0.0, 0.0, 0.0, rgb);
    erase_area(xr, yr - 0.03, 0.4, 0.12, rgb);
    
    glColor3f(1.0, 1.0, 1.0);
    if (npart > 1) sprintf(message, "%i particles", npart);
    else sprintf(message, "%i particle", npart);
    write_text(xr, yr, message);
}

void print_collision_number(int ncollisions, double x, double y)
{
    char message[50];
    double rgb[3];
            
    hsl_to_rgb(0.0, 0.0, 0.0, rgb);
    erase_area(x, y, 0.5, 0.1, rgb);
    
    glColor3f(1.0, 1.0, 1.0);
    sprintf(message, "%i collisions", ncollisions);
    write_text(x, y, message);
    
}

void animation()
{
    double time, dt, alpha, r, rgb[3], alphamax;
    double *configs[NPARTMAX];
    int i, j, resamp = 1, s, i1, i2, c, lengthmax;
    int *color, *newcolor, *active, *heatmap_number, *heatmap_total;
    short int *heatmap_visited;
    char message[100];
    t_exit *exits; 
//     t_circle *circles;      /* experimental */
    
    /* Since NPARTMAX can be big, it seemed wiser to use some memory allocation here */
    color = malloc(sizeof(int)*(NPARTMAX));
    newcolor = malloc(sizeof(int)*(NPARTMAX));
    active = malloc(sizeof(int)*(NPARTMAX));
//     circles = malloc(sizeof(t_circle)*(NMAXCIRCLES));      /* experimental */
    
    if (HEATMAP) 
    {
        heatmap_number = malloc(sizeof(int)*(NXMAZE*NYMAZE+1));
        heatmap_total = malloc(sizeof(int)*(NXMAZE*NYMAZE+1));
        heatmap_visited = malloc(sizeof(short int)*(NXMAZE*NYMAZE+1));
        for (i=0; i<NXMAZE*NYMAZE+1; i++) heatmap_number[i] = 0;
        for (i=0; i<NXMAZE*NYMAZE+1; i++) heatmap_total[i] = 0;
        for (i=0; i<NXMAZE*NYMAZE+1; i++) heatmap_visited[i] = 0;
    }
    
    if (PRINT_LEFT_RIGHT_PARTICLE_NUMBER)
    {
        exits = malloc(sizeof(t_exit)*(NPARTMAX));
        for (i=0; i<NPARTMAX; i++)
        {
            exits[i].left = 0;
            exits[i].right = 0;
        }
    }
    
    for (i=0; i<NPARTMAX; i++)
        configs[i] = (double *)malloc(8*sizeof(double));

    /* init circle configuration if the domain is D_CIRCLES */
    if ((B_DOMAIN == D_CIRCLES)||(B_DOMAIN == D_CIRCLES_IN_RECT)||(B_DOMAIN == D_CIRCLES_IN_GENUSN)
        ||(B_DOMAIN == D_CIRCLES_IN_TORUS)) init_circles(circles);
    
    else if ((B_DOMAIN == D_POLYLINE)||(B_DOMAIN == D_POLYLINE_ARCS)) init_polyline(polyline, circles, arcs);

    if (POLYLINE_PATTERN == P_TOKA_PRIME)
    {
        x_shooter = -polyline[84].x1;
        y_shooter = polyline[84].y1;
        x_target = polyline[84].x1;
        y_target = polyline[84].y1;
    }

    /* initialize system by putting particles in a given point with a range of velocities */
    r = cos(PI/(double)NPOLY)/cos(DPI/(double)NPOLY);
    
//     init_partial_drop_config(LAMBDA, 0.0, 0.0, DPI, 0, NPART/4, 0, configs, color, newcolor);
//     init_partial_drop_config(-LAMBDA, 0.0, 0.0, DPI, NPART/4, NPART/2, 0, configs, color, newcolor);
//     init_partial_drop_config(0.0, LAMBDA, 0.0, DPI, NPART/2, 3*NPART/4, 0, configs, color, newcolor);
//     init_partial_drop_config(0.0, -LAMBDA, 0.0, DPI, 3*NPART/4, NPART, 0, configs, color, newcolor);
    
//     init_drop_config(-1.0 + 0.3*sqrt(2.0), -1.0 + 0.5*sqrt(2.0), 0.0, DPI, configs);

//     init_line_config(0.0, 0.0, 0.0, 0.9, 0.0, configs);   
//     init_drop_config(-0.95, 0.0, -0.103 + DPI/15.0, -0.1 + DPI/15.0, configs);
    
    /* find long trajectory */
//     alphamax = 0.0;
//     lengthmax = 1;
//     for (alpha = 0.0; alpha < DPI; alpha += 0.00001)
//     {
//         init_drop_config(x_shooter, y_shooter, alpha, alpha + DPI, configs);
//         i = 0;
//         c = 1;
//         while ((c >= 0)&&(i<=1000))
//         {    
//             c = vbilliard(configs[0]);
//             i++;
//         }
//         if (i > 100) printf("Angle %.6lg, length %i\n", alpha, i);
//         if (i > lengthmax)
//         {
//             lengthmax = i;
//             alphamax = alpha;
//         }
//     }
//     printf("Angle %.6lg, max length %i\n", alphamax, lengthmax);
    
//     alphamax = 2.50949;
//     init_drop_config(x_shooter, y_shooter, alphamax, alphamax + DPI, configs);
    
    init_drop_config(0.05, 0.05, 0.0, DPI, configs);   
//     init_drop_config(-0.95, 0.95, 0.0, DPI, configs);   
    
//     init_sym_drop_config(-1.0, 0.5, -PID, PID, configs);
//     init_drop_config(-0.999, 0.0, -alpha, alpha, configs);

//  other possible initial conditions :
//     init_line_config(-1.3, -0.3, -1.2, -0.3, PID, configs);
//     init_line_config(0.0, 0.0, 0.5, 0.0, PID, configs);
//     init_line_config(0.0, 0.0, 0.0, -0.5, PI, configs);
//     init_line_config(-1.25, -0.5, -1.25, 0.5, 0.0*PID, configs);
//     init_line_config(-1.0, -0.3, -1.0, 0.3, 0.0, configs);
//     init_line_config(-0.7, -0.45, -0.7, 0.45, 0.0, configs);
//     init_line_config(-1.5, 0.1, -0.1, 1.0, -0.5*PID, configs);
    
    if (!SHOWTRAILS) blank();
    glColor3f(0.0, 0.0, 0.0);
    if (DRAW_BILLIARD) draw_billiard();
    if (PRINT_PARTICLE_NUMBER) print_part_number(configs, active, XMIN + 0.1, YMIN + 0.1);
    else if (PRINT_LEFT_RIGHT_PARTICLE_NUMBER) 
        print_left_right_part_number(configs, active, XMIN + 0.05, YMIN + 0.05, XMAX - 0.35, YMIN + 0.05, exits);
    else if (PRINT_CIRCLE_PARTICLE_NUMBER) print_circle_part_number(configs, active, XMAX - 0.45, YMIN + 0.05);
    else if (PRINT_COLLISION_NUMBER) print_collision_number(ncollisions, XMIN + 0.1, YMIN + 0.1);
    
    glutSwapBuffers();   
  
  
    for (i=0; i<NPARTMAX; i++) 
    {
        color[i] = 0;
        newcolor[i] = 0;
        active[i] = 1;
    }
    
    if (FLOWER_COLOR)   /* adapt color scheme to flower configuration (beta implementation) */
    {
//         i1 = (int)((double)NPART*0.2538);     /* the 0.27 is just a trial-and-error guess, to be improved */
//         i1 = (int)((double)NPART*0.1971);     /* the 0.27 is just a trial-and-error guess, to be improved */
        i1 = (int)((double)NPART*0.3015);     /* the 0.27 is just a trial-and-error guess, to be improved */
        i2 = NPART-i1;
        for (i=i1; i<i2; i++) 
        {
            color[i] += NCOLORS/3;
            newcolor[i] = NCOLORS/3;
        }
        for (i=i2; i<NPART; i++) 
        {
            color[i] += 2*NCOLORS/3;
            newcolor[i] = 2*NCOLORS/3;
        }
    }
  
    if (RAINBOW_COLOR)      /* rainbow color scheme */
        for (i=0; i<NPART; i++) 
        {
            color[i] = (i*NCOLORS)/NPART;
            newcolor[i] = (i*NCOLORS)/NPART;  
        }
        
    if (TEST_INITIAL_COND) nparticles = test_initial_condition(configs, active, newcolor);

    sleep(SLEEP1);
    
    
    /* initialize drops in different colors */
//     init_partial_drop_config(0.0, 0.0, 0.0, DPI, 0, 2*NPART/5, 0, configs, color, newcolor);
//     init_partial_drop_config(0.0, 0.8, 0.0, DPI, 2*NPART/5, 4*NPART/5, 30, configs, color, newcolor);
//     init_partial_drop_config(LAMBDA - 0.05, 0.1, 0.0, DPI, 4*NPART/5, NPART, 60, configs, color, newcolor);
  
    for (i=0; i<=NSTEPS; i++)
    {
        graph_movie(TIME, newcolor, configs, active);
        
        if (SHOWTRAILS) draw_config_showtrails(newcolor, configs, active);
        else if (HEATMAP) 
        {
            draw_config_heatmap(configs, active, heatmap_number, heatmap_total, heatmap_visited, DRAW_HEATMAP_PARTICLES);
            if (DRAW_HEATMAP_HISTOGRAM) draw_heatmap_histogram(heatmap_number, heatmap_total, i);
//             draw_config(newcolor, configs, active);
        }
        else draw_config(newcolor, configs, active);
//         draw_config(newcolor, configs, active);
        if (DRAW_BILLIARD) draw_billiard();
        if (PRINT_PARTICLE_NUMBER) print_part_number(configs, active, XMIN + 0.1, YMIN + 0.1);
        else if (PRINT_LEFT_RIGHT_PARTICLE_NUMBER) 
            print_left_right_part_number(configs, active, XMIN + 0.05, YMIN + 0.05, XMAX - 0.35, YMIN + 0.05, exits);
        else if (PRINT_CIRCLE_PARTICLE_NUMBER) print_circle_part_number(configs, active, XMAX - 0.45, YMIN + 0.05);
//             print_left_right_part_number(configs, XMIN + 0.1, YMIN + 0.1, XMAX - 0.45, YMIN + 0.1, YMIN + MAZE_XSHIFT, YMAX + MAZE_XSHIFT);
        else if (PRINT_COLLISION_NUMBER) print_collision_number(ncollisions, XMIN + 0.1, YMIN + 0.1);
    
        for (j=0; j<NPARTMAX; j++) color[j] = newcolor[j];
        
        /* draw initial points */
//         draw_initial_condition_circle(0.0, 0.0, 0.02, 0);
//         draw_initial_condition_circle(0.0, 0.8, 0.02, 10);
//         draw_initial_condition_circle(1.2, 0.1, 0.02, 36);
        
	if (MOVIE) 
        {
            if (INVERT_COUNTER) save_frame_counter(NSTEPS+1-i);
            else save_frame();
            
            /* it seems that saving too many files too fast can cause trouble with the file system */
            /* so this is to make a pause from time to time - parameter PAUSE may need adjusting   */
            if (i % PAUSE == PAUSE - 1) 
            {
                printf("Making a short pause\n");
                sleep(PSLEEP);
                s = system("mv part*.tif tif_part/");
            }
        }
        else printf("Frame %i\n", i);
    }
 
    if (MOVIE)
    {
        if (DRAW_FINAL_HEATMAP) 
            draw_config_heatmap(configs, active, heatmap_number, heatmap_total, heatmap_visited, 0);
        if (DRAW_HEATMAP_HISTOGRAM) draw_heatmap_histogram(heatmap_number, heatmap_total, NSTEPS);
        for (i=0; i<END_FRAMES; i++) 
        {
            if (INVERT_COUNTER) 
            {
                if (i == 0)
                {
                    sprintf(message, "mv part.%05i.tif tif_part/", NSTEPS+1);
                    s = system(message);
                }
                sprintf(message, "cp tif_part/part.%05i.tif tif_part/part.%05i.tif", NSTEPS+1, NSTEPS+i+2);
                s = system(message); 
            }
            else save_frame();
        }
        s = system("mv part*.tif tif_part/");
    }
    
    free(color);
    free(newcolor);
    if (HEATMAP) 
    {
        free(heatmap_number);
        free(heatmap_total);
        free(heatmap_visited);
    }
    if (PRINT_LEFT_RIGHT_PARTICLE_NUMBER) free(exits);
    for (i=0; i<NPARTMAX; i++) free(configs[i]);
 
}


void display(void)
{
    time_t rawtime;
    struct tm * timeinfo;

    time(&rawtime);
    timeinfo = localtime(&rawtime);
    
    glPushMatrix();

    blank();
    
    if (!SHOWTRAILS)
    {
        glutSwapBuffers();
        blank();
        glutSwapBuffers();
    }

    animation();        

    sleep(SLEEP2); 

    glPopMatrix();
    
    glutDestroyWindow(glutGetWindow());

    printf("Start local time and date: %s", asctime(timeinfo));
    time(&rawtime);
    timeinfo = localtime(&rawtime);
    printf("Current local time and date: %s", asctime(timeinfo));
}


int main(int argc, char** argv)
{
    glutInit(&argc, argv);
    if (SHOWTRAILS) glutInitDisplayMode(GLUT_RGB | GLUT_SINGLE | GLUT_DEPTH);
    else glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH);
//     glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH);
    glutInitWindowSize(WINWIDTH,WINHEIGHT);
    glutCreateWindow("Billiard animation");
       
    init();

    glutDisplayFunc(display);

    glutMainLoop();
  
    return 0;
}

