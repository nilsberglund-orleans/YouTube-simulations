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
#include <time.h>

#define MOVIE 0         /* set to 1 to generate movie */

#define WINWIDTH 	1760  /* window width */
#define WINHEIGHT 	990   /* window height */

#define SAVE_MEMORY 1   /* set to 1 to save memory when writing tiff images */
#define NO_EXTRA_BUFFER_SWAP 1    /* some OS require one less buffer swap when recording images */

#define XMIN -1.2
#define XMAX 2.8	/* x interval */
#define YMIN -1.125
#define YMAX 1.125	/* y interval for 9/16 aspect ratio */

#define SCALING_FACTOR 1.0       /* scaling factor of drawing, needed for flower billiards, otherwise set to 1.0 */

/* Choice of the billiard table, see global_particles.c */

#define B_DOMAIN 30     /* choice of domain shape */

#define CIRCLE_PATTERN 1   /* pattern of circles */
#define POLYLINE_PATTERN 8  /* pattern of polyline */

#define ABSORBING_CIRCLES 1 /* set to 1 for circular scatterers to be absorbing */
#define NABSCIRCLES 10       /* number of absorbing circles */

#define NMAXCIRCLES 50000        /* total number of circles (must be at least NCX*NCY for square grid) */
#define NMAXPOLY 100004        /* total number of sides of polygonal line */   
#define NCX 9             /* number of circles in x direction */
#define NCY 20            /* number of circles in y direction */
#define NPOISSON 500        /* number of points for Poisson C_RAND_POISSON arrangement */
#define NGOLDENSPIRAL 2000  /* max number of points for C_GOLDEN_SPIRAL arrandement */
#define SDEPTH 2            /* Sierpinski gastket depth */

#define LAMBDA 1.0 	/* parameter controlling shape of domain */
#define MU 0.0075       /* second parameter controlling shape of billiard */
#define FOCI 1          /* set to 1 to draw focal points of ellipse */
#define NPOLY 5             /* number of sides of polygon */
#define APOLY 0.2           /* angle by which to turn polygon, in units of Pi/2 */ 
#define LAMBDA_B 1.0  /* parameter controlling shape of domain (for P_POLYRING) */
#define NPOLY_B 100000           /* number of sides of second polygon */
#define APOLY_B 1.0         /* angle by which to turn second polygon, in units of Pi/2 */ 
#define DRAW_BILLIARD 1     /* set to 1 to draw billiard */
#define DRAW_CONSTRUCTION_LINES 0   /* set to 1 to draw additional construction lines for billiard */
#define PERIODIC_BC 0       /* set to 1 to enforce periodic boundary conditions when drawing particles */
#define PENROSE_RATIO 2.5    /* parameter controlling the shape of small ellipses in Penrose room */

#define RESAMPLE 0      /* set to 1 if particles should be added when dispersion too large */
#define DEBUG 0         /* draw trajectories, for debugging purposes */

/* Simulation parameters */

#define NPART 1         /* number of particles */
#define NPARTMAX 100000	/* maximal number of particles after resampling */
#define TRAJ_LENGTH 120 /* length of trajectory */
#define PLOT_NMAX 100   /* max length on length plot */
#define LMAX 0.01       /* minimal segment length triggering resampling */ 
#define DMIN 0.02       /* minimal distance to boundary for triggering resampling */ 
#define CYCLE 0         /* set to 1 for closed curve (start in all directions) */
#define SHOWTRAILS 0    /* set to 1 to keep trails of the particles */
#define HEATMAP 0       /* set to 1 to show heat map of particles */
#define PLOT_HEATMAP_AVERAGE 0      /* set to 1 to plot average number of particles in heat map */
#define SHOWZOOM 0      /* set to 1 to show zoom on specific area */
#define TEST_ACTIVE 0   /* set to 1 to test whether particle is in billiard */
#define PRINT_TRAJECTORY_LENGTH 1   /* set to 1 to print length of trajectory 0 */
#define PRINT_TRAJECTORY_ANGLE 1    /* set to 1 to print angle of trajectory 0 */
#define PRINT_TRAJECTORY_SLOPE 0    /* set to 1 to print slope of trajectory 0 */
#define PRINT_TRAJECTORY_PERIOD 0   /* set to 1 to print period of trajectory 0 */
#define DRAW_LENGTHS_PLOT 1         /* set to 1 to plot trajectory lengths */
#define LENGTHS_LOG_SCALE 0         /* set to 1 to use log scale for plot of lengths */
#define LENGTH_PLOT_POLAR 0         /* set to 1 to plot lengths in polar coordinates */
#define MIN_ANGLE 0.0           /* range of angles of trajectory */
#define MAX_ANGLE 135.0         /* range of angles of trajectory */
#define TEST_INITIAL_COND 0     /* set to 1 to allow only initial conditions that pass a test */

#define SLOW_AT_LONG_TRAJ 0     /* set to 1 to slow down movie for long trajectories */
#define ADD_SUCCESS_GALLERY 0   /* set to 1 to add gallery of successful trajectories at end of movie */
#define SUCCESS_GALLERY_FRAMES 25   /* number of frames per success */
#define EXIT_BOTH_WAYS 1        /* set to 1 to add exits to he left to succesful trajectories */

#define NSTEPS 6400       /* number of frames of movie */
#define TIME 2500         /* time between movie frames, for fluidity of real-time simulation */ 
#define DPHI 0.00001     /* integration step */
#define NVID 150         /* number of iterations between images displayed on screen */

/* Decreasing TIME accelerates the animation and the movie                               */
/* For constant speed of movie, TIME*DPHI should be kept constant                        */
/* However, increasing DPHI too much deterioriates quality of simulation                 */
/* NVID tells how often a picture is drawn in the animation, increase it for faster anim */
/* For a good quality movie, take for instance TIME = 400, DPHI = 0.00005, NVID = 100    */

/* Colors and other graphical parameters */

#define COLOR_PALETTE 11    /* Color palette, see list in global_pdes.c  */

#define NCOLORS 64       /* number of colors */
#define COLORSHIFT 0     /* hue of initial color */ 
#define COLOR_HUEMIN 0  /* minimal color hue */
#define COLOR_HUEMAX 360 /* maximal color hue */
#define RAINBOW_COLOR 0  /* set to 1 to use different colors for all particles */
#define FLOWER_COLOR 0   /* set to 1 to adapt initial colors to flower billiard (tracks vs core) */
#define NSEG 100         /* number of segments of boundary */
#define LENGTH 0.03       /* length of velocity vectors */
#define BILLIARD_WIDTH 2    /* width of billiard */
#define PARTICLE_WIDTH 2    /* width of particles */
#define FRONT_WIDTH 3       /* width of wave front */
#define GRAPH_HUE 200.0      /* color hue of graph */
#define SUCCESS_HUE 30.0   /* color hue of success */

#define BLACK 1             /* set to 1 for black background */
#define COLOR_OUTSIDE 0     /* set to 1 for colored outside */ 
#define OUTER_COLOR 270.0   /* color outside billiard */
#define PAINT_INT 0         /* set to 1 to paint interior in other color (for polygon/Reuleaux) */
#define PAINT_EXT 1         /* set to 1 to paint exterior */

#define PAUSE 1000       /* number of frames after which to pause */
#define PSLEEP 2         /* sleep time during pause */
#define SLEEP1  1        /* initial sleeping time */
#define SLEEP2  1        /* final sleeping time */
#define END_FRAMES  250   /* number of frames at end of movie */

#define NXMAZE 25      /* width of maze */
#define NYMAZE 25      /* height of maze */
#define MAZE_MAX_NGBH 8     /* max number of neighbours of maze cell */
#define RAND_SHIFT 10       /* seed of random number generator */
#define MAZE_XSHIFT -0.7     /* horizontal shift of maze */
#define MAZE_RANDOM_FACTOR 0.1     /* randomization factor for S_MAZE_RANDOM */
#define MAZE_CORNER_RADIUS 0.4     /* radius of tounded corners in maze */
#define CLOSE_MAZE 0        /* set to 1 to close maze exits */

#include "global_particles.c"
#include "sub_maze.c"
#include "sub_part_billiard.c"


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


void draw_zoom(int color[NPARTMAX], double *trajectory[NPARTMAX], int traj_length[NPARTMAX], double x_target, double y_target, double width, double shiftx, double shifty, double zoomwidth, int shooter)
/* draw zoom around target (for laser in room of mirrors) */
{
    int i, j, col, xw, yw, ww;
    double x1, y1, x2, y2, xb, yb, cosphi, sinphi, rgb[3], tradius, phi, lum;
    
    x1 = shiftx - zoomwidth;
    y1 = shifty - zoomwidth;
    
    /* convert to window coordinates for use of GL_SCISSOR_TEST - should be improvable */
    xw = (int)((x1 - XMIN)/(XMAX - XMIN)*(double)WINWIDTH);
    yw = (int)((y1 - YMIN)/(YMAX - YMIN)*(double)WINHEIGHT);
    ww = (int)(2.0*zoomwidth/(XMAX - XMIN)*(double)WINWIDTH);
    
    /* draw the trajectories */
    glEnable(GL_SCISSOR_TEST);
    glScissor(xw, yw, ww, ww);
    
    glLineWidth(BILLIARD_WIDTH*2);
    for (j=TRAJ_LENGTH-1; j>=1; j--)
    {
        lum = 0.5*(1.0 - (double)j/(double)TRAJ_LENGTH);
                
        for (i=0; i<NPART; i++) 
        {
            x1 = (trajectory[i][2*j] - x_target)/width;
            y1 = (trajectory[i][2*j+1] - y_target)/width;
            x2 = (trajectory[i][2*j-2] - x_target)/width;
            y2 = (trajectory[i][2*j-1] - y_target)/width;
            
            if (j < traj_length[i])
            {
                if (RAINBOW_COLOR) col = color[i]; 
                else col = (color[i] + j)%NCOLORS;
                rgb_color_scheme_lum(col, lum, rgb);
                glColor3f(rgb[0], rgb[1], rgb[2]);
                glBegin(GL_LINE_STRIP);
                glVertex2d(shiftx + zoomwidth*SCALING_FACTOR*x1, shifty + zoomwidth*SCALING_FACTOR*y1);
                glVertex2d(shiftx + zoomwidth*SCALING_FACTOR*x2, shifty + zoomwidth*SCALING_FACTOR*y2);
                glEnd ();
            }
        }
    }
    glDisable(GL_SCISSOR_TEST);
    
    /* erase area outside zoomed area */
//     x1 = shiftx - zoomwidth;
//     y1 = shifty - zoomwidth;
//     x2 = shiftx + zoomwidth;
//     y2 = shifty + zoomwidth;
//     rgb[0] = 0.0; rgb[1] = 0.0; rgb[2] = 0.0;
//     erase_rectangle(XMIN, YMIN, x1, YMAX, rgb);
//     erase_rectangle(x2, YMIN, XMAX, YMAX, rgb);
//     erase_rectangle(XMIN, YMIN, XMAX, y1, rgb);
//     erase_rectangle(XMIN, y2, XMAX, YMAX, rgb);
    
    /* draw zoom boundary */
    x1 = shiftx - zoomwidth;
    y1 = shifty - zoomwidth;
    x2 = shiftx + zoomwidth;
    y2 = shifty + zoomwidth;
    glLineWidth(BILLIARD_WIDTH*2);
    glColor3f(1.0, 1.0, 1.0);
    glBegin(GL_LINE_LOOP);    
    glVertex2d(x1, y1);
    glVertex2d(x2, y1);
    glVertex2d(x2, y2);
    glVertex2d(x1, y2);
    glEnd();
    

    /* draw billiard boundaries in zoom */
    glLineWidth(BILLIARD_WIDTH*2);
    
    if (POLYLINE_PATTERN == P_ISOCELES_TRIANGLE)
    {
        draw_line(shiftx, shifty, x1, shifty);
        draw_line(shiftx, shifty, x1, y2);
    }
    else if (y_target + width > 1.0)
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
    
    glBegin(GL_LINE_LOOP);
    tradius = zoomwidth*MU/width;
    for (i=0; i<=NSEG; i++)
    {
        phi = (double)i*DPI/(double)NSEG;
        x1 = shiftx + tradius*cos(phi);
        y1 = shifty + tradius*sin(phi);
        glVertex2d(x1, y1);
    }
    glEnd ();
    
//     glLineWidth(PARTICLE_WIDTH*2);
}

void draw_zoom_area(double x_target, double y_target, double width)
/* draw rectangle marking zoomed area */
{
    double x1, x2, y1, y2;
    
    glEnable(GL_LINE_SMOOTH);
    glLineWidth(BILLIARD_WIDTH/2);
    
    x1 = x_target - width;
    y1 = y_target - width;
    x2 = x_target + width;
    y2 = y_target + width;

    glColor3f(1.0, 1.0, 1.0);
    glBegin(GL_LINE_LOOP);    
    glVertex2d(x1, y1);
    glVertex2d(x2, y1);
    glVertex2d(x2, y2);
    glVertex2d(x1, y2);
    glEnd();
}

void draw_trajectories_lshape(int color[NPARTMAX], double *trajectory[NPARTMAX], int traj_length[NPARTMAX])
/* draw the trajectories */
{
    int i, j;
    double rgb[3], lum, x1, y1, x2, y2, x3, y3, x4, y4;
    
    for (j=TRAJ_LENGTH-2; j>=1; j--)
    {
        lum = 0.5*(1.0 - (double)j/(double)TRAJ_LENGTH);
        
        for (i=0; i<NPART; i++) 
        {
            x1 = trajectory[i][2*j];
            y1 = trajectory[i][2*j+1];
            x2 = trajectory[i][2*j-2];
            y2 = trajectory[i][2*j-1];
            
//             printf("j = %i, z1 = (%.2f, %.2f), z2 = (%.2f, %.2f)\n", j, x1, y1, x2, y2); 
            
            /* take care of boundary conditions */
            if (x1 == -1.0) 
            {
                if (y1 >= 0.0) x3 = 0.0;
                else x3 = 1.0;
                y3 = y1;
            }
            else if ((x1 == 0.0)||(x1 == 1.0))
            {
                x3 = -1.0;
                y3 = y1;
            }
            if (y1 == -1.0) 
            {
                if (x1 >= 0.0) y3 = 0.0;
                else y3 = 1.0;
                x3 = x1;
            }
            else if ((y1 == 0.0)||(y1 == 1.0))
            {
                y3 = -1.0;
                x3 = x1;
            }
            
            if (j < traj_length[i])
            {
                rgb_color_scheme_lum(color[i], lum, rgb);
                glColor3f(rgb[0], rgb[1], rgb[2]);
//                 glBegin(GL_LINE_STRIP);
//                 glVertex2d(x1, y1);
//                 glVertex2d(x2, y2);
//                 glEnd ();
                glBegin(GL_LINE_STRIP);
                glVertex2d(x2, y2);
                glVertex2d(x3, y3);
                glEnd ();
            }
        }
    }
 }


void draw_trajectories(int color[NPARTMAX], double *trajectory[NPARTMAX], int traj_length[NPARTMAX])
/* draw the trajectories */
{
    int i, j, col;
    double rgb[3], lum, x1, y1, x2, y2;
    
    glutSwapBuffers(); 
    if (!SHOWTRAILS) blank();
    if (PAINT_INT) paint_billiard_interior();
    
    if (SHOWZOOM) switch (POLYLINE_PATTERN) {
        case (P_TOKA_PRIME):
        {
            draw_zoom(color, trajectory, traj_length, x_target, y_target, 0.05, 1.65, 0.75, 0.3, 0);
            draw_zoom(color, trajectory, traj_length, x_shooter, y_shooter, 0.05, -1.65, 0.75, 0.3, 1);
            break;
        }
        case (P_TOKA_NONSELF):
        {
            draw_zoom(color, trajectory, traj_length, 0.0, 0.0, 0.1, 1.65, 0.75, 0.3, 1);
            break;
        }
        case (P_ISOCELES_TRIANGLE):
        {
            draw_zoom(color, trajectory, traj_length, 1.0, -1.0, 0.2, 0.6, 0.6, 0.4, 0);
            break;
        }
    }
      
    glLineWidth(PARTICLE_WIDTH);
    
    glEnable(GL_LINE_SMOOTH);
    
    /* if domain is L shape with periodic boundary conditions, we have to take care of those */
    if (B_DOMAIN == D_LSHAPE) draw_trajectories_lshape(color, trajectory, traj_length);

    else for (j=TRAJ_LENGTH-1; j>=1; j--)
    {
        lum = 0.5*(1.0 - 0.5*(double)j/(double)TRAJ_LENGTH);
        
        for (i=0; i<NPART; i++) 
        {
            x1 = trajectory[i][2*j];
            y1 = trajectory[i][2*j+1];
            x2 = trajectory[i][2*j-2];
            y2 = trajectory[i][2*j-1];
            
//             printf("j = %i, z1 = (%.2f, %.2f), z2 = (%.2f, %.2f)\n", j, x1, y1, x2, y2); 
            if (j < traj_length[i])
            {
                if (RAINBOW_COLOR) col = color[i]; 
                else col = (color[i] + j)%NCOLORS;
//                 rgb_color_scheme_minmax_lum(col, lum, rgb);
//                 rgb_color_scheme_lum(col, lum, rgb);
                rgb_color_scheme_minmax_lum(col, lum, rgb); 
                glColor3f(rgb[0], rgb[1], rgb[2]);
                glBegin(GL_LINE_STRIP);
                glVertex2d(x1, y1);
                glVertex2d(x2, y2);
                glEnd ();
            }
        }
        
    }
    
    if (SHOWZOOM) switch (POLYLINE_PATTERN) {
        case (P_TOKA_PRIME):
        {
            draw_zoom_area(x_target, y_target, 0.05);
            draw_zoom_area(x_shooter, y_shooter, 0.05);
            break;
        }
        case (P_TOKA_NONSELF):
        {
            draw_zoom_area(0.0, 0.0, 0.05);
            break;
        }
    }
}


void draw_config_heatmap(double *trajectory[NPARTMAX], int traj_length[NPARTMAX], double *configs[NPARTMAX], int active[NPARTMAX], int heatmap_number[NXMAZE*NYMAZE+1], int heatmap_total[NXMAZE*NYMAZE+1], int draw_particles)
/* draw a heat map of particle distribution (for mazes) */
{
    int i, j, k, n, part_number;
    double x, y, cosphi, sinphi, rgb[3], len, padding = 0.02, x1, y1, x2, y2, l, pathstep = 0.001;
    double *xtable, *ytable;
    short int *drawtable;
    static int time, first = 1;
    static double minprop;
    
    if (first)
    {
        time = 0;
        first = 0;
        minprop = 1000.0;
//         if (PLOT_HEATMAP_AVERAGE) minprop = 0.005;
//         else minprop = 0.01;
    }
    time++;
    
    drawtable = malloc(sizeof(short int)*(NPARTMAX));
    xtable = malloc(sizeof(double)*(NPARTMAX));
    ytable = malloc(sizeof(double)*(NPARTMAX));
    
//     glutSwapBuffers(); 
//     blank();
    if (PAINT_INT) paint_billiard_interior();
    
    for (i=0; i<NXMAZE*NYMAZE+1; i++) heatmap_number[i] = 0;
    
    for (i=0; i<NPART; i++) 
    {
        for (j=TRAJ_LENGTH-1; j>=1; j--) if (j < traj_length[i])
        {
            x1 = trajectory[i][2*j];
            y1 = trajectory[i][2*j+1];
            x2 = trajectory[i][2*j-2];
            y2 = trajectory[i][2*j-1];
            
            len = module2(x2 - x1, y2 - y1);
            
            for (k=0; k<(int)(len/pathstep); k++)
            {
                x = x1 + (x2 - x1)*(double)k*pathstep/len;
                y = y1 + (y2 - y1)*(double)k*pathstep/len;
                
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
                    }
            
//             if (HEATMAP_MAX_PART_BY_CELL > 0)
//             {
//                 drawtable[i] = ((n == -2)||((n > 0)&&(heatmap_number[n] <= HEATMAP_MAX_PART_BY_CELL)));
//             }
//             else 
                    drawtable[i] = (n > -1);
                }
            }
        }
        
    }
    
    for (n=0; n<NXMAZE*NYMAZE+1; n++) 
    {
        if (PLOT_HEATMAP_AVERAGE) part_number = (int)((double)heatmap_total[n]/(double)time);
        else part_number = heatmap_number[n];
        draw_maze_cell(n, part_number/10, minprop);
//         printf("%i particles in maze cell %i\n", color, n);
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

double angle_schedule(double time)
{
    return((MIN_ANGLE + time*(MAX_ANGLE - MIN_ANGLE))*DPI/360.0);
//     return((1.0 + time)*0.25*PID);
//     return(time*0.5*PID);
}

double log_scale_lplot(int length)
{
    return(0.15*log((double)length));
}

void plot_lengths_linear(int traj_length_table[NSTEPS+1], short int traj_test_table[NSTEPS+1], int i)
/* draw a plot of thajectory lengths */
{
    int j, k, nmax = PLOT_NMAX, window = 5, imin, l;
    double x, y, x0, y0, x1, y1, xmin = 1.2, xmax = 2.7, ymin = -1.0, ymax = 0.8, ymax_test = 0.5, dx, dy, angle_diff, angle0, shiftx; 
    char message[50];
    static int iw, first = 1, jmin, jmax;
    static double c, rgb[3], rgb_success[3], rgb1[3], rgb1_success[3], gradx_shift;
    
    ymin += (YMIN + 1.125);
    ymax += (YMIN + 1.125);
    
    if (first)
    {
        angle0 = angle_schedule(0.0)*360.0/DPI;
        angle_diff = (angle_schedule(1.0) - angle_schedule(0.0))*360.0/DPI;
        gradx_shift = angle0 - (double)((int)angle0);
        iw = (int)((double)NSTEPS*(double)window/angle_diff);
        c = ((angle_diff -(double)window)/(double)window)*(xmax - xmin)/(double)(NSTEPS - iw);
//         iw = (int)((double)NSTEPS*(double)window/MAX_ANGLE);
//         c = ((MAX_ANGLE -(double)window)/(double)window)*(xmax - xmin)/(double)(NSTEPS - iw);
        rgb_color_scheme(GRAPH_HUE, rgb);
        rgb_color_scheme(SUCCESS_HUE, rgb_success);
        for (k=0; k<3; k++)
        {
            rgb1[k] = rgb[k]*0.75;
            rgb1_success[k] = rgb_success[k]*0.75;
        }
        jmin = (int)(angle_schedule(0.0)*360.0/DPI);
        jmax = (int)(angle_schedule(1.0)*360.0/DPI);
        first = 0;
        printf("jmin = %i, jmax = %i\n", jmin, jmax);
//         sleep(1);
    }
    if (i < iw) imin = 0;
    else imin = i - iw;
    
    dx = (xmax - xmin)/((double)iw);
    dy = (ymax - ymin)/((double)nmax);
    x0 = xmin;
    y0 = ymin;
    
    glColor3f(1.0, 1.0, 1.0);
    glLineWidth(1);
    
    draw_line(xmin, ymin, xmax, ymin);
    draw_line(xmin, ymin, xmin, ymax + 0.0);
    
    sprintf(message, "angle");
    write_text_fixedwidth(xmax - 0.12, ymin - 0.05, message);
    
    sprintf(message, "length");
    write_text_fixedwidth(xmin - 0.1, ymax + 0.05, message);
    
    /* vertical axis graduation */
    x1 = xmin - 0.025;
    x = xmin + 0.025;
    l = 1;
    if (LENGTHS_LOG_SCALE) for (j=0; j<5; j++)
    {
        y = ymin + log_scale_lplot(l);
        draw_line(x1, y, x, y);
        sprintf(message, "%i", l);
        write_text_fixedwidth(xmin - 0.1 - 0.02*(double)j, y - 0.01, message);
        
        for (k=2; k<10; k++)
        {
            y = ymin + log_scale_lplot(k*l);
            if (k*l < 10000) draw_line(xmin - 0.015, y, xmin + 0.015, y);
        }
        l *= 10;
    }
    else for (j=0; j<nmax; j+=10)
    {
        y = ymin + (double)j*dy;
        draw_line(x1, y, x, y);
        sprintf(message, "%i", j);
        if (j >= 100) shiftx = 0.1;
        else if (j >= 10) shiftx = 0.0775;
        else shiftx = 0.055;
        write_text_fixedwidth(xmin - shiftx, y - 0.01, message);
    }
    
    /* horizontal axis graduation */
    for (j=0; j<jmax - jmin; j+=1) 
    {
        x = xmin + ((double)j-gradx_shift)*(xmax - xmin)/(double)window;
        if (i > iw) x -= c*(double)(i - iw);
        if ((x > xmin)&&(x < xmax - 0.05))
        {
            y = ymin - 0.025;
            y1 = ymin + 0.025;
            draw_line(x, y, x, y1);
            if (x < xmax - 0.15)
            {
                sprintf(message, "%i", jmin + j);
                if (j < 100) x += 0.01;
                if (j < 10) x += 0.015;
                write_text_fixedwidth(x - 0.04, ymin - 0.075, message);
            }
        }
    }
    
    /* to reduce aliasing, plot a broader graph at smaller luminosity */
    glLineWidth(3);
    if (TEST_INITIAL_COND)
    {
        glColor3f(rgb1_success[0], rgb1_success[1], rgb1_success[2]);
        for (j=imin; j<i; j++) if (traj_test_table[j])
        {
//             printf("Success path for j = %i\n", j);
            x = xmin + (double)(j-imin)*dx;
            draw_line(x, ymin, x, ymax_test);
        }
    }

    /* plot of path lengths */
    glColor3f(rgb1[0], rgb1[1], rgb1[2]);
    for (j=imin; j<i; j++)
    {
        x = xmin + (double)(j-imin)*dx;
        if (LENGTHS_LOG_SCALE) y = ymin +log_scale_lplot(traj_length_table[j]);
        else y = ymin + (double)traj_length_table[j]*dy;
        
        if (j>imin)
        {
            draw_line(x0, y0, x, y0);
            draw_line(x, y0, x, y);
        }
        
        x0 = x;
        y0 = y;
    }
    
    glLineWidth(2);
    x0 = xmin;
    y0 = ymin;
    
    if (TEST_INITIAL_COND)
    {
        glColor3f(rgb_success[0], rgb_success[1], rgb_success[2]);
        for (j=imin; j<i; j++) if (traj_test_table[j])
        {
//             printf("Success path for j = %i\n", j);
            x = xmin + (double)(j-imin)*dx;
            draw_line(x, ymin, x, ymax_test);
        }
    }

    /* plot of path lengths */
    glColor3f(rgb[0], rgb[1], rgb[2]);
    for (j=imin; j<i; j++)
    {
        x = xmin + (double)(j-imin)*dx;
        if (LENGTHS_LOG_SCALE) y = ymin +log_scale_lplot(traj_length_table[j]);
        else y = ymin + (double)traj_length_table[j]*dy;
        
        if (j>imin)
        {
            draw_line(x0, y0, x, y0);
            draw_line(x, y0, x, y);
        }
        
        x0 = x;
        y0 = y;
    }
    
    draw_initial_condition_circle(x0, y0, 0.02, 30);
    
//     draw_initial_condition_circle(x0, y0, 0.02, traj_length_table[i]);
    
}

void plot_lengths_polar(int traj_length_table[NSTEPS+1], short int traj_test_table[NSTEPS+1], int i)
/* draw a plot of thajectory lengths */
{
    int j, k;
    double xc, yc, x0, y0, x1, y1, xmin = 0.5, xmax = 1.9, ymin = -0.7, ymax = 0.3, phi, dphi, r; 
    char message[50];
    static int first = 1;
    static double rgb[3], rgb_success[3], logfactor;

    if (first)
    {
        rgb_color_scheme(GRAPH_HUE, rgb);
        rgb_color_scheme(SUCCESS_HUE, rgb_success);
        logfactor = 0.07;
    }
    
    dphi = DPI/(double)NSTEPS;
    xc = 0.5*(xmin + xmax);
    yc = 0.5*(ymin + ymax);
    
    glColor3f(1.0, 1.0, 1.0);
    glLineWidth(1);
    
    /* draw grid */
    r = 0.5*(xmax - xmin);
    for (j=0; j<12; j++) 
    {
        phi = (double)j*DPI/12.0;
        draw_line(xc, yc, xc + r*cos(phi), yc + r*sin(phi));
    }
    for (j=1; j<5; j++)
    {
        r = logfactor*(double)j*log(10.0);
        draw_circle(xc, yc, r, NSEG);
    }
            
    glLineWidth(2);
    /* draw plot of succesful paths */
//     if (TEST_INITIAL_COND)
//     {
//         r = 0.5*(xmax - xmin);
//         glColor3f(rgb_success[0], rgb_success[1], rgb_success[2]);
//         for (j=0; j<i; j++) if (traj_test_table[j])
//         {
//             phi = (MIN_ANGLE + (MAX_ANGLE - MIN_ANGLE)*(double)j/(double)NSTEPS)*DPI/360.0;
//             draw_line(xc, yc, xc + r*cos(phi), yc + r*sin(phi));
//         }
//     }
    
    /* draw plot of path lengths */
    glColor3f(rgb[0], rgb[1], rgb[2]);
    r = logfactor*log((double)traj_length_table[0]);
    phi = DPI*MIN_ANGLE/360.0;
    x0 = xc + r*cos(phi);
    y0 = yc + r*sin(phi);
    for (j=0; j<i; j++)
    {
        phi = (MIN_ANGLE + (MAX_ANGLE - MIN_ANGLE)*(double)j/(double)NSTEPS)*DPI/360.0;
        r = logfactor*log((double)traj_length_table[j]);
        x1 = xc + r*cos(phi);
        y1 = yc + r*sin(phi);
//         printf("phi = %.3f, r = %.3f, plot point (%.2f, %.2f)\n", phi, r, x1, y1);
        draw_line(x0, y0, x1, y1);
        x0 = x1;
        y0 = y1;
    }
    
    glColor3f(1.0, 1.0, 1.0);
    for (j=1; j<5; j++)
    {
        r = logfactor*(double)j*log(10.0);
        sprintf(message, "%i", (int)ipow(10.0, j));
        write_text_fixedwidth(xc + 0.01, yc + 0.02 + r, message);
    }
            
    
}

void plot_lengths(int traj_length_table[NSTEPS+1], short int traj_test_table[NSTEPS+1], int i)
/* draw a plot of thajectory lengths */
{
    if (LENGTH_PLOT_POLAR) plot_lengths_polar(traj_length_table, traj_test_table, i);
    else plot_lengths_linear(traj_length_table, traj_test_table, i);
}

int compute_trajectories_xy(double x, double y, double angle1, double angle2, double *configs[NPARTMAX], 
                             double *trajectory[NPARTMAX], int traj_length[NPARTMAX], double *xmax)
/* compute trajectories when starting in (x,y), with angles between angle1 and angle2 */
/* returns period of first trajectory */
{
    int i, j, c, c0, period = 0, loop = 0;
    double dalpha, alpha;
    double conf[2], pos[2], s0, u0;
  
    while (angle2 < angle1) angle2 += DPI;
    if (NPART > 1) dalpha = (angle2 - angle1)/((double)(NPART));
    else dalpha = 0.0;
    for (i=0; i<NPART; i++) 
    {
        alpha = angle1 + dalpha*((double)i);  
        period = 0;
//         printf("Angle = %.3lg\n", alpha*360.0/DPI); 
        
        trajectory[i][0] = x;
        trajectory[i][1] = y;
        traj_length[i] = 1;
      
        pos[0] = x;
        pos[1] = y;
        
        c = vbilliard_xy(configs[i], alpha, pos);
        
        s0 = configs[0][0];
        u0 = configs[0][1];
        
//         while ((c != DUMMY_SIDE_ABS)&&(c != -1)&&(traj_length[i] < TRAJ_LENGTH - 1)) 
        while ((c != DUMMY_SIDE_ABS)&&(c >= 0.0)&&(traj_length[i] < TRAJ_LENGTH - 1)) 
        {
//             c = vbilliard(configs[i]);
            pos_billiard(configs[i], pos, &alpha);
            c = vbilliard_xy(configs[i], alpha, pos);
            j = traj_length[i];
            trajectory[i][2*j] = configs[i][4];
            trajectory[i][2*j+1] = configs[i][5];
            traj_length[i]++;
            
            if ((!loop)&&(i == 0))
            {
                if (module2(configs[i][0] - s0, configs[i][1] - u0) > 0.001) period++;
                else loop = 1;
            }
        }
        if (i==0) c0 = c;
        j = traj_length[i];
        trajectory[i][2*j] = configs[i][6];
        trajectory[i][2*j+1] = configs[i][7];
        traj_length[i]++;
        
    }
    
    *xmax = trajectory[0][2*(traj_length[0]-1)];
    
//     printf("Period = %i\n", period);
    if ((c0 == DUMMY_SIDE_ABS)||(c0 < 0)) return(0);
    else return(period + 1);
//     if ((c != DUMMY_SIDE_ABS)&&(c >= 0.0)) return (period);
//     else return(0);
}

void compute_trajectories_ellipse(double x, double y, double *configs[NPARTMAX], 
                             double *trajectory[NPARTMAX], int traj_length[NPARTMAX])
/* compute trajectories when starting in (x,y), with angles between angle1 and angle2 */
{
    int i, j, c;
    double dalpha, alpha, a, b, cc, dx, d;
    double conf[2], pos[2];
  
    cc = sqrt(LAMBDA*LAMBDA - 1.0);
    dx = 4.0*cc/(double)NPART;
    
    for (i=0; i<NPART; i++) 
    {
        d = dx*(double)i;
        
        if (i%2 == 0) alpha = argument(d - x, -y);  
        else alpha = argument(-d + dx - x, -y);
        
        trajectory[i][0] = x;
        trajectory[i][1] = y;
        traj_length[i] = 1;
      
        pos[0] = x;
        pos[1] = y;
        
        c = vbilliard_xy(configs[i], alpha, pos);
        
        while ((c != DUMMY_SIDE_ABS)&&(traj_length[i] < TRAJ_LENGTH - 1)) 
        {
            c = vbilliard(configs[i]);
            j = traj_length[i];
            trajectory[i][2*j] = configs[i][4];
            trajectory[i][2*j+1] = configs[i][5];
            traj_length[i]++;
        }
        j = traj_length[i];
        trajectory[i][2*j] = configs[i][6];
        trajectory[i][2*j+1] = configs[i][7];
        traj_length[i]++;
        
    }
    
}

void compute_trajectory(double *configs[NPARTMAX], double *trajectory[NPARTMAX], int traj_length[NPARTMAX])
/* compute (x,y) coordinates of trajectory for all particles and fixed number of collisions */
{
    int i, j, k, c = 1;
    double config[8];
    
    for (i=0; i<nparticles; i++)
    {
//         printf("Computing trajectory %i\n", i);
        traj_length[i] = 1;
        
        for (k=0; k<8; k++) config[k] = configs[i][k];
        if (config[0] == DUMMY_ABSORBING) c = DUMMY_SIDE_ABS;
        else c = 1;
        
        for (j=0; j<TRAJ_LENGTH; j++) 
        {
            trajectory[i][2*j] = config[4];
            trajectory[i][2*j+1] = config[5];
            
            /* deal with absorbing circles */
            if (c != DUMMY_SIDE_ABS) c = vbilliard(config);
//             printf("Particle %i, c = %i\n", i, c);
            if (c != DUMMY_SIDE_ABS) traj_length[i]++;            
        }
        printf("Particle %i, trajectory length %i\n", i, traj_length[i]);
    }
    
    
}


void print_parameters(int i, int *traj_length, int *traj_length_table, int period, double alpha)
{
    char message[50];
    
    if (PRINT_TRAJECTORY_LENGTH)
    {
        printf("Length %i\n", traj_length[0] - 1);
        traj_length_table[i] = traj_length[0] - 1;
        sprintf(message, "Length %i\n", traj_length[0] - 1);
        glColor3f(1.0, 1.0, 1.0);
        write_text(XMAX - 0.5, YMAX - 0.125, message);
    }
    else if (PRINT_TRAJECTORY_PERIOD)
    {
        printf("Period %i\n", period);
        if (period == 0) sprintf(message, "Not periodic\n");
        else sprintf(message, "Period %i\n", period);
        glColor3f(1.0, 1.0, 1.0);
        write_text(XMIN + 0.1, YMIN + 0.04, message);
    }
    if (PRINT_TRAJECTORY_ANGLE)
    {
        printf("Angle %.3f\n", alpha*360.0/DPI);
        sprintf(message, "Angle %.3f\n", alpha*360.0/DPI);
        glColor3f(1.0, 1.0, 1.0);
        write_text(XMAX - 0.5, YMAX - 0.225, message);
    }
    else if (PRINT_TRAJECTORY_SLOPE)
    {
        printf("SLOPE %.3f\n", vabs(tan(alpha)));
        sprintf(message, "Slope %.3f\n", vabs(tan(alpha)));
        glColor3f(1.0, 1.0, 1.0);
        write_text(XMAX - 0.5, 0.9, message);
    }
}

void animation()
{
    double time, dt, alpha, r, rgb[3], x, y, a, b, c, xmax, xmin = 0.0;
    double *configs[NPARTMAX], *trajectory[NPARTMAX];
    int i, j, resamp = 1, s, i1, i2, period, npause, nsuccess = 0, success_counter = 0, gallery_counter, nend = END_FRAMES, firstswap = 1, frame;
    int *color, *newcolor, *active, *traj_length, *traj_length_table, *heatmap_number, *heatmap_total;
    short int *traj_test_table, first = 0, second = 0;
    char message[50];
    
    /* Since NPARTMAX can be big, it seemed wiser to use some memory allocation here */
    color = malloc(sizeof(int)*(NPARTMAX));
    newcolor = malloc(sizeof(int)*(NPARTMAX));
    active = malloc(sizeof(int)*(NPARTMAX));
    traj_length = malloc(sizeof(int)*(NPARTMAX));
    traj_length_table = malloc(sizeof(int)*(NSTEPS+1));
    traj_test_table = malloc(sizeof(short int)*(NSTEPS+1));
    for (i=0; i<NPARTMAX; i++)
    {
        configs[i] = (double *)malloc(8*sizeof(double));
        trajectory[i] = (double *)malloc(2*TRAJ_LENGTH*sizeof(double));
    }
    if (HEATMAP) 
    {
        heatmap_number = malloc(sizeof(int)*(NXMAZE*NYMAZE+1));
        heatmap_total = malloc(sizeof(int)*(NXMAZE*NYMAZE+1));
        for (i=0; i<NXMAZE*NYMAZE+1; i++) heatmap_number[i] = 0;
        for (i=0; i<NXMAZE*NYMAZE+1; i++) heatmap_total[i] = 0;
    }
    
    /* init circle configuration if the domain is D_CIRCLES */
//     if ((B_DOMAIN == D_CIRCLES)||(B_DOMAIN == D_CIRCLES_IN_RECT)||(B_DOMAIN == D_CIRCLES_IN_GENUSN)) init_circle_config();
    
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
    
    if (ABSORBING_CIRCLES)
    {
        if (B_DOMAIN == D_ELLIPSE)
        {
            ncircles = 1;
            circles[0].xc = -LAMBDA;
            circles[0].yc = 0.0;
            circles[0].radius = MU;
            circles[0].active = 1;
        }
        
    }
 
    /* initialize system by putting particles in a given point with a range of velocities */
    r = cos(PI/(double)NPOLY)/cos(DPI/(double)NPOLY);

//     init_line_config(-1.25, -0.5, -1.25, 0.5, 0.0, configs);   
//     init_drop_config(0.5, 0.7, 0.0, DPI, configs);    
//     init_drop_config(-1.3, -0.1, 0.0, DPI, configs);    
    init_drop_config(1.0, -1.0, 0.0, DPI, configs);    
//     init_drop_config(0.5, 0.5, -1.0, 1.0, configs);    
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
  
    glutSwapBuffers();   
  
  
    for (i=0; i<NPARTMAX; i++) 
    {
        color[i] = 0;
        newcolor[i] = 0;
        active[i] = 1;
    }
      
    if (RAINBOW_COLOR)      /* rainbow color scheme */
        for (i=0; i<NPART; i++) 
        {
            color[i] = (i*NCOLORS)/NPART;
            newcolor[i] = (i*NCOLORS)/NPART;  
        }

        sleep(SLEEP1);
    
  
    for (i=0; i<=NSTEPS; i++)
    {
        time = (double)i/(double)(NSTEPS-1);
        
//         x = 0.0;
//         if (time < 0.5) y = 0.999 - 0.798*2.0*time;
//         else y = 0.999 - 0.798 + 0.798*2.0*(time - 0.5);
//         alpha = 0.0;
//         y = 0.0;
//         alpha = (DPI/(double)NPOLY)*time;
//         alpha = PI*(0.5 + time)/(double)NPOLY;
//         alpha = time*DPI;
//         x =  MAZE_XSHIFT - 0.025;
//         y = 0.0;
//         x =  MAZE_XSHIFT - 1.2;
//         y = 0.05;
        
        x = -0.4999*LAMBDA;
        y = -0.4999*LAMBDA;
        
        alpha = angle_schedule(time);
        printf("Angle %.4lg\n", alpha*360.0/DPI); 

        period = compute_trajectories_xy(x, y, alpha, alpha + DPI, configs, trajectory, traj_length, &xmax); 
        
        if (TEST_INITIAL_COND) 
        {
            if ((xmax > XMAX)||((EXIT_BOTH_WAYS)&&(xmax < XMIN)))
            {
                traj_test_table[i] = 1;
                nsuccess++;
            }
            else traj_test_table[i] = 0;
//             traj_test_table[i] = test_initial_condition(configs, active, newcolor);
//             if (traj_test_table[i] >= 1) nsuccess++;
        }
        
//         printf("Period = %i\n", period);
        
        draw_trajectories(color, trajectory, traj_length);
        if (HEATMAP) draw_config_heatmap(trajectory, traj_length, configs, active, heatmap_number, heatmap_total, 0);
//         else draw_trajectories(color, trajectory, traj_length);
        
        print_parameters(i, traj_length, traj_length_table, period, alpha);
        if (DRAW_BILLIARD) draw_billiard();
        if (DRAW_LENGTHS_PLOT) plot_lengths(traj_length_table, traj_test_table, i+1);
        
        /* draw the shooter in red */
        rgb[0] = 1.0;   rgb[1] = 0.0;   rgb[2] = 0.0;
        draw_colored_circle(x, y, 0.01, NSEG, rgb);
//         draw_colored_circle(x_shooter, y_shooter, 0.01, NSEG, rgb);
        
        for (j=0; j<NPARTMAX; j++) color[j] = newcolor[j];
        
        if (!((NO_EXTRA_BUFFER_SWAP)&&(MOVIE))) glutSwapBuffers(); 
//         if (i == 0)  glutSwapBuffers();
        
	if (MOVIE) 
        {
            save_frame();
            
            if ((SLOW_AT_LONG_TRAJ)&&(traj_length[0] >= 10)) 
            {
                npause = (int)(0.5*sqrt((double)traj_length[0]));
                for (j=0; j<npause; j++) save_frame();
            }
            else if ((ADD_SUCCESS_GALLERY)&&(i>0)&&(traj_test_table[i]))
            {
//                 sleep(0.5);
                /* additional buffer swap seems to be necessary on some OS */
                glutSwapBuffers();   
                if (first == 0) 
                {
//                     gallery_counter = NSTEPS + END_FRAMES + success_counter*SUCCESS_GALLERY_FRAMES;
                    glutSwapBuffers(); 
                    first = 1;
                }
                else 
                {
                    for (j=0; j<SUCCESS_GALLERY_FRAMES; j++)
                    {
                        frame = NSTEPS + END_FRAMES + (success_counter-1)*SUCCESS_GALLERY_FRAMES + j;
                        if (HEATMAP) frame += (success_counter)*SUCCESS_GALLERY_FRAMES;
                        save_frame_counter(frame);
                    }
                    if (HEATMAP)
                    {
                        draw_trajectories(color, trajectory, traj_length);
                        print_parameters(i, traj_length, traj_length_table, period, alpha);
                        if (DRAW_BILLIARD) draw_billiard();
                        if (DRAW_LENGTHS_PLOT) plot_lengths(traj_length_table, traj_test_table, i+1);
                        rgb[0] = 1.0;   rgb[1] = 0.0;   rgb[2] = 0.0;
                        draw_colored_circle(x, y, 0.01, NSEG, rgb);
                        glutSwapBuffers();   
                        if (second == 0) 
                        {
                            glutSwapBuffers(); 
                            second = 1;
                        }
                        for (j=0; j<SUCCESS_GALLERY_FRAMES; j++)
                        {
                            frame = NSTEPS + END_FRAMES + (2*success_counter)*SUCCESS_GALLERY_FRAMES + j;
                            save_frame_counter(frame);
                        }
                    }
                    success_counter++;
                }
                
            }
            
            /* it seems that saving too many files too fast can cause trouble with the file system */
            /* so this is to make a pause from time to time - parameter PAUSE may need adjusting   */
            if (i % PAUSE == PAUSE - 1) 
            {
                printf("Making a short pause\n");
                sleep(PSLEEP);
                s = system("mv part*.tif tif_part/");
            }
        }
    }
 
    printf("1\n");
    
    if (MOVIE) 
    {
        printf("gallery_counter == %i, nend = %i\n", gallery_counter, nend); 
        if (gallery_counter - NSTEPS >= nend) nend = gallery_counter - NSTEPS;
        for (i=0; i<nend; i++) save_frame();
        s = system("mv part*.tif tif_part/");
    }
    
    printf("2\n");
    
    free(color);
    printf("3\n");
    
    free(newcolor);
    printf("4\n");
    
    free(active);
    printf("5\n");
    
    free(traj_length);
    printf("6\n");
    
    free(traj_length_table);
    printf("7\n");
    
    free(traj_test_table);
    printf("8\n");
    
    if (HEATMAP) 
    {
        free(heatmap_number);
        free(heatmap_total);
    }
    printf("9\n");
    for (i=0; i<NPARTMAX; i++) 
    {
        free(configs[i]);
        free(trajectory[i]);
    }
 
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


