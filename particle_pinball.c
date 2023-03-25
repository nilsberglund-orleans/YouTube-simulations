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

#define MOVIE 0         /* set to 1 to generate movie */
#define SAVE_MEMORY 1           /* set to 1 to save memory when writing tiff images */

#define WINWIDTH 	1280  /* window width */
#define WINHEIGHT 	720   /* window height */

#define XMIN -4.0
#define XMAX 4.0	/* x interval */
#define YMIN -1.25
#define YMAX 3.25	/* y interval for 9/16 aspect ratio */
// #define XMIN -2.0
// #define XMAX 2.0	/* x interval */
// #define YMIN -1.125
// #define YMAX 1.125	/* y interval for 9/16 aspect ratio */

#define BOXYMIN -1.0
#define BOXYMAX 1.0     /* y dimensions of box (for circles in rectangle) */

#define SCALING_FACTOR 1.0       /* scaling factor of drawing, needed for flower billiards, otherwise set to 1.0 */

/* Choice of the billiard table, see global_particles.c */

#define B_DOMAIN 23      /* choice of domain shape */

#define CIRCLE_PATTERN 7    /* pattern of circles */
#define POLYLINE_PATTERN 1  /* pattern of polyline */
// #define CIRCLE_PATTERN 21    /* pattern of circles */

#define ABSORBING_CIRCLES 0 /* set to 1 for circular scatterers to be absorbing */

#define NMAXCIRCLES 5000        /* total number of circles (must be at least NCX*NCY for square grid) */
#define NMAXPOLY 1000        /* total number of sides of polygonal line */   
#define NCX 44            /* number of circles in x direction */
#define NCY 10            /* number of circles in y direction */
#define NPOISSON 350        /* number of points for Poisson C_RAND_POISSON arrangement */
#define NGOLDENSPIRAL 2000  /* max number of points for C_GOLDEN_SPIRAL arrandement */
#define SDEPTH 1            /* Sierpinski gastket depth */

#define LAMBDA 3.8     /* parameter controlling shape of domain */
// #define MU 0.1         /* second parameter controlling shape of billiard */
#define MU 0.07         /* second parameter controlling shape of billiard */
// #define MU 0.085         /* second parameter controlling shape of billiard */
// #define MU 0.09         /* second parameter controlling shape of billiard */
// #define MU 0.034         /* second parameter controlling shape of billiard */
#define FOCI 1          /* set to 1 to draw focal points of ellipse */
#define NPOLY 4             /* number of sides of polygon */
#define APOLY 0.5          /* angle by which to turn polygon, in units of Pi/2 */ 
#define DRAW_BILLIARD 1     /* set to 1 to draw billiard */
#define DRAW_CONSTRUCTION_LINES 0   /* set to 1 to draw additional construction lines for billiard */
#define PERIODIC_BC 0       /* set to 1 to enforce periodic boundary conditions when drawing particles */
#define PENROSE_RATIO 2.5    /* parameter controlling the shape of small ellipses in Penrose room */

#define RESAMPLE 0      /* set to 1 if particles should be added when dispersion too large */
#define DEBUG 0         /* draw trajectories, for debugging purposes */

/* Simulation parameters */

#define NPART 1         /* number of particles */
#define NPARTMAX 100000	/* maximal number of particles after resampling */
#define LMAX 0.01       /* minimal segment length triggering resampling */ 
#define DMIN 0.02       /* minimal distance to boundary for triggering resampling */ 
#define CYCLE 1         /* set to 1 for closed curve (start in all directions) */
#define SHOWTRAILS 0    /* set to 1 to keep trails of the particles */
#define TEST_ACTIVE 0   /* set to 1 to test whether particle is in billiard */

#define NSTEPS 10400     /* number of frames of movie */
// #define NSTEPS 1000     /* number of frames of movie */
// #define TIME 1500         /* time between movie frames, for fluidity of real-time simulation */ 
#define TIME 4000         /* time between movie frames, for fluidity of real-time simulation */ 
// #define DPHI 0.0001     /* integration step */
// #define DPHI 0.00002    /* integration step */
#define DPHI 0.00007     /* integration step */
// #define DPHI 0.000035     /* integration step */
#define NVID 150         /* number of iterations between images displayed on screen */

/* Decreasing TIME accelerates the animation and the movie                               */
/* For constant speed of movie, TIME*DPHI should be kept constant                        */
/* However, increasing DPHI too much deterioriates quality of simulation                 */
/* NVID tells how often a picture is drawn in the animation, increase it for faster anim */
/* For a good quality movie, take for instance TIME = 400, DPHI = 0.00005, NVID = 100    */

/* Colors and other graphical parameters */

#define COLOR_PALETTE 1     /* Color palette, see list in global_pdes.c  */

// #define NCOLORS 256       /* number of colors */
#define NCOLORS 48       /* number of colors */
#define COLORSHIFT 2     /* hue of initial color */ 
#define COLOR_HUEMIN 0   /* minimal color hue */
#define COLOR_HUEMAX 360 /* maximal color hue */
#define RAINBOW_COLOR 1  /* set to 1 to use different colors for all particles */
#define SINGLE_COLOR 1   /* set to 1 to make all particles a single color */
#define FLOWER_COLOR 0   /* set to 1 to adapt initial colors to flower billiard (tracks vs core) */
#define NSEG 100         /* number of segments of boundary */
#define LENGTH 0.1      /* length of velocity vectors */
#define BILLIARD_WIDTH 2    /* width of billiard */
#define PARTICLE_WIDTH 4    /* width of particles */
#define FRONT_WIDTH 3       /* width of wave front */
// #define COLOR_TRAJECTORY 0    /* hue for single color */
#define COLOR_TRAJECTORY 8    /* hue for single color */

#define BLACK 1             /* set to 1 for black background */
#define COLOR_OUTSIDE 0     /* set to 1 for colored outside */ 
#define OUTER_COLOR 270.0   /* color outside billiard */
#define PAINT_INT 0         /* set to 1 to paint interior in other color (for polygon/Reuleaux) */
#define PAINT_EXT 0         /* set to 1 to paint exterior in other color */
#define ERASE_OUTSIDE 1     /* set to 1 to erase outside of rectangular billiard (beta) */


#define PAUSE 500        /* number of frames after which to pause */
#define PSLEEP 5         /* sleep time during pause */
#define SLEEP1  1        /* initial sleeping time */
#define SLEEP2  1000     /* final sleeping time */
#define END_FRAMES 100   /* number of still frames at end of movie */

#define NXMAZE 8      /* width of maze */
#define NYMAZE 8      /* height of maze */
#define MAZE_MAX_NGBH 4     /* max number of neighbours of maze cell */
#define RAND_SHIFT 58       /* seed of random number generator */
#define MAZE_XSHIFT 0.0     /* horizontal shift of maze */
#define MAZE_RANDOM_FACTOR 0.1     /* randomization factor for S_MAZE_RANDOM */
#define MAZE_CORNER_RADIUS 0.5     /* radius of tounded corners in maze */
#define CLOSE_MAZE 0        /* set to 1 to close maze exits */

#define NPATHBINS 200     /* number of bins for path length histogramm */
#define PATHLMAX 1.8     /* max free path on graph */

#include "global_particles.c"
#include "sub_maze.c"
#include "sub_part_billiard.c"
#include "sub_part_pinball.c"

int ncol = 0, nobst = 0, nmaxpeg = 0;
int npath[NPATHBINS];
double max_free_path = 0.0;

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
            rgb_color_scheme(color[i], rgb);
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
}



void draw_config(int color[NPARTMAX], double *configs[NPARTMAX], int active[NPARTMAX])
/* draw the particles */
{
    int i;
    double x0, y0, x1, y1, x2, y2, cosphi, sinphi, rgb[3];

    glutSwapBuffers(); 
    if (!SHOWTRAILS) blank();
    if (PAINT_INT) paint_billiard_interior();
      
    glLineWidth(PARTICLE_WIDTH);
    
    glEnable(GL_LINE_SMOOTH);

    for (i=0; i<nparticles; i++)
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
        x1 = configs[i][4] + configs[i][2]*cosphi;
        y1 = configs[i][5] + configs[i][2]*sinphi;
        x2 = configs[i][4] + (configs[i][2] + LENGTH)*cosphi;
        y2 = configs[i][5] + (configs[i][2] + LENGTH)*sinphi;
        
        /* test whether particle does not escape billiard */
        if ((TEST_ACTIVE)&&(active[i])) active[i] = xy_in_billiard(x1, y1);
        
        if (active[i])  
        {
            rgb_color_scheme(color[i], rgb);
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
}


void graph_movie(int time, int color[NPARTMAX], double *configs[NPARTMAX], int active[NPARTMAX])
/* compute next movie frame */
{
    int i, j, k, c;
    double rgb[3];
    static double total_pathlength = 0.0;

    for (j=0; j<time; j++)
    {
        for (i=0; i<nparticles; i++)
        {      
            if (configs[i][2]<0.0) 
            {    
//                 printf("reflecting particle %i\n", i);
                if ((SHOWTRAILS)&&(active[i]))  
                {
                    glLineWidth(PARTICLE_WIDTH);
                    rgb_color_scheme(color[i], rgb);
                    glColor3f(rgb[0], rgb[1], rgb[2]);
        
                    glBegin(GL_LINE_STRIP);
                    glVertex2d(SCALING_FACTOR*configs[i][4], SCALING_FACTOR*configs[i][5]);
                    glVertex2d(SCALING_FACTOR*configs[i][6], SCALING_FACTOR*configs[i][7]);
                    glEnd ();
                }

                total_pathlength += configs[i][3];
                
                c = vbilliard(configs[i]);
                
                /* update number of collisions, not counting boundary for periodic b.c. */
                if (B_DOMAIN != D_CIRCLES_IN_TORUS) ncol++;
                else if (c >= 0) ncol++;
                
                if ((c >= 0)&&(circles[c].color == 0)) nobst++;
                circles[c].color++;
                
                /* take care of circles doubled because of periodic boundary conditions */ 
                if ((circles[c].active)&&(B_DOMAIN == D_CIRCLES_IN_TORUS)&&(circles[c].partner != c))                     circles[circles[c].partner].color++;
                    
                circles[c].new = 10;
                
                /* update free path statistics */
                if (ncol > 1) /* disregard very first collision */
                {
                    if (B_DOMAIN != D_CIRCLES_IN_TORUS)
                    {
                        k = (int)((double)NPATHBINS*configs[i][3]/PATHLMAX);
                        if (k < NPATHBINS) npath[k]++;
                        if (total_pathlength > max_free_path) max_free_path = total_pathlength;
                        total_pathlength = 0.0;
                    }
                    else    /* case with periodic boundary conditions */
                    {
                        if (c >= 0)     /* a circle is hit, update histogram */
                        {
//                             printf("total path length %.3lg\n", total_pathlength);
                            k = (int)((double)NPATHBINS*total_pathlength/PATHLMAX);
                            if (k < NPATHBINS) npath[k]++;
                            if (total_pathlength > max_free_path) max_free_path = total_pathlength;
                            total_pathlength = 0.0;
                        }
                    }
                }
                
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


void print_particle_numbers(double *configs[NPARTMAX])
{
    char message[50], message1[50];
    double cosphi, x1;
    static double rgb[3], xleft, xright;
    static short int first = 1;
    int i, nleft = 0, nmid = 0, nright = 0;
    
    rgb[0] = 0.0; rgb[1] = 0.0; rgb[2] = 0.0;

    if (first) /* compute box limits */
    {
        /* find leftmost and rightmost circle */
        for (i=0; i<ncircles; i++) 
            if ((circles[i].active)&&(circles[i].xc - circles[i].radius < xleft)) xleft = circles[i].xc - circles[i].radius; 
        for (i=0; i<ncircles; i++) 
            if ((circles[i].active)&&(circles[i].xc + circles[i].radius > xright)) xright = circles[i].xc + circles[i].radius; 
        
        first = 0;
        
        printf("xleft = %.3lg, xright = %.3lg", xleft, xright);
    }
    
    for (i=0; i<nparticles; i++)
    {
        cosphi = (configs[i][6] - configs[i][4])/configs[i][3];
        x1 = configs[i][4] + configs[i][2]*cosphi;
        
        if (x1 > xright) nright++;
        else if (x1 < xleft) nleft++;
        else nmid++; 
        
//         if (i == nparticles-1) printf("x1 = %.3lg, nleft = %i, nright = %i\n", x1, nleft, nright);
    }
    
    erase_area(XMIN + 0.31, YMIN + 0.07, 0.24, 0.05, rgb);
    sprintf(message, "%4d particles", nleft);
    glColor3f(1.0, 1.0, 1.0);
    write_text_fixedwidth(XMIN + 0.1, YMIN + 0.04, message);
    
    erase_area(0.0, YMIN + 0.07, 0.24, 0.05, rgb);
    sprintf(message, "%4d particles", nmid);
    glColor3f(1.0, 1.0, 1.0);
    write_text_fixedwidth(-0.21, YMIN + 0.04, message);    
    
    erase_area(XMAX - 0.29, YMIN + 0.07, 0.24, 0.05, rgb);
    sprintf(message, "%4d particles", nright);
    glColor3f(1.0, 1.0, 1.0);
    write_text_fixedwidth(XMAX - 0.5, YMIN + 0.04, message);
}

void draw_statistics()
{
    int i, n, colmax = 55, pegcollisions[90], nypegs = 70, meanpegs = 0, meansquarepegs = 0, total_coll = 0, ymax = 0, 
        meanbins = 0, meansquarebins, total_bin = 0, stephits = 10, tickstephits = 1;
    double x, y, yscale = 110.0, y0, dx, rgb[3], xshift, coll_mean, coll_stdv, path_mean, path_stdv, ebin, len_over_bin, 
        steppath = 0.5;
    char message[50];
    
    glLineWidth(1);
        
    y0 = 0.5*(YMAX + YMIN) + 0.2;
    dx = (XMAX-0.6)/(double)colmax;
    xshift = XMIN + 0.3;
    rgb[0] = 0.0; rgb[1] = 0.0; rgb[2] = 1.0;
    
    /* histogram of number of collisions per peg */
    for (i=0; i<colmax; i++) pegcollisions[i] = 0;
    
    for (i=0; i<ncircles; i++) if ((circles[i].active)&&(!circles[i].double_circle))
    {
        n = circles[i].color;
        if (n < colmax) pegcollisions[n]++;
    }
    for (i=1; i<colmax; i++) 
    {
        total_coll += pegcollisions[i];
        meanpegs += i*pegcollisions[i];
        meansquarepegs += i*i*pegcollisions[i];
    }
    
    for (i=1; i<colmax; i++)
    {
        x = xshift + (double)i*dx;
        y = y0 + (double)pegcollisions[i]*YMAX/yscale;
        
        rgb_color_scheme(i, rgb);
        erase_rectangle(x, y0, x+dx, y, rgb);
    }

    glColor3f(1.0, 1.0, 1.0);
    glBegin(GL_LINE_STRIP);
    /* histogram */
    for (i=0; i<colmax; i++)
    {
        x = xshift + (double)i*dx;
        y = y0 + (double)pegcollisions[i]*YMAX/yscale;
        glVertex2d(x, y0);
        glVertex2d(x, y);
        glVertex2d(x+dx, y);
        glVertex2d(x+dx, y0);
        glVertex2d(x, y0);
    }
    glVertex2d(xshift, y0);
    glVertex2d(xshift, y0 + (double)nypegs*YMAX/yscale);
    glEnd ();
    
    /* graduation and labels */
    for (i=10; i<nypegs; i+=10)
    {
        glBegin(GL_LINE_STRIP);
        glVertex2d(xshift - 0.025, y0 + (double)i*YMAX/yscale);
        glVertex2d(xshift + 0.025, y0 + (double)i*YMAX/yscale);
        glEnd ();
    }
    
    for (i=tickstephits; i<colmax; i+=tickstephits)
    {
        glBegin(GL_LINE_STRIP);
        glVertex2d(xshift + (double)i*dx, y0 - 0.025);
        glVertex2d(xshift + (double)i*dx, y0 + 0.025);
        glEnd ();
    }

    for (i=stephits; i<nypegs - stephits; i+=stephits)
    {
        sprintf(message, "%4d", i);
        write_text_fixedwidth(xshift + (double)i*dx - 0.15, y0 - 0.12, message);
    }

    for (i=10; i<nypegs; i+=10)
    {
        sprintf(message, "%4d", i);
        write_text_fixedwidth(xshift - 0.3, y0 - 0.025 + (double)i*YMAX/yscale, message);
    }
    
    coll_mean = (double)meanpegs/(double)total_coll;
    coll_stdv = sqrt((double)meansquarepegs/(double)total_coll - coll_mean*coll_mean);
    
    
    sprintf(message, "hits");
    write_text_fixedwidth(xshift + (double)(colmax-3)*dx, y0 - 0.12, message);
    sprintf(message, "pegs");
    write_text_fixedwidth(xshift - 0.25, y0 - 0.025 + (double)(nypegs - 3)*YMAX/yscale, message);
//     sprintf(message, "Mean %.4lg", (double)meanpegs/(double)total_coll);
    
    sprintf(message, "Max hits/peg %d", nmaxpeg);
    write_text(-1.5, YMAX - 0.3, message);
    
    sprintf(message, "Mean hits/peg %.4lg", coll_mean);
    write_text(-1.5, YMAX - 0.5, message);
    
    sprintf(message, "Stdv hits/peg %.4lg", coll_stdv);
    write_text(-1.5, YMAX - 0.7, message);
    
    
    /* histogram of path lengths */
    
    for (i=1; i<NPATHBINS; i++)
    {
        if (npath[i] > ymax) ymax = npath[i];
        total_bin += npath[i];
        meanbins += i*npath[i];
        meansquarebins += i*i*npath[i];
    }
    
    yscale = 0.9*(YMAX-y0)*(double)ymax;
    dx = (XMAX-0.6)/(double)NPATHBINS;
    xshift = 0.3;
    rgb[0] = 1.0; rgb[1] = 0.0; rgb[2] = 0.0;

    for (i=1; i<NPATHBINS; i++)
    {
        x = xshift + (double)i*dx;
        y = y0 + (double)npath[i]*YMAX/yscale;
        
        if (y > y0) erase_rectangle(x, y0, x+dx, y, rgb);
    }
    
    glColor3f(1.0, 1.0, 1.0);
    glBegin(GL_LINE_STRIP);
    /* histogram */
    for (i=1; i<NPATHBINS; i++) if (npath[i] > 0)
    {
        x = xshift + (double)i*dx;
        y = y0 + (double)npath[i]*YMAX/yscale;
        glVertex2d(x, y0);
        glVertex2d(x, y);
        glVertex2d(x+dx, y);
        glVertex2d(x+dx, y0);
        glVertex2d(x, y0);
    }
    glEnd ();
    
    glBegin(GL_LINE_STRIP);
    glVertex2d(xshift, YMAX - 0.1);
    glVertex2d(xshift, y0);
    glVertex2d(xshift + (double)NPATHBINS*dx, y0);
    glEnd ();
    
    
    for (x = steppath; x < PATHLMAX; x+=steppath)
    {
        i = (int)(x*(double)NPATHBINS/PATHLMAX);
        sprintf(message, "%.2f", x);
        write_text_fixedwidth(xshift + (double)i*dx - 0.1, y0 - 0.12, message);
    }

    for (x = steppath; x < PATHLMAX; x+=steppath)
    {
        i = (int)(x*(double)NPATHBINS/PATHLMAX);
        glBegin(GL_LINE_STRIP);
        glVertex2d(xshift + (double)i*dx, y0 - 0.025);
        glVertex2d(xshift + (double)i*dx, y0 + 0.025);
        glEnd ();
    }
    
    
    ebin = (double)meanbins/(double)total_bin;      /* mean bin */
    len_over_bin = PATHLMAX/(double)NPATHBINS;      /* conversion from bin to path length */
    path_mean = ebin*len_over_bin;                  /* mean free path */
    path_stdv = sqrt((double)meansquarebins/(double)total_bin - ebin*ebin)*len_over_bin;
    
    sprintf(message, "free path");
    write_text_fixedwidth(XMAX - 0.6, y0 - 0.12, message);
    
    sprintf(message, "Max free path %.4lg", max_free_path);
    write_text(2.2, YMAX - 0.3, message);
    
    sprintf(message, "Mean free path %.4lg", path_mean);
    write_text(2.2, YMAX - 0.5, message);
    
    sprintf(message, "Stdv free path %.4lg", path_stdv);
    write_text(2.2, YMAX - 0.7, message);
}


void animation(int circle_config)
{
    double time, dt, alpha, r;
    double *configs[NPARTMAX];
    int i, j, resamp = 1, s, i1, i2;
    int *color, *newcolor, *active;
    char message[50];
    
    /* Since NPARTMAX can be big, it seemed wiser to use some memory allocation here */
    color = malloc(sizeof(int)*(NPARTMAX));
    newcolor = malloc(sizeof(int)*(NPARTMAX));
    active = malloc(sizeof(int)*(NPARTMAX));
    for (i=0; i<NPARTMAX; i++)
        configs[i] = (double *)malloc(8*sizeof(double));
    
    /* init circle configuration if the domain is D_CIRCLES */
//     if ((B_DOMAIN == D_CIRCLES)||(B_DOMAIN == D_CIRCLES_IN_RECT)||(B_DOMAIN == D_CIRCLES_IN_GENUSN)
//         ||(B_DOMAIN == D_CIRCLES_IN_TORUS)) 
//         init_circle_config_pinball(circle_config);
    
    /* init circle configuration if the domain is D_CIRCLES */
    if ((B_DOMAIN == D_CIRCLES)||(B_DOMAIN == D_CIRCLES_IN_RECT)||(B_DOMAIN == D_CIRCLES_IN_GENUSN)
        ||(B_DOMAIN == D_CIRCLES_IN_TORUS)) init_circles_pinball(circle_config, circles);
    
    
    /* remove discs that are not in domain */
    if ((B_DOMAIN == D_CIRCLES_IN_RECT)||(B_DOMAIN == D_CIRCLES_IN_GENUSN))
//         ||(B_DOMAIN == D_CIRCLES_IN_TORUS)) 
        for (i=0; i<ncircles; i++)
        {
            if (vabs(circles[i].xc) + circles[i].radius > 0.99) circles[i].active = 0;
            if (vabs(circles[i].xc) + circles[i].radius > 0.99*LAMBDA) circles[i].active = 0;
        }
    else if (B_DOMAIN == D_CIRCLES_IN_TORUS)
        for (i=0; i<ncircles; i++)
        {
            if (vabs(circles[i].yc) - circles[i].radius > 1.0) circles[i].active = 0;
            if (vabs(circles[i].xc) - circles[i].radius > LAMBDA) circles[i].active = 0;
        }
        
//     for (i=0; i<ncircles; i++) 
//         printf("Circle %i at (%.2f, %.2f), double %i, partner %i\n", i, circlex[i], circley[i], double_circle[i], partner_circle[i]);
        
//             if (vabs(circley[i]) > 1.0) circleactive[i] = 0;
      
    /* initialize system by putting particles in a given point with a range of velocities */
    r = cos(PI/(double)NPOLY)/cos(DPI/(double)NPOLY);

//     init_drop_config(0.4, 0.0, PID, 0.5*PID, configs);    
    init_drop_config(0.0, 0.0, -0.45*PID, 0.45*PID, configs);    

//     init_line_config(-1.25, -0.5, -1.25, 0.5, 0.0, configs);   
//     init_drop_config(0.5, 0.1, -0.5*PID, 0.5*PID, configs);    
//     init_drop_config(-1.4, 0.0, -0.5*PID, 0.5*PID, configs);    
//     init_drop_config(0.5, 0.5, -1.0, 1.0, configs);    
//     init_sym_drop_config(-1.0, 0.5, -PID, PID, configs);
//     init_drop_config(-0.999, 0.0, -alpha, alpha, configs);

//  other possible initial conditions :
//     init_line_config(0.0, -0.5, 0.0, 0.5, 0.0, configs);
//     init_line_config(-1.25, -0.5, -1.25, 0.5, 0.0*PID, configs);
//     init_line_config(-1.0, -0.3, -1.0, 0.3, 0.0, configs);
//     init_line_config(-0.7, -0.45, -0.7, 0.45, 0.0, configs);
//     init_line_config(-1.5, 0.1, -0.1, 1.0, -0.5*PID, configs);
  
//     if (!SHOWTRAILS) blank();
    blank();
    glColor3f(0.0, 0.0, 0.0);
    if (DRAW_BILLIARD) draw_billiard();
    if (ERASE_OUTSIDE) erase_rectangle_outside(270.0, 0.1, 0.15);
//     print_particle_numbers(configs);
  
    glutSwapBuffers();   
    
//     if (MOVIE) 
//     {
//         for (i=0; i<20; i++) save_frame();
//         s = system("mv part*.tif tif_part/");
//     }
  
  
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
        
    if (SINGLE_COLOR)
        for (i=0; i<NPART; i++) 
        {
            color[i] = COLOR_TRAJECTORY;
            newcolor[i] = COLOR_TRAJECTORY;  
        }
        
    for (i=0; i<NPATHBINS; i++) npath[i] = 0;

    sleep(SLEEP1);
  
    for (i=0; i<=NSTEPS; i++)
    {
        graph_movie(TIME, newcolor, configs, active);
        
        if (SHOWTRAILS) draw_config_showtrails(newcolor, configs, active);
        else draw_config(newcolor, configs, active);
        if (DRAW_BILLIARD) draw_billiard();
        if (ERASE_OUTSIDE) erase_rectangle_outside(270.0, 0.1, 0.15);
//         print_particle_numbers(configs);
 
        sprintf(message, "%d collisions", ncol);
        glColor3f(1.0, 1.0, 1.0);
        write_text(XMIN + 0.6, YMIN + 0.08, message);
//         write_text(XMIN + 0.3, YMIN + 0.04, message);
  
        sprintf(message, "%d pegs hit", nobst);
        glColor3f(1.0, 1.0, 1.0);
        write_text(XMAX-1.4, YMIN + 0.08, message);
//         write_text(XMAX-0.7, YMIN + 0.04, message);
        
        /* count max number a peg is hit */
        nmaxpeg = 0;
        for (j=0; j<ncircles; j++)
            if (circles[j].color > nmaxpeg) nmaxpeg = circles[j].color;

//         sprintf(message, "max hits per peg: %d", nmaxpeg);
//         glColor3f(1.0, 1.0, 1.0);
//         write_text(-0.6, YMIN + 0.08, message);
//         write_text(-0.3, YMIN + 0.04, message);
        
        draw_statistics();

        for (j=0; j<NPARTMAX; j++) color[j] = newcolor[j];

        
	if (MOVIE) 
        {
            save_frame();
            
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
 
    if (MOVIE) 
    {
        for (i=0; i<END_FRAMES; i++) save_frame();
        printf("Making a short pause\n");
        sleep(PSLEEP);
        s = system("mv part*.tif tif_part/");
    }
    
    free(color);
    free(newcolor);
    for (i=0; i<NPARTMAX; i++) free(configs[i]);
 
}


void display(void)
{
    glPushMatrix();

    blank();
    
    if (!SHOWTRAILS)
    {
        glutSwapBuffers();
        blank();
        glutSwapBuffers();
    }

    animation(CIRCLE_PATTERN);
//     animation(C_TRI);        

//     animation(C_GOLDEN_SPIRAL);        
//     animation(C_SQUARE);        
//     animation(C_HEX);        
//     animation(C_GOLDEN_MEAN);        
//     animation(C_RAND_DISPLACED);        
//     animation(C_POISSON_DISC);        
//     animation(C_RAND_POISSON);
    
    sleep(SLEEP2); 

    glPopMatrix();
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


