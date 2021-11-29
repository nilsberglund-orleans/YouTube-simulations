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

#define MOVIE 1         /* set to 1 to generate movie */

// #define WINWIDTH 	1280  /* window width */
// #define WINHEIGHT 	720   /* window height */

#define WINWIDTH 	1920  /* window width */
#define WINHEIGHT 	1080   /* window height */

#define XMIN -2.0
#define XMAX 2.0	/* x interval */
#define YMIN -1.125
#define YMAX 1.125	/* y interval for 9/16 aspect ratio */

#define SCALING_FACTOR 1.0       /* scaling factor of drawing, needed for flower billiards, otherwise set to 1.0 */

/* Choice of the billiard table, see global_particles.c */

#define B_DOMAIN 30     /* choice of domain shape */

#define CIRCLE_PATTERN 1    /* pattern of circles */
#define POLYLINE_PATTERN 1  /* pattern of polyline */

#define ABSORBING_CIRCLES 1 /* set to 1 for circular scatterers to be absorbing */

#define NMAXCIRCLES 100000     /* total number of circles (must be at least NCX*NCY for square grid) */
#define NMAXPOLY 100000        /* total number of sides of polygonal line */   
// #define NCX 10            /* number of circles in x direction */
// #define NCY 10            /* number of circles in y direction */
#define NCX 30            /* number of circles in x direction */
#define NCY 20            /* number of circles in y direction */
#define NPOISSON 500        /* number of points for Poisson C_RAND_POISSON arrangement */
#define NGOLDENSPIRAL 2000  /* max number of points for C_GOLDEN_SPIRAL arrandement */
#define SDEPTH 1            /* Sierpinski gastket depth */

#define LAMBDA 1.8	/* parameter controlling shape of domain */
#define MU 0.01          /* second parameter controlling shape of billiard */
// #define LAMBDA 0.3	/* parameter controlling shape of domain */
// #define MU 0.7          /* second parameter controlling shape of billiard */
#define FOCI 1          /* set to 1 to draw focal points of ellipse */
#define NPOLY 6             /* number of sides of polygon */
#define APOLY 0.0           /* angle by which to turn polygon, in units of Pi/2 */ 
#define DRAW_BILLIARD 1     /* set to 1 to draw billiard */
#define DRAW_CONSTRUCTION_LINES 0   /* set to 1 to draw additional construction lines for billiard */
#define PERIODIC_BC 0       /* set to 1 to enforce periodic boundary conditions when drawing particles */

#define RESAMPLE 0      /* set to 1 if particles should be added when dispersion too large */
#define DEBUG 0         /* draw trajectories, for debugging purposes */

/* Simulation parameters */

#define NPART 5000        /* number of particles */
#define NPARTMAX 100000	/* maximal number of particles after resampling */
#define LMAX 0.01       /* minimal segment length triggering resampling */ 
#define DMIN 0.02       /* minimal distance to boundary for triggering resampling */ 
#define CYCLE 1         /* set to 1 for closed curve (start in all directions) */
#define SHOWTRAILS 0    /* set to 1 to keep trails of the particles */
#define SHOWZOOM 1      /* set to 1 to show zoom on specific area */
#define PRINT_PARTICLE_NUMBER 0 /* set to 1 to print number of particles */
#define TEST_ACTIVE 1   /* set to 1 to test whether particle is in billiard */

#define NSTEPS 5000      /* number of frames of movie */
#define TIME 400         /* time between movie frames, for fluidity of real-time simulation */ 
// #define DPHI 0.00001    /* integration step */
#define DPHI 0.00005    /* integration step */
#define NVID 150         /* number of iterations between images displayed on screen */

/* Decreasing TIME accelerates the animation and the movie                               */
/* For constant speed of movie, TIME*DPHI should be kept constant                        */
/* However, increasing DPHI too much deterioriates quality of simulation                 */
/* NVID tells how often a picture is drawn in the animation, increase it for faster anim */
/* For a good quality movie, take for instance TIME = 400, DPHI = 0.00005, NVID = 100    */

/* Colors and other graphical parameters */

#define COLOR_PALETTE 10     /* Color palette, see list in global_pdes.c  */

#define NCOLORS 16       /* number of colors */
#define COLORSHIFT 0     /* hue of initial color */ 
#define RAINBOW_COLOR 1  /* set to 1 to use different colors for all particles */
#define FLOWER_COLOR 0   /* set to 1 to adapt initial colors to flower billiard (tracks vs core) */
#define NSEG 100         /* number of segments of boundary */
#define LENGTH 0.025        /* length of velocity vectors */
#define BILLIARD_WIDTH 3    /* width of billiard */
#define PARTICLE_WIDTH 2    /* width of particles */
#define FRONT_WIDTH 3       /* width of wave front */

#define BLACK 1             /* set to 1 for black background */
#define COLOR_OUTSIDE 0     /* set to 1 for colored outside */ 
#define OUTER_COLOR 270.0   /* color outside billiard */
#define PAINT_INT 0         /* set to 1 to paint interior in other color (for polygon/Reuleaux) */
#define PAINT_EXT 0         /* set to 1 to paint exterior */

#define PAUSE 1000       /* number of frames after which to pause */
#define PSLEEP 1         /* sleep time during pause */
#define SLEEP1  1        /* initial sleeping time */
#define SLEEP2  1000      /* final sleeping time */


#include "global_particles.c"
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

void draw_zoom(int color[NPARTMAX], double *configs[NPARTMAX], int active[NPARTMAX], double x_target, double y_target, double width)
/* draw zoom around target (for laser in room of mirrors) */
{
    int i;
    double x1, y1, x2, y2, xb, yb, cosphi, sinphi, rgb[3], shiftx = 0.0, shifty = 0.65, tradius, phi, zoomwidth = 0.4;
    
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
    
    glColor3f(0.0, 0.8, 0.2);
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
            rgb_color_scheme(color[i], rgb);
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
    
    if (SHOWZOOM) draw_zoom(color, configs, active, 0.95, 0.0, 0.1);
}


void draw_config(int color[NPARTMAX], double *configs[NPARTMAX], int active[NPARTMAX])
/* draw the particles */
{
    int i, c;
    double x1, y1, x2, y2, cosphi, sinphi, rgb[3];

    glutSwapBuffers(); 
    if (!SHOWTRAILS) blank();
    if (PAINT_INT) paint_billiard_interior();
      
    glLineWidth(PARTICLE_WIDTH);
    
    glEnable(GL_LINE_SMOOTH);

    for (i=0; i<nparticles; i++)
    {
        if (configs[i][2]<0.0) 
        {    
            c = vbilliard(configs[i]);
            if (!RAINBOW_COLOR)
            {
                color[i]++;
                if (color[i] >= NCOLORS) color[i] -= NCOLORS;
            }
            if ((ABSORBING_CIRCLES)&&(c < 0)) active[i] = 0;
        }

        configs[i][2] += DPHI; 
        
        cosphi = (configs[i][6] - configs[i][4])/configs[i][3];
        sinphi = (configs[i][7] - configs[i][5])/configs[i][3];
        x1 = configs[i][4] + configs[i][2]*cosphi;
        y1 = configs[i][5] + configs[i][2]*sinphi;
        x2 = configs[i][4] + (configs[i][2] + LENGTH)*cosphi;
        y2 = configs[i][5] + (configs[i][2] + LENGTH)*sinphi;
        
        /* test whether particle does not escape billiard */
        if (active[i]) active[i] = xy_in_billiard(x1, y1);
        
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
    
        if (configs[i][2] > configs[i][3] - DPHI) configs[i][2] -= configs[i][3];
    }
    if (DRAW_BILLIARD) draw_billiard();    
    
    if (SHOWZOOM) draw_zoom(color, configs, active, 0.95, 0.0, 0.1);
}


void graph_movie(int time, int color[NPARTMAX], double *configs[NPARTMAX], int active[NPARTMAX])
/* compute next movie frame */
{
    int i, j, c;

    for (j=0; j<time; j++)
    {
        for (i=0; i<nparticles; i++)
        {      
            if (configs[i][2]<0.0) 
            {    
//                 printf("reflecting particle %i\n", i);
                c = vbilliard(configs[i]);
                if ((ABSORBING_CIRCLES)&&(c < 0)) active[i] = 0;
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



void animation()
{
    double time, dt, alpha, r, rgb[3];
    double *configs[NPARTMAX];
    int i, j, resamp = 1, s, i1, i2;
    int *color, *newcolor, *active;
//     t_circle *circles;      /* experimental */
    
    /* Since NPARTMAX can be big, it seemed wiser to use some memory allocation here */
    color = malloc(sizeof(int)*(NPARTMAX));
    newcolor = malloc(sizeof(int)*(NPARTMAX));
    active = malloc(sizeof(int)*(NPARTMAX));
//     circles = malloc(sizeof(t_circle)*(NMAXCIRCLES));      /* experimental */
    
    for (i=0; i<NPARTMAX; i++)
        configs[i] = (double *)malloc(8*sizeof(double));

    /* init circle configuration if the domain is D_CIRCLES */
    if ((B_DOMAIN == D_CIRCLES)||(B_DOMAIN == D_CIRCLES_IN_RECT)||(B_DOMAIN == D_CIRCLES_IN_GENUSN)
        ||(B_DOMAIN == D_CIRCLES_IN_TORUS)) init_circles(circles);
    
    else if (B_DOMAIN == D_POLYLINE) init_polyline(polyline, circles);
      
    /* initialize system by putting particles in a given point with a range of velocities */
    r = cos(PI/(double)NPOLY)/cos(DPI/(double)NPOLY);
    
//     init_partial_drop_config(LAMBDA, 0.0, 0.0, DPI, 0, NPART/4, 0, configs, color, newcolor);
//     init_partial_drop_config(-LAMBDA, 0.0, 0.0, DPI, NPART/4, NPART/2, 0, configs, color, newcolor);
//     init_partial_drop_config(0.0, LAMBDA, 0.0, DPI, NPART/2, 3*NPART/4, 0, configs, color, newcolor);
//     init_partial_drop_config(0.0, -LAMBDA, 0.0, DPI, 3*NPART/4, NPART, 0, configs, color, newcolor);
    
//     init_drop_config(-1.0 + 0.3*sqrt(2.0), -1.0 + 0.5*sqrt(2.0), 0.0, DPI, configs);

//     init_line_config(-1.25, -0.5, -1.25, 0.5, 0.0, configs);   
//     init_drop_config(0.0, 0.0, -PI, PI, configs);    
    init_drop_config(-0.95, 0.0, PI, 3.0*PI, configs);
//     init_drop_config(-1.3, -0.1, 0.0, DPI, configs);    
//     init_drop_config(1.4, 0.1, 0.0, DPI, configs);    
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
    if (PRINT_PARTICLE_NUMBER) print_part_number(configs, active, XMIN + 0.1, YMIN + 0.1);
  
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

    sleep(SLEEP1);
    
    
    /* initialize drops in different colors */
//     init_partial_drop_config(0.0, 0.0, 0.0, DPI, 0, 2*NPART/5, 0, configs, color, newcolor);
//     init_partial_drop_config(0.0, 0.8, 0.0, DPI, 2*NPART/5, 4*NPART/5, 10, configs, color, newcolor);
//     init_partial_drop_config(1.2, 0.1, 0.0, DPI, 4*NPART/5, NPART, 36, configs, color, newcolor);
  
    for (i=0; i<=NSTEPS; i++)
    {
        graph_movie(TIME, newcolor, configs, active);
        
        if (SHOWTRAILS) draw_config_showtrails(newcolor, configs, active);
        else draw_config(newcolor, configs, active);
//         draw_config(newcolor, configs, active);
        if (DRAW_BILLIARD) draw_billiard();
        if (PRINT_PARTICLE_NUMBER) print_part_number(configs, active, XMIN + 0.1, YMIN + 0.1);
        for (j=0; j<NPARTMAX; j++) color[j] = newcolor[j];
        
        /* draw initial points */
//         draw_initial_condition_circle(0.0, 0.0, 0.02, 0);
//         draw_initial_condition_circle(0.0, 0.8, 0.02, 10);
//         draw_initial_condition_circle(1.2, 0.1, 0.02, 36);
        
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
        for (i=0; i<20; i++) save_frame();
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

    animation();        

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


