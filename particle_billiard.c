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

#define WINWIDTH 	1280  /* window width */
#define WINHEIGHT 	720   /* window height */

#define XMIN -2.0
#define XMAX 2.0	/* x interval */
#define YMIN -1.125
#define YMAX 1.125	/* y interval for 9/16 aspect ratio */

#define SCALING_FACTOR 1.0       /* scaling factor of drawing, needed for flower billiards, otherwise set to 1.0 */

// #define XMIN -1.8
// #define XMAX 1.8	/* x interval */
// #define YMIN -0.91
// #define YMAX 1.115	/* y interval for 9/16 aspect ratio */

/* Choice of the billiard table, see global_particles.c */

#define B_DOMAIN 14      /* choice of domain shape */

#define CIRCLE_PATTERN 0    /* pattern of circles */

#define NMAXCIRCLES 1000        /* total number of circles (must be at least NCX*NCY for square grid) */
// #define NCX 10            /* number of circles in x direction */
// #define NCY 15            /* number of circles in y direction */
#define NCX 15            /* number of circles in x direction */
#define NCY 20            /* number of circles in y direction */

#define LAMBDA 0.75	/* parameter controlling shape of billiard */
// #define LAMBDA -3.346065215	/* sin(60°)/sin(15°) for Reuleaux-type triangle with 90° angles */
// #define LAMBDA 3.0	/* parameter controlling shape of billiard */
// #define LAMBDA 0.6	/* parameter controlling shape of billiard */
// #define LAMBDA 0.4175295	/* sin(20°)/sin(55°) for 9-star shape with 30° angles */
// #define LAMBDA -1.949855824	/* 7-sided Reuleaux triangle */
// #define LAMBDA 3.75738973	/* sin(36°)/sin(9°) for 5-star shape with 90° angles */
// #define LAMBDA -1.73205080756888	/* -sqrt(3) for Reuleaux triangle */
// #define LAMBDA 1.73205080756888	/* sqrt(3) for triangle tiling plane */
#define MU 0.035          /* second parameter controlling shape of billiard */
#define FOCI 1          /* set to 1 to draw focal points of ellipse */
#define NPOLY 8             /* number of sides of polygon */
#define APOLY 0.25           /* angle by which to turn polygon, in units of Pi/2 */ 
#define DRAW_BILLIARD 1     /* set to 1 to draw billiard */
#define DRAW_CONSTRUCTION_LINES 1   /* set to 1 to draw additional construction lines for billiard */
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
#define TIME 1000       /* time between movie frames, for fluidity of real-time simulation */ 
// #define DPHI 0.000002    /* integration step */
// #define DPHI 0.00002    /* integration step */
#define DPHI 0.000005    /* integration step */
#define NVID 150         /* number of iterations between images displayed on screen */

/* Decreasing TIME accelerates the animation and the movie                               */
/* For constant speed of movie, TIME*DPHI should be kept constant                        */
/* However, increasing DPHI too much deterioriates quality of simulation                 */
/* NVID tells how often a picture is drawn in the animation, increase it for faster anim */
/* For a good quality movie, take for instance TIME = 400, DPHI = 0.00005, NVID = 100    */

/* Colors and other graphical parameters */

#define NCOLORS 16      /* number of colors */
#define COLORSHIFT 0     /* hue of initial color */ 
#define RAINBOW_COLOR 1  /* set to 1 to use different colors for all particles */
#define FLOWER_COLOR 0   /* set to 1 to adapt initial colors to flower billiard (tracks vs core) */
#define NSEG 100         /* number of segments of boundary */
#define LENGTH 0.02      /* length of velocity vectors */
#define BILLIARD_WIDTH 2    /* width of billiard */
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


#include "global_particles.c"
#include "sub_part_billiard.c"


/*********************/
/* animation part    */
/*********************/

void init_boundary_config(smin, smax, anglemin, anglemax, configs)
/* initialize configuration: drop on the boundary, beta version */
/* WORKS FOR ELLIPSE, HAS TO BE ADAPTED TO GENERAL BILLIARD */
double smin, smax, anglemin, anglemax;
double *configs[NPARTMAX];
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

void init_drop_config(x0, y0, angle1, angle2, configs)   /* initialize configuration: drop at (x0,y0) */
double x0, y0, angle1, angle2;
double *configs[NPARTMAX];
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
 
void init_sym_drop_config(x0, y0, angle1, angle2, configs)   
/* initialize configuration with two symmetric partial drops */
double x0, y0, angle1, angle2;
double *configs[NPARTMAX];
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
 
void init_line_config(x0, y0, x1, y1, angle, configs)   /* initialize configuration: line (x0,y0)-(x1,y1) in direction alpha */
double x0, y0, x1, y1, angle;
double *configs[NPARTMAX];
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


void draw_config(color, configs, active)
/* draw the particles */
int color[NPARTMAX], active[NPARTMAX];
double *configs[NPARTMAX];
{
    int i;
    double x1, y1, x2, y2, cosphi, sinphi, rgb[3];

    glutSwapBuffers(); 
    blank();
    if (PAINT_INT) paint_billiard_interior();
      
    glLineWidth(PARTICLE_WIDTH);
    
    glEnable(GL_LINE_SMOOTH);

    for (i=0; i<nparticles; i++)
    {
        if (configs[i][2]<0.0) 
        {    
            vbilliard(configs[i]);
            if (!RAINBOW_COLOR)
            {
                color[i]++;
                if (color[i] >= NCOLORS) color[i] -= NCOLORS;
            }
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
}


void graph_movie(time, color, configs, active)
/* compute next movie frame */
int time, color[NPARTMAX], active[NPARTMAX];
double *configs[NPARTMAX];
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
//                 if (c>=0) color[i]++;
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


void init_circle_config()
{
    int i, j, n; 
    double dx, dy;
    
    switch (CIRCLE_PATTERN) {
        case (C_FOUR_CIRCLES):
        {
            ncircles = 4;
            
            circlex[0] = 1.0;
            circley[0] = 0.0;
            circlerad[0] = 0.8;
            
            circlex[1] = -1.0;
            circley[1] = 0.0;
            circlerad[1] = 0.8;
            
            circlex[2] = 0.0;
            circley[2] = 0.8;
            circlerad[2] = 0.4;
            
            circlex[3] = 0.0;
            circley[3] = -0.8;
            circlerad[3] = 0.4;
            
            for (i=0; i<4; i++) circleactive[i] = 1;

            break;
        }
        case (C_SQUARE):
        {
            ncircles = NCX*NCY;
            dy = (YMAX - YMIN)/((double)NCY);
            for (i = 0; i < NCX; i++)
                for (j = 0; j < NCY; j++)
                {
                    n = NCY*i + j;
                    circlex[n] = ((double)(i-NCX/2) + 0.5)*dy;
                    circley[n] = YMIN + ((double)j + 0.5)*dy;
                    circlerad[n] = MU;
                    circleactive[n] = 1;
                }
            break;
        }
        case (C_HEX):
        {
            ncircles = NCX*(NCY+1);
            dy = (YMAX - YMIN)/((double)NCY);
            dx = dy*0.5*sqrt(3.0);
            for (i = 0; i < NCX; i++)
                for (j = 0; j < NCY+1; j++)
                {
                    n = (NCY+1)*i + j;
                    circlex[n] = ((double)(i-NCX/2) + 0.5)*dy;
                    circley[n] = YMIN + ((double)j - 0.5)*dy;
                    if ((i+NCX)%2 == 1) circley[n] += 0.5*dy;
                    circlerad[n] = MU;
                    circleactive[n] = 1;
                }
            break;
        }
        default: 
        {
            printf("Function init_circle_config not defined for this pattern \n");
        }
    }
}


void animation()
{
    double time, dt, alpha, r;
    double *configs[NPARTMAX];
    int i, j, resamp = 1, s, i1, i2;
    int *color, *newcolor, *active;
    
    /* Since NPARTMAX can be big, it seemed wiser to use some memory allocation here */
    color = malloc(sizeof(int)*(NPARTMAX));
    newcolor = malloc(sizeof(int)*(NPARTMAX));
    active = malloc(sizeof(int)*(NPARTMAX));
    for (i=0; i<NPARTMAX; i++)
        configs[i] = (double *)malloc(8*sizeof(double));
    
    /* init circle configuration if the domain is D_CIRCLES */
    if (B_DOMAIN == D_CIRCLES) init_circle_config();
    
      
    /* initialize system by putting particles in a given point with a range of velocities */
    r = cos(PI/(double)NPOLY)/cos(DPI/(double)NPOLY);

//     init_drop_config(0.0, 0.0, 0.0, PI, configs);    
//     init_drop_config(0.5, 0.5, -1.0, 1.0, configs);    
//     init_sym_drop_config(-1.0, 0.5, -PID, PID, configs);
//     init_drop_config(-0.999, 0.0, -alpha, alpha, configs);

//  other possible initial conditions :
//     init_line_config(-1.25, -0.5, -1.25, 0.5, 0.0, configs);
    init_line_config(0.0, -1.0, -1.0, 1.0, 0.25*PID, configs);
//     init_line_config(-0.7, -0.45, -0.7, 0.45, 0.0, configs);
//     init_line_config(-1.5, 0.1, -0.1, 1.0, -0.5*PID, configs);
  
    blank();  
    glColor3f(0.0, 0.0, 0.0);
    if (DRAW_BILLIARD) draw_billiard();
  
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
  
    for (i=0; i<=NSTEPS; i++)
    {
        graph_movie(TIME, newcolor, configs, active);
        
        draw_config(newcolor, configs, active);
        if (DRAW_BILLIARD) draw_billiard();
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
    glutSwapBuffers();
    blank();
    glutSwapBuffers();

    animation();        

    sleep(SLEEP2); 

    glPopMatrix();
}


int main(int argc, char** argv)
{
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH);
    glutInitWindowSize(WINWIDTH,WINHEIGHT);
    glutCreateWindow("Billiard animation");
       
    init();

    glutDisplayFunc(display);

    glutMainLoop();
  
    return 0;
}


