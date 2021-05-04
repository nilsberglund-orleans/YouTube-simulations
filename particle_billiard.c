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

/* Choice of the billiard table */

#define B_DOMAIN 3      /* choice of domain shape */

#define D_RECTANGLE 0   /* rectangular domain */
#define D_ELLIPSE 1     /* elliptical domain */
#define D_STADIUM 2     /* stadium-shaped domain */
#define D_SINAI 3       /* Sinai billiard */
#define D_DIAMOND 4     /* diamond-shaped billiard */
#define D_TRIANGLE 5    /* triangular billiard */

#define LAMBDA 0.5	/* parameter controlling shape of billiard */
// #define LAMBDA 1.73205080756888	/* sqrt(3) for triangle tiling plane */
#define FOCI 1          /* set to 1 to draw focal points of ellipse */

#define RESAMPLE 0      /* set to 1 if particles should be added when dispersion too large */

#define NPART 5000      /* number of particles */
#define NPARTMAX 50000	/* maximal number of particles after resampling */
#define LMAX 0.01       /* minimal segment length triggering resampling */ 
#define DMIN 0.02       /* minimal distance to boundary for triggering resampling */ 
#define CYCLE 1         /* set to 1 for closed curve (start in all directions) */

#define NSTEPS 5000     /* number of frames of movie */
#define TIME 1500       /* time between movie frames, for fluidity of real-time simulation */ 
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
/*        
        vbilliard(configs[i], alpha, pos);*/
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
    dalpha = (angle2 - angle1)/((double)(NPART-1));
    for (i=0; i<NPART; i++) 
    {
        alpha = angle1 + dalpha*((double)i);  
        
//         printf("alpha=%.5lg\n", alpha);
      
        pos[0] = x0;
        pos[1] = y0;
        vbilliard_xy(configs[i], alpha, pos);
    }
}
 
void init_line_config(x0, y0, x1, y1, angle, configs)   /* initialize configuration: line (x0,y0)-(x1,y1) in direction alpha */
double x0, y0, x1, y1, angle;
double *configs[NPARTMAX];
{
    int i;
    double dx, dy;
    double conf[2], pos[2];
  
    dx = (x1-x0)/((double)(NPART-1));
    dy = (y1-y0)/((double)(NPART-1));
    for (i=0; i<NPART; i++) 
    {
        pos[0] = x0 + ((double)i)*dx;
        pos[1] = y0 + ((double)i)*dy;
        vbilliard_xy(configs[i], angle, pos);
    }
}


void draw_config(color, configs)
/* draw the particles */
int color[NPARTMAX];
double *configs[NPARTMAX];
{
    int i;
    double x1, y1, x2, y2, cosphi, sinphi, rgb[3];

    glutSwapBuffers(); 
    blank();
      
    glLineWidth(3.0);
    
    glEnable(GL_LINE_SMOOTH);

    for (i=0; i<nparticles; i++)
    {
        if (configs[i][2]<0.0) 
        {    
            vbilliard(configs[i]);
            color[i]++;
            if (color[i] >= NCOLORS) color[i] -= NCOLORS;
        }

        configs[i][2] += DPHI; 
        
        cosphi = (configs[i][6] - configs[i][4])/configs[i][3];
        sinphi = (configs[i][7] - configs[i][5])/configs[i][3];
        x1 = configs[i][4] + configs[i][2]*cosphi;
        y1 = configs[i][5] + configs[i][2]*sinphi;
        x2 = configs[i][4] + (configs[i][2] + LENGTH)*cosphi;
        y2 = configs[i][5] + (configs[i][2] + LENGTH)*sinphi;
        
        rgb_color_scheme(color[i], rgb);
        glColor3f(rgb[0], rgb[1], rgb[2]);
        
        glBegin(GL_LINE_STRIP);
        glVertex2d(x1, y1);
        glVertex2d(x2, y2);
        glEnd ();
        
        /* taking care of boundary conditions - only needed for periodic boundary conditions */
        if (x2 > XMAX)
        {
            glBegin(GL_LINE_STRIP);
            glVertex2d(x1+XMIN-XMAX, y1);
            glVertex2d(x2+XMIN-XMAX, y2);
            glEnd ();
        }
        
        if (x2 < XMIN)
        {
            glBegin(GL_LINE_STRIP);
            glVertex2d(x1-XMIN+XMAX, y1);
            glVertex2d(x2-XMIN+XMAX, y2);
            glEnd ();
        }

        if (y2 > YMAX)
        {
            glBegin(GL_LINE_STRIP);
            glVertex2d(x1, y1+YMIN-YMAX);
            glVertex2d(x2, y2+YMIN-YMAX);
            glEnd ();
        }

        if (y2 < YMIN)
        {
            glBegin(GL_LINE_STRIP);
            glVertex2d(x1, y1+YMAX-YMIN);
            glVertex2d(x2, y2+YMAX-YMIN);
            glEnd ();
        }

        
        /* for debugging purpose */
//         glLineWidth(1.0);
//         glBegin(GL_LINES);
//         glVertex2d(configs[i][4], configs[i][5]);
//         glVertex2d(configs[i][6], configs[i][7]);
//         glEnd ();
//         glLineWidth(3.0);
    
        if (configs[i][2] > configs[i][3] - DPHI) configs[i][2] -= configs[i][3];
    }
    draw_billiard(LAMBDA);    
}


void graph_movie(time, color, configs)
/* compute next movie frame */
int time, color[NPARTMAX];
double *configs[NPARTMAX];
{
    int i, j, c;

    for (j=0; j<time; j++)
    {
        for (i=0; i<nparticles; i++)
        {
//      print_config(configs[i]);
      
            if (configs[i][2]<0.0) 
            {    
                c = vbilliard(configs[i]);
                if (c>=0) color[i]++;
                if (color[i] >= NCOLORS) color[i] -= NCOLORS;
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

void graph_no_movie(time, color, configs)
/* plot next image without making a movie */
int time, color[NPARTMAX];
double *configs[NPARTMAX];
{
    int i, j, c;
    
    for (j=0; j<time; j++)
    {        
        for (i=0; i<nparticles; i++)
        {
//      print_config(configs[i]);
      
            if (configs[i][2]<0.0) 
            {    
                c = vbilliard(configs[i]);
                if (c>=0) color[i]++;
                if (color[i] >= NCOLORS) color[i] -= NCOLORS;
            }

            configs[i][2] += DPHI; 
        
    
            if (configs[i][2] > configs[i][3] - DPHI) configs[i][2] -= configs[i][3];
        }
    }
}





void animation()
{
    double time, dt;
    double *configs[NPARTMAX];
    int i, j, resamp = 1, s;
    int *color, *newcolor;
    
    /* Since NPARTMAX can be big, it seemed wiser to use some memory allocation here */
    color = malloc(sizeof(int)*(NPARTMAX));
    newcolor = malloc(sizeof(int)*(NPARTMAX));
    for (i=0; i<NPARTMAX; i++)
        configs[i] = (double *)malloc(8*sizeof(double));
  
    /* initialize system by putting particles in a given point with a range of velocities */
    init_drop_config(-1.5, 0.3, -0.1, 0.1, configs);

//  other possible initial conditions :
//     init_line_config(0.0, -0.05, 0.0, 0.05, 0.0, configs);
  
    blank();  
    glColor3f(0.0, 0.0, 0.0);
    draw_billiard(LAMBDA);
  
    glutSwapBuffers();   
  
  
    for (i=0; i<NPARTMAX; i++) 
    {
        color[i] = 0;
        newcolor[i] = 0;
    }
  
    sleep(SLEEP1);
  
    for (i=0; i<=NSTEPS; i++)
    {
	if (MOVIE) graph_movie(TIME, newcolor, configs);
        else graph_no_movie(NVID, newcolor, configs);
        
        draw_config(color, configs);
        draw_billiard(color, configs);
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


