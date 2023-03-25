/*********************************************************************************/
/*                                                                               */
/*  Animation of wave front in billiard                                          */
/*                                                                               */
/*  N. Berglund, december 2012, april 2021                                       */
/*                                                                               */
/*  Feel free to reuse, but if doing so it would be nice to drop a               */
/*  line to nils.berglund@univ-orleans.fr - Thanks!                              */
/*                                                                               */
/*  compile with                                                                 */
/*  gcc -o drop_billiard drop_billiard.c                                         */
/*  -O3 -L/usr/X11R6/lib -ltiff -lm -lGL -lGLU -lX11 -lXmu -lglut                */
/*                                                                               */
/*                                                                               */
/*  To make a video, set MOVIE to 1 and create subfolder tif_drop                */
/*  It may be possible to increase parameter PAUSE                               */
/*                                                                               */
/*  create movie using                                                           */
/*  ffmpeg -i part.%05d.tif -vcodec libx264 drop.mp4                             */
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

#define XMIN -2.0
#define XMAX 2.0	/* x interval */
#define YMIN -1.125
#define YMAX 1.125	/* y interval for 9/16 aspect ratio */

#define SCALING_FACTOR 1.0       /* scaling factor of drawing, needed for flower billiards, otherwise set to 1.0 */

/* Choice of the billiard table, see global_particles.c  */

#define B_DOMAIN 1      /* choice of domain shape */

#define CIRCLE_PATTERN 0    /* pattern of circles */
#define POLYLINE_PATTERN 1  /* pattern of polyline */

#define ABSORBING_CIRCLES 0 /* set to 1 for circular scatterers to be absorbing */

#define NMAXCIRCLES 1000        /* total number of circles (must be at least NCX*NCY for square grid) */
#define NMAXPOLY 1000        /* total number of sides of polygonal line */   
// #define NCX 10            /* number of circles in x direction */
// #define NCY 15            /* number of circles in y direction */
#define NCX 15            /* number of circles in x direction */
#define NCY 20            /* number of circles in y direction */
#define NPOISSON 500        /* number of points for Poisson C_RAND_POISSON arrangement */
#define NGOLDENSPIRAL 2000  /* max number of points for C_GOLDEN_SPIRAL arrandement */
#define SDEPTH 1            /* Sierpinski gastket depth */

#define LAMBDA 1.5	    /* parameter controlling the dimensions of domain */
// #define LAMBDA 1.124950941	/* sin(36°)/sin(31.5°) for 5-star shape with 45° angles */
// #define LAMBDA 1.445124904	/* sin(36°)/sin(24°) for 5-star shape with 60° angles */
// #define LAMBDA 3.75738973	/* sin(36°)/sin(9°) for 5-star shape with 90° angles */
// #define LAMBDA -1.73205080756888	/* -sqrt(3) for Reuleaux triangle */
// #define LAMBDA 1.73205080756888	/* sqrt(3) for triangle tiling plane */
#define MU 0.7              /* parameter controlling the dimensions of domain */
#define FOCI 1          /* set to 1 to draw focal points of ellipse */
#define NPOLY 4             /* number of sides of polygon */
#define APOLY 0.0           /* angle by which to turn polygon, in units of Pi/2 */ 
#define DRAW_BILLIARD 1     /* set to 1 to draw billiard */
#define DRAW_CONSTRUCTION_LINES 1   /* set to 1 to draw additional construction lines for billiard */
#define PERIODIC_BC 0       /* set to 1 to enforce periodic boundary conditions when drawing particles */
#define PENROSE_RATIO 2.5    /* parameter controlling the shape of small ellipses in Penrose room */

#define RESAMPLE 1      /* set to 1 if particles should be added when dispersion too large */

#define NPART 20000	/* number of particles */
#define NPARTMAX 100000	/* maximal number of particles after resampling */

#define NSTEPS 5250         /* number of frames of movie */
#define TIME 150             /* time between movie frames, for fluidity of real-time simulation */ 
#define DPHI 0.0001         /* integration step */
#define NVID 75             /* number of iterations between images displayed on screen */

/* Decreasing TIME accelerates the animation and the movie               */
/* For constant speed of movie, TIME*DPHI should be kept constant        */
/* However, increasing DPHI too much deterioriates quality of simulation */
/* For a good quality movie, take for instance TIME = 50, DPHI = 0.0002  */

/* simulation parameters */

#define LMAX 0.1       /* minimal segment length triggering resampling */ 
#define LPERIODIC 0.5   /* lines longer than this are not drawn (useful for Sinai billiard) */
#define LCUT 10000000.0        /* controls the max size of segments not considered as being cut */
#define DMIN 0.1       /* minimal distance to boundary for triggering resampling */ 
#define CYCLE 0         /* set to 1 for closed curve (start in all directions) */
#define ORDER_COLORS 0  /* set to 1 if colors should be drawn in order */ 

/* color and other graphical parameters */

#define COLOR_PALETTE 14     /* Color palette, see list in global_pdes.c  */

#define NCOLORS 16           /* number of colors */
#define COLORSHIFT 3        /* hue of initial color */ 
#define COLOR_HUEMIN 0   /* minimal color hue */
#define COLOR_HUEMAX 360 /* maximal color hue */
#define RAINBOW_COLOR 0     /* set to 1 to use different colors for all particles */
#define NSEG 100            /* number of segments of boundary */
#define BILLIARD_WIDTH 4    /* width of billiard */
#define FRONT_WIDTH 4       /* width of wave front */

#define BLACK 1             /* set to 1 for black background */
#define COLOR_OUTSIDE 0     /* set to 1 for colored outside */ 
#define OUTER_COLOR 300.0   /* color outside billiard */
#define PAINT_INT 0         /* set to 1 to paint interior in other color (for polygon) */
#define PAINT_EXT 1         /* set to 1 to paint exterior of billiard */


#define PAUSE 1000          /* number of frames after which to pause */
#define PSLEEP 1            /* sleep time during pause */
#define SLEEP1  1           /* initial sleeping time */
#define SLEEP2  100         /* final sleeping time */
#define END_FRAMES 25       /* number of frames at end of movie */

#define NXMAZE 8      /* width of maze */
#define NYMAZE 8      /* height of maze */
#define MAZE_MAX_NGBH 4     /* max number of neighbours of maze cell */
#define RAND_SHIFT 58       /* seed of random number generator */
#define MAZE_XSHIFT 0.0     /* horizontal shift of maze */
#define MAZE_RANDOM_FACTOR 0.1     /* randomization factor for S_MAZE_RANDOM */
#define MAZE_CORNER_RADIUS 0.5     /* radius of tounded corners in maze */
#define CLOSE_MAZE 0        /* set to 1 to close maze exits */



#define PI 	3.141592654
#define DPI 	6.283185307
#define PID 	1.570796327


#include "global_particles.c"
#include "sub_maze.c"
#include "sub_part_billiard.c"


/*********************/
/* animation part    */
/*********************/

void init_boundary_config(double smin, double smax, double anglemin, double anglemax, double *configs[NPARTMAX])
/* initialize configuration: drop on the boundary, beta version */
{
    int i;
    double ds, da, s, angle, alpha, pos[2], conf[2];
  
    if (anglemin <= 0.0) anglemin = PI/((double)NPART);
    if (anglemax >= PI) anglemax = PI*(1.0 - 1.0/((double)NPART));
    ds = (smax - smin)/((double)NPART);
    da = (anglemax - anglemin)/((double)NPART);
    for (i=0; i<NPART; i++) 
    {
        s = smin + ds*((double)i);
        angle = anglemin + da*((double)i);
        conf[0] = s;
        conf[1] = angle;
        
        pos_billiard(conf, pos, &alpha);

        vbilliard_xy(configs[i], alpha, pos);
    }
}

void init_drop_config(double x0, double y0, double angle1, double angle2, double *configs[NPARTMAX])   
/* initialize configuration: drop at (x0,y0) */
{
    int i;
    double dalpha, alpha, pos[2];
  
    while (angle2 < angle1) angle2 += DPI;
    dalpha = (angle2 - angle1)/((double)(NPART));
//     dalpha = (angle2 - angle1)/((double)(NPART-1));
    for (i=0; i<NPART; i++) 
    {
        alpha = angle1 + dalpha*((double)i);  
      
        pos[0] = x0;
        pos[1] = y0;
        vbilliard_xy(configs[i], alpha, pos);
    }
}
 
void init_partial_drop_config(int i1, int i2, double x0, double y0, double angle1, double angle2, double *configs[NPARTMAX])   
/* initialize configuration: drop at (x0,y0) */
{
    int i;
    double dalpha, alpha, pos[2];
  
    while (angle2 < angle1) angle2 += DPI;
    dalpha = (angle2 - angle1)/((double)(i2 - i1));
    if (i2 >= NPART) i2 = NPART;
    
    for (i=i1; i<i2; i++) 
    {
        alpha = angle1 + dalpha*((double)i);  
      
        pos[0] = x0;
        pos[1] = y0;
        vbilliard_xy(configs[i], alpha, pos);
    }
}
 


int resample(int color[NPARTMAX], double *configs[NPARTMAX])     
/* add particles where the front is stretched too thin */
{
    int len, i, j, k, iplus, newnparticles=nparticles, *newcolor;
    double dx, dy, pos[2], s1, s2, s, x, y, x1, y1, theta, alpha, beta, length2; 
    double *newconfigs[NPARTMAX];
    
    /* Since NPARTMAX can be big, it seemed wiser to use some memory allocation here */
    newcolor = malloc(sizeof(int)*(NPARTMAX));
    for (i=0; i<NPARTMAX; i++)
        newconfigs[i] = (double *)malloc(8*sizeof(double));

    
    printf("resampling, %i particles\n", nparticles);
    newnparticles=nparticles;
    j = 0;
    for (i=0; i<nparticles; i++)
    {
        iplus = i+1;
        if (iplus==nparticles) 
            if (CYCLE) iplus = 0;
            else iplus = nparticles - 1;
        for (k=0; k<8; k++) newconfigs[j][k] = configs[i][k];
        newcolor[j] = color[i];
        dx = configs[iplus][4] - configs[i][4];
        dy = configs[iplus][5] - configs[i][5];
        length2 = dx*dx + dy*dy;
        /* add a particle if length too big, but only if particles are of the same color, 
         and not too close to the boundary, to avoid problems due to roundoff */
        if ((color[i]==color[iplus])&&(length2 > LMAX*LMAX)&&(configs[i][2] < configs[i][3] - DMIN))
        {
//             print_config(configs[i]);
            if (newnparticles < NPARTMAX)
            {
                j++;
                newnparticles++;
//                 printf("Adding one point at %i, %i particles \n", j, newnparticles);
                newcolor[j] = color[i];
                s1 = configs[i][0];
                s2 = configs[iplus][0];
                s = 0.5*(s1 + s2);
                if (vabs(s-s1) > PID) s += PI;  /* needed if s1, s2 close to 0 and Pi */
                while (s<0) s += DPI;
                while (s>DPI) s -= DPI;
                x1 = LAMBDA*cos(s);
                y1 = sin(s);
                x = 0.5*(configs[i][4] + configs[iplus][4]);
                y = 0.5*(configs[i][5] + configs[iplus][5]);
                theta = argument(-LAMBDA*y1,x1/LAMBDA);
                alpha = argument(x1-x,y1-y);
                beta = theta-alpha;
                while (beta<0) beta += PI;
                while (beta>PI) beta -= PI;
                
                newconfigs[j][0] = s;
                newconfigs[j][1] = theta - alpha;
                newconfigs[j][2] = 0.5*(configs[i][2] + configs[iplus][2]);
                newconfigs[j][3] = module2(x-x1,y-y1);
                newconfigs[j][4] = x;
                newconfigs[j][5] = y;
                newconfigs[j][6] = x1;
                newconfigs[j][7] = y1;
//                 print_config(newconfigs[j]);
            }
        }
        j++;
    }
    
    if ((newnparticles > nparticles)&&(newnparticles < NPARTMAX))
    {
        for (i=0; i<newnparticles; i++)
        {
            for (k=0; k<8; k++) 
                configs[i][k] = newconfigs[i][k];
            color[i] = newcolor[i]; 
        }
    }
    
//     if (newnparticles == NPARTMAX) printf("Warning: Cannot add more particles\n");
    nparticles = newnparticles;
    
    free(newcolor);
    for (i=0; i<NPARTMAX; i++) free(newconfigs[i]);
    
    if (newnparticles == NPARTMAX) return(0); 
    else return(1);
}


void draw_config(int color[NPARTMAX], double *configs[NPARTMAX])
/* draw the particles */
{
    int i;
    double x1, y1, x2, y2, cosphi, sinphi, rgb[3], dist, dmax;

    glutSwapBuffers(); 
    blank();
    if (PAINT_INT) paint_billiard_interior();
      
    glLineWidth(FRONT_WIDTH);
    
    glEnable(GL_LINE_SMOOTH);
    if (CYCLE) glBegin (GL_LINE_LOOP);
    else glBegin(GL_LINE_STRIP);
    
    for (i=0; i<nparticles; i++)
    {        
        cosphi = (configs[i][6] - configs[i][4])/configs[i][3];
        sinphi = (configs[i][7] - configs[i][5])/configs[i][3];
        x2 = configs[i][4] + configs[i][2]*cosphi;
        y2 = configs[i][5] + configs[i][2]*sinphi;
        
        /* determine length of segment to avoid drawing too long segments */
        if (i>0) dist = module2(x2-x1,y2-y1);
        else dist = 0.0;
        
        dmax = DPI*((double)global_time)*DPHI/((double)nparticles);
        /* expected maximal distance between points for growing circle */
    
        rgb_color_scheme(color[i], rgb);
        glColor3d(rgb[0], rgb[1], rgb[2]);
        
        /* draw line only if it does not exceed LPERIODIC and 2*dmax */
        if ((xy_in_billiard(x2, y2))&&(dist < LPERIODIC)&&(dist < LCUT*dmax)) glVertex2d(x2, y2);
        else 
        {
            glEnd();
            glBegin (GL_LINE_STRIP);
        }
        
        if (configs[i][2] > configs[i][3] - DPHI) configs[i][2] -= configs[i][3];

        /* keep track of previous point to determine segment length */
        x1 = x2;
        y1 = y2;
    }
    glEnd ();
    draw_billiard(LAMBDA);    
}


void draw_ordered_config(int color[NPARTMAX], double *configs[NPARTMAX])
/* draw the wave front, one color after the other */
{
    int i, col;
    double x1, y1, x2, y2, cosphi, sinphi, rgb[3], dist, dmax;

    glutSwapBuffers(); 
    blank();
    if (PAINT_INT) paint_billiard_interior();
      
    glLineWidth(FRONT_WIDTH);
    
    glEnable(GL_LINE_SMOOTH);
    
    
    for (col=0; col<NCOLORS; col++)
    {
        glBegin(GL_LINE_STRIP);
        for (i=0; i<nparticles; i++)
        {
            if (color[i] == col)
            {        
                cosphi = (configs[i][6] - configs[i][4])/configs[i][3];
                sinphi = (configs[i][7] - configs[i][5])/configs[i][3];
                x2 = configs[i][4] + configs[i][2]*cosphi;
                y2 = configs[i][5] + configs[i][2]*sinphi;
        
                /* determine length of segment to avoid drawing too long segments */
                if (i>0) dist = module2(x2-x1,y2-y1);
                else dist = 0.0;
        
                dmax = DPI*((double)global_time)*DPHI/((double)nparticles);
                /* expected maximal distance between points for growing circle */
    
                rgb_color_scheme(color[i], rgb);
                glColor3d(rgb[0], rgb[1], rgb[2]);
        
                /* draw line only if it does not exceed LPERIODIC and 2*dmax */
                if ((i>0)&&(xy_in_billiard(x2, y2))&&(dist < LPERIODIC)&&(dist < LCUT*dmax)) glVertex2d(x2, y2);
                else 
                {
                    glEnd();
                    glBegin (GL_LINE_STRIP);
                }
        
                if (configs[i][2] > configs[i][3] - DPHI) configs[i][2] -= configs[i][3];

                /* keep track of previous point to determine segment length */
                x1 = x2;
                y1 = y2;
            }
        }
        glEnd ();
    }
    
    draw_billiard(LAMBDA);    
}


void graph_movie(int time, int color[NPARTMAX], double *configs[NPARTMAX])
/* compute next movie frame */
{
    int i, j;

    for (j=0; j<time; j++)
    {
        global_time++;
        for (i=0; i<nparticles; i++)
        {
//      print_config(configs[i]);
      
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
        
            if (configs[i][2] > configs[i][3] - DPHI) configs[i][2] -= configs[i][3];
        }
    }
}

void graph_no_movie(int time, int color[NPARTMAX], double *configs[NPARTMAX])
/* plot next image without making a movie */
{
    int i, j;

    for (j=0; j<time; j++)
    {
        global_time++;
        for (i=0; i<nparticles; i++)
        {
            if (configs[i][2]<0.0) 
            {    
                vbilliard(configs[i]);
//                 print_config(configs[i]);
                color[i]++;
                if (color[i] >= NCOLORS) color[i] -= NCOLORS;
            }

            configs[i][2] += DPHI; 
        
    
            if (configs[i][2] > configs[i][3] - DPHI) configs[i][2] -= configs[i][3];
        }
    }
}





void animation()
{
//     double time, dt;
    double *configs[NPARTMAX];
    int i, resamp = 1, s;
    int *color;
    
    /* Since NPARTMAX can be big, it seemed wiser to use some memory allocation here */
    color = malloc(sizeof(int)*(NPARTMAX));
    for (i=0; i<NPARTMAX; i++)
        configs[i] = (double *)malloc(8*sizeof(double));
  
    init_drop_config(0.0, 0.0, 0.0, DPI, configs);
//     init_partial_drop_config(0, NPART/4, LAMBDA, 0.0, 0.0, DPI, configs);
//     init_partial_drop_config(NPART/4, NPART/2, -LAMBDA, 0.0, 0.0, DPI, configs);
//     init_partial_drop_config(NPART/2, 3*NPART/4, 0.0, LAMBDA, 0.0, DPI, configs);
//     init_partial_drop_config(3*NPART/4, NPART, 0.0, -LAMBDA, 0.0, DPI, configs);
//     init_drop_config(-1.0 + 0.3*sqrt(2.0), -1.0 + 0.5*sqrt(2.0), 0.0, DPI, configs);
//     init_drop_config(-0.5, -0.5, 0.0, DPI, configs);
//     init_boundary_config(1.5, 1.5, 0.0, PI, configs);

  
    blank();  
    glColor3d(0.0, 0.0, 0.0);
    draw_billiard(LAMBDA);
    if (PAINT_INT) paint_billiard_interior();

  
    glutSwapBuffers();   
  
  
    for (i=0; i<NPARTMAX; i++) color[i] = 0;
    
    if (RAINBOW_COLOR)      /* rainbow color scheme */
        for (i=0; i<NPART; i++) color[i] = (i*NCOLORS)/NPART;
  
    sleep(SLEEP1);
  
    for (i=0; i<=NSTEPS; i++)
    {
	if (MOVIE) graph_movie(TIME, color, configs);
        else graph_no_movie(NVID, color, configs);
        
        if (ORDER_COLORS) draw_ordered_config(color, configs);
        else draw_config(color, configs);
        draw_billiard();

        
        /* for the ellipse, paths passing close to the foci are stronly divergent 
         * and the configurations may need to be resampled be adding extra points */ 
        if ((RESAMPLE)&&(i % 5 == 0)&&(nparticles < NPARTMAX)) resamp = resample(color, configs);
        if (!resamp) printf("Warning: Cannot add more particles\n");
        
	if (MOVIE) 
        {
            save_frame();
            
            /* it seems that saving too many files too fast can cause trouble with the file system */
            /* so this is to make a pause from time to time - parameter PAUSE may need adjusting   */
            if (i % PAUSE == PAUSE - 1) 
            {
                printf("Making a short pause\n");
                sleep(PSLEEP);
                s = system("mv part*.tif tif_drop/");
            }
        }
    }
 
    if (MOVIE) 
    {
        for (i=0; i<END_FRAMES; i++) save_frame();
        s = system("mv part*.tif tif_drop/");
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
    glutCreateWindow("Wave front in billiard");
       
    init();

    glutDisplayFunc(display);

    glutMainLoop();
  
    return 0;
}


