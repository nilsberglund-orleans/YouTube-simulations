/*********************************************************************************/
/*                                                                               */
/*  Animation of particles in billiard                                           */
/*                                                                               */
/*  N. Berglund, august 2022                                                     */
/*                                                                               */
/*  Feel free to reuse, but if doing so it would be nice to drop a               */
/*  line to nils.berglund@univ-orleans.fr - Thanks!                              */
/*                                                                               */
/*  compile with                                                                 */
/*  gcc -o billiard_phase_space billiard_phase_space.c                           */
/*  -O3 -L/usr/X11R6/lib -ltiff -lm -lGL -lGLU -lX11 -lXmu -lglut                */
/*                                                                               */
/*  OMP acceleration may be more effective after executing                       */
/*  export OMP_NUM_THREADS=2 in the shell before running the program             */
/*                                                                               */
/*  To make a video, set MOVIE to 1 and create subfolder tif_part                */
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
#define SAVE_MEMORY 1       /* set to 1 to save memory when writing tiff images */

#define WINWIDTH 	1280  /* window width */
#define WINHEIGHT 	720   /* window height */

#define XMIN -2.0
#define XMAX 2.0	/* x interval */
#define YMIN -1.125
#define YMAX 1.125	/* y interval for 9/16 aspect ratio */

#define XPHASEMAX 0.0   /* max x coordinate of phase portrait */

#define PIXELIZE 1      /* set to 1 to pixelize phase portrait (beta) */
#define NGRID 200       /* size of grid to draw orbits with positive Lyapunov exponent */

#define SYMMETRIZE_S 0  /* set to 1 to symmetrize orbits wrt s */

#define SCALING_FACTOR 1.0       /* scaling factor of drawing, needed for flower billiards, otherwise set to 1.0 */

/* Choice of the billiard table, see global_particles.c */

#define B_DOMAIN 1      /* choice of domain shape */

#define CIRCLE_PATTERN 6   /* pattern of circles */
#define POLYLINE_PATTERN 4  /* pattern of polyline */

#define ABSORBING_CIRCLES 1 /* set to 1 for circular scatterers to be absorbing */

#define NMAXCIRCLES 50000        /* total number of circles (must be at least NCX*NCY for square grid) */
#define NMAXPOLY 50000        /* total number of sides of polygonal line */   
#define NCX 9             /* number of circles in x direction */
#define NCY 20            /* number of circles in y direction */
#define NPOISSON 500        /* number of points for Poisson C_RAND_POISSON arrangement */
#define NGOLDENSPIRAL 2000  /* max number of points for C_GOLDEN_SPIRAL arrandement */
#define SDEPTH 2            /* Sierpinski gastket depth */

#define LAMBDAMIN 0.0	/* parameter controlling shape of domain */
#define LAMBDA 1.5	/* parameter controlling shape of domain */
#define MU 1.0          /* second parameter controlling shape of billiard */
#define FOCI 1          /* set to 1 to draw focal points of ellipse */
#define NPOLY 6             /* number of sides of polygon */
// #define NPOLY 3             /* number of sides of polygon */
#define APOLY -1.0           /* angle by which to turn polygon, in units of Pi/2 */ 
#define DRAW_BILLIARD 1     /* set to 1 to draw billiard */
#define DRAW_CONSTRUCTION_LINES 1   /* set to 1 to draw additional construction lines for billiard */
#define PERIODIC_BC 0       /* set to 1 to enforce periodic boundary conditions when drawing particles */
#define PENROSE_RATIO 1.0    /* parameter controlling the shape of small ellipses in Penrose room */

#define RESAMPLE 0      /* set to 1 if particles should be added when dispersion too large */
#define DEBUG 0         /* draw trajectories, for debugging purposes */

/* Simulation parameters */

#define NPART 1         /* number of particles */
#define NPARTMAX 100000	/* maximal number of particles after resampling */
#define TRAJ_LENGTH 10000 /* length of trajectory */
#define LMAX 0.01       /* minimal segment length triggering resampling */ 
#define DMIN 0.02       /* minimal distance to boundary for triggering resampling */ 
#define CYCLE 1         /* set to 1 for closed curve (start in all directions) */
#define SHOWTRAILS 0    /* set to 1 to keep trails of the particles */
#define SHOWZOOM 0      /* set to 1 to show zoom on specific area */
#define TEST_ACTIVE 1   /* set to 1 to test whether particle is in billiard */
#define PRINT_TRAJECTORY_LENGTH 0   /* set to 1 to print length of trajectory 0 */
#define PRINT_TRAJECTORY_PERIOD 0   /* set to 1 to print period of trajectory 0 */
#define DRAW_LENGTHS_PLOT 0         /* set to 1 to plot trajectory lengths */
#define LENGTHS_LOG_SCALE 0         /* set to 1 to use log scale for plot of lengths */
#define MAX_ANGLE 90.0         /* range of angles of trajectory */

#define NSTEPS 4000      /* number of frames of movie */
// #define NSTEPS 500       /* number of frames of movie */
#define TIME 2500        /* time between movie frames, for fluidity of real-time simulation */ 
#define DPHI 0.00001     /* integration step */
#define NVID 150         /* number of iterations between images displayed on screen */
#define SYMMETRIC_PARAMETER 0    /* set to 1 if parameters depend symmetrically on time */ 

/* Decreasing TIME accelerates the animation and the movie                               */
/* For constant speed of movie, TIME*DPHI should be kept constant                        */
/* However, increasing DPHI too much deterioriates quality of simulation                 */
/* NVID tells how often a picture is drawn in the animation, increase it for faster anim */
/* For a good quality movie, take for instance TIME = 400, DPHI = 0.00005, NVID = 100    */

/* Colors and other graphical parameters */

#define COLOR_PALETTE 11     /* Color palette, see list in global_pdes.c  */

#define NCOLORS 64       /* number of colors */
#define CFACTOR 3        /* color step */
#define COLORSHIFT 0     /* hue of initial color */ 
#define LYAP_PLOT_COLOR 100.0   /* color hue of Lyapunov exponent plot */
#define RAINBOW_COLOR 0  /* set to 1 to use different colors for all particles */
#define FLOWER_COLOR 0   /* set to 1 to adapt initial colors to flower billiard (tracks vs core) */
#define NSEG 100         /* number of segments of boundary */
#define LENGTH 0.03       /* length of velocity vectors */
#define BILLIARD_WIDTH 2    /* width of billiard */
#define PARTICLE_WIDTH 2    /* width of particles */
#define FRONT_WIDTH 3       /* width of wave front */

#define BLACK 1             /* set to 1 for black background */
#define COLOR_OUTSIDE 0     /* set to 1 for colored outside */ 
#define OUTER_COLOR 270.0   /* color outside billiard */
#define PAINT_INT 0         /* set to 1 to paint interior in other color (for polygon/Reuleaux) */
#define PAINT_EXT 1         /* set to 1 to paint exterior */

#define PRINT_LAMBDA 1      /* set to 1 to print value of lambda */
#define PRINT_MU 0          /* set to 1 to print value of mu */
#define PRINT_PENROSE_RATIO 0      /* set to 1 to print value of the Penrose billiard ratio */
#define PRINT_LYAPUNOV 1    /* set to 1 to print mean Lyapunov exponent */
#define PLOT_LYAPUNOV 1     /* set to 1 to add plot of Lyapunov exponents */
#define LOGSCALEX_LYAP 0    /* set to 1 to use log scale on parameter axis of Lyapunov exponents */
#define LYAP_MAX 1.0        /* maximal Lyapunov exponent */
#define ADAPT_TO_SYMMETRY 0 /* set to 1 to show only one symmetric part of phase space */
#define SYMMETRY_FACTOR 3   /* proportion of phase space to be shown */

#define PAUSE 1000       /* number of frames after which to pause */
#define PSLEEP 1         /* sleep time during pause */
#define SLEEP1  1        /* initial sleeping time */
#define SLEEP2  1        /* final sleeping time */
#define END_FRAMES  100   /* number of frames at end of movie */

#include "global_particles.c"
#include "sub_part_billiard_phasespace.c"


/*********************/
/* animation part    */
/*********************/

double print_parameters()
/* print billiard parameters */
{
    double x, xtext;
    char message[50];
    
    if (PRINT_LAMBDA) switch(B_DOMAIN) {
         case (D_ELLIPSE): 
        {
            x = sqrt(1.0 - 1.0/(lambda*lambda));
            sprintf(message, "Eccentricity %.3f", x);
            printf("Eccentricity %.3f\n", x);
            xtext = 0.75;
            break;
        }
        case (D_STADIUM): 
        {
            x = lambda;
            sprintf(message, "Linear part %.3f", x);
            printf("Linear part %.3f\n", x);
            xtext = 0.75;
            break;
        }
        case (D_REULEAUX): 
        {
            x = -lambda;
            sprintf(message, "Radius %.3f", x);
            printf("Radius %.3f\n", x);
            xtext = 0.75;
            break;
        }
        case (D_ALT_REU): 
        {
            x = lambda;
            sprintf(message, "Radius %.3f", x);
            printf("Radius %.3f\n", x);
            xtext = 0.75;
            break;
        }
        case (D_PARABOLAS): 
        {
            x = lambda;
            sprintf(message, "Focal distance %.3f", x);
            printf("Focal distance %.3f\n", x);
            xtext = 0.65;
            break;
        }
        case (D_PENROSE): 
        {
            x = lambda;
            sprintf(message, "Aspect ratio %.3f", x);
            printf("Aspect ratio %.3f\n", x);
            xtext = 0.65;
            break;
        }
        case (D_FLOWER): 
        {
            x = lambda;
            sprintf(message, "Parameter %.3f", x);
            printf("Aspect ratio %.3f\n", x);
            xtext = 0.65;
            break;
        }
        default: sprintf(message, " ");
    }
    
    if (PRINT_MU) switch(B_DOMAIN) {
        case (D_ANNULUS): 
        {
            x = mu;
            sprintf(message, "Distance to center %.3f", x);
            printf("Distance to center %.3f\n", x);
            xtext = 0.65;
            break;
        }
//         case (D_PARABOLAS): 
//         {
//             x = mu;
//             sprintf(message, "Radius %.3f", x);
//             printf("Radius %.3f\n", x);
//             xtext = 0.65;
//             break;
//         }
        default: sprintf(message, " ");
    }
    
    if (PRINT_PENROSE_RATIO) switch(B_DOMAIN) {
        case (D_PENROSE): 
        {
            x = penrose_ratio;
            sprintf(message, "Aspect ratio %.3f", x);
            printf("Aspect ratio %.3f\n", x);
            xtext = 0.65;
            break;
        }
        default: sprintf(message, " ");
    }
    
    glColor3f(1.0, 1.0, 1.0);
    glLineWidth(1);
    write_text(xtext, -0.25, message); 
    
}

double print_lyap_exponent(double lyap)
/* print Lyapunov exponent */
{
    char message[50];
    
    sprintf(message, "Average Lyapunov exponent %.3f", lyap);
    glColor3f(1.0, 1.0, 1.0);
    glLineWidth(1);
    write_text(0.45, -0.35, message); 
}

double s_range(double lambda, double mu)
/* return range of s values (boundary parametrization of billiard) */
{
    double psi0, omega, omega2, beta2, aa, bb, cc, ymax, width, l1, l2, co, so, axis1, axis2, phimax;
    
    switch(B_DOMAIN) {
        case (D_ELLIPSE): return(DPI);
        case (D_STADIUM): 
        {
            if (lambda > 0.0) return (DPI + 2.0*lambda);
            else 
            {
                psi0 = asin(-lambda/2);
                return(DPI-4.0*psi0);
            }
        }
        case (D_ANNULUS): return(2.0*DPI);
        case (D_REULEAUX):
        {
            omega2 = PI/((double)NPOLY);
            beta2 = asin(sin(omega2)/vabs(lambda));
            return(beta2*2.0*(double)NPOLY);
        }
        case (D_ALT_REU):
        {
            omega2 = PI/((double)NPOLY);
            beta2 = asin(sin(omega2)/vabs(lambda));
            return(beta2*2.0*(double)NPOLY);
        }
        case (D_PARABOLAS):
        {
            omega2 = PI/((double)NPOLY);
            aa = 0.25/mu;
            bb = 1.0/tan(omega2);
            cc = lambda + mu;
            ymax = ( - bb + sqrt(bb*bb + 4.0*aa*cc))/(2.0*aa);
            
            return((double)NPOLY*2.0*ymax);
        }
        case (D_PENROSE):
        {
            cc = sqrt(lambda*lambda - (1.0-mu)*(1.0-mu));
            width = 0.1*mu;
            l1 = mu - width;
            l2 = lambda - cc;
//             printf("l2 = %.3lg\n", l2);
            if (SYMMETRIZE_S) return(2.0*(PI + 2.0*l1 + l2));
            else return(4.0*(PI + 2.0*l1 + l2));
        }
        case (D_FLOWER):
        {
            compute_flower_parameters(&omega, &co, &so, &axis1, &axis2, &phimax);
            return(2.0*(double)NPOLY*phimax);
        }
        default: return(1.0);
    }
    
}

double lyapunov_exponent(double s, double u, double srange)
/* estimate Lyapunov exponent of orbit starting in (s, u) */
{
    double s1, u1, ds, du, nlyap = 0.0, delta = 0.01, config[8], config2[8], delta2;
    int i, n = 1000;
    
    config[0] = s;
    config[1] = u;
    vbilliard(config);
    s1 = s + delta;
    if (s1 > srange) s1 -= srange;
    u1 = u;
    config2[0] = s1;
    config2[1] = u1;
    vbilliard(config2);
    
    for (i=0; i<n; i++)
    {
        vbilliard(config);
        vbilliard(config2);
        ds = config2[0] - config[0];
        du = config2[1] - config[1];
        delta2 = module2(ds, du);
        
        nlyap += log(delta2/delta);
        config2[0] = config[0] + ds*delta/delta2;
        config2[1] = config[1] + du*delta/delta2;
    }
    return(nlyap/(double)n);
}

void draw_one_trajectory(double s, double u, int length, int draw)
/* draw a trajectory with initial condition (s, u) and length iterations */
{
    int i;
    double config[8], x1, y1, x2, y2, xw, yw;
    
    /* convert to window coordinates for use of GL_SCISSOR_TEST - should be improvable */
    xw = (int)((XPHASEMAX - XMIN)/(XMAX - XMIN)*(double)WINWIDTH);
    yw = (int)((-0.1 - YMIN)/(YMAX - YMIN)*(double)WINHEIGHT);
    
    /* draw the trajectories */
    glEnable(GL_SCISSOR_TEST);
    glScissor(xw, yw, WINWIDTH, WINHEIGHT);
    
    config[0] = s;
    config[1] = u;
    vbilliard(config);
    
    glPushMatrix();
    glTranslatef(1.0, 0.5, 0.0);
    glScalef(0.5, 0.5, 1.0);
    glBegin(GL_LINE_STRIP);
    for (i=0; i<length; i++)
    {
        vbilliard(config);
        x1 = scaling_factor*config[4];
        y1 = scaling_factor*config[5];
        if (x1 > XMIN) glVertex2d(x1, y1);
    }
    glEnd();
    
    if (draw) draw_billiard();
    
    glPopMatrix();
    glDisable(GL_SCISSOR_TEST);
}

void draw_trajectory(double s, double u, int length, double range)
{
    int i;
    
    if (ADAPT_TO_SYMMETRY) 
    {
        for (i=0; i<SYMMETRY_FACTOR; i++) 
            draw_one_trajectory(s + (double)i*range/(double)SYMMETRY_FACTOR, u, length, i == SYMMETRY_FACTOR-1);
    }
    else draw_one_trajectory(s, u, length, 1);
    
}

double adapt_to_symmetry(double x, double range)
{
    static double factor = (double)SYMMETRY_FACTOR;
//     printf("x = %.2lg\t", x);
    x *= factor;
    while (x > range) x -= range;
//     printf("x = %.2lg\n", x);
    return(x);
}

void draw_orbit(double s, double u, int length)
/* draw an orbit with initial condition (s, u) and length iterations */
{
    int i;
    double config[8], r, range, x;
    
    config[0] = s;
    config[1] = u;
    vbilliard(config);
    r = 0.001;
    range = s_range(lambda, mu);
    
    glEnable(GL_SCISSOR_TEST);
    glScissor(0, 0, WINWIDTH/2, WINHEIGHT);
    
    for (i=0; i<length; i++)
    {
        vbilliard(config);
        
        x = config[0];
//         while (x > range) x -= range;
//         while (x < 0.0) x += range;
        
        if (ADAPT_TO_SYMMETRY) x = adapt_to_symmetry(x, range);
        
//         x = x*3.0 - 3.0*range*(double)((int)x/range);
        
        draw_circle(XMIN + x*(XPHASEMAX-XMIN)/range, cos(config[1])*YMAX, r, NSEG);
    }
    
    glDisable(GL_SCISSOR_TEST);
}

void update_grid(double s, double u, int length, int *grid)
/* update pixelization grid */
{
    int i, is, iu, n;
    double config[8], range, rs, ru;
    
    config[0] = s;
    config[1] = u;
    vbilliard(config);
    range = s_range(lambda, mu);
    rs = (double)NGRID/range;
    ru = (double)NGRID/2.0;
    
    for (i=0; i<length; i++)
    {
        vbilliard(config);
        
        is = (int)(config[0]*rs);
        iu = (int)((cos(config[1]) + 1.0)*ru);
        
        n = iu*NGRID + is;
//         printf("conf = (%.2lg, %.2lg), is = %i, iu = %i, n = %i\n", config[0], config[1], is, iu, n);
        if ((n >= 0)&&(n < NGRID*NGRID)) grid[n]++;
    }
    
}

double draw_phase_portrait(int ns, int nu, int length, int trajlength, double lyapmax)
/* draw several orbits */
{
    int i, j, k, color, nlyap = 0;
    double s, u, range, ds, du, rgb[3], lyap, x, y, dx, dy, cratio, total_lyap = 0.0, mean_lyap, x1;
    int *grid; 
    double *lyap_exp;
    
    grid = (int *)malloc(NGRID*NGRID*sizeof(int));
    lyap_exp = (double *)malloc(ns*nu*sizeof(double));
    
    for (i=0; i<NGRID*NGRID; i++) grid[i] = 0;
    
    range = s_range(lambda, mu);
    
    ds = 1.0/(double)ns;
    du = 1.0/(double)nu;
    
    /* TODO */
    if (ADAPT_TO_SYMMETRY) ds *= 1.0/(double)SYMMETRY_FACTOR;
    
    print_parameters();
    
    /* compute Lyapunov exponents */
    #pragma omp parallel for private(i,j)
    for (i=0; i<ns; i++)
        for (j=0; j<nu; j++)
        {
            s = ((double)i + 0.5)*ds*range;
            u = acos(((double)j + 0.5)*du*2.0 - 1.0);
            lyap = lyapunov_exponent(s, u, range);
            if ((lyap > -1.0e20)&&(lyap < 1.0e20)) 
            {
                total_lyap += lyap;
                nlyap++;
            }
            lyap_exp[i + ns*j] = lyap;
        }

    /* draw orbits */
    #pragma omp parallel for private(i,j)
    for (i=0; i<ns; i++)
        for (j=0; j<nu; j++)
        {
            s = ((double)i + 0.5)*ds*range;
            u = acos(((double)j + 0.5)*du*2.0 - 1.0);
            lyap = lyap_exp[i + ns*j];
            color = (j*CFACTOR)%NCOLORS;
            rgb_color_scheme_lum(color, 0.5, rgb);
            if (lyap > lyapmax) 
            {
                cratio = lyapmax/lyap;
                for (k=0; k<3; k++) rgb[k]*=cratio;
            }
            glColor3f(rgb[0], rgb[1], rgb[2]);
                
            if ((PIXELIZE)&&(lyap >= lyapmax)) update_grid(s, u, 150*length, grid);

//             if ((i%1 == 0)&&(j%4 == 0)) draw_trajectory(s, u, trajlength);
            draw_orbit(s, u, length);
        }
        
    /* draw trajectories, large Lyapunov exponent */
    #pragma omp parallel for private(i,j)
    for (i=0; i<ns; i++)
        for (j=0; j<nu; j++) if (lyap_exp[i + ns*j] > lyapmax)
        {
            s = ((double)i + 0.5)*ds*range;
            u = acos(((double)j + 0.5)*du*2.0 - 1.0);
            color = (j*CFACTOR)%NCOLORS;
            rgb_color_scheme_lum(color, 0.5, rgb);
            cratio = lyapmax/lyap_exp[i + ns*j];
            if (cratio < 0.1) cratio = 0.1;
            for (k=0; k<3; k++) rgb[k]*=cratio;
            glColor3f(rgb[0], rgb[1], rgb[2]);
                
            if ((i%2 == 0)&&(j%4 == 0)) 
            {
                draw_trajectory(s, u, trajlength, range);
//                 if (SYMMETRIZE_S) 
//                 {
//                     glColor3f(rgb[0], rgb[1], rgb[2]);
//                     draw_trajectory(range + s, u, trajlength);
//                 }
            }
        }
    /* draw trajectories, small Lyapunov exponent */
    #pragma omp parallel for private(i,j)
    for (i=0; i<ns; i++)
        for (j=0; j<nu; j++) if (lyap_exp[i + ns*j] <= lyapmax)
        {
            s = ((double)i + 0.5)*ds*range;
            u = acos(((double)j + 0.5)*du*2.0 - 1.0);
            color = (j*CFACTOR)%NCOLORS;
            rgb_color_scheme_lum(color, 0.5, rgb);
            glColor3f(rgb[0], rgb[1], rgb[2]);
                
            if ((i%2 == 0)&&(j%4 == 0)) 
            {
                draw_trajectory(s, u, trajlength, range);
//                 if (SYMMETRIZE_S) 
//                 {
//                     glColor3f(rgb[0], rgb[1], rgb[2]);
//                     draw_trajectory(range + s, u, trajlength);
//                 }
            }
        }

    /* draw pixelized background */
    if (PIXELIZE)
    {
        dx = (XPHASEMAX -  XMIN)/(double)NGRID;
        dy = (YMAX - YMIN)/(double)NGRID;
    
        glBegin(GL_QUADS);
        #pragma omp parallel for private(i,j)
        for (i=0; i<NGRID; i++)
            for (j=0; j<NGRID; j++) if (grid[j*NGRID + i] > 0)
            {
                cratio = (double)grid[j*NGRID + i]/200.0;
                if (cratio > 1.0) cratio = 1.0;
                glColor3f(cratio, cratio, cratio);
        
                x = XMIN + (double)i*dx;
                x1 = x + dx;
                if (ADAPT_TO_SYMMETRY) 
                {
                    x = adapt_to_symmetry(x, range);
                    x1 = adapt_to_symmetry(x1, range);
                }
                y = YMIN + (double)j*dy;
                glVertex2d(x, y);
                glVertex2d(x1, y);
                glVertex2d(x1, y+dy);
                glVertex2d(x, y+dy);
            }
        glEnd();
    }
    
    /* redraw some orbits */
    #pragma omp parallel for private(i,j)
    for (i=0; i<ns; i++)
        for (j=0; j<nu; j++)
        {
            s = ((double)i + 0.5)*ds*range;
            u = acos(((double)j + 0.5)*du*2.0 - 1.0);
            lyap = lyap_exp[i + ns*j];
            color = (j*CFACTOR)%NCOLORS;
            rgb_color_scheme_lum(color, 0.5, rgb);
            if (lyap > lyapmax) 
            {
                cratio = lyapmax/lyap;
                for (k=0; k<3; k++) rgb[k]*=cratio;
            }
            glColor3f(rgb[0], rgb[1], rgb[2]);
                
            if ((!PIXELIZE)||(lyap < 2.0*lyapmax)) draw_orbit(s, u, length);
        }
        
    glColor3f(1.0, 1.0, 1.0);
    glLineWidth(1);
    draw_line(XPHASEMAX, YMIN, XPHASEMAX, YMAX);
    
    free(grid);
    free(lyap_exp);
    
    if (nlyap > 0) mean_lyap = total_lyap/(double)(nlyap);
    else mean_lyap = 0.0;
    
    if (PRINT_LYAPUNOV) print_lyap_exponent(mean_lyap);
    
    return(mean_lyap);
}


double lambda_schedule(double time)
{
    switch (B_DOMAIN){
        case (D_ELLIPSE): return(1.0 + 0.5*(LAMBDA-1.0)*(1.0 - cos(time*DPI)));
        case (D_STADIUM): 
        {
            if (time < 0.5) return(0.5*LAMBDA*(1.0 - cos(2.0*time*DPI)));
            else return(-0.5*(1.0 - cos(2.0*(time-0.5)*DPI)));
        }
//         case (D_REULEAUX): return(-0.0001 - exp(0.5*log(vabs(LAMBDA))*(1.0 - cos(time*DPI))));
//         case (D_REULEAUX): return(-1.00001 - 0.5*(LAMBDA - 1.0)*(1.0 - cos(time*DPI)));
        case (D_REULEAUX): 
        {
            if (LOGSCALEX_LYAP) return(-0.0001 - exp(0.5*log(vabs(LAMBDA))*(1.0 - cos(time*DPI))));
            else return(-1.00001 - 0.5*(LAMBDA - 1.0)*(1.0 - cos(time*PI)));
        }
        case (D_ALT_REU): 
        {
            if (LOGSCALEX_LYAP) return(0.0001 + exp(0.5*log(vabs(LAMBDA))*(1.0 - cos(time*DPI))));
            else return(1.00001 + 0.5*(LAMBDA - 1.0)*(1.0 - cos(time*PI)));
        }
        case (D_PARABOLAS): 
        {
            if (SYMMETRIC_PARAMETER) return((-1.0*cos(time*DPI)));
            else return((-1.0*cos(time*PI)));
        }
        case (D_PENROSE): 
        {
            if (PRINT_LAMBDA) return(0.5*(LAMBDA + LAMBDAMIN - (LAMBDA - LAMBDAMIN)*cos(time*DPI)));
            else return(LAMBDA);
        }
        case (D_FLOWER):
        {
            return(LAMBDAMIN + 0.5*(LAMBDA-LAMBDAMIN)*(1.0 - cos(time*DPI)));
        }
        default: return(LAMBDA);
    }
}

double mu_schedule(double time)
{
    switch (B_DOMAIN){
        case (D_ANNULUS): return(0.5*MU*(1.0 - cos(time*DPI)));
        case (D_PARABOLAS): 
        {
            if (SYMMETRIC_PARAMETER) return(MU*(1.0 + 1.0*cos(time*DPI)));
            else return(MU*(1.0 + 1.0*cos(time*PI)));
        }
        default: return(MU);
    }
}

double penrose_ratio_schedule(double time)
{
    switch (B_DOMAIN){
        case (D_PENROSE): 
        {
            if (PRINT_PENROSE_RATIO) return(PENROSE_RATIO*0.5*(1.0 - 0.999*cos(time*DPI)));
            else return(PENROSE_RATIO);
        }
        default: return(MU);
    }
}

double plot_coord(double x, double xmin, double xmax)
{
    return(xmin + x*(xmax - xmin));
}

double plot_coord_log(double x, double min, double xmin, double xmax)
{
    return(xmin + (1.0 - log(x)/log(min))*(xmax - xmin));
}

void plot_lyapunov_exponents_linscale(int i, double *lyap_exponents)
/* add plot of lyapunov exponents */
{
    int j, k, l, n1, n2, n3, n4, jmin, jmax;
    char message[100];
    static double xmin, xmax, ymin, ymax, xmid, ymid, dx, dy, plotxmin, plotxmax, plotymin, plotymax, lmin, lmax;
    double pos[2], x1, y1, x2, y2, rgb[3], x, y, lambda0, t, xshift;
    static int first = 1;
    
    if (first)
    {
        xmin = 0.1;
        xmax = XMAX - 0.1;
        ymin = YMIN + 0.05;
        ymax = YMIN + 0.75;
                
        xmid = 0.5*(xmin + xmax);
        ymid = 0.5*(ymin + ymax);
        
        dx = 0.5*(xmax - xmin);
        dy = 0.5*(ymax - ymin);
        
        plotxmin = xmin + 0.1;
        plotxmax = xmax - 0.1;
        plotymin = ymin + 0.07;
        plotymax = ymax - 0.15;
        
        if (PRINT_LAMBDA)
        {
            lmin = vabs(lambda_schedule(0.0));
//             lmax = lambda_schedule(0.5);
            lmax = vabs(lambda_schedule(1.0));
        }
        else if (PRINT_PENROSE_RATIO)
        {
            lmin = penrose_ratio_schedule(0.0);
            lmax = penrose_ratio_schedule(0.5);
        }
                
        printf("lmin = %.3lg, lmax = %.3lg\n", lmin, lmax);
        
        first = 0;
    }
    
    glColor3f(1.0, 1.0, 1.0);
    glLineWidth(2);
    
    /* axes and labels */
    draw_line(plotxmin, plotymin, plotxmax + 0.1, plotymin);
    draw_line(plotxmin, plotymin, plotxmin, plotymax + 0.1);
    
    switch (B_DOMAIN) {
        case (D_PARABOLAS) :
        {
            sprintf(message, "Focal dist");
            jmin = -10;
            jmax = 11;
            break;
        }
        case (D_REULEAUX) :
        {
            sprintf(message, "Radius");
            jmin = 1;
            jmax = (int)(LAMBDA) + 1;
            break;
        }
        case (D_ALT_REU) :
        {
            sprintf(message, "Radius");
            jmin = 1;
            jmax = (int)(LAMBDA) + 1;
            break;
        }
        case (D_PENROSE) : 
        {
            sprintf(message, "Aspect ratio");
            jmin = (int)(LAMBDAMIN*10.0);
            jmax = (int)(LAMBDA*10.0) + 1;
//             jmin = 0;
//             jmax = 10;
            break;
        }
        case (D_FLOWER) :
        {
            sprintf(message, "Parameter");
            jmin = (int)(LAMBDAMIN*10.0);
            jmax = (int)(LAMBDA*10.0) + 1;
//             jmin = 0;
//             jmax = 10;
            break;
        }
        default:
        {
            sprintf(message, "Parameter");
            jmin = 0;
            jmax = 10;
        }
    }
    
    write_text_fixedwidth(plotxmax - 0.15, plotymin + 0.075, message);
    
    /* graduations */
    for (j=jmin; j<jmax; j++) 
    {
//         lambda0 = 0.1*(double)j;
        lambda0 = (double)j;
        x = (lambda0 - lmin)/(lmax - lmin);
//         printf("Graduation %.1f x = %.3f\n", lambda0, x);
        x1 = plot_coord(x, plotxmin, plotxmax);      /* TO FIX */
//         printf("Graduation %.1f at %.3f\n", lambda0, x1);
        if (x1 < 0.0) x1 = -x1;
        if (lambda0 >= 0.0) xshift = -0.04;
        else xshift = -0.07;
        if (j%10 == 0)
        {
            draw_line(x1, plotymin - 0.02, x1, plotymin + 0.02);
            sprintf(message, "%.1f", lambda0);
            write_text_fixedwidth(x1 + xshift, plotymin - 0.075, message);
        }
        else
        {
            draw_line(x1, plotymin - 0.01, x1, plotymin + 0.01);
            if ((j+10)%2 == 0)
            {
                sprintf(message, "%.1f", lambda0);
                write_text_fixedwidth(x1 + xshift, plotymin - 0.075, message);
            }
        }
    }
    
    sprintf(message, "Avrg Lyap");
    write_text_fixedwidth(plotxmin - 0.05, plotymax + 0.15, message);
    
    for (j=1; j<=(int)(10.0*LYAP_MAX); j++)
    {
        y = (double)j/(10.0*LYAP_MAX);
        y1 = plot_coord(y, plotymin, plotymax);
        draw_line(plotxmin - 0.025, y1, plotxmin + 0.025, y1);
        sprintf(message, "%.1f", 0.1*(double)j);
        write_text_fixedwidth(plotxmin - 0.14, y1 - 0.015, message);
    }
    
    /* plot */
    hsl_to_rgb_palette(LYAP_PLOT_COLOR, 0.9, 0.5, rgb, COLOR_PALETTE);
    glColor3f(rgb[0], rgb[1], rgb[2]);
    x1 = plotxmin;
    y1 = plotymin;
    for (j=0; j<i; j++)
    {
        t = (double)j/(double)(NSTEPS-1);
        
        if (PRINT_LAMBDA) lambda0 = vabs(lambda_schedule(t));
        else if (PRINT_PENROSE_RATIO) lambda0 = penrose_ratio_schedule(t);
        
//         switch (B_DOMAIN){
//             case (D_PARABOLAS):
//             {
//                 lambda0 = lambda_schedule(t);
//                 break;
//             }
//             case (D_PENROSE):
//             {
//                 lambda0 = penrose_ratio_schedule(t);
//                 break;
//             }
//         }
        
//         x = (-lambda0 - 1.0)/(LAMBDA - 1.0);
        x = (lambda0 - lmin)/(lmax - lmin);
        y = lyap_exponents[j]/LYAP_MAX;
        
//         printf("lambda = %.3lg, x = %.3lg\n", lambda0, x);
                
        x2 = plot_coord(x, plotxmin, plotxmax);
        y2 = plot_coord(y, plotymin, plotymax);
        
        draw_line(x1, y1, x2, y2);
        x1 = x2;
        y1 = y2;
    }
}

void plot_lyapunov_exponents_logscale(int i, double *lyap_exponents)
/* add plot of lyapunov exponents */
{
    int j, k, l, n1, n2, n3, n4;
    char message[100];
    static double xmin, xmax, ymin, ymax, xmid, ymid, dx, dy, plotxmin, plotxmax, plotymin, plotymax;
    double pos[2], x1, y1, x2, y2, rgb[3], x, y, lambda0, t;
    static int first = 1;
    
    if (first)
    {
        xmin = 0.1;
        xmax = XMAX - 0.1;
        ymin = YMIN + 0.05;
        ymax = YMIN + 0.75;
                
        xmid = 0.5*(xmin + xmax);
        ymid = 0.5*(ymin + ymax);
        
        dx = 0.5*(xmax - xmin);
        dy = 0.5*(ymax - ymin);
        
        plotxmin = xmin + 0.1;
        plotxmax = xmax - 0.1;
        plotymin = ymin + 0.07;
        plotymax = ymax - 0.15;
        
        first = 0;
    }
    
    glColor3f(1.0, 1.0, 1.0);
    glLineWidth(2);
    
    /* axes and labels */
    draw_line(plotxmin, plotymin, plotxmax + 0.05, plotymin);
    draw_line(plotxmin, plotymin, plotxmin, plotymax + 0.1);
    
    sprintf(message, "Radius");
    write_text_fixedwidth(plotxmax, plotymin - 0.075, message);
    
//     n1 = (int)(log(LAMBDA)/log(10.0));
    n1 = (int)(log(LAMBDA)/log(10.0)) + 1;
    for (k=1; k<n1; k++)
    {
        n1 = (int)ipow(10, k);
        j = n1;
        lambda0 = (double)j;
        
        if (lambda0 < LAMBDA)
        {
            x = lambda0/LAMBDA;
            x1 = plot_coord_log(x, 1.0/LAMBDA, plotxmin, plotxmax);
        
            draw_line(x1, plotymin - 0.02, x1, plotymin + 0.02);
            sprintf(message, "%i", (int)lambda0);
            write_text_fixedwidth(x1 - 0.015 - 0.015*(double)(k-1), plotymin - 0.075, message);
        }
        
        for (l=1; l<10; l++)
        {
            j = (n1/10)*l;
            lambda0 = (double)j;
            if (lambda0 < LAMBDA)
            {
                x = lambda0/LAMBDA;
                x1 = plot_coord_log(x, 1.0/LAMBDA, plotxmin, plotxmax);
                draw_line(x1, plotymin - 0.01, x1, plotymin + 0.01);
                
                if ((l == 2)||(l == 5))
                {
                    sprintf(message, "%i", (int)lambda0);
                    write_text_fixedwidth(x1 - 0.015 - 0.015*(double)(k-1), plotymin - 0.075, message);
                }
            }
        }
    }
    
    sprintf(message, "Avrg Lyap");
    write_text_fixedwidth(plotxmin - 0.05, plotymax + 0.15, message);
    
//     for (j=1; j<=5; j++)
//     {
//         y = 0.2*(double)j;
    for (j=1; j<=(int)(10.0*LYAP_MAX); j++)
    {
        y = (double)j/(10.0*LYAP_MAX);
        y1 = plot_coord(y, plotymin, plotymax);
        draw_line(plotxmin - 0.025, y1, plotxmin + 0.025, y1);
        sprintf(message, "%.1f", 0.1*(double)j);
        write_text_fixedwidth(plotxmin - 0.14, y1 - 0.015, message);
    }
    
    /* plot */
    hsl_to_rgb_palette(LYAP_PLOT_COLOR, 0.9, 0.5, rgb, COLOR_PALETTE);
    glColor3f(rgb[0], rgb[1], rgb[2]);
    x1 = plotxmin;
    y1 = plotymin;
    for (j=0; j<i; j++)
    {
        t = (double)j/(double)(NSTEPS-1);
        lambda0 = lambda_schedule(t);
//         x = (-lambda0 - 1.0)/(LAMBDA - 1.0);
        x = -lambda0/LAMBDA;
        y = lyap_exponents[j]/LYAP_MAX;
                
        x2 = plot_coord_log(x, 1.0/LAMBDA, plotxmin, plotxmax);
        y2 = plot_coord(y, plotymin, plotymax);
        
        draw_line(x1, y1, x2, y2);
        x1 = x2;
        y1 = y2;
    }
}

void plot_lyapunov_exponents(int i, double *lyap_exponents)
{
    
    if (LOGSCALEX_LYAP) plot_lyapunov_exponents_logscale(i, lyap_exponents);
    else plot_lyapunov_exponents_linscale(i, lyap_exponents);
}

void animation()
{
    double time, dt, alpha, r, rgb[3], x, y, a, b, c, nsteps;
    double *lyap_exponents;
    int i, j, resamp = 1, s, i1, i2, period;
    char message[50];
    
    if (PLOT_LYAPUNOV) lyap_exponents = (double *)malloc(NSTEPS*sizeof(double));
    
    if (SYMMETRIC_PARAMETER) nsteps = NSTEPS/2;
    else nsteps = NSTEPS;
  
    for (i=0; i<=nsteps; i++)
    {
        time = (double)i/(double)(NSTEPS-1);
        
        blank();
        
        lambda = lambda_schedule(time);
        mu = mu_schedule(time);
        penrose_ratio = penrose_ratio_schedule(time);
        if (B_DOMAIN == D_FLOWER) scaling_factor = 0.7/lambda;
        else scaling_factor = 1.0;
        
//         lyap_exponents[i] = draw_phase_portrait(14, 20, 200, 20, 0.1);
//         lyap_exponents[i] = draw_phase_portrait(10, 20, 100, 20, 0.2);
//         lyap_exponents[i] = draw_phase_portrait(16, 20, 100, 20, 0.2);
        lyap_exponents[i] = draw_phase_portrait(6, 20, 200, 20, 0.3);
        
        if (PLOT_LYAPUNOV) plot_lyapunov_exponents(i, lyap_exponents);
        
        glutSwapBuffers();
        

        
        
	if (MOVIE) 
        {
            save_frame();
            
            if (SYMMETRIC_PARAMETER) save_frame_counter (NSTEPS + 1 - i);
            
            /* it seems that saving too many files too fast can cause trouble with the file system */
            /* so this is to make a pause from time to time - parameter PAUSE may need adjusting   */
            if (i % PAUSE == PAUSE - 1) 
            {
                printf("Making a short pause\n");
                sleep(PSLEEP);
                s = system("mv part*.tif tif_part/");
            }
        }
        
        if ((i==0)&&(SYMMETRIC_PARAMETER)) 
            for (j=0; j<END_FRAMES; j++) save_frame_counter(NSTEPS + j + 1);
    }
 
    if ((MOVIE)&&(!SYMMETRIC_PARAMETER)) 
        for (i=0; i<END_FRAMES; i++) save_frame_counter(NSTEPS + i + 1);
    
    if (MOVIE) s = system("mv part*.tif tif_part/");
    
    if (PLOT_LYAPUNOV) free(lyap_exponents);
     
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


