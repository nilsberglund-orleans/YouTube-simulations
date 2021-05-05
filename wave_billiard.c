/*********************************************************************************/
/*                                                                               */
/*  Animation of wave equation in a planar domain                                */
/*                                                                               */
/*  N. Berglund, december 2012, april 2021                                       */
/*                                                                               */
/*  UPDATE 24/04: distinction between damping and "elasticity" parameters        */
/*  UPDATE 27/04: new billiard shapes, bug in color scheme fixed                 */
/*  UPDATE 28/04: code made more efficient, with help of Marco Mancini           */
/*                                                                               */
/*  Feel free to reuse, but if doing so it would be nice to drop a               */
/*  line to nils.berglund@univ-orleans.fr - Thanks!                              */
/*                                                                               */
/*  compile with                                                                 */
/*  gcc -o wave_billiard wave_billiard.c                                         */
/* -L/usr/X11R6/lib -ltiff -lm -lGL -lGLU -lX11 -lXmu -lglut -O3 -fopenmp        */
/*                                                                               */
/*  To make a video, set MOVIE to 1 and create subfolder tif_wave                */
/*  It may be possible to increase parameter PAUSE                               */
/*                                                                               */
/*  create movie using                                                           */
/*  ffmpeg -i wave.%05d.tif -vcodec libx264 wave.mp4                             */
/*                                                                               */
/*********************************************************************************/

/*********************************************************************************/
/*                                                                               */
/* NB: The algorithm used to simulate the wave equation is highly paralellizable */
/* One could make it much faster by using a GPU                                  */
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

#define B_DOMAIN 8      /* choice of domain shape */

#define D_RECTANGLE 0   /* rectangular domain */
#define D_ELLIPSE 1     /* elliptical domain */
#define D_STADIUM 2     /* stadium-shaped domain */
#define D_SINAI 3       /* Sinai billiard */
#define D_DIAMOND 4     /* diamond-shaped billiard */
#define D_TRIANGLE 5    /* triangular billiard */
#define D_FLAT 6        /* flat interface */
#define D_ANNULUS 7     /* annulus */
#define D_POLYGON 8     /* polygon */

#define LAMBDA 0.3	    /* parameter controlling the dimensions of domain */
#define NPOLY 6             /* number of sides of polygon */
#define APOLY 1.0           /* angle by which to turn polygon, in units of Pi/2 */ 
#define FOCI 1              /* set to 1 to draw focal points of ellipse */

/* You can add more billiard tables by adapting the functions */
/* xy_in_billiard and draw_billiard below */

/* Physical patameters of wave equation */

#define COURANT 0.01       /* Courant number */
#define GAMMA 0.0      /* damping factor in wave equation */
// #define GAMMA 5.0e-10      /* damping factor in wave equation */
#define KAPPA 5.0e-7       /* "elasticity" term enforcing oscillations */
// #define KAPPA 5.0e-9       /* "elasticity" term enforcing oscillations */
// #define KAPPA 5.0e-8       /* "elasticity" term enforcing oscillations */
/* The Courant number is given by c*DT/DX, where DT is the time step and DX the lattice spacing */
/* The physical damping coefficient is given by GAMMA/(DT)^2 */
/* Increasing COURANT speeds up the simulation, but decreases accuracy */
/* For similar wave forms, COURANT^2*GAMMA should be kept constant */

/* For debugging purposes only */
#define FLOOR 0         /* set to 1 to limit wave amplitude to VMAX */
#define VMAX 10.0       /* max value of wave amplitude */

/* Parameters for length and speed of simulation */

#define NSTEPS 4575     //7500      /* number of frames of movie */
#define NVID 50          /* number of iterations between images displayed on screen */
#define NSEG 100         /* number of segments of boundary */

#define PAUSE 1000         /* number of frames after which to pause */
#define PSLEEP 1         /* sleep time during pause */
#define SLEEP1  1        /* initial sleeping time */
#define SLEEP2  1   /* final sleeping time */

/* Color schemes */

#define COLOR_SCHEME 1   /* choice of color scheme */

#define C_LUM 0          /* color scheme modifies luminosity (with slow drift of hue) */
#define C_HUE 1          /* color scheme modifies hue */

#define SCALE 1          /* set to 1 to adjust color scheme to variance of field */
#define SLOPE 1.0        /* sensitivity of color on wave amplitude */
#define ATTENUATION 0.0  /* exponential attenuation coefficient of contrast with time */

#define COLORHUE 260     /* initial hue of water color for scheme C_LUM */
#define COLORDRIFT 0.0   /* how much the color hue drifts during the whole simulation */
#define LUMMEAN 0.5      /* amplitude of luminosity variation for scheme C_LUM */
#define LUMAMP 0.3       /* amplitude of luminosity variation for scheme C_LUM */
#define HUEMEAN 235.0    /* mean value of hue for color scheme C_HUE */
#define HUEAMP 60.0      /* amplitude of variation of hue for color scheme C_HUE */
// #define HUEMEAN 320.0    /* mean value of hue for color scheme C_HUE */
// #define HUEAMP 100.0      /* amplitude of variation of hue for color scheme C_HUE */

/* Basic math */

#define PI 	3.141592654
#define DPI 	6.283185307
#define PID 	1.570796327

double courant2;  /* Courant parameter squared */


/*********************/
/* Graphics routines */
/*********************/

int writetiff(char *filename, char *description, int x, int y, int width, int height, int compression)
{
  TIFF *file;
  GLubyte *image, *p;
  int i;

  file = TIFFOpen(filename, "w");
  if (file == NULL)
  {
    return 1;
  }

  image = (GLubyte *) malloc(width * height * sizeof(GLubyte) * 3);

  /* OpenGL's default 4 byte pack alignment would leave extra bytes at the
     end of each image row so that each full row contained a number of bytes
     divisible by 4.  Ie, an RGB row with 3 pixels and 8-bit componets would
     be laid out like "RGBRGBRGBxxx" where the last three "xxx" bytes exist
     just to pad the row out to 12 bytes (12 is divisible by 4). To make sure
     the rows are packed as tight as possible (no row padding), set the pack
     alignment to 1. */

  glPixelStorei(GL_PACK_ALIGNMENT, 1);

  glReadPixels(x, y, width, height, GL_RGB, GL_UNSIGNED_BYTE, image);
  TIFFSetField(file, TIFFTAG_IMAGEWIDTH, (uint32) width);
  TIFFSetField(file, TIFFTAG_IMAGELENGTH, (uint32) height);
  TIFFSetField(file, TIFFTAG_BITSPERSAMPLE, 8);
  TIFFSetField(file, TIFFTAG_COMPRESSION, compression);
  TIFFSetField(file, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_RGB);
  TIFFSetField(file, TIFFTAG_ORIENTATION, ORIENTATION_BOTLEFT);
  TIFFSetField(file, TIFFTAG_SAMPLESPERPIXEL, 3);
  TIFFSetField(file, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
  TIFFSetField(file, TIFFTAG_ROWSPERSTRIP, 1);
  TIFFSetField(file, TIFFTAG_IMAGEDESCRIPTION, description);
  p = image;
  for (i = height - 1; i >= 0; i--)
  {
//     if (TIFFWriteScanline(file, p, height - i - 1, 0) < 0)
    if (TIFFWriteScanline(file, p, i, 0) < 0)
    {
      free(image);
      TIFFClose(file);
      return 1;
    }
    p += width * sizeof(GLubyte) * 3;
  }
  TIFFClose(file);
  return 0;
}


void init()		/* initialisation of window */
{
    glLineWidth(3);

    glClearColor(0.0, 0.0, 0.0, 1.0);
    glClear(GL_COLOR_BUFFER_BIT);

//     glOrtho(XMIN, XMAX, YMIN, YMAX , -1.0, 1.0);
    glOrtho(0.0, NX, 0.0, NY, -1.0, 1.0);
}


void hsl_to_rgb(h, s, l, rgb)       /* color conversion from HSL to RGB */
/* h = hue, s = saturation, l = luminosity */
double h, s, l, rgb[3];
{
    double c = 0.0, m = 0.0, x = 0.0;

    c = (1.0 - fabs(2.0 * l - 1.0)) * s;
    m = 1.0 * (l - 0.5 * c);
    x = c * (1.0 - fabs(fmod(h / 60.0, 2) - 1.0));

    if (h >= 0.0 && h < 60.0)
    {
        rgb[0] = c+m; rgb[1] = x+m; rgb[2] = m;
    }
    else if (h < 120.0)
    {
        rgb[0] = x+m; rgb[1] = c+m; rgb[2] = m;
    }
    else if (h < 180.0)
    {
        rgb[0] = m; rgb[1] = c+m; rgb[2] = x+m;
    }
    else if (h < 240.0)
    {
        rgb[0] = m; rgb[1] = x+m; rgb[2] = c+m;
    }
    else if (h < 300.0)
    {
        rgb[0] = x+m; rgb[1] = m; rgb[2] = c+m;
    }
    else if (h < 360.0)
    {
        rgb[0] = c+m; rgb[1] = m; rgb[2] = x+m;
    }
    else
    {
        rgb[0] = m; rgb[1] = m; rgb[2] = m;
    }
}


double color_amplitude(value, scale, time)
/* transforms the wave amplitude into a double in [-1,1] to feed into color scheme */
double value, scale;
int time;
{
    return(tanh(SLOPE*value/scale)*exp(-((double)time*ATTENUATION)));
}


void color_scheme(scheme, value, scale, time, rgb) /* color scheme */
double value, scale;
int scheme, time;
double rgb[3];
{
    double hue, y, r, amplitude;
    int intpart;

    /* saturation = r, luminosity = y */
    switch (scheme) {
        case C_LUM:
        {
            hue = COLORHUE + (double)time*COLORDRIFT/(double)NSTEPS;
            if (hue < 0.0) hue += 360.0;
            if (hue >= 360.0) hue -= 360.0;
            r = 0.9;
            amplitude = color_amplitude(value, scale, time);
            y = LUMMEAN + amplitude*LUMAMP;
            intpart = (int)y;
            y -= (double)intpart;
            hsl_to_rgb(hue, r, y, rgb);
            break;
        }
        case C_HUE:
        {
            r = 0.9;
            amplitude = color_amplitude(value, scale, time);
            y = 0.5;
            hue = HUEMEAN + amplitude*HUEAMP;
            if (hue < 0.0) hue += 360.0;
            if (hue >= 360.0) hue -= 360.0;
            hsl_to_rgb(hue, r, y, rgb);
            break;
        }
    }
}


void blank()
{
    glClearColor(1.0, 1.0, 1.0, 1.0);
    glClear(GL_COLOR_BUFFER_BIT);
}


void save_frame()
{
  static int counter = 0;
  char *name="wave.", n2[100];
  char format[6]=".%05i";

    counter++;
//     printf (" p2 counter = %d \n",counter);
    strcpy(n2, name);
    sprintf(strstr(n2,"."), format, counter);
    strcat(n2, ".tif");
    printf(" saving frame %s \n",n2);
    writetiff(n2, "Wave equation in a planar domain", 0, 0,
         WINWIDTH, WINHEIGHT, COMPRESSION_LZW);

}

/*********************/
/* some basic math   */
/*********************/

 double vabs(x)     /* absolute value */
 double x;
 {
	double res;

	if (x<0.0) res = -x;
	else res = x;
	return(res);
 }

 double module2(x, y)   /* Euclidean norm */
 double x, y;

 {
	double m;

	m = sqrt(x*x + y*y);
	return(m);
 }

 double argument(x, y)
 double x, y;

 {
	double alph;

	if (x!=0.0)
	{
		alph = atan(y/x);
		if (x<0.0)
			alph += PI;
	}
	else
	{
		alph = PID;
		if (y<0.0)
			alph = PI*1.5;
	}
	return(alph);
 }

 /*********************/
/* drawing routines  */
/*********************/

/* The billiard boundary is drawn in (x,y) coordinates               */
/* However for the grid points, we use integer coordinates (i,j)     */
/* GL would allow to always work in (x,y) coordinates but using both */
/* sets of coordinates decreases number of double computations when  */
/* drawing the field                                                 */

void xy_to_ij(x, y, ij)
/* convert (x,y) position to (i,j) in table representing wave */
double x, y;
int ij[2];
{
    double x1, y1;

    x1 = (x - XMIN)/(XMAX - XMIN);
    y1 = (y - YMIN)/(YMAX - YMIN);

    ij[0] = (int)(x1 * (double)NX);
    ij[1] = (int)(y1 * (double)NY);
}


void xy_to_pos(x, y, pos)
/* convert (x,y) position to double-valued position in table representing wave */
double x, y, pos[2];
{
    double x1, y1;

    x1 = (x - XMIN)/(XMAX - XMIN);
    y1 = (y - YMIN)/(YMAX - YMIN);

    pos[0] = x1 * (double)NX;
    pos[1] = y1 * (double)NY;
}


void ij_to_xy(i, j, xy)
/* convert (i,j) position in table representing wave to (x,y) */
int i, j;
double xy[2];
{
    double x1, y1;

    xy[0] = XMIN + ((double)i)*(XMAX-XMIN)/((double)NX);
    xy[1] = YMIN + ((double)j)*(YMAX-YMIN)/((double)NY);
}




int xy_in_billiard(x, y)
/* returns 1 if (x,y) represents a point in the billiard */
double x, y;
{
    double l2, r2, omega, c, angle;
    int k, condition;

    switch (B_DOMAIN) {
        case D_RECTANGLE:
        {
            if ((vabs(x) <LAMBDA)&&(vabs(y) < 1.0)) return(1);
            else return(0);
            break;
        }
        case D_ELLIPSE:
        {
            if (x*x/(LAMBDA*LAMBDA) + y*y < 1.0) return(1);
            else return(0);
            break;
        }
        case D_STADIUM:
        {
            if ((x > -0.5*LAMBDA)&&(x < 0.5*LAMBDA)&&(y > -1.0)&&(y < 1.0)) return(1);
            else if (module2(x+0.5*LAMBDA, y) < 1.0) return(1);
            else if (module2(x-0.5*LAMBDA, y) < 1.0) return(1);
            else return(0);
            break;
        }
        case D_SINAI:
        {
            if (x*x + y*y > LAMBDA*LAMBDA) return(1);
            else return(0);
            break;
        }
        case D_DIAMOND:
        {
            l2 = LAMBDA*LAMBDA;
            r2 = l2 + (LAMBDA-1.0)*(LAMBDA-1.0);
            if ((x*x + y*y < 1.0)&&((x-LAMBDA)*(x-LAMBDA) + (y-LAMBDA)*(y-LAMBDA) > r2)
                &&((x-LAMBDA)*(x-LAMBDA) + (y+LAMBDA)*(y+LAMBDA) > r2)
                &&((x+LAMBDA)*(x+LAMBDA) + (y-LAMBDA)*(y-LAMBDA) > r2)
                &&((x+LAMBDA)*(x+LAMBDA) + (y+LAMBDA)*(y+LAMBDA) > r2)) return(1);
            else return(0);
            break;
        }
        case D_TRIANGLE:
        {
            if ((x>-LAMBDA)&&(y>-1.0)&&(LAMBDA*y+x<0.0)) return(1);
            else return(0);
            break;
        }
        case D_FLAT:
        {
            if (y > -LAMBDA) return(1);
            else return(0);
            break;
        }
        case D_ANNULUS:
        {
            l2 = LAMBDA*LAMBDA;
            r2 = x*x + y*y;
            if ((r2 > l2)&&(r2 < 1.0)) return(1);
            else return(0);
        }
        case D_POLYGON:
        {
            condition = 1;
            omega = DPI/((double)NPOLY);
            c = cos(omega*0.5);
            for (k=0; k<NPOLY; k++)  
            {
                angle = APOLY*PID + (k+0.5)*omega;
                condition = condition*(x*cos(angle) + y*sin(angle) < c);
            }
//             for (k=0; k<NPOLY; k++)  condition = condition*(-x*sin((k+0.5)*omega) + y*cos((k+0.5)*omega) < c);
            return(condition);
        }   
        default:
        {
            printf("Function ij_in_billiard not defined for this billiard \n");
            return(0);
        }
    }
}

void init_wave(x, y, phi, psi, xy_in)
/* initialise field with drop at (x,y) - phi is wave height, psi is phi at time t-1 */
    double x, y, *phi[NX], *psi[NX]; short int * xy_in[NX];

{
    int i, j;
    double xy[2], dist2;

    for (i=0; i<NX; i++)
        for (j=0; j<NY; j++)
        {
            ij_to_xy(i, j, xy);
            dist2 = (xy[0]-x)*(xy[0]-x) + (xy[1]-y)*(xy[1]-y);
	    xy_in[i][j] = xy_in_billiard(xy[0],xy[1]);
	    phi[i][j] = 0.2*exp(-dist2/0.001)*cos(-sqrt(dist2)/0.01);
            psi[i][j] = 0.0;
        }
}

void add_drop_to_wave(factor, x, y, phi, psi)
/* add drop at (x,y) to the field with given prefactor */
double factor, x, y, *phi[NX], *psi[NX];
{
    int i, j;
    double xy[2], dist2;

    for (i=0; i<NX; i++)
        for (j=0; j<NY; j++)
        {
            ij_to_xy(i, j, xy);
            dist2 = (xy[0]-x)*(xy[0]-x) + (xy[1]-y)*(xy[1]-y);
            phi[i][j] += 0.2*factor*exp(-dist2/0.001)*cos(-sqrt(dist2)/0.01);
        }
}




int ij_in_billiard(i, j)
/* returns 1 if (i,j) represents a point in the billiard */
int i, j;
{
    double xy[2];

    ij_to_xy(i, j, xy);

    return(xy_in_billiard(xy[0], xy[1]));
}


void draw_billiard()      /* draws the billiard boundary */
{
    double x0, x, y, phi, r = 0.01, pos[2], pos1[2], alpha, dphi, omega;
    int i;

    glColor3f(0.0, 0.0, 0.0);
    glLineWidth(5);

    glEnable(GL_LINE_SMOOTH);

    switch (B_DOMAIN) {
        case (D_RECTANGLE):
        {
            glBegin(GL_LINE_LOOP);
            xy_to_pos(LAMBDA, -1.0, pos);
            glVertex2d(pos[0], pos[1]);
            xy_to_pos(LAMBDA, 1.0, pos);
            glVertex2d(pos[0], pos[1]);
            xy_to_pos(-LAMBDA, 1.0, pos);
            glVertex2d(pos[0], pos[1]);
            xy_to_pos(-LAMBDA, -1.0, pos);
            glVertex2d(pos[0], pos[1]);
            glEnd();
            break;
        }
        case (D_ELLIPSE):
        {
            glBegin(GL_LINE_LOOP);
            for (i=0; i<=NSEG; i++)
            {
                phi = (double)i*DPI/(double)NSEG;
                x = LAMBDA*cos(phi);
                y = sin(phi);
                xy_to_pos(x, y, pos);
                glVertex2d(pos[0], pos[1]);
            }
            glEnd ();

            /* draw foci */
            if (FOCI)
            {
                glColor3f(0.3, 0.3, 0.3);
                x0 = sqrt(LAMBDA*LAMBDA-1.0);

                glLineWidth(2);
                glEnable(GL_LINE_SMOOTH);
                glBegin(GL_LINE_LOOP);
                for (i=0; i<=NSEG; i++)
                {
                    phi = (double)i*DPI/(double)NSEG;
                    x = x0 + r*cos(phi);
                    y = r*sin(phi);
                    xy_to_pos(x, y, pos);
                    glVertex2d(pos[0], pos[1]);
                }
                glEnd();

                glBegin(GL_LINE_LOOP);
                for (i=0; i<=NSEG; i++)
                {
                    phi = (double)i*DPI/(double)NSEG;
                    x = -x0 + r*cos(phi);
                    y = r*sin(phi);
                    xy_to_pos(x, y, pos);
                    glVertex2d(pos[0], pos[1]);
                }
                glEnd();
            }
            break;
        }
        case D_STADIUM:
        {
            glBegin(GL_LINE_LOOP);
            for (i=0; i<=NSEG; i++)
            {
                phi = -PID + (double)i*PI/(double)NSEG;
                x = 0.5*LAMBDA + cos(phi);
                y = sin(phi);
                xy_to_pos(x, y, pos);
                glVertex2d(pos[0], pos[1]);
            }
            for (i=0; i<=NSEG; i++)
            {
                phi = PID + (double)i*PI/(double)NSEG;
                x = -0.5*LAMBDA + cos(phi);
                y = sin(phi);
                xy_to_pos(x, y, pos);
                glVertex2d(pos[0], pos[1]);
            }
            glEnd();
            break;
        }
        case D_SINAI:
        {
            glBegin(GL_LINE_LOOP);
            for (i=0; i<=NSEG; i++)
            {
                phi = (double)i*DPI/(double)NSEG;
                x = LAMBDA*cos(phi);
                y = LAMBDA*sin(phi);
                xy_to_pos(x, y, pos);
                glVertex2d(pos[0], pos[1]);
            }
            glEnd();
            break;
        }
        case D_DIAMOND:
        {
            alpha = atan(1.0 - 1.0/LAMBDA);
            dphi = (PID - 2.0*alpha)/(double)NSEG;
            r = sqrt(LAMBDA*LAMBDA + (LAMBDA-1.0)*(LAMBDA-1.0));
            glBegin(GL_LINE_LOOP);
            for (i=0; i<=NSEG; i++)
            {
                phi = alpha + (double)i*dphi;
                x = -LAMBDA + r*cos(phi);
                y = -LAMBDA + r*sin(phi);
                xy_to_pos(x, y, pos);
                glVertex2d(pos[0], pos[1]);
            }
            for (i=0; i<=NSEG; i++)
            {
                phi = alpha - PID + (double)i*dphi;
                x = -LAMBDA + r*cos(phi);
                y = LAMBDA + r*sin(phi);
                xy_to_pos(x, y, pos);
                glVertex2d(pos[0], pos[1]);
            }
            for (i=0; i<=NSEG; i++)
            {
                phi = alpha + PI + (double)i*dphi;
                x = LAMBDA + r*cos(phi);
                y = LAMBDA + r*sin(phi);
                xy_to_pos(x, y, pos);
                glVertex2d(pos[0], pos[1]);
            }
            for (i=0; i<=NSEG; i++)
            {
                phi = alpha + PID + (double)i*dphi;
                x = LAMBDA + r*cos(phi);
                y = -LAMBDA + r*sin(phi);
                xy_to_pos(x, y, pos);
                glVertex2d(pos[0], pos[1]);
            }
            glEnd();
            break;
        }
        case (D_TRIANGLE):
        {
            glBegin(GL_LINE_LOOP);
            xy_to_pos(-LAMBDA, -1.0, pos);
            glVertex2d(pos[0], pos[1]);
            xy_to_pos(LAMBDA, -1.0, pos);
            glVertex2d(pos[0], pos[1]);
            xy_to_pos(-LAMBDA, 1.0, pos);
            glVertex2d(pos[0], pos[1]);
            glEnd();
            break;
        }
        case (D_FLAT):
        {
            glBegin(GL_LINE_LOOP);
            xy_to_pos(XMIN, -LAMBDA, pos);
            glVertex2d(pos[0], pos[1]);
            xy_to_pos(XMAX, -LAMBDA, pos);
            glVertex2d(pos[0], pos[1]);
            glEnd();
            break;
        }
        case (D_ANNULUS):
        {
            glBegin(GL_LINE_LOOP);
            for (i=0; i<=NSEG; i++)
            {
                phi = (double)i*DPI/(double)NSEG;
                x = LAMBDA*cos(phi);
                y = LAMBDA*sin(phi);
                xy_to_pos(x, y, pos);
                glVertex2d(pos[0], pos[1]);
            }
            glEnd ();
            glBegin(GL_LINE_LOOP);
            for (i=0; i<=NSEG; i++)
            {
                phi = (double)i*DPI/(double)NSEG;
                x = cos(phi);
                y = sin(phi);
                xy_to_pos(x, y, pos);
                glVertex2d(pos[0], pos[1]);
            }
            glEnd ();      
            break;
        }
        case (D_POLYGON):
        {
            omega = DPI/((double)NPOLY);
            glBegin(GL_LINE_LOOP);
            for (i=0; i<=NPOLY; i++)
            {
                x = cos(i*omega + APOLY*PID);
                y = sin(i*omega + APOLY*PID);
                xy_to_pos(x, y, pos);
                glVertex2d(pos[0], pos[1]);
            }
            glEnd ();
            break;
        }
        default:
        {
            printf("Function draw_billiard not defined for this billiard \n");
        }
    }
}


/*********************/
/* animation part    */
/*********************/


void draw_wave(phi, psi, xy_in, scale, time)
/* draw the field */
double *phi[NX], *psi[NX], scale;
short int *xy_in[NX];
int time;
{
    int i, j;
    double rgb[3], xy[2], x1, y1, x2, y2;

    glBegin(GL_QUADS);

    for (i=0; i<NX; i++)
        for (j=0; j<NY; j++)
        {
            if (xy_in[i][j])
            {
                color_scheme(COLOR_SCHEME, phi[i][j], scale, time, rgb);
                glColor3f(rgb[0], rgb[1], rgb[2]);

                glVertex2i(i, j);
                glVertex2i(i+1, j);
                glVertex2i(i+1, j+1);
                glVertex2i(i, j+1);
            }
        }

    glEnd ();
}

void evolve_wave(phi, psi, xy_in)
/* time step of field evolution */
/* phi is value of field at time t, psi at time t-1 */
    double *phi[NX], *psi[NX]; short int *xy_in[NX];
{
    int i, j, iplus, iminus, jplus, jminus;
    double delta, x, y;

    #pragma omp parallel for private(i,j,iplus,iminus,delta,x,y)
    for (i=0; i<NX; i++){
        for (j=0; j<NY; j++){
            if (xy_in[i][j]){
                /* discretized Laplacian */
		iplus = (i+1) % NX;
		iminus = (i-1) % NX;
		if (iminus < 0) iminus += NX;
                jplus = (j+1) % NY;
                jminus = (j-1) % NY;
                if (jminus < 0) jminus += NY;
                delta = phi[iplus][j] + phi[iminus][j] + phi[i][jplus] + phi[i][jminus] - 4.0*phi[i][j];

                x = phi[i][j];
		y = psi[i][j];

                /* evolve phi */
                phi[i][j] = -y + 2*x + courant2*delta - KAPPA*x - GAMMA*(x-y);
//                 phi[i][j] = -psi[i][j] + 2*x + courant2*delta - GAMMA*x;

                /* Old versions of the simulation used this: */
//                 phi[i][j] = (-psi[i][j] + 2*phi[i][j] + courant2*delta)*damping;
//                 where damping = 1.0 - 0.0001;

                psi[i][j] = x;

                if (FLOOR)
                {
                    if (phi[i][j] > VMAX) phi[i][j] = VMAX;
                    if (phi[i][j] < -VMAX) phi[i][j] = -VMAX;
                    if (psi[i][j] > VMAX) psi[i][j] = VMAX;
                    if (psi[i][j] < -VMAX) psi[i][j] = -VMAX;
                }
            }
        }
    }
//     printf("phi(0,0) = %.3lg, psi(0,0) = %.3lg\n", phi[NX/2][NY/2], psi[NX/2][NY/2]);
}


double compute_variance(phi, psi, xy_in)
/* compute the variance of the field, to adjust color scheme */
    double *phi[NX], *psi[NX]; short int * xy_in[NX];
{
    int i, j, n = 0;
    double variance = 0.0;

    for (i=1; i<NX; i++)
        for (j=1; j<NY; j++)
        {
            if (xy_in[i][j])
            {
                n++;
                variance += phi[i][j]*phi[i][j];
            }
        }
    if (n==0) n=1;
    return(variance/(double)n);
}


void animation()
{
    double time, scale;
    double *phi[NX], *psi[NX];
    short int *xy_in[NX];
    int i, j, s;

    /* Since NX and NY are big, it seemed wiser to use some memory allocation here */
    for (i=0; i<NX; i++)
    {
        phi[i] = (double *)malloc(NY*sizeof(double));
        psi[i] = (double *)malloc(NY*sizeof(double));
        xy_in[i] = (short int *)malloc(NY*sizeof(short int));
    }

    courant2 = COURANT*COURANT;

    /* initialize wave with a drop at one point, zero elsewhere */
    init_wave(0.0, 0.0, phi, psi, xy_in);

    /* add a drop at another point */
//     add_drop_to_wave(1.0, 0.7, 0.0, phi, psi);
//     add_drop_to_wave(1.0, -0.7, 0.0, phi, psi);
//     add_drop_to_wave(1.0, 0.0, -0.7, phi, psi);

    blank();
    glColor3f(0.0, 0.0, 0.0);
    draw_wave(phi, psi, xy_in, 1.0, 0);
    draw_billiard();

    glutSwapBuffers();



    sleep(SLEEP1);

    for (i=0; i<=NSTEPS; i++)
    {
	//printf("%d\n",i);
        /* compute the variance of the field to adjust color scheme */
        /* the color depends on the field divided by sqrt(1 + variance) */
        if (SCALE)
        {
            scale = sqrt(1.0 + compute_variance(phi,psi, xy_in));
//             printf("Scaling factor: %5lg\n", scale);
        }

        else scale = 1.0;


        draw_wave(phi, psi, xy_in, scale, i);
        for (j=0; j<NVID; j++) evolve_wave(phi, psi, xy_in);
        draw_billiard();


	glutSwapBuffers();

	if (MOVIE)
        {
            save_frame();

            /* it seems that saving too many files too fast can cause trouble with the file system */
            /* so this is to make a pause from time to time - parameter PAUSE may need adjusting   */
            if (i % PAUSE == PAUSE - 1)
            {
                printf("Making a short pause\n");
                sleep(PSLEEP);
                s = system("mv wave*.tif tif_wave/");
            }
        }

    }

    if (MOVIE) for (i=0; i<20; i++){
		       save_frame();
		       s = system("mv wave*.tif tif_wave/");
		   }
    for (i=0; i<NX; i++)
    {
        free(phi[i]);
        free(psi[i]);
    }

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

    glutDestroyWindow(glutGetWindow());

}


int main(int argc, char** argv)
{
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH);
    glutInitWindowSize(WINWIDTH,WINHEIGHT);
    glutCreateWindow("Wave equation in a planar domain");

    init();

    glutDisplayFunc(display);

    glutMainLoop();

    return 0;
}

