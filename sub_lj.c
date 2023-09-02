/*********************/
/* Graphics routines */
/*********************/

#include "colors_waves.c"

// #define HUE_TYPE0 260.0     /* hue of particles of type 0 */
// #define HUE_TYPE0 300.0     /* hue of particles of type 0 */
// #define HUE_TYPE1 90.0      /* hue of particles of type 1 */


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
  if (SAVE_MEMORY) free(image);
  return 0;
}


void init()		/* initialisation of window */
{
    glLineWidth(3);

    glClearColor(0.0, 0.0, 0.0, 1.0);
    glClear(GL_COLOR_BUFFER_BIT);

    glOrtho(XMIN, XMAX, YMIN, YMAX , -1.0, 1.0);
//     glOrtho(0.0, NX, 0.0, NY, -1.0, 1.0);
}


void blank()
{
    if (BLACK) glClearColor(0.0, 0.0, 0.0, 1.0);
    else glClearColor(1.0, 1.0, 1.0, 1.0);
    glClear(GL_COLOR_BUFFER_BIT);
}

void save_frame_lj()
{
  static int counter = 0;
  char *name="lj.", n2[100];
  char format[6]=".%05i";

    counter++;
//     printf (" p2 counter = %d \n",counter);
    strcpy(n2, name);
    sprintf(strstr(n2,"."), format, counter);
    strcat(n2, ".tif");
    printf(" saving frame %s \n",n2);
    writetiff(n2, "Particles with Lennard-Jones interaction in a planar domain", 0, 0,
         WINWIDTH, WINHEIGHT, COMPRESSION_LZW);

}

void save_frame_lj_counter(int counter)
{
  char *name="lj.", n2[100];
  char format[6]=".%05i";

    strcpy(n2, name);
    sprintf(strstr(n2,"."), format, counter);
    strcat(n2, ".tif");
    printf(" saving frame %s \n",n2);
    writetiff(n2, "Particles with Lennard-Jones interaction in a planar domain", 0, 0,
         WINWIDTH, WINHEIGHT, COMPRESSION_LZW);

}


void write_text_fixedwidth( double x, double y, char *st)
{
    int l, i;

    l=strlen( st ); // see how many characters are in text string.
    glRasterPos2d( x, y); // location to start printing text
    for( i=0; i < l; i++) // loop until i is greater then l
    {
//         glutBitmapCharacter(GLUT_BITMAP_TIMES_ROMAN_24, st[i]); // Print a character on the screen
//    glutBitmapCharacter(GLUT_BITMAP_8_BY_13, st[i]); // Print a character on the screen
   glutBitmapCharacter(GLUT_BITMAP_9_BY_15, st[i]); // Print a character on the screen
    }
} 


void write_text( double x, double y, char *st)
{
    int l,i;

    l=strlen( st ); // see how many characters are in text string.
    glRasterPos2d( x, y); // location to start printing text
    for( i=0; i < l; i++) // loop until i is greater then l
    {
        glutBitmapCharacter(GLUT_BITMAP_TIMES_ROMAN_24, st[i]); // Print a character on the screen
//    glutBitmapCharacter(GLUT_BITMAP_8_BY_13, st[i]); // Print a character on the screen
    }
} 




/*********************/
/* some basic math   */
/*********************/

 double vabs(double x)     /* absolute value */
 {
	double res;

	if (x<0.0) res = -x;
	else res = x;
	return(res);
 }

 double module2(double x, double y)   /* Euclidean norm */
 {
	double m;

	m = sqrt(x*x + y*y);
	return(m);
 }

 double argument(double x, double y)
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
 
 double ipow(double x, int n)
 {
    double y;
    int i;
    
    y = x;
    for (i=1; i<n; i++) y *= x;
    
    return(y);
 }
 

double gaussian()
/* returns standard normal random variable, using Box-Mueller algorithm */
{
    static double V1, V2, S;
    static int phase = 0;
    double X;

    if (phase == 0) 
    {
        do 
        {
        double U1 = (double)rand() / RAND_MAX;
        double U2 = (double)rand() / RAND_MAX;
        V1 = 2 * U1 - 1;
        V2 = 2 * U2 - 1;
        S = V1 * V1 + V2 * V2;
        } 
        while(S >= 1 || S == 0);
        X = V1 * sqrt(-2 * log(S) / S);
    } 
    else X = V2 * sqrt(-2 * log(S) / S);

    phase = 1 - phase;

    return X;
}


/*********************/
/* drawing routines  */
/*********************/



void erase_area(double x, double y, double dx, double dy)
{
    double pos[2], rgb[3];
    
    hsl_to_rgb(220.0, 0.8, 0.7, rgb);
    
    glColor3f(rgb[0], rgb[1], rgb[2]);
    glBegin(GL_QUADS);
    glVertex2d(x - dx, y - dy);
    glVertex2d(x + dx, y - dy);
    glVertex2d(x + dx, y + dy);
    glVertex2d(x - dx, y + dy);
    glEnd();
}


void erase_area_rgb(double x, double y, double dx, double dy, double rgb[3])
{
    double pos[2];
    
    glColor3f(rgb[0], rgb[1], rgb[2]);
    glBegin(GL_QUADS);
    glVertex2d(x - dx, y - dy);
    glVertex2d(x + dx, y - dy);
    glVertex2d(x + dx, y + dy);
    glVertex2d(x - dx, y + dy);
    glEnd();
}


void erase_area_hsl(double x, double y, double dx, double dy, double h, double s, double l)
{
    double pos[2], rgb[3];
    
    hsl_to_rgb(h, s, l, rgb);
    erase_area_rgb(x, y, dx, dy, rgb);
}

void erase_area_hsl_turbo(double x, double y, double dx, double dy, double h, double s, double l)
{
    double pos[2], rgb[3];
    
    hsl_to_rgb_turbo(h, s, l, rgb);
    erase_area_rgb(x, y, dx, dy, rgb);
}

void draw_line(double x1, double y1, double x2, double y2)
{
    glBegin(GL_LINE_STRIP);
    glVertex2d(x1, y1);
    glVertex2d(x2, y2);
    glEnd();    
}

void draw_arrow(double x1, double y1, double x2, double y2, double angle, double length)
{
    double alpha, beta, x3, y3, x4, y4, x5, y5;
    
    alpha = argument(x2 - x1, y2 - y1);
    beta = angle*PI/180.0;
    x3 = x2 - length*cos(alpha - beta);
    y3 = y2 - length*sin(alpha - beta);
    x4 = x2 - length*cos(alpha + beta);
    y4 = y2 - length*sin(alpha + beta);
    x5 = x2 - 0.5*length*cos(alpha);
    y5 = y2 - 0.5*length*sin(alpha);
    
    glBegin(GL_LINE_STRIP);
    glVertex2d(x1, y1);
    glVertex2d(x5, y5);
    glEnd();  
    
    glBegin(GL_TRIANGLE_FAN);
    glVertex2d(x2, y2);
    glVertex2d(x3, y3);
    glVertex2d(x4, y4);
    glEnd(); 
}

void draw_rectangle(double x1, double y1, double x2, double y2)
{    
    glBegin(GL_LINE_LOOP);
    glVertex2d(x1, y1);
    glVertex2d(x2, y1);
    glVertex2d(x2, y2);
    glVertex2d(x1, y2);
    glEnd();    
}

void draw_colored_rectangle(double x1, double y1, double x2, double y2, double rgb[3])
{    
    glColor3f(rgb[0], rgb[1], rgb[2]);
    glBegin(GL_TRIANGLE_FAN);
    glVertex2d(x1, y1);
    glVertex2d(x2, y1);
    glVertex2d(x2, y2);
    glVertex2d(x1, y2);
    glEnd();    
}

void draw_triangle(double x1, double y1, double x2, double y2, double x3, double y3)
{    
    glBegin(GL_LINE_LOOP);
    glVertex2d(x1, y1);
    glVertex2d(x2, y2);
    glVertex2d(x3, y3);
    glEnd();    
}

void draw_colored_triangle(double x1, double y1, double x2, double y2, double x3, double y3, double rgb[3])
{    
    glColor3f(rgb[0], rgb[1], rgb[2]);
    glBegin(GL_TRIANGLE_FAN);
    glVertex2d(x1, y1);
    glVertex2d(x2, y2);
    glVertex2d(x3, y3);
    glEnd();    
}


void draw_circle(double x, double y, double r, int nseg)
{
    int i;
    double pos[2], alpha, dalpha;
    
    dalpha = DPI/(double)nseg;
    
    glBegin(GL_LINE_LOOP);
    for (i=0; i<=nseg; i++)
    {
        alpha = (double)i*dalpha;
        glVertex2d(x + r*cos(alpha), y + r*sin(alpha));
    }
    glEnd();
}

void draw_colored_circle(double x, double y, double r, int nseg, double rgb[3])
{
    int i;
    double pos[2], alpha, dalpha;
    
    dalpha = DPI/(double)nseg;
    
    glColor3f(rgb[0], rgb[1], rgb[2]);
    glBegin(GL_TRIANGLE_FAN);
    glVertex2d(x, y);
    for (i=0; i<=nseg; i++)
    {
        alpha = (double)i*dalpha;
        glVertex2d(x + r*cos(alpha), y + r*sin(alpha));
    }
    
    glEnd();
}

void draw_circle_precomp(double x, double y, double r)
{
    int i;
        
    glBegin(GL_LINE_LOOP);
    for (i=0; i<=NSEG; i++) glVertex2d(x + r*cosangle[i], y + r*sinangle[i]);
    glEnd();
}

void draw_colored_circle_precomp(double x, double y, double r, double rgb[3])
{
    int i;
    
    glColor3f(rgb[0], rgb[1], rgb[2]);
    glBegin(GL_TRIANGLE_FAN);
    glVertex2d(x, y);
    for (i=0; i<=NSEG; i++) glVertex2d(x + r*cosangle[i], y + r*sinangle[i]);
    glEnd();
}

void draw_polygon(double x, double y, double r, int nsides, double angle)
{
    int i;
    double pos[2], alpha, dalpha;
    
    dalpha = DPI/(double)nsides;
    
    glColor3f(1.0, 1.0, 1.0);
    glBegin(GL_LINE_LOOP);
    for (i=0; i<=nsides; i++)
    {
        alpha = angle + (double)i*dalpha;
        glVertex2d(x + r*cos(alpha), y + r*sin(alpha));
    }
    glEnd();
}

void draw_colored_polygon(double x, double y, double r, int nsides, double angle, double rgb[3])
{
    int i;
    double pos[2], alpha, dalpha;
    
    dalpha = DPI/(double)nsides;
    
    glColor3f(rgb[0], rgb[1], rgb[2]);
    glBegin(GL_TRIANGLE_FAN);
    glVertex2d(x, y);
    for (i=0; i<=nsides; i++)
    {
        alpha = angle + (double)i*dalpha;
        glVertex2d(x + r*cos(alpha), y + r*sin(alpha));
    }
    
    glEnd();
}

void draw_rhombus(double x, double y, double r, double angle)
{
    int i;
    static int first = 1;
    static double ratio;
    
    if (first)
    {
        ratio = tan(0.1*PI);
        first = 0;
    }
    
    glBegin(GL_LINE_LOOP);
    glVertex2d(x + r*cos(angle), y + r*sin(angle));
    glVertex2d(x - ratio*r*sin(angle), y + ratio*r*cos(angle));
    glVertex2d(x - r*cos(angle), y - r*sin(angle));
    glVertex2d(x + ratio*r*sin(angle), y - ratio*r*cos(angle));
    glEnd();
}

void draw_colored_rhombus(double x, double y, double r, double angle, double rgb[3])
{
    int i;
    static int first = 1;
    static double ratio;
    
    if (first)
    {
        ratio = tan(0.1*PI);
        first = 0;
    }
    
    glColor3f(rgb[0], rgb[1], rgb[2]);
    glBegin(GL_TRIANGLE_FAN);
    glVertex2d(x, y);
    glVertex2d(x + r*cos(angle), y + r*sin(angle));
    glVertex2d(x - ratio*r*sin(angle), y + ratio*r*cos(angle));
    glVertex2d(x - r*cos(angle), y - r*sin(angle));
    glVertex2d(x + ratio*r*sin(angle), y - ratio*r*cos(angle));
    glVertex2d(x + r*cos(angle), y + r*sin(angle));
    glEnd();
}

void draw_colored_sector(double xc, double yc, double r1, double r2, double angle1, double angle2, double rgb[3], int nsides)
{    
    int i;
    double angle, dangle;
    
    dangle = (angle2 - angle1)/(double)nsides;
    
    glColor3f(rgb[0], rgb[1], rgb[2]);
    glBegin(GL_TRIANGLE_FAN);
    glVertex2d(xc + r1*cos(angle1), yc + r1*sin(angle1));
    for (i = 0; i < nsides+1; i++)
    {
        angle = angle1 + dangle*(double)i;
        glVertex2d(xc + r2*cos(angle), yc + r2*sin(angle));
    }
    glEnd();  
    glBegin(GL_TRIANGLE_FAN);
    glVertex2d(xc + r2*cos(angle2), yc + r2*sin(angle2));
    for (i = 0; i < nsides+1; i++)
    {
        angle = angle1 + dangle*(double)i;
        glVertex2d(xc + r1*cos(angle), yc + r1*sin(angle));
    }
    glEnd();   
}

double type_hue(int type)
{
    int hue;
    double t2;
    static double b, hmax;
    static int first = 1;
    
    if (first)
    {
        hmax = 360.0;
        b = 16.0*(hmax - HUE_TYPE3);
        first = 0;
    }
    
    if ((RD_REACTION == CHEM_CATALYTIC_A2D)&&(type == 4)) return(HUE_TYPE3); 
    
    if ((RD_REACTION == CHEM_ABDACBE)&&(type == 4)) return(HUE_TYPE3); 
    if ((RD_REACTION == CHEM_ABDACBE)&&(type == 5)) return(280.0); 
    
    switch (type) {
        case (0): return(HUE_TYPE0);
        case (1): return(HUE_TYPE0);
        case (2): return(HUE_TYPE1);
        case (3): return(HUE_TYPE2);
        default:
        {
            if (RD_REACTION == CHEM_BZ)
            {
                if (type == 7) return(HUE_TYPE2);
                if (type == 8) type = 5;
            }
            else if ((RD_REACTION == CHEM_BRUSSELATOR)&&(type >= 5)) return(70.0);
            t2 = (double)(type*type);
            hue = (hmax*t2 - b)/t2;
            return(hue);
        }
    }
}


void set_type_color(int type, double lum, double rgb[3])
{
    int hue;
    
    hue = type_hue(type);
    hsl_to_rgb_palette(hue, 0.9, 0.5, rgb, COLOR_PALETTE);
    glColor3f(lum*rgb[0], lum*rgb[1], lum*rgb[2]);
}

double distance_to_segment(double x, double y, double x1, double y1, double x2, double y2)
/* distance of (x,y) to segment from (x1,y1) to (x2,y2) */
{
    double xp, yp, angle, length, ca, sa;
    
    angle = argument(x2 - x1, y2 - y1);
    length = module2(x2 - x1, y2 - y1);
    
    ca = cos(angle);
    sa = sin(angle);
    
    xp = ca*(x - x1) + sa*(y - y1);
    yp = -sa*(x - x1) + ca*(y - y1);
    
    if ((xp >= 0)&&(xp <= length)) return(vabs(yp));
    else if (xp < 0) return(module2(xp, yp));
    else return(module2(xp-length, yp));
}

int in_polygon(double x, double y, double r, int npoly, double apoly)
/* test whether (x,y) is in regular polygon of npoly sides inscribed in circle of radius r, turned by apoly Pi/2 */
{
    int condition = 1, k;
    double omega, cw, angle; 
    
    omega = DPI/((double)npoly);
    cw = cos(omega*0.5);
    for (k=0; k<npoly; k++)  
    {
        angle = -apoly*PID + ((double)k+0.5)*omega;
        condition = condition*(x*cos(angle) + y*sin(angle) < r*cw);
    }
    return(condition);
}


void init_angles()
/* initialise cos and sin of angles to save computing time */
{
    int i;
    double alpha, dalpha;
    
    dalpha = DPI/(double)NSEG;
    
    for (i=0; i<=NSEG; i++)
    {
        alpha = (double)i*dalpha;
        cosangle[i] = cos(alpha);
        sinangle[i] = sin(alpha);        
    }
}

void init_particle_config(t_particle particles[NMAXCIRCLES])
/* initialise particle configuration */
{
    int i, j, k, n, ncirc0, n_p_active, ncandidates = PDISC_CANDIDATES, naccepted; 
    double dx, dy, p, phi, r, r0, ra[5], sa[5], height, x, y = 0.0, gamma, dpoisson = PDISC_DISTANCE*MU, xx[4], yy[4];
    short int active_poisson[NMAXCIRCLES], far;
    
    switch (CIRCLE_PATTERN) {
        case (C_SQUARE):
        {
            ncircles = NGRIDX*NGRIDY;
            dy = (YMAX - YMIN)/((double)NGRIDY);
            for (i = 0; i < NGRIDX; i++)
                for (j = 0; j < NGRIDY; j++)
                {
                    n = NGRIDY*i + j;
                    particles[n].xc = ((double)(i-NGRIDX/2) + 0.5)*dy;
                    particles[n].yc = YMIN + ((double)j + 0.5)*dy;
                    particles[n].radius = MU;
                    particles[n].active = 1;
                }
            break;
        }
        case (C_HEX):
        {
            ncircles = NGRIDX*(NGRIDY+1);
            dx = (INITXMAX - INITXMIN)/((double)NGRIDX);
            dy = (INITYMAX - INITYMIN)/((double)NGRIDY);
//             dx = dy*0.5*sqrt(3.0);
            for (i = 0; i < NGRIDX; i++)
                for (j = 0; j < NGRIDY+1; j++)
                {
                    n = (NGRIDY+1)*i + j;
//                     particles[n].xc = ((double)(i-NGRIDX/2) + 0.5)*dx;   /* is +0.5 needed? */
                    particles[n].xc = INITXMIN + ((double)i - 0.5)*dx;   
                    particles[n].yc = INITYMIN + ((double)j - 0.5)*dy;
                    if ((i+NGRIDX)%2 == 1) particles[n].yc += 0.5*dy;
                    if (particles[n].yc > YMAX) particles[n].yc += YMIN - YMAX;
//                     else if (particles[n].yc < YMIN) particles[n].yc += YMAX - YMIN;
                    particles[n].radius = MU;
                    /* activate only circles that intersect the domain */
                    if ((particles[n].yc < INITYMAX + MU)&&(particles[n].yc > INITYMIN - MU)&&(particles[n].xc < INITXMAX + MU)&&(particles[n].xc > INITXMIN - MU)) particles[n].active = 1;
                    else particles[n].active = 0;
                }
            break;
        }
        case (C_RAND_DISPLACED):
        {
            ncircles = NGRIDX*NGRIDY;
            dy = (YMAX - YMIN)/((double)NGRIDY);
            for (i = 0; i < NGRIDX; i++)
                for (j = 0; j < NGRIDY; j++)
                {
                    n = NGRIDY*i + j;
                    particles[n].xc = ((double)(i-NGRIDX/2) + 0.5*((double)rand()/RAND_MAX - 0.5))*dy;
                    particles[n].yc = YMIN + ((double)j + 0.5 + 0.5*((double)rand()/RAND_MAX - 0.5))*dy;
                    particles[n].radius = MU;
//                     particles[n].radius = MU*sqrt(1.0 + 0.8*((double)rand()/RAND_MAX - 0.5));
                    particles[n].active = 1;
                }
            break;
        }
        case (C_RAND_PERCOL):
        {
            ncircles = NGRIDX*NGRIDY;
            dy = (YMAX - YMIN)/((double)NGRIDY);
            for (i = 0; i < NGRIDX; i++)
                for (j = 0; j < NGRIDY; j++)
                {
                    n = NGRIDY*i + j;
                    particles[n].xc = ((double)(i-NGRIDX/2) + 0.5)*dy;
                    particles[n].yc = YMIN + ((double)j + 0.5)*dy;
                    particles[n].radius = MU;
                    p = (double)rand()/RAND_MAX;
                    if (p < P_PERCOL) particles[n].active = 1;
                    else particles[n].active = 0;
                }
            break;
        }
        case (C_RAND_POISSON):
        {
            ncircles = NPOISSON;
            for (n = 0; n < NPOISSON; n++)
            {
//                 particles[n].xc = LAMBDA*(2.0*(double)rand()/RAND_MAX - 1.0);
                particles[n].xc = (XMAX - XMIN)*(double)rand()/RAND_MAX + XMIN;
                particles[n].yc = (YMAX - YMIN)*(double)rand()/RAND_MAX + YMIN;
                particles[n].radius = MU;
                particles[n].active = 1;
            }
            break;
        }
        case (C_CLOAK):
        {
            ncircles = 200;
            for (i = 0; i < 40; i++)
                for (j = 0; j < 5; j++)
                {
                    n = 5*i + j;
                    phi = (double)i*DPI/40.0;
                    r = LAMBDA*0.5*(1.0 + (double)j/5.0);
                    particles[n].xc = r*cos(phi);
                    particles[n].yc = r*sin(phi);
                    particles[n].radius = MU;
                    particles[n].active = 1;
                }
            break;
        }
        case (C_CLOAK_A):   /* optimized model A1 by C. Jo et al */
        {
            ncircles = 200;
            ra[0] = 0.0731;     sa[0] = 1.115;
            ra[1] = 0.0768;     sa[1] = 1.292;
            ra[2] = 0.0652;     sa[2] = 1.464;
            ra[3] = 0.056;      sa[3] = 1.633;
            ra[4] = 0.0375;     sa[4] = 1.794;
            for (i = 0; i < 40; i++)
                for (j = 0; j < 5; j++)
                {
                    n = 5*i + j;
                    phi = (double)i*DPI/40.0;
                    r = LAMBDA*sa[j];
                    particles[n].xc = r*cos(phi);
                    particles[n].yc = r*sin(phi);
                    particles[n].radius = LAMBDA*ra[j];
                    particles[n].active = 1;
                }
            break;
        }
        case (C_LASER):
        {
            ncircles = 17;
            
            xx[0] = 0.5*(X_SHOOTER + X_TARGET);
            xx[1] = LAMBDA - 0.5*(X_TARGET - X_SHOOTER);
            xx[2] = -xx[0];
            xx[3] = -xx[1];
            
            yy[0] = 0.5*(Y_SHOOTER + Y_TARGET);
            yy[1] = 1.0 - 0.5*(Y_TARGET - Y_SHOOTER);
            yy[2] = -yy[0];
            yy[3] = -yy[1];

            for (i = 0; i < 4; i++)
                for (j = 0; j < 4; j++)
                {
                    particles[4*i + j].xc = xx[i];
                    particles[4*i + j].yc = yy[j];
                    
                }
                
            particles[ncircles - 1].xc = X_TARGET;
            particles[ncircles - 1].yc = Y_TARGET;
            
            for (i=0; i<ncircles - 1; i++)
            {
                particles[i].radius = MU;
                particles[i].active = 1;
            }
            
            particles[ncircles - 1].radius = 0.5*MU;
            particles[ncircles - 1].active = 2;
            
            break;
        }        
        case (C_POISSON_DISC):
        {
            printf("Generating Poisson disc sample\n");
            /* generate first circle */
//             particles[0].xc = LAMBDA*(2.0*(double)rand()/RAND_MAX - 1.0);
            particles[0].xc = (INITXMAX - INITXMIN)*(double)rand()/RAND_MAX + INITXMIN;
            particles[0].yc = (INITYMAX - INITYMIN)*(double)rand()/RAND_MAX + INITYMIN;
            active_poisson[0] = 1;
//             particles[0].active = 1;
            n_p_active = 1;
            ncircles = 1;
            
            while ((n_p_active > 0)&&(ncircles < NMAXCIRCLES))
            {
                /* randomly select an active circle */
                i = rand()%(ncircles);
                while (!active_poisson[i]) i = rand()%(ncircles);                 
//                 printf("Starting from circle %i at (%.3f,%.3f)\n", i, particles[i].xc, particles[i].yc);
                /* generate new candidates */
                naccepted = 0;
                for (j=0; j<ncandidates; j++)
                {
                    r = dpoisson*(2.0*(double)rand()/RAND_MAX + 1.0);
                    phi = DPI*(double)rand()/RAND_MAX;
                    x = particles[i].xc + r*cos(phi);
                    y = particles[i].yc + r*sin(phi);
//                        printf("Testing new circle at (%.3f,%.3f)\t", x, y);
                    far = 1;
                    for (k=0; k<ncircles; k++) if ((k!=i))
                    {
                        /* new circle is far away from circle k */
                        far = far*((x - particles[k].xc)*(x - particles[k].xc) + (y - particles[k].yc)*(y - particles[k].yc) >=     dpoisson*dpoisson);
                        /* new circle is in domain */
                        far = far*(x < INITXMAX)*(x > INITXMIN)*(y < INITYMAX)*(y > INITYMIN);
//                         far = far*(vabs(x) < LAMBDA)*(y < INITYMAX)*(y > INITYMIN);
                    }
                    if (far)    /* accept new circle */
                    {
                        printf("New circle at (%.3f,%.3f) accepted\n", x, y);
                        particles[ncircles].xc = x;
                        particles[ncircles].yc = y;
                        particles[ncircles].radius = MU;
                        particles[ncircles].active = 1;
                        active_poisson[ncircles] = 1;
                        ncircles++;
                        n_p_active++;
                        naccepted++;
                    }
//                        else printf("Rejected\n");
                }
                if (naccepted == 0)    /* inactivate circle i */ 
                {
//                     printf("No candidates work, inactivate circle %i\n", i);
                    active_poisson[i] = 0;
                    n_p_active--;
                }
                printf("%i active circles\n", n_p_active);
            }
            
            printf("Generated %i circles\n", ncircles);
            break;
        }
        case (C_GOLDEN_MEAN):
        {
            ncircles = 300;
            gamma = (sqrt(5.0) - 1.0)*0.5;    /* golden mean */
            height = YMAX - YMIN;
            dx = 2.0*LAMBDA/((double)ncircles);
            for (n = 0; n < ncircles; n++)
            {
                particles[n].xc = -LAMBDA + n*dx;
                particles[n].yc = y;
                y += height*gamma; 
                if (y > YMAX) y -= height;
                particles[n].radius = MU;
                particles[n].active = 1;
            }
            
            /* test for circles that overlap top or bottom boundary */
            ncirc0 = ncircles;
            for (n=0; n < ncirc0; n++)
            {
                if (particles[n].yc + particles[n].radius > YMAX)
                {
                    particles[ncircles].xc = particles[n].xc;
                    particles[ncircles].yc = particles[n].yc - height;
                    particles[ncircles].radius = MU;
                    particles[ncircles].active = 1;
                    ncircles ++;
                }
                else if (particles[n].yc - particles[n].radius < YMIN)
                {
                    particles[ncircles].xc = particles[n].xc;
                    particles[ncircles].yc = particles[n].yc + height;
                    particles[ncircles].radius = MU;
                    particles[ncircles].active = 1;
                    ncircles ++;
                }
            }
            break;
        }
        case (C_GOLDEN_SPIRAL):
        {
            ncircles = 1;
            particles[0].xc = 0.0;
            particles[0].yc = 0.0;
            
            gamma = (sqrt(5.0) - 1.0)*PI;    /* golden mean times 2Pi */
            phi = 0.0;
            r0 = 2.0*MU;
            r = r0 + MU;
            
            for (i=0; i<1000; i++) 
            {
                x = r*cos(phi);
                y = r*sin(phi);
                
                phi += gamma;
                r += MU*r0/r;
                
                if ((vabs(x) < LAMBDA)&&(vabs(y) < YMAX + MU))
                {
                    particles[ncircles].xc = x;
                    particles[ncircles].yc = y;
                    ncircles++;
                }
            }
            
            for (i=0; i<ncircles; i++)
            {
                particles[i].radius = MU;
                /* inactivate circles outside the domain */
                if ((particles[i].yc < YMAX + MU)&&(particles[i].yc > YMIN - MU)) particles[i].active = 1;
//                 printf("i = %i, circlex = %.3lg, circley = %.3lg\n", i, particles[i].xc, particles[i].yc);
            }
        break;
        }
        case (C_SQUARE_HEX):
        {
            ncircles = NGRIDX*(NGRIDY+1);
            dy = (YMAX - YMIN)/((double)NGRIDY);
            dx = dy*0.5*sqrt(3.0);
            for (i = 0; i < NGRIDX; i++)
                for (j = 0; j < NGRIDY+1; j++)
                {
                    n = (NGRIDY+1)*i + j;
                    particles[n].xc = ((double)(i-NGRIDX/2) + 0.5)*dy;   /* is +0.5 needed? */
                    particles[n].yc = YMIN + ((double)j - 0.5)*dy;
                    if (((i+NGRIDX)%4 == 2)||((i+NGRIDX)%4 == 3)) particles[n].yc += 0.5*dy;
                    particles[n].radius = MU;
                    /* activate only circles that intersect the domain */
                    if ((particles[n].yc < YMAX + MU)&&(particles[n].yc > YMIN - MU)) particles[n].active = 1;
                    else particles[n].active = 0;
                }
            break;
        }
        case (C_POOL_TABLE):
        {
            for (i=1; i<6; i++) for (j=0; j<i; j++)
            {
                particles[ncircles].xc = INITXMIN + (double)i*0.25*(INITXMAX - INITXMIN);
                particles[ncircles].yc = 0.5*(INITYMIN + INITYMAX) + ((double)j - 0.5*(double)(i-1))*0.25*(INITYMAX - INITYMIN);
                particles[ncircles].radius = MU;
                particles[ncircles].active = 1;
                ncircles += 1;
            }
            break;
        }
        case (C_ONE):
        {
            particles[ncircles].xc = 0.0;
            particles[ncircles].yc = 0.0;
            particles[ncircles].radius = MU;
            particles[ncircles].active = 1;
            ncircles += 1;
            break;
        }
        case (C_TWO):   /* used for comparison with cloak */
        {
            particles[ncircles].xc = 0.0;
            particles[ncircles].yc = 0.0;
            particles[ncircles].radius = MU;
            particles[ncircles].active = 2;
            ncircles += 1;

            particles[ncircles].xc = 0.0;
            particles[ncircles].yc = 0.0;
            particles[ncircles].radius = 2.0*MU;
            particles[ncircles].active = 1;
            ncircles += 1;
            break;
        }
        case (C_NOTHING):
        {
            ncircles += 0;
            break;
        }
        default: 
        {
            printf("Function init_circle_config not defined for this pattern \n");
        }
    }
}

void add_particle_config(t_particle particles[NMAXCIRCLES], double xmin, double xmax, double ymin, double ymax, double radius)
/* add particles to configuration */
{
    int i, j, k, n, ncirc0, n_p_active, ncandidates = PDISC_CANDIDATES, naccepted, newcircles; 
    double dx, dy, p, phi, r, r0, ra[5], sa[5], height, x, y = 0.0, gamma, dpoisson, xx[4], yy[4];
    short int active_poisson[NMAXCIRCLES], far;
    
    dpoisson = PDISC_DISTANCE*radius;
    
    switch (CIRCLE_PATTERN_B) {
        case (C_SQUARE):
        {
            ncircles = NGRIDX*NGRIDY;
            dy = (YMAX - YMIN)/((double)NGRIDY);
            for (i = 0; i < NGRIDX; i++)
                for (j = 0; j < NGRIDY; j++)
                {
                    n = NGRIDY*i + j;
                    particles[n].xc = ((double)(i-NGRIDX/2) + 0.5)*dy;
                    particles[n].yc = YMIN + ((double)j + 0.5)*dy;
                    particles[n].radius = MU;
                    particles[n].active = 1;
                }
            break;
        }
        case (C_HEX):
        {
            dx = (xmax - xmin)/((double)NGRIDX);
            dy = (ymax - ymin)/((double)NGRIDY);
//             dx = dy*0.5*sqrt(3.0);
            for (i = 0; i < NGRIDX; i++)
                for (j = 0; j < NGRIDY+1; j++)
                {
                    n = ncircles + (NGRIDY+1)*i + j;
//                     particles[n].xc = ((double)(i-NGRIDX/2) + 0.5)*dx;   /* is +0.5 needed? */
                    particles[n].xc = xmin + ((double)i - 0.5)*dx;   
                    particles[n].yc = ymin + ((double)j - 0.5)*dy;
                    if ((i+NGRIDX)%2 == 1) particles[n].yc += 0.5*dy;
                    particles[n].radius = radius;
                    /* activate only circles that intersect the domain */
                    if ((particles[n].yc < ymax + radius)&&(particles[n].yc > ymin - radius)&&(particles[n].xc < xmax + radius)&&(particles[n].xc > xmin - radius)) particles[n].active = 1;
                    else particles[n].active = 0;
                }
            ncircles += NGRIDX*(NGRIDY+1);
            break;
        }
        case (C_POISSON_DISC):
        {
            ncirc0 = ncircles;
            printf("Generating new Poisson disc sample\n");
            /* generate first circle */
            particles[ncirc0].xc = (xmax - xmin)*(double)rand()/RAND_MAX + xmin;
            particles[ncirc0].yc = (ymax - ymin)*(double)rand()/RAND_MAX + ymin;
            active_poisson[0] = 1;
// //             particles[0].active = 1;
            n_p_active = 1;
            newcircles = 1;
            
            while ((n_p_active > 0)&&(ncircles < NMAXCIRCLES))
            {
                /* randomly select an active circle */
                i = rand()%(newcircles);
                while (!active_poisson[i]) i = rand()%(ncircles);                 
//                 printf("Starting from circle %i at (%.3f,%.3f)\n", i, particles[i].xc, particles[i].yc);
                /* generate new candidates */
                naccepted = 0;
                for (j=0; j<ncandidates; j++)
                {
                    r = dpoisson*(2.0*(double)rand()/RAND_MAX + 1.0);
                    phi = DPI*(double)rand()/RAND_MAX;
                    x = particles[ncirc0 + i].xc + r*cos(phi);
                    y = particles[ncirc0 + i].yc + r*sin(phi);
//                        printf("Testing new circle at (%.3f,%.3f)\t", x, y);
                    far = 1;
                    for (k=0; k<ncircles; k++) if ((k!=i))
                    {
                        /* new circle is far away from circle k */
                        far = far*((x - particles[k].xc)*(x - particles[k].xc) + (y - particles[k].yc)*(y - particles[k].yc) >= dpoisson*dpoisson);
                        /* new circle is in domain */
                        far = far*(x < xmax)*(x > xmin)*(y < ymax)*(y > ymin);
//                         far = far*(vabs(x) < LAMBDA)*(y < INITYMAX)*(y > INITYMIN);
                    }
                    if (far)    /* accept new circle */
                    {
                        printf("New circle at (%.3f,%.3f) accepted\n", x, y);
                        particles[ncircles].xc = x;
                        particles[ncircles].yc = y;
                        particles[ncircles].radius = radius;
                        particles[ncircles].active = 1;
                        active_poisson[ncircles] = 1;
                        ncircles++;
                        ncirc0++;
                        n_p_active++;
                        naccepted++;
                    }
//                        else printf("Rejected\n");
                }
                if (naccepted == 0)    /* inactivate circle i */ 
                {
//                     printf("No candidates work, inactivate circle %i\n", i);
                    active_poisson[i] = 0;
                    n_p_active--;
                }
                printf("%i active circles\n", n_p_active);
            }
            
            printf("Generated %i circles\n", ncircles);
            break;
        }
        default: 
        {
            printf("Function init_circle_config not defined for this pattern \n");
        }
    }
}

void init_people_config(t_person people[NMAXCIRCLES])
/* initialise particle configuration */
{
    t_particle particles[NMAXCIRCLES];
    int n;
    
    init_particle_config(particles);
    
    for (n=0; n<ncircles; n++)
    {
        people[n].xc = particles[n].xc;
        people[n].yc = particles[n].yc;
        people[n].radius = particles[n].radius;
        people[n].active = particles[n].active;
    }
    
}


void add_obstacle(double x, double y, double radius, t_obstacle obstacle[NMAXOBSTACLES])
/* add a circular obstacle to obstacle configuration */
{
    if (nobstacles + 1 < NMAXOBSTACLES)
    {
        obstacle[nobstacles].xc = x;
        obstacle[nobstacles].yc = y;
        obstacle[nobstacles].radius = radius;
        obstacle[nobstacles].active = 1;
        
        nobstacles++;
    }
    else printf("Warning: NMAXOBSTACLES should be increased\n");
}


void init_obstacle_config(t_obstacle obstacle[NMAXOBSTACLES])
/* initialise circular obstacle configuration */
{
    int i, j, n; 
    double x, y, dx, dy, width, lpocket, xmid = 0.5*(BCXMIN + BCXMAX), radius;
    
    switch (OBSTACLE_PATTERN) {
        case (O_CORNERS):
        {
            n = 0;
            for (i = 0; i < 2; i++)
                for (j = 0; j < 2; j++)
                {
                    obstacle[n].xc = BCXMIN + ((double)i)*(BCXMAX - BCXMIN);
                    obstacle[n].yc = BCYMIN + ((double)j)*(BCYMAX - BCYMIN);
                    obstacle[n].radius = OBSTACLE_RADIUS;
                    obstacle[n].active = 1;
                    n++;
                }
            nobstacles = n;
            break;
        }
        case (O_GALTON_BOARD):
        {
            dy = (YMAX - YMIN)/((double)NGRIDX + 3);
            dx = dy/cos(PI/6.0);
            n = 0;
            for (i = 0; i < NGRIDX + 1; i++)
                for (j = 0; j < i; j++)
                {
                    obstacle[n].yc = YMAX - ((double)i)*dy;
                    obstacle[n].xc = ((double)j - 0.5*(double)i + 0.5)*dx;
                    obstacle[n].radius = OBSTACLE_RADIUS;
                    obstacle[n].active = 1;
                    n++;
                }
            nobstacles = n;
            break;
        }
        case (O_GENUS_TWO):
        {
            n = 0;
            for (i = 0; i < 3; i++)
                for (j = 0; j < 3; j++)
                {
                    obstacle[n].xc = BCXMIN + 0.5*((double)i)*(BCXMAX - BCXMIN);
                    obstacle[n].yc = BCYMIN + 0.5*((double)j)*(BCYMAX - BCYMIN);
                    obstacle[n].radius = OBSTACLE_RADIUS;
                    obstacle[n].active = 1;
                    n++;
                }
            nobstacles = n;
            break;
        }
        case (O_POOL_TABLE):
        {
            lpocket = 0.1;
            width = 0.5*MU;
            radius = 2.0*width;
            
            add_obstacle(BCXMIN + lpocket, BCYMIN - width, radius, obstacle);
            add_obstacle(xmid - lpocket, BCYMIN - width, radius, obstacle);

            add_obstacle(xmid + lpocket, BCYMIN - width, radius, obstacle);
            add_obstacle(BCXMAX - lpocket, BCYMIN - width, radius, obstacle);
            
            add_obstacle(BCXMAX + width, BCYMIN + lpocket, radius, obstacle);
            add_obstacle(BCXMAX + width, BCYMAX - lpocket, radius, obstacle);
            
            add_obstacle(BCXMIN + lpocket, BCYMAX + width, radius, obstacle);
            add_obstacle(xmid - lpocket, BCYMAX + width, radius, obstacle);

            add_obstacle(xmid + lpocket, BCYMAX + width, radius, obstacle);
            add_obstacle(BCXMAX - lpocket, BCYMAX + width, radius, obstacle);
            
            add_obstacle(BCXMIN - width, BCYMIN + lpocket, radius, obstacle);
            add_obstacle(BCXMIN - width, BCYMAX - lpocket, radius, obstacle);
            
            break;
        }
        case (O_HLINE_HOLE_SPOKES):
        {
            radius = 2.0*MU;
            
            for (i=-3; i<4; i+=2)
                add_obstacle(0.5*(double)i, YMIN + 0.3, radius, obstacle);
            break;
        }
        case (O_CIRCLE):
        {
            n = 0;
            obstacle[n].xc = 0.0;
            obstacle[n].yc = 0.0;
            obstacle[n].radius = OBSTACLE_RADIUS;
            obstacle[n].active = 1;
            nobstacles = 1;
            break;
        }case (O_FOUR_CIRCLES):
        {
            n = 0;
            
            obstacle[n].xc = -1.5;
            obstacle[n].yc = -0.5;
            obstacle[n].radius = OBSTACLE_RADIUS;
            obstacle[n].active = 1;
            n++;

            obstacle[n].xc = -0.5;
            obstacle[n].yc = 0.5;
            obstacle[n].radius = OBSTACLE_RADIUS;
            obstacle[n].active = 1;
            n++;

            obstacle[n].xc = 0.5;
            obstacle[n].yc = -0.5;
            obstacle[n].radius = OBSTACLE_RADIUS;
            obstacle[n].active = 1;
            n++;

            obstacle[n].xc = 1.5;
            obstacle[n].yc = 0.5;
            obstacle[n].radius = OBSTACLE_RADIUS;
            obstacle[n].active = 1;
            n++;
            
            nobstacles = 4;
            break;
        }
        
        default: 
        {
            printf("Function init_obstacle_config not defined for this pattern \n");
        }
    }
}

void add_rotated_angle_to_segments(double x1, double y1, double x2, double y2, double width, int center, t_segment segment[NMAXSEGMENTS], int group)
/* add four segments forming a rectangle, specified by two adjacent corners and width */
{
    double tx, ty, ux, uy, norm, x3, y3, x4, y4;
    int i, n = nsegments; 
    
    tx = x2 - x1;
    ty = y2 - y1;
    norm = module2(tx, ty);
    tx = tx/norm;
    ty = ty/norm;
    if (center)
    {
        x2 -= 0.5*width*ty;
        y2 += 0.5*width*tx;
        x1 -= 0.5*width*ty;
        y1 += 0.5*width*tx;
    }
    x3 = x2 + width*ty;
    y3 = y2 - width*tx;
    x4 = x1 + width*ty;
    y4 = y1 - width*tx;
    
    if (nsegments + 4 < NMAXSEGMENTS)
    {
        segment[n].x1 = x1;
        segment[n].y1 = y1;
        segment[n].x2 = x2;
        segment[n].y2 = y2;
        
        segment[n+1].x1 = x2;
        segment[n+1].y1 = y2;
        segment[n+1].x2 = x3;
        segment[n+1].y2 = y3;
        
        segment[n+2].x1 = x3;
        segment[n+2].y1 = y3;
        segment[n+2].x2 = x4;
        segment[n+2].y2 = y4;
        
        segment[n+3].x1 = x4;
        segment[n+3].y1 = y4;
        segment[n+3].x2 = x1;
        segment[n+3].y2 = y1;
        
        for (i=0; i<4; i++) 
        {
            segment[n+i].concave = 1;
            segment[n+i].group = group;
        }
        nsegments += 4;
    }
    else printf("Warning: NMAXSEGMENTS too small\n");
}
 
void add_rectangle_to_segments(double x1, double y1, double x2, double y2, t_segment segment[NMAXSEGMENTS], int group)
/* add four segements forming a rectangle to linear obstacle configuration */
{
    int i, n = nsegments, nplus, nminus; 
    
    if (nsegments + 4 < NMAXSEGMENTS)
    {
        segment[n].x1 = x1;
        segment[n].y1 = y1;
        segment[n].x2 = x2;
        segment[n].y2 = y1;
        
        segment[n+1].x1 = x2;
        segment[n+1].y1 = y1;
        segment[n+1].x2 = x2;
        segment[n+1].y2 = y2;
        
        segment[n+2].x1 = x2;
        segment[n+2].y1 = y2;
        segment[n+2].x2 = x1;
        segment[n+2].y2 = y2;
        
        segment[n+3].x1 = x1;
        segment[n+3].y1 = y2;
        segment[n+3].x2 = x1;
        segment[n+3].y2 = y1;
        
        segment[n].angle1 = -PID;
        segment[n].angle2 = 0.0;
        
        segment[n+1].angle1 = PI;
        segment[n+1].angle2 = 1.5*PI;
        
        segment[n+2].angle1 = PID;
        segment[n+2].angle2 = PI;
        
        segment[n+3].angle1 = 0.0;
        segment[n+3].angle2 = PID;
        
        for (i=0; i<4; i++) 
        {
            segment[n+i].concave = 1;
            segment[n+i].group = group;
        }
        
        nsegments += 4;
    }
    else printf("Warning: NMAXSEGMENTS too small\n");
}


void add_circle_to_segments(double x, double y, double r, int nsegs, double angle0, t_segment segment[NMAXSEGMENTS], int group)
/* add segments forming a circle/polygon to linear obstacle configuration */
{
    int i, n = nsegments, nplus, nminus; 
    double angle;
    
    angle = DPI/(double)nsegs;
    
    if (nsegments + nsegs < NMAXSEGMENTS) 
    {
        for (i=0; i<nsegs; i++)
        {
            segment[n+i].x1 = x + r*cos(((double)i)*angle + angle0*PID);
            segment[n+i].y1 = y - r*sin(((double)i)*angle + angle0*PID);
            segment[n+i].x2 = x + r*cos(((double)(i+1))*angle + angle0*PID);
            segment[n+i].y2 = y - r*sin(((double)(i+1))*angle + angle0*PID);
            segment[n+i].angle1 = -((double)i + 0.5)*angle - angle0*PID;
            segment[n+i].angle2 = -((double)i - 0.5)*angle - angle0*PID;
            while (segment[n+i].angle1 < 0.0) segment[n+i].angle1 += DPI;
            while (segment[n+i].angle2 < segment[n+i].angle1) segment[n+i].angle2 += DPI;
            segment[n+i].concave = 1;
            segment[n+i].group = group;
        }
    
        nsegments += nsegs;
    }
    else printf("Warning: NMAXSEGMENTS too small\n");
}
        
        


double nozzle_width(double x, double width, int nozzle_shape)
/* width of bell-shaped nozzle */
{
    double lam  = 0.5*LAMBDA, a, b;
    
    if (x >= 0.0) return(width);
    else switch (nozzle_shape) {
        case (NZ_STRAIGHT): return(width);
        case (NZ_BELL): return(sqrt(width*width - 0.5*x));
        case (NZ_GLAS): return(sqrt(width*width - 1.2*x) + 1.0*x);
        case (NZ_CONE): return(width - (sqrt(width*width + 0.5) - width)*x);
        case (NZ_TRUMPET): return(width + (sqrt(width*width + LAMBDA)-width)*x*x);
        case (NZ_BROAD): 
        {
            if (-x < 0.1) return(width - (0.5 - width)*x/0.1);
            else return(0.5);
        }
        case (NZ_DELAVAL): 
        {
            a = 1.5;
            b = 0.05;
            return(sqrt(width*width - 0.5*x) + a*b*x*(1.0 + x)/(b + x*x));
//             a = (sqrt(width*width+0.5) - width)/sqrt(0.5);
//             c = (a*sqrt(0.5*h) - width)/(h*h);
//             if (-x < h) return(width + c*x*x);
//             else return(width + a * sqrt(-0.5*x));
        }
        default: return(0.0);
    }
}


void add_rocket_to_segments(t_segment segment[NMAXSEGMENTS], double x0, double y0, int rocket_shape, int nozzle_shape, int nsides, int group)
/* add one or several rocket_shaped set of segments */
{
    int i, j, cycle = 0, nsegments0; 
    double angle, dx, x1, y1, x2, y2, nozx, nozy, a, b, c, w;
    
    nsegments0 = nsegments;

    /* compute intersection point of nozzle and ellipse */
    angle = -PID + DPI/(double)NPOLY;
    nozx = 0.7*LAMBDA*cos(angle);
    nozy = y0 + YMIN + LAMBDA*(1.7 + 0.7*sin(angle));
            
    /* form of combustion chamber */
    switch (rocket_shape) {
        case (RCK_DISC):    /* circular chamber */
        {
            for (i=1; i<NPOLY-1; i++)
            {
                angle = -PID + (double)i*DPI/(double)NPOLY;
                x1 = x0 + 0.7*LAMBDA*cos(angle);
                y1 = y0 + YMIN + LAMBDA*(1.7 + 0.7*sin(angle));
                angle = -PID + (double)(i+1)*DPI/(double)NPOLY;
                x2 = x0 + 0.7*LAMBDA*cos(angle);
                y2 = y0 + YMIN + LAMBDA*(1.7 + 0.7*sin(angle));
                add_rotated_angle_to_segments(x1, y1, x2, y2, 0.02, 0, segment, 0);
            }
            break;
        }
        case (RCK_RECT):    /* rectangular chamber */
        {
            /* dimensions chosen to have same area as circular chamber */
            a = 0.5*LAMBDA;
            b = 0.49*PI*LAMBDA;
            add_rotated_angle_to_segments(x0+nozx, nozy, x0+a, nozy, 0.02, 0, segment, 0);
            add_rotated_angle_to_segments(x0+a, nozy, x0+a, nozy+b, 0.02, 0, segment, 0);
            add_rotated_angle_to_segments(x0+a, nozy+b, x0-a, nozy+b, 0.02, 0, segment, 0);
            add_rotated_angle_to_segments(x0-a, nozy+b, x0-a, nozy, 0.02, 0, segment, 0);
            add_rotated_angle_to_segments(x0-a, nozy, x0-nozx, nozy, 0.02, 0, segment, 0);
            break;
        }
        case (RCK_RECT_HAT):    /* rectangular chamber with a hat */
        {
            a = 0.5*LAMBDA;
            b = (0.49*PI-0.25)*LAMBDA;
            add_rotated_angle_to_segments(x0+nozx, nozy, x0+a, nozy, 0.02, 0, segment, 0);
            add_rotated_angle_to_segments(x0+a, nozy, x0+a, nozy+b, 0.02, 0, segment, 0);
            add_rotated_angle_to_segments(x0+a, nozy+b, x0, nozy+b+a, 0.02, 0, segment, 0);
            add_rotated_angle_to_segments(x0, nozy+b+a, x0-a, nozy+b, 0.02, 0, segment, 0);
            add_rotated_angle_to_segments(x0-a, nozy+b, x0-a, nozy, 0.02, 0, segment, 0);
            add_rotated_angle_to_segments(x0-a, nozy, x0-nozx, nozy, 0.02, 0, segment, 0);
            break;
        }
        case (RCK_RECT_BAR):    /* rectangular chamber with a hat and separating bar */
        {
            a = 0.5*LAMBDA;
            b = (0.49*PI-0.25)*LAMBDA;
            c = 0.5*a;
            w = 0.025;
            add_rotated_angle_to_segments(x0+nozx, nozy, x0+nozx, nozy+0.5*c, 0.02, 0, segment, 0);
            add_rotated_angle_to_segments(x0+nozx, nozy+0.5*c, x0+nozx+0.5*c, nozy+c, 0.02, 0, segment, 0);
            add_rotated_angle_to_segments(x0+nozx+0.5*c, nozy+c, x0+a, nozy+c, 0.02, 0, segment, 0);
            add_rotated_angle_to_segments(x0+a, nozy+c, x0+a, nozy+b+c, 0.02, 0, segment, 0);
            add_rotated_angle_to_segments(x0+a, nozy+b+c, x0, nozy+b+a+c, 0.02, 0, segment, 0);
            add_rotated_angle_to_segments(x0-w, nozy+a+c, x0-w, nozy+b+a+c, 2.0*w, 0, segment, 0);
            add_rotated_angle_to_segments(x0, nozy+b+a+c, x0-a, nozy+b+c, 0.02, 0, segment, 0);
            add_rotated_angle_to_segments(x0-a, nozy+b+c, x0-a, nozy+c, 0.02, 0, segment, 0);
            add_rotated_angle_to_segments(x0-a, nozy+c, x0-nozx-0.5*c, nozy+c, 0.02, 0, segment, 0);
            add_rotated_angle_to_segments(x0-nozx-0.5*c, nozy+c, x0-nozx+w, nozy+0.5*c, 0.02, 0, segment, 0);
            add_rotated_angle_to_segments(x0-nozx+w, nozy+0.5*c, x0-nozx+w, nozy, 0.02, 0, segment, 0);
            break;
        }
    }
            
    dx = LAMBDA/(double)(nsides);
    
    /* nozzle */
    if (nozzle_shape != NZ_NONE)
    {
        for (i=0; i<nsides; i++)
        {
            y1 = y0 - LAMBDA + dx*(double)(i-1);
            x1 = x0 + nozzle_width(y1 - y0, nozx, nozzle_shape);
            y2 = y1 + dx;
            x2 = x0 + nozzle_width(y2 - y0, nozx, nozzle_shape);
            add_rotated_angle_to_segments(x1, y1 + YMIN + LAMBDA, x2, y2 + YMIN + LAMBDA, 0.02, 0, segment, 0);
        }
        add_rotated_angle_to_segments(x2, y2 + YMIN + LAMBDA, x0 + nozx, nozy, 0.02, 0, segment, 0);
        for (i=0; i<nsides; i++)
        {
            y1 = y0 - LAMBDA + dx*(double)(i-1);
            x1 = x0 - nozzle_width(y1 - y0, nozx, nozzle_shape);
            y2 = y1 + dx;
            x2 = x0 - nozzle_width(y2 - y0, nozx, nozzle_shape);
            add_rotated_angle_to_segments(x1, y1 + YMIN + LAMBDA, x2, y2 + YMIN + LAMBDA, 0.02, 0, segment, 0);
        }
        add_rotated_angle_to_segments(x2, y2 + YMIN + LAMBDA, x0 - nozx, nozy, 0.02, 0, segment, 0);
    }
    
    for (i=nsegments0; i<nsegments; i++) segment[i].inactivate = 0;
    
    /* closing segment */
    segment[nsegments].x1 = x0 - nozx;
    segment[nsegments].y1 = nozy;
    segment[nsegments].x2 = x0 + nozx;
    segment[nsegments].y2 = nozy;
    segment[nsegments].inactivate = 1;
    nsegments++;
    
    /* set group of segments */
    for (i=nsegments0; i<nsegments; i++) segment[i].group = group;
}
 
int init_maze_segments(t_segment segment[NMAXSEGMENTS], int diag)
/* init segments forming a maze */
{
    t_maze* maze;
    int i, j, n;
    double x1, y1, x2, y2, dx, dy, padding = 0.02, width = MAZE_WIDTH;
    
    maze = (t_maze *)malloc(NXMAZE*NYMAZE*sizeof(t_maze));
    
    init_maze(maze);
        
    /* build walls of maze */
    dx = (YMAX - YMIN - 2.0*padding)/(double)(NXMAZE);
    dy = (YMAX - YMIN - 2.0*padding)/(double)(NYMAZE);
    
    for (i=0; i<NXMAZE; i++)
        for (j=0; j<NYMAZE; j++)
        {
            n = nmaze(i, j);
            x1 = YMIN + padding + (double)i*dx + MAZE_XSHIFT;
            y1 = YMIN + padding + (double)j*dy;
            
            if (diag)
            {
                if (((i>0)||(j<NYMAZE-1))&&(maze[n].west)) add_rectangle_to_segments(x1, y1, x1 - width, y1 + dy, segment, 0);
            }
            else if (((i>0)||(j!=NYMAZE/2))&&(maze[n].west)) add_rectangle_to_segments(x1, y1, x1 - width, y1 + dy, segment, 0);
            if (maze[n].south) add_rectangle_to_segments(x1, y1, x1 + dx, y1 - width, segment, 0);            
        }
    
    /* top side of maze */
    add_rectangle_to_segments(YMIN + padding + MAZE_XSHIFT, YMAX - padding, YMAX - padding + MAZE_XSHIFT, YMAX - padding - width, segment, 0);
    
    /* right side of maze */
    x1 = YMAX - padding + MAZE_XSHIFT;
    if (diag)
    {
        y1 = YMIN + padding + dy;
        add_rectangle_to_segments(x1, y1, x1 - width, YMAX + 10.0, segment, 0);
    }
    else
    {
        y1 = YMIN + padding + dy*((double)NYMAZE/2);
        add_rectangle_to_segments(x1, YMIN - 1.0, x1 - width, y1 - dy, segment, 0);
        add_rectangle_to_segments(x1, y1, x1 - width, YMAX + 1.0, segment, 0);
    }
    
    /* left side of maze */
    x1 = YMIN + padding + MAZE_XSHIFT;
    add_rectangle_to_segments(x1, YMIN - 1.0, x1 - width, YMIN + padding, segment, 0);
    add_rectangle_to_segments(x1, YMAX - padding, x1 - width, YMAX + 10.0, segment, 0);
    
    if (diag) 
    {
        add_rotated_angle_to_segments(XMIN, YMAX - 0.5*dy, x1, YMAX - dy - 2.0*width, width, 0, segment, 0);
        add_rectangle_to_segments(XMIN, YMAX - 0.5*dy, XMIN - width, YMAX + 10.0, segment, 0);
    }
    
    free(maze);
}


void init_segment_config(t_segment segment[NMAXSEGMENTS])
/* initialise linear obstacle configuration */
{
    int i, j, cycle = 0, iminus, iplus, nsides, n, concave = 1; 
    double angle, angle2, dangle, dx, width, height, a, b, length, xmid = 0.5*(BCXMIN + BCXMAX), lpocket, r, x, x1, y1, x2, y2, nozx, nozy, y, dy, ca, sa;
    
    switch (SEGMENT_PATTERN) {
        case (S_RECTANGLE):
        {
            segment[0].x1 = BCXMIN;
            segment[0].y1 = BCYMAX;
            
            segment[1].x1 = BCXMIN;
            segment[1].y1 = BCYMIN;

            segment[2].x1 = BCXMAX;
            segment[2].y1 = BCYMIN;

            segment[3].x1 = BCXMAX;
            segment[3].y1 = BCYMAX;
            
            cycle = 1;
            nsegments = 4;
            
            for (i=0; i<nsegments; i++) 
            {
                segment[i].concave = 0;
                segment[i].group = 0;
                segment[i].inactivate = 0;
            }
            break;
        }
        case (S_CUP):
        {
            angle = APOLY*PID;
            dx = (BCYMAX - BCYMIN)/tan(angle);
            
            segment[0].x1 = BCXMIN;
            segment[0].y1 = BCYMAX;
            
            segment[1].x1 = BCXMIN + dx;
            segment[1].y1 = BCYMIN;

            segment[2].x1 = BCXMAX - dx;
            segment[2].y1 = BCYMIN;

            segment[3].x1 = BCXMAX;
            segment[3].y1 = BCYMAX;
            
            cycle = 1;
            nsegments = 4;
            
            for (i=0; i<nsegments; i++) 
            {
                segment[i].concave = 0;
                segment[i].group = 0;
                segment[i].inactivate = 0;
            }
            
            break;
        }
        case (S_HOURGLASS):
        {
            angle = APOLY*PID;
            width = 2.5*MU;
            height = 2.5*MU;
            
            segment[0].x1 = BCXMIN;
            segment[0].y1 = BCYMAX;
            segment[0].concave = 0;
            
            segment[1].x1 = -width;
            segment[1].y1 = height;
            segment[1].concave = 1;
            
            segment[2].x1 = -width;
            segment[2].y1 = -height;
            segment[2].concave = 1;
            
            segment[3].x1 = BCXMIN;
            segment[3].y1 = BCYMIN;
            segment[3].concave = 0;
            
            segment[4].x1 = BCXMAX;
            segment[4].y1 = BCYMIN;
            segment[4].concave = 0;
            
            segment[5].x1 = width;
            segment[5].y1 = -height;
            segment[5].concave = 1;
            
            segment[6].x1 = width;
            segment[6].y1 = height;
            segment[6].concave = 1;
            
            segment[7].x1 = BCXMAX;
            segment[7].y1 = BCYMAX;
            segment[7].concave = 0;
            
            cycle = 1;
            nsegments = 8;
            for (i=0; i<nsegments; i++) 
            {
                segment[i].group = 0;
                segment[i].inactivate = 0;
            }
            
            break;
        }
        case (S_PENTA):
        {
            height = 0.5*(BCYMAX - BCYMIN);
            width = height/sqrt(3.0);
            
            segment[0].x1 = BCXMIN;
            segment[0].y1 = 0.5*(BCYMIN + BCYMAX);
            
            segment[1].x1 = BCXMIN + width;
            segment[1].y1 = BCYMIN;

            segment[2].x1 = BCXMAX;
            segment[2].y1 = BCYMIN;

            segment[3].x1 = BCXMAX;
            segment[3].y1 = BCYMAX;
            
            segment[4].x1 = BCXMIN + width;
            segment[4].y1 = BCYMAX;

            cycle = 1;
            nsegments = 5;
            
            for (i=0; i<nsegments; i++) 
            {
                segment[i].concave = 0;
                segment[i].group = 0;
                segment[i].inactivate = 0;
            }
            
            break;
        }
        case (S_CENTRIFUGE):
        {
            angle = DPI/(double)NPOLY;
            
            for (i=0; i<NPOLY; i++)
            {
                segment[i*4].x1 = LAMBDA*cos(((double)i + 0.02)*angle);
                segment[i*4].y1 = LAMBDA*sin(((double)i + 0.02)*angle);
                segment[i*4].concave = 1;
                
                segment[i*4 + 1].x1 = cos(((double)i + 0.05)*angle);
                segment[i*4 + 1].y1 = sin(((double)i + 0.05)*angle);
                segment[i*4 + 1].concave = 0;
                
                segment[i*4 + 2].x1 = cos(((double)i + 0.95)*angle);
                segment[i*4 + 2].y1 = sin(((double)i + 0.95)*angle);
                segment[i*4 + 2].concave = 0;
                
                segment[i*4 + 3].x1 = LAMBDA*cos(((double)i + 0.98)*angle);
                segment[i*4 + 3].y1 = LAMBDA*sin(((double)i + 0.98)*angle);
                segment[i*4 + 3].concave = 1;
            }
            
            cycle = 1;
            nsegments = 4*NPOLY;
            
            for (i=0; i<nsegments; i++) 
            {
                segment[i].group = 0;
                segment[i].inactivate = 0;
            }
            
            break;
        }
        case (S_POLY_ELLIPSE):
        {
            angle = DPI/(double)NPOLY;
            
            for (i=0; i<NPOLY; i++)
            {
                segment[i].x1 = -cos(((double)i + 0.5)*angle);
                segment[i].y1 = -LAMBDA*sin(((double)i + 0.5)*angle);
                segment[i].concave = 0;
            }
            segment[0].concave = 1;
            segment[NPOLY-1].concave = 1;
            for (i=0; i<NPOLY; i++)
            {
                segment[NPOLY+i].x1 = 1.05*segment[NPOLY-1-i].x1;
                segment[NPOLY+i].y1 = 1.05*segment[NPOLY-1-i].y1;
                segment[NPOLY+i].concave = 1;
            }
            
            cycle = 1;
            nsegments = 2*NPOLY;
            
            for (i=0; i<nsegments; i++) segment[i].inactivate = 0;
            break;
        }
        case (S_POOL_TABLE):
        {
            width = MU;
            lpocket = 0.1;
            
            add_rectangle_to_segments(BCXMIN + lpocket, BCYMIN, xmid - lpocket, BCYMIN - width, segment, 0);
            add_rectangle_to_segments(xmid + lpocket, BCYMIN, BCXMAX - lpocket, BCYMIN - width, segment, 0);
            add_rectangle_to_segments(BCXMAX + width, BCYMIN + lpocket, BCXMAX, BCYMAX - lpocket, segment, 0); 
            
            add_rectangle_to_segments(BCXMAX - lpocket, BCYMAX, xmid + lpocket, BCYMAX + width, segment, 0);
            add_rectangle_to_segments(xmid - lpocket, BCYMAX, BCXMIN + lpocket, BCYMAX + width, segment, 0);
            add_rectangle_to_segments(BCXMIN - width, BCYMAX - lpocket, BCXMIN, BCYMIN + lpocket, segment, 0); 
            
            cycle = 0;
            for (i=0; i<nsegments; i++) 
            {
                segment[i].concave = 0;
                segment[i].group = 0;
                segment[i].inactivate = 0;
            }
            
            break;
        }   
        case (S_CENTRIFUGE_RND):
        {
            angle = DPI/(double)NPOLY;
            nsides = 24;
            if ((nsides+2)*NPOLY > NMAXSEGMENTS)
            {
                printf("Error: NMAXSEGMENTS is too small\n");
                exit(1);
            }
            
            for (i=0; i<NPOLY; i++)
            {
                segment[i*(nsides+2)].x1 = LAMBDA*cos(((double)i + 0.02)*angle);
                segment[i*(nsides+2)].y1 = LAMBDA*sin(((double)i + 0.02)*angle);
                segment[i*(nsides+2)].concave = 1;
                
                for (j=1; j<=nsides; j++)
                {
                    x = (double)j/(double)(nsides+1);
                    r = 0.5 + sqrt(x*(1.0-x));
                    angle2 = (double)i*angle + angle*(double)j/(double)(nsides+1);
                    segment[i*(nsides+2) + j].x1 = r*cos(angle2);
                    segment[i*(nsides+2) + j].y1 = r*sin(angle2);
                    segment[i*(nsides+2) + j].concave = 0;
                }
                
                segment[i*(nsides+2) + nsides + 1].x1 = LAMBDA*cos(((double)i + 0.98)*angle);
                segment[i*(nsides+2) + nsides + 1].y1 = LAMBDA*sin(((double)i + 0.98)*angle);
                segment[i*(nsides+2) + nsides + 1].concave = 1;
            }
            
            cycle = 1;
            nsegments = (nsides+2)*NPOLY;
            for (i=0; i<nsegments; i++) 
            {
                segment[i].group = 0;
                segment[i].inactivate = 0;
            }
            
            break;
        }
        case (S_CENTRIFUGE_LEAKY):
        {
            angle = DPI/(double)NPOLY;
            nsides = 20;
            if (2*(nsides+2)*NPOLY > NMAXSEGMENTS)
            {
                printf("Error: NMAXSEGMENTS is too small\n");
                exit(1);
            }
            
            for (i=0; i<NPOLY; i++)
            {
                angle2 = (double)i*angle;
                x1 = LAMBDA*cos(angle2);
                y1 = LAMBDA*sin(angle2);
                x2 = 0.7*cos(angle2);
                y2 = 0.7*sin(angle2);
                add_rotated_angle_to_segments(x1, y1, x2, y2, MU, 0, segment, 0);
                
                for (j=0; j<nsides; j++) if (j!=nsides/2)
                {
                    x = (double)j/(double)(nsides);
                    r = 0.5 + sqrt(x*(1.0-x) + 0.04);
                    if (j < nsides/2) angle2 = (double)i*angle + angle*(double)j/((double)(nsides) - 0.15);
                    else angle2 = (double)i*angle + angle*(double)j/((double)(nsides) + 0.15);
                    x1 = r*cos(angle2);
                    y1 = r*sin(angle2);
                    
                    x = (double)(j+1)/(double)(nsides);
                    r = 0.5 + sqrt(x*(1.0-x) + 0.04);
                    if (j < nsides/2) angle2 = (double)i*angle + angle*(double)(j+1)/((double)(nsides) - 0.15);
                    else angle2 = (double)i*angle + angle*(double)(j+1)/((double)(nsides) + 0.15);
                    x2 = r*cos(angle2);
                    y2 = r*sin(angle2);
                    
                    add_rotated_angle_to_segments(x1, y1, x2, y2, 0.5*MU, 0, segment, 0);
                }
            }
            
            cycle = 0;
            
            for (i=0; i<nsegments; i++) 
            {
                segment[i].concave = 0;
                segment[i].group = 0;
                segment[i].inactivate = 0;
            }
            
            break;
        }
        case (S_CIRCLE_EXT):
        {
            angle = DPI/(double)NPOLY;
            
            for (i=0; i<NPOLY; i++)
            {
                segment[i].x1 = SEGMENTS_X0 + LAMBDA*cos(((double)i)*angle);
                segment[i].y1 = SEGMENTS_Y0 - LAMBDA*sin(((double)i)*angle);
                segment[i].concave = 1;
            }
            
            cycle = 1;
            nsegments = NPOLY;
            ngroups = 2;
            
            for (i=0; i<nsegments; i++) 
            {
                segment[i].group = 1;
                segment[i].inactivate = 0;
            }
            
            break;
        }
        case (S_ROCKET_NOZZLE):
        {
            /* ellipse */
            for (i=1; i<NPOLY-1; i++)
            {
                angle = -PI + (double)i*DPI/(double)NPOLY;
                x1 = 0.7*LAMBDA*(1.0 + cos(angle));
                y1 = 0.5*LAMBDA*sin(angle);
                angle = -PI + (double)(i+1)*DPI/(double)NPOLY;
                x2 = 0.7*LAMBDA*(1.0 + cos(angle));
                y2 = 0.5*LAMBDA*sin(angle);                
                add_rotated_angle_to_segments(x1, y1, x2, y2, 0.02, 0, segment, 0);
            }
             
            /* compute intersection point of nozzle and ellipse */
            angle = PI - DPI/(double)NPOLY;
            nozx = 0.7*LAMBDA*(1.0 + cos(angle));
            nozy = 0.5*LAMBDA*sin(angle);
            
            nsides = 10;
            dx = LAMBDA/(double)(nsides);
            
            /* nozzle */
            for (i=0; i<nsides; i++)
            {
                x1 = -LAMBDA + dx*(double)(i-1);
                y1 = nozzle_width(x1, nozy, NOZZLE_SHAPE);
                x2 = x1 + dx;
                y2 = nozzle_width(x2, nozy, NOZZLE_SHAPE);
                add_rotated_angle_to_segments(x1, y1, x2, y2, 0.02, 0, segment, 0);
            }
            add_rotated_angle_to_segments(x2, y2, nozx, nozy, 0.02, 0, segment, 0);
            for (i=0; i<nsides; i++)
            {
                x1 = -LAMBDA + dx*(double)(i-1);
                y1 = -nozzle_width(x1, nozy, NOZZLE_SHAPE);
                x2 = x1 + dx;
                y2 = -nozzle_width(x2, nozy, NOZZLE_SHAPE);
                add_rotated_angle_to_segments(x1, y1, x2, y2, 0.02, 0, segment, 0);
            }
            add_rotated_angle_to_segments(x2, y2, nozx, -nozy, 0.02, 0, segment, 0);
            
            /* closing segment */
            segment[nsegments].x1 = nozx;
            segment[nsegments].y1 = nozy;
            segment[nsegments].x2 = nozx;
            segment[nsegments].y2 = -nozy;
            nsegments++;
        
            cycle = 0;
            for (i=0; i<nsegments; i++) 
            {
                segment[i].group = 0;
                segment[i].inactivate = 0;
            }
            segment[nsegments-1].inactivate = 1;
            
            break;
        }
        case (S_ROCKET_NOZZLE_ROTATED):
        {
            add_rocket_to_segments(segment, 0.0, SEGMENTS_Y0, ROCKET_SHAPE, NOZZLE_SHAPE, 10, 1);
            /* segments from group 0 are immobile by convention */
            ngroups = 2;
            cycle = 0; 
//             for (i=0; i<nsegments; i++) segment[i].group = 1;
            break;
        }
        case (S_TWO_ROCKETS):
        {
            add_rocket_to_segments(segment, -SEGMENTS_X0, SEGMENTS_Y0, ROCKET_SHAPE, NOZZLE_SHAPE, 10, 1);
            add_rocket_to_segments(segment, SEGMENTS_X0, SEGMENTS_Y0, ROCKET_SHAPE_B, NOZZLE_SHAPE_B, 10, 2);
            ngroups = 3;
            cycle = 0;
            break;
        }
        case (S_TWO_CIRCLES_EXT):
        {
            angle = DPI/(double)NPOLY;
            
            for (i=0; i<NPOLY; i++)
            {
                segment[i].x1 = SEGMENTS_X0 + LAMBDA*cos(((double)i)*angle);
                segment[i].y1 = SEGMENTS_Y0 - LAMBDA*sin(((double)i)*angle);
                segment[i].x2 = SEGMENTS_X0 + LAMBDA*cos(((double)(i+1))*angle);
                segment[i].y2 = SEGMENTS_Y0 - LAMBDA*sin(((double)(i+1))*angle);
                segment[i].concave = 1;
                segment[i].group = 0;
                segment[i].inactivate = 0;
            }
            for (i=NPOLY; i<2*NPOLY; i++)
            {
                segment[i].x1 = -SEGMENTS_X0 + TWO_CIRCLES_RADIUS_RATIO*LAMBDA*cos(((double)i)*angle);
                segment[i].y1 = SEGMENTS_Y0 - TWO_CIRCLES_RADIUS_RATIO*LAMBDA*sin(((double)i)*angle);
                segment[i].x2 = -SEGMENTS_X0 + TWO_CIRCLES_RADIUS_RATIO*LAMBDA*cos(((double)(i+1))*angle);
                segment[i].y2 = SEGMENTS_Y0 - TWO_CIRCLES_RADIUS_RATIO*LAMBDA*sin(((double)(i+1))*angle);
                segment[i].concave = 1;
                segment[i].group = 1;
                segment[i].inactivate = 0;
            }
            
            cycle = 0;
            nsegments = 2*NPOLY;
            
            break;
        }
        case (S_DAM):
        {
            add_rectangle_to_segments(DAM_WIDTH, BCYMIN - 0.5, -DAM_WIDTH, LAMBDA, segment, 0);
            
            cycle = 0;
            for (i=0; i<nsegments; i++) 
            {
                segment[i].group = 0;
                segment[i].inactivate = 1;
            }
            break;
        }
        case (S_DAM_WITH_HOLE):
        {
            add_rectangle_to_segments(DAM_WIDTH, BCYMIN - 0.5, -DAM_WIDTH, BCYMIN + 0.1, segment, 0);
            add_rectangle_to_segments(DAM_WIDTH, BCYMIN + 0.3, -DAM_WIDTH, LAMBDA, segment, 0);
            add_rectangle_to_segments(DAM_WIDTH, BCYMIN + 0.1, -DAM_WIDTH, BCYMIN + 0.3, segment, 0);
            
            cycle = 0;
            concave = 0;    /* add_rectangle_to_segments already deals with concave corners */
            
            for (i=0; i<nsegments; i++) 
            {
                segment[i].group = 0;
                if (i > 7) segment[i].inactivate = 1;
                else segment[i].inactivate = 0;
            }
            break;
        }
        case (S_DAM_WITH_HOLE_AND_RAMP):
        {
            add_rectangle_to_segments(DAM_WIDTH, BCYMIN - 0.5, -DAM_WIDTH, BCYMIN + 0.2, segment, 0);
            add_rectangle_to_segments(DAM_WIDTH, BCYMIN + 0.3, -DAM_WIDTH, LAMBDA, segment, 0);
            add_rectangle_to_segments(DAM_WIDTH, BCYMIN + 0.2, -DAM_WIDTH, BCYMIN + 0.3, segment, 0);
            
            r = 1.0;
            for (i=0; i<10; i++)
            {
                angle = 0.1*PID*(double)i;
                dangle = 0.1*PID;
                x1 = XMAX - r + (r + MU)*cos(angle);
                y1 = YMIN + r - (r + MU)*sin(angle);
                x2 = XMAX - r + (r + MU)*cos(angle + dangle);
                y2 = YMIN + r - (r + MU)*sin(angle + dangle);
                add_rotated_angle_to_segments(x1, y1, x2, y2, MU, 0, segment, 0);
            }
            
            cycle = 0;
            concave = 0;    /* add_rectangle_to_segments already deals with concave corners */
            
            for (i=0; i<nsegments; i++) 
            {
                segment[i].group = 0;
                if ((i > 7)&&(i < 12)) segment[i].inactivate = 1;
                else segment[i].inactivate = 0;
            }
            break;
        }
        case (S_MAZE):
        {
            init_maze_segments(segment, 0);
            
            cycle = 0;
            for (i=0; i<nsegments; i++) 
            {
                segment[i].group = 0;
                segment[i].inactivate = 0;
            }
            
            break;
        }
        case (S_MAZE_DIAG):
        {
            init_maze_segments(segment, 1);
            
            cycle = 0;
            for (i=0; i<nsegments; i++) 
            {
                segment[i].group = 0;
                segment[i].inactivate = 0;
            }
            
            break;
        }
        case (S_EXT_RECTANGLE):
        {
            width = 0.1*LAMBDA;
            
            segment[0].x1 = SEGMENTS_X0 - LAMBDA;
            segment[0].y1 = SEGMENTS_Y0 - width;
            
            segment[1].x1 = SEGMENTS_X0 - LAMBDA;
            segment[1].y1 = SEGMENTS_Y0 + width;

            segment[2].x1 = SEGMENTS_X0 + LAMBDA;
            segment[2].y1 = SEGMENTS_Y0 + width;

            segment[3].x1 = SEGMENTS_X0 + LAMBDA;
            segment[3].y1 = SEGMENTS_Y0 - width;
            
            cycle = 1;
            nsegments = 4;
            ngroups = 2;
            
            for (i=0; i<nsegments; i++) 
            {
                segment[i].concave = 1;
                segment[i].group = 1;
                segment[i].inactivate = 0;
            }
            break;
        }
        case (S_DAM_BRICKS):
        {
            add_rectangle_to_segments(DAM_WIDTH, BCYMIN - 0.5, -DAM_WIDTH, BCYMIN, segment, 0);
            dy = 0.1*(LAMBDA - BCYMIN);
            
            for (i=1; i<11; i++)
            {
                y = BCYMIN + (double)i*dy;
                add_rectangle_to_segments(DAM_WIDTH, y-dy+MU, -DAM_WIDTH, y, segment, i);
            }
            
            ngroups = 11;
            cycle = 0;
            concave = 0;    /* add_rectangle_to_segments already deals with concave corners */
            
            for (i=0; i<nsegments; i++) 
            {
                segment[i].inactivate = 0;
            }
            break;
            
        }
        case (S_HLINE_HOLE):
        {
            x = 0.15;
            x1 = XMAX + 1.0;
            y1 = 0.0;
            width = 0.05;
            
            add_rectangle_to_segments(x1, y1 - width, x, y1, segment, 0);
            add_rectangle_to_segments(-x, y1 - width, -x1, y1, segment, 0);
                        
            /* closing segment */
            segment[nsegments].x1 = -x;
            segment[nsegments].y1 = y1;
            segment[nsegments].x2 = x;
            segment[nsegments].y2 = y1;
            nsegments++;
        
            cycle = 0;
            concave = 0;    /* add_rectangle_to_segments already deals with concave corners */
            
            for (i=0; i<nsegments; i++) 
            {
                segment[i].inactivate = 0;
            }
            segment[nsegments-1].inactivate = 1;
            
            break;
        }
        case (S_HLINE_HOLE_SPOKES):
        {
            x = 0.15;
            x1 = XMAX + 1.0;
            y1 = 0.0;
            width = 0.05;
            
            add_rectangle_to_segments(x1, y1 - width, x, y1, segment, 0);
            add_rectangle_to_segments(-x, y1 - width, -x1, y1, segment, 0);
            
            /* closing segment */
            segment[nsegments].x1 = -x;
            segment[nsegments].y1 = y1 - 0.5*width;
            segment[nsegments].x2 = x;
            segment[nsegments].y2 = y1 - 0.5*width;
            nsegments++;
                        
            /* spokes */
            for (i=-3; i<4; i+=2)
            {
                x = 0.5*(double)i;
                segment[nsegments].x1 = x - 0.1;
                segment[nsegments].y1 = BCYMIN - 0.1;
                segment[nsegments].x2 = x;
                segment[nsegments].y2 = YMIN + 0.3;
                nsegments++;
                segment[nsegments].x1 = x;
                segment[nsegments].y1 = YMIN + 0.3;
                segment[nsegments].x2 = x + 0.1;
                segment[nsegments].y2 = BCYMIN - 0.1;
                nsegments++;
            }
        
            
            cycle = 0;
            concave = 0;    /* add_rectangle_to_segments already deals with concave corners */
            
            for (i=0; i<nsegments; i++) 
            {
                segment[i].inactivate = 0;
            }
            segment[8].inactivate = 1;
            
            break;
        }
        case (S_EXT_CIRCLE_RECT):
        {
            width = 0.1*LAMBDA;
            
            add_rectangle_to_segments(SEGMENTS_X0 + LAMBDA, SEGMENTS_Y0 - width, SEGMENTS_X0 - LAMBDA, SEGMENTS_Y0 + width, segment, 1);
            
            add_circle_to_segments(-SEGMENTS_X0, SEGMENTS_Y0, 0.5*LAMBDA, NPOLY, APOLY, segment, 2);
                        
            cycle = 0;
            concave = 0;
            nsegments = NPOLY + 4;
            ngroups = 3;
            
            for (i=0; i<nsegments; i++) 
            {
                segment[i].concave = 1;
                segment[i].inactivate = 0;
            }
            
            break;
        }
        case (S_BIN_OPENING):
        {
            add_rectangle_to_segments(LAMBDA, 1.0 - LAMBDA, LAMBDA - MU, YMAX + 1.0, segment, 0);
            add_rectangle_to_segments(-LAMBDA + MU, 1.0 - LAMBDA, -LAMBDA, YMAX + 1.0, segment, 0);
            add_rectangle_to_segments(LAMBDA, 1.0 - LAMBDA, -LAMBDA, 1.0 - LAMBDA + MU, segment, 0);
            
            cycle = 0;
            concave = 0;    /* add_rectangle_to_segments already deals with concave corners */
            
            for (i=0; i<nsegments; i++) 
            {
                segment[i].group = 0;
                segment[i].inactivate = 1;
            }
            break;
        }
        case (S_POLYGON_EXT):
        {
            add_circle_to_segments(0.0, 0.0, LAMBDA, NPOLY, APOLY, segment, 0);
                        
            cycle = 0;
            concave = 0;
            nsegments = NPOLY;
            ngroups = 1;
            
            for (i=0; i<nsegments; i++) 
            {
                segment[i].concave = 1;
                segment[i].inactivate = 0;
            }
            
            break;
        }
        case (S_WEDGE_EXT):
        {
            angle = AWEDGE*PID;
            
            segment[0].x1 = LAMBDA;
            segment[0].y1 = 0.0;
            
            segment[1].x1 = -LAMBDA*cos(angle);
            segment[1].y1 = -LAMBDA*sin(angle);

            segment[2].x1 = 0.0;
            segment[2].y1 = 0.0;

            segment[3].x1 = -LAMBDA*cos(angle);
            segment[3].y1 = LAMBDA*sin(angle);
            
            cycle = 1;
            nsegments = 4;
            ngroups = 2;
            
            ca = cos(APOLY*PID);
            sa = sin(APOLY*PID);
            
            for (i=0; i<nsegments; i++) 
            {
                x = segment[i].x1;
                y = segment[i].y1;
                segment[i].x1 = x*ca + y*sa;
                segment[i].y1 = -x*sa + y*ca;
                segment[i].concave = 1;
                segment[i].group = 1;
                segment[i].inactivate = 0;
            }
            break;
        }
        case (S_MIXER):
        {
            for (i=0; i<NPOLY; i++)
            {
                angle = (double)i*DPI/(double)NPOLY;
                add_rotated_angle_to_segments(0.0, 0.0, LAMBDA*cos(angle), LAMBDA*sin(angle), 0.05, 1, segment, 0);
            }
            
            cycle = 0;
            concave = 1;    /* add_rectangle_to_segments already deals with concave corners */
            
            for (i=0; i<nsegments; i++) 
            {
                segment[i].group = 0;
                segment[i].inactivate = 0;
            }
            
            break;
        }
        case (S_AIRFOIL):
        {
            dangle = DPI/(double)NPOLY;
            angle = APOLY*PID;
            ca = cos(angle);
            sa = sin(angle);
            for (i=0; i<NPOLY; i++)
            {
                angle = (double)i*dangle;
                x = LAMBDA*cos(angle);
                y1 = -0.2*LAMBDA*sin(angle);
                y1 -= 0.5*x*x;
                segment[i].x1 = x*ca + y1*sa;
                segment[i].y1 = -x*sa + y1*ca;
            }
            
            cycle = 1;
            concave = 1;
            nsegments = NPOLY;
            ngroups = 1;
            
            for (i=0; i<nsegments; i++) 
            {
                segment[i].concave = 1;
                segment[i].inactivate = 0;
            }
            
            break;
        }
        default: 
        {
            printf("Function init_segment_config not defined for this pattern \n");
        }
    }
    
    if (cycle) for (i=0; i<nsegments; i++)
    {
        segment[i].x2 = segment[(i+1)%(nsegments)].x1;
        segment[i].y2 = segment[(i+1)%(nsegments)].y1;
    }
    else if (SEGMENT_PATTERN != S_TWO_CIRCLES_EXT) for (i=0; i<nsegments; i++) if (segment[i].cycle)
    {
        segment[i].x2 = segment[(i+1)%(nsegments)].x1;
        segment[i].y2 = segment[(i+1)%(nsegments)].y1;
    }
    
    /* add one segment for S_POLY_ELLIPSE configuration */
    if (SEGMENT_PATTERN == S_POLY_ELLIPSE)
    {
        segment[nsegments].x1 = -cos(((double)i + 0.5)*angle);
        segment[nsegments].y1 = LAMBDA*sin(((double)i + 0.5)*angle);
        segment[nsegments].x2 = -cos(((double)i + 0.5)*angle);
        segment[nsegments].y2 = -LAMBDA*sin(((double)i + 0.5)*angle);
        segment[nsegments].inactivate = 1;
        nsegments++;
    }
    
    /* activate all segments */
    for (i=0; i<nsegments; i++) segment[i].active = 1;
    
    /* inactivate some segments of leaky centrifuge */
//     if (SEGMENT_PATTERN == S_CENTRIFUGE_LEAKY)
//         for (i=0; i<NPOLY; i++) segment[(nsides+2)*i + nsides/2].active = 0;
    
    /* compute parameters for slope and normal of segments */
    for (i=0; i<nsegments; i++)
    {
        a = segment[i].y1 - segment[i].y2;
        b = segment[i].x2 - segment[i].x1;
        length = module2(a, b);
        segment[i].nx = a/length;
        segment[i].ny = b/length;
        segment[i].c = segment[i].nx*segment[i].x1 + segment[i].ny*segment[i].y1;
        segment[i].length = length;
        segment[i].fx = 0.0;
        segment[i].fy = 0.0;
        segment[i].torque = 0.0;
        
        segment[i].xc = 0.5*(segment[i].x1 + segment[i].x2);
        segment[i].yc = 0.5*(segment[i].y1 + segment[i].y2);
    }
    
    /* deal with concave corners */
    if (concave) for (i=0; i<nsegments; i++) if (segment[i].concave)
        {
            iminus = i-1;  
            iplus = i+1;  
            if (SEGMENT_PATTERN == S_TWO_CIRCLES_EXT)
            {
                if (iminus == -1) iminus += nsegments/2;
                else if (iminus == nsegments/2 - 1) iminus += nsegments/2;
                if (iplus == nsegments/2) iplus = 0;
                else if (iplus == nsegments) iplus = nsegments/2;
            }
            else
            {
                if (iminus < 0) iminus = nsegments - 1;
                if (iplus > nsegments - 1) iplus = 0;
            }
            angle = argument(segment[iplus].x1 - segment[i].x1, segment[iplus].y1 - segment[i].y1) + PID;
            angle2 = argument(segment[i].x1 - segment[iminus].x1, segment[i].y1 - segment[iminus].y1) + PID;
            if (angle2 < angle) angle2 += DPI;
            segment[i].angle1 = angle;
            segment[i].angle2 = angle2;
            
            printf("i = %i, iplus = %i, iminus = %i, angle1 = %.0f, angle2 = %.0f\n", i, iplus, iminus, angle*360.0/DPI, angle2*360.0/DPI);
        }
    
    /* make copy of initial values in case of rotation/translation */
    if ((ROTATE_BOUNDARY)||(MOVE_BOUNDARY)||(MOVE_SEGMENT_GROUPS)) for (i=0; i<nsegments; i++) 
    {
        segment[i].x01 = segment[i].x1;
        segment[i].x02 = segment[i].x2;
        segment[i].y01 = segment[i].y1;
        segment[i].y02 = segment[i].y2;
        segment[i].nx0 = segment[i].nx;
        segment[i].ny0 = segment[i].ny;
        segment[i].angle01 = segment[i].angle1;
        segment[i].angle02 = segment[i].angle2;
    }
    
//     for (i=0; i<nsegments; i++) 
//     {
//         printf("Segment %i: (x1, y1) = (%.3lg,%.3lg), (x2, y2) = (%.3lg,%.3lg)\n (nx, ny) = (%.3lg,%.3lg), c = %.3lg, length = %.3lg\n", i, segment[i].x1, segment[i].y1, segment[i].x2, segment[i].y2, segment[i].nx, segment[i].ny, segment[i].c, segment[i].length);
//         if (segment[i].concave) printf("Concave with angles %.3lg Pi, %.3lg Pi\n", segment[i].angle1/PI, segment[i].angle2/PI);
//     }
//     sleep(4);
}

int in_rocket(double x, double y, int rocket_shape)
/* returns 1 if (x,y) is in rocket chamber, with translated coordinates */
{
    double l, y1, a, b, c;
    
    switch (rocket_shape) {
        case (RCK_DISC) :
        {
            l = 0.7*LAMBDA;
            y1 = y - YMIN - 1.7*LAMBDA;
            return ((x*x + y1*y1)/(l*l) + MU*MU < 0.875);
//             return ((x*x + y1*y1)/(l*l) + MU*MU < 0.925*LAMBDA*LAMBDA);
        }
        case (RCK_RECT) :
        {
            a = 0.5*LAMBDA;
            b = 0.49*PI*LAMBDA;
            y1 = y - YMIN - LAMBDA;
            return ((vabs(x) < 0.95*a)&&(y1 > 0.05)&&(y1 < b - 0.05));
        }
        case (RCK_RECT_HAT) :
        {
//             printf("(%.2lg,%.2lg) in rocket?\n", x, y);
            a = 0.5*LAMBDA;
            b = (0.49*PI-0.25)*LAMBDA;
            y1 = y - YMIN - LAMBDA;
            if (vabs(x) > 0.95*a) return(0);
            if (y1 < 0.05) return(0);
            if (y1 < b - 0.05) return(1);
            return(y1 < a + b - 0.05 - vabs(x));
//             return(1);
        }
        case (RCK_RECT_BAR) :
        {
            a = 0.5*LAMBDA;
            b = (0.49*PI-0.25)*LAMBDA;
            c = 0.5*a;
            y1 = y - YMIN - LAMBDA;
            if (vabs(x) > 0.95*a) return(0);
            if (vabs(x) < 0.1*a) return(0);
            if (y1 < 0.05 + c) return(0);
            if (y1 < b - 0.05 + c) return(1);
            return(y1 < a + b - 0.05 + c - vabs(x));
        }
    }
}

int in_segment_region(double x, double y, t_segment segment[NMAXSEGMENTS])
/* returns 1 if (x,y) is inside region delimited by obstacle segments */
{
    int i;
    double angle, dx, height, width, theta, lx, ly, x1, y1, x2, y2, padding, ca, sa, r;
    
    if (x >= BCXMAX) return(0);
    if (x <= BCXMIN) return(0);
    if (y >= BCYMAX) return(0);
    if (y <= BCYMIN) return(0);
            
    switch (SEGMENT_PATTERN) {
        case (S_CUP):
        {
            angle = APOLY*PID;
            dx = (BCYMAX - BCYMIN)/tan(angle);
            
            if (y < BCYMAX - (BCYMAX - BCYMIN)*(x - BCXMIN)/dx) return(0);
            if (y < BCYMAX - (BCYMAX - BCYMIN)*(BCXMAX - x)/dx) return(0);
        }
        case (S_HOURGLASS):
        {
            angle = APOLY*PID;
            width = 2.5*MU;
            height = 2.5*MU;
            
            x = vabs(x);
            y = vabs(y);
            
            if ((x >= width)&&(x - width >= (y - height)*(BCXMAX - width)/(BCYMAX - height))) return(0);
            return(1);
        }
        case (S_PENTA):
        {
            height = 0.5*(BCYMAX - BCYMIN);
            width = height/sqrt(3.0);
            
            if (y < BCYMIN + height*(1.0 - (x - BCXMIN)/width)) return(0);
            if (y > BCYMAX - height*(1.0 - (x - BCXMIN)/width)) return(0);
            
            return(1);
        }
        case (S_CENTRIFUGE):
        {
            angle = argument(x,y);
            theta = DPI/(double)NPOLY;
            while (angle > theta) angle -= theta;
            while (angle < 0.0) angle += theta;
            if (angle < 0.1) return(0);
            if (angle > 0.9) return(0);
            return(1);
        }
        case (S_POLY_ELLIPSE):
        {
            if (x*x + y*y/(LAMBDA*LAMBDA) < 0.95) return(1);
            else return(0);
        }
        case (S_CENTRIFUGE_RND):
        {
            if (module2(x,y) > 0.75) return(0);
            angle = argument(x,y);
            theta = DPI/(double)NPOLY;
            while (angle > theta) angle -= theta;
            while (angle < 0.0) angle += theta;
            if (angle < 0.1) return(0);
            if (angle > 0.9) return(0);
            return(1);
        }
        case (S_CENTRIFUGE_LEAKY):
        {
            if (module2(x,y) > 0.75) return(0);
            angle = argument(x,y);
            theta = DPI/(double)NPOLY;
            while (angle > theta) angle -= theta;
            while (angle < 0.0) angle += theta;
            if (angle < 0.1) return(0);
            if (angle > 0.9) return(0);
            return(1);
        }
        case (S_CIRCLE_EXT):
        {
            if (module2(x - SEGMENTS_X0, y - SEGMENTS_Y0) > LAMBDA + 2.0*MU) return(1);
            else return(0);
        }
        case (S_TWO_CIRCLES_EXT):
        {
            if ((module2(x - SEGMENTS_X0, y - SEGMENTS_Y0) > LAMBDA + 2.0*MU)&&(module2(x + SEGMENTS_X0, y - SEGMENTS_Y0) > TWO_CIRCLES_RADIUS_RATIO*LAMBDA + 2.0*MU)) return(1);
            else return(0);
        }
        case (S_ROCKET_NOZZLE):
        {
            if (x < 0.0) return(0);
            else if (x > 1.4*LAMBDA) return(0);
            else 
            {
                lx = 0.7*LAMBDA;
                ly = 0.7*LAMBDA;
                x -= lx;
                if (x*x/(lx*lx) + y*y/(ly*ly) < 0.95) return(1);
            }
            return(0);
        }
        case (S_ROCKET_NOZZLE_ROTATED):
        {
            y1 = y - ysegments[0];
            if (y1 < YMIN + LAMBDA) return(0);
            else if (y1 > YMIN + 2.4*LAMBDA) return(0);
            else if (in_rocket(x - xsegments[0], y1, ROCKET_SHAPE)) return(1);
//             {
//                 ly = 0.7*LAMBDA;
//                 lx = 0.7*LAMBDA;
//                 x1 = x - xsegments[0];
//                 y1 -= YMIN + 1.7*LAMBDA;
//                 if (x1*x1/(lx*lx) + y1*y1/(ly*ly) + MU*MU < 0.925) return(1);
//             }
            return(0);
        }
        case (S_TWO_ROCKETS):
        {
            y1 = y - ysegments[0];
            y2 = y - ysegments[1];
            if ((y1 < YMIN + LAMBDA)&&(y2 < YMIN + LAMBDA)) return(0);
            else if ((y1 > YMIN + 3.5*LAMBDA)&&(y2 > YMIN + 3.5*LAMBDA)) return(0);
            else 
            {
                if (in_rocket(x - xsegments[0], y1, ROCKET_SHAPE_B)) return(1);
                if (in_rocket(x - xsegments[1], y2, ROCKET_SHAPE)) return(1);
            }
            return(0);
        }
        case (S_DAM):
        {
            if (vabs(x) > DAM_WIDTH) return(1);
            else if (y > LAMBDA) return(1);
            else return(0);
        }
        case (S_EXT_RECTANGLE):
        {
            padding = 0.1;
            if (vabs(x) > LAMBDA + padding) return(1);
            else if (vabs(y) > 0.1*LAMBDA + padding) return(1);
            else return(0);
        }
        case (S_EXT_CIRCLE_RECT):
        {
            padding = 0.1;
            if ((vabs(x - SEGMENTS_X0) < LAMBDA + padding)&&(vabs(y - SEGMENTS_Y0) < 0.1*LAMBDA + padding)) return(0);
            else if (module2(x + SEGMENTS_X0, y - SEGMENTS_Y0) < 0.5*LAMBDA + padding) return(0);
            else return(1);
        }
        case (S_BIN_OPENING):
        {
            padding = 3.0*MU;
            if (y < 0.0) return(1);
            if ((y > 1.0 - LAMBDA + padding)&&(vabs(x) < LAMBDA - padding)) return(1);
            return(0);
        }
        case (S_POLYGON_EXT):
        {
            padding = 3.0*MU;
            if (in_polygon(x, y, LAMBDA + padding, NPOLY, APOLY)) return(0);
            else return(1);
        }
        case (S_WEDGE_EXT):
        {
            padding = 3.0*MU;
            angle = AWEDGE*PID;
            ca = cos(APOLY*PID);
            sa = sin(APOLY*PID);
            x1 = x*ca - y*sa;
            y1 = x*sa + y*ca;
            if (vabs(y1) - padding > (LAMBDA-x1)*sin(angle)/(1.0+cos(angle))) return(1);
            if (vabs(y1) + padding < -x1*tan(angle)) return(1);
            return(0);
        }
        case (S_MIXER):
        {
            padding = 1.5*MU;
            r = module2(x,y);
            if (r > LAMBDA + padding) return(1);
            if (r < 2.0*padding) return(0);
            for (i=0; i<NPOLY; i++)
            {
                angle = (double)i*DPI/(double)NPOLY;
                ca = cos(angle);
                sa = sin(angle);
                x1 = x*ca - y*sa;
                y1 = x*sa + y*ca;
                if ((x1 > 0.0)&&(vabs(y1) < 0.025 + padding)) return(0); 
            }
            return(1);
        }
        case (S_AIRFOIL):
        {
            if (vabs(x) > LAMBDA + 4.0*MU) return(1);
            padding = 0.35;
            angle = APOLY*PID;
            ca = cos(angle);
            sa = sin(angle);
            x1 = x*ca - y*sa;
            y1 = x*sa + y*ca;
            y1 += 0.5*x1*x1;
            return(x1*x1 + 25.0*y1*y1 > LAMBDA*LAMBDA*(1.0 + padding));
        }
        case (S_MAZE):
        {
            for (i=0; i<nsegments; i++)
            {
                if (distance_to_segment(x, y, segment[i].x1, segment[i].y1, segment[i].x2, segment[i].y2) < 5.0*MAZE_WIDTH) return(0);
            }
            return(1);
        }
        case (S_MAZE_DIAG):
        {
            for (i=0; i<nsegments; i++)
            {
                if (distance_to_segment(x, y, segment[i].x1, segment[i].y1, segment[i].x2, segment[i].y2) < 5.0*MAZE_WIDTH) return(0);
            }
            return(1);
        }
        default: return(1);
    }
}


void rotate_segments(t_segment segment[NMAXSEGMENTS], double angle)
/* rotates the repelling segments by given angle */
{
    int i;
    double ca, sa;
    
    ca = cos(angle);
    sa = sin(angle);
    
    for (i=0; i<nsegments; i++) 
    {
        segment[i].x1 = ca*segment[i].x01 - sa*segment[i].y01;
        segment[i].y1 = sa*segment[i].x01 + ca*segment[i].y01;
        segment[i].x2 = ca*segment[i].x02 - sa*segment[i].y02;
        segment[i].y2 = sa*segment[i].x02 + ca*segment[i].y02;
        segment[i].nx = ca*segment[i].nx0 - sa*segment[i].ny0;
        segment[i].ny = sa*segment[i].nx0 + ca*segment[i].ny0;
        
        if (segment[i].concave)
        {
            segment[i].angle1 = segment[i].angle01 + angle;
            segment[i].angle2 = segment[i].angle02 + angle;
            while (segment[i].angle1 > DPI)
            {
                segment[i].angle1 -= DPI;
                segment[i].angle2 -= DPI;
            }
        }
    }
    
}
 
void translate_segments(t_segment segment[NMAXSEGMENTS], double deltax[2], double deltay[2])
/* rotates the repelling segments by given angle */
{
    int i, group;
    
    for (i=0; i<nsegments; i++) 
    {
        group = segment[i].group;
        if (group == 0)
        {
            segment[i].x1 = segment[i].x01 + deltax[group] - SEGMENTS_X0;
            segment[i].x2 = segment[i].x02 + deltax[group] - SEGMENTS_X0;
        }
        else
        {
            segment[i].x1 = segment[i].x01 + deltax[group] + SEGMENTS_X0;
            segment[i].x2 = segment[i].x02 + deltax[group] + SEGMENTS_X0;
        }
            
        segment[i].y1 = segment[i].y01 + deltay[group] - SEGMENTS_Y0;
        segment[i].y2 = segment[i].y02 + deltay[group] - SEGMENTS_Y0;
        segment[i].c = segment[i].nx*segment[i].x1 + segment[i].ny*segment[i].y1;
    }
}

void translate_one_segment(t_segment segment[NMAXSEGMENTS], int i, double deltax, double deltay)
/* translates the repelling segment by given vector */
{
    segment[i].x1 += deltax;
    segment[i].x2 += deltax;
        
    segment[i].y1 += deltay;
    segment[i].y2 += deltay;
    segment[i].c = segment[i].nx*segment[i].x1 + segment[i].ny*segment[i].y1;
    
    segment[i].xc += deltax;
    segment[i].yc += deltay;
    
}

void rotate_one_segment(t_segment segment[NMAXSEGMENTS], int i, double dalpha, double xc, double yc)
/* rotates the repelling segment by given angle around (xc, yc) */
{
    double ca, sa, x, y, nx, ny;
    
    ca = cos(dalpha);
    sa = sin(dalpha);
    
    x = segment[i].x1 - xc;
    y = segment[i].y1 - yc;
    
    segment[i].x1 = xc + x*ca - y*sa;
    segment[i].y1 = yc + x*sa + y*ca;
    
    x = segment[i].x2 - xc;
    y = segment[i].y2 - yc;
    
    segment[i].x2 = xc + x*ca - y*sa;
    segment[i].y2 = yc + x*sa + y*ca;
    
    segment[i].xc = 0.5*(segment[i].x1 + segment[i].x2);
    segment[i].yc = 0.5*(segment[i].y1 + segment[i].y2);
    
    nx = segment[i].nx;
    ny = segment[i].ny;
    
    segment[i].nx = ca*nx - sa*ny;
    segment[i].ny = sa*nx + ca*ny;
    
    segment[i].c = segment[i].nx*segment[i].x1 + segment[i].ny*segment[i].y1;
    
    if (segment[i].concave)
    {
        segment[i].angle1 += dalpha;
        segment[i].angle2 += dalpha;
        while (segment[i].angle1 > DPI)
        {
            segment[i].angle1 -= DPI;
            segment[i].angle2 -= DPI;
        }
        while (segment[i].angle2 < 0.0)
        {
            segment[i].angle1 += DPI;
            segment[i].angle2 += DPI;
        }
    }
}

/* Computation of interaction force */

double lennard_jones_force(double r, t_particle particle)
{
    int i;
    double rmin = 0.01, rplus, ratio = 1.0;
    
    if (r > REPEL_RADIUS*particle.radius) return(0.0);
    else
    {
        if (r > rmin) rplus = r;
        else rplus = rmin;
    
//         ratio = ipow(particle.eq_dist*particle.radius/rplus, 6);
        ratio = ipow(particle.eq_dist*particle.radius/rplus, 3);
        ratio = ratio*ratio;
        
        return((ratio - 2.0*ratio*ratio)/rplus);
    }
}

double harmonic_force(double r, t_particle particle)
{
    int i;
    double rmin = 0.01, rplus, ratio = 1.0;
    
    if (r > REPEL_RADIUS*particle.radius) return(0.0);
    else
    {
//         if (r > rmin) rplus = r;
//         else rplus = rmin;
//         ratio = rplus/particle.eq_dist*particle.radius;
    
        ratio = r/particle.eq_dist*particle.radius;
    
        if (ratio < 1.0) return(ratio - 1.0);
        else return(0.0);
    }
}

double coulomb_force(double r, t_particle particle)
{
    int i;
    double rmin = 0.01, rplus, ratio = 1.0;
    
    if (r > REPEL_RADIUS*particle.radius) return(0.0);
    else
    {
//         if (r > rmin) rplus = r;
//         else rplus = rmin;
//         ratio = rplus/particle.eq_dist*particle.radius;
    
        ratio = r/particle.eq_dist*particle.radius;
        
        return(-1.0/(ratio*ratio + rmin*rmin));    
//         if (ratio < 1.0) return(ratio - 1.0);
//         else return(0.0);
    }
}

void aniso_lj_force(double r, double ca, double sa, double ca_rel, double sa_rel, double force[2], t_particle particle)
{
    int i;
    double rmin = 0.01, rplus, ratio = 1.0, c2, s2, c4, s4, a, aprime, f1, f2;
    
    if (r > REPEL_RADIUS*particle.radius)
    {
        force[0] = 0.0;
        force[1] = 0.0;
    }
    else
    {
        if (r > rmin) rplus = r;
        else rplus = rmin;
    
//         ratio = ipow(particle.eq_dist*particle.radius/rplus, 6);
        ratio = ipow(particle.eq_dist*particle.radius/rplus, 3);
        ratio = ratio*ratio;

        /* cos(2phi) and sin(2phi) */
        c2 = ca_rel*ca_rel - sa_rel*sa_rel;
        s2 = 2.0*ca_rel*sa_rel;

        /* cos(4phi) and sin(4phi) */
        c4 = c2*c2 - s2*s2;
        s4 = 2.0*c2*s2;
        
        a = 0.5*(9.0 - 7.0*c4);
        aprime = 14.0*s4;
        
        f1 = ratio*(a - ratio)/rplus;
        f2 = ratio*aprime/rplus;
        
        force[0] = f1*ca - f2*sa;
        force[1] = f1*sa + f2*ca;
    }
}

void penta_lj_force(double r, double ca, double sa, double ca_rel, double sa_rel, double force[2], t_particle particle)
{
    int i;
    double rmin = 0.01, rplus, ratio = 1.0, c2, s2, c4, s4, c5, s5, a, aprime, f1, f2;
    static double a0, b0;
    static int first = 1;
    
    if (first)
    {
        a0 = cos(0.1*PI) + 0.5;
        b0 = a0 - 1.0;
        first = 0;
    }
    
    if (r > REPEL_RADIUS*particle.radius)
    {
        force[0] = 0.0;
        force[1] = 0.0;
    }
    else
    {
        if (r > rmin) rplus = r;
        else rplus = rmin;
    
//         ratio = ipow(particle.eq_dist*particle.radius/rplus, 6);
        ratio = ipow(particle.eq_dist*particle.radius/rplus, 3);
        ratio = ratio*ratio;

    
        /* cos(2phi) and sin(2phi) */
        c2 = ca_rel*ca_rel - sa_rel*sa_rel;
        s2 = 2.0*ca_rel*sa_rel;

        /* cos(4phi) and sin(4phi) */
        c4 = c2*c2 - s2*s2;
        s4 = 2.0*c2*s2;
        
        /* cos(5phi) and sin(5phi) */
        c5 = ca_rel*c4 - sa_rel*s4;
        s5 = sa_rel*c4 + ca_rel*s4;
        
        a = a0 - b0*c5;
        aprime = 5.0*b0*s5;
        
        f1 = ratio*(a - ratio)/rplus;
        f2 = ratio*aprime/rplus;
        
        force[0] = f1*ca - f2*sa;
        force[1] = f1*sa + f2*ca;
    }
}

double old_golden_ratio_force(double r, t_particle particle)
/* potential with two minima at distances whose ratio is the golden ratio Phi */
/* old version that does not work very well */
{
    int i;
    double x, y, z, rplus, ratio = 1.0, phi, a, phi3;
    static int first = 1;
    static double rmin, b, c, d;
    
    if (first) 
    {
        rmin = 0.5*particle.radius;
        phi = 0.5*(1.0 + sqrt(5.0));
        phi3 = 1.0/(phi*phi*phi);
        a = 0.66;
        b = 1.0 + phi3 + a;
        d = phi3*a;
        c = phi3 + a + d;
//         b = 7.04;
//         c = 13.66;
//         d = 6.7;
        first = 0;
        printf("a = %.4lg, b = %.4lg, c = %.4lg, d = %.4lg\n", a, b, c, d);
    }
    
    if (r > REPEL_RADIUS*particle.radius) return(0.0);
    else
    {
        if (r > rmin) rplus = r;
        else rplus = rmin;
    
        x = particle.eq_dist*particle.radius/rplus;
        y = x*x*x;
        z = d - c*y + b*y*y - y*y*y;
        return(x*z/rplus);
    }
}

double golden_ratio_force(double r, t_particle particle)
/* potential with two minima at distances whose ratio is the golden ratio Phi */
/* piecewise polynomial/LJ version */
{
    int i;
    double x, rplus, xm6, y1;
    static int first = 1;
    static double rmin, phi, a, h1, h2, phi6;
    
    if (first) 
    {
        rmin = 0.5*particle.radius;
        phi = 0.5*(1.0 + sqrt(5.0));
        a = 1.2;
        
        h1 = 1.0;       /* inner potential well depth */
        h2 = 10.0;      /* outer potential well depth */
        phi6 = ipow(phi, 6);
        
        first = 0;
    }
    
    if (r > REPEL_RADIUS*particle.radius) return(0.0);
    else
    {
        if (r > rmin) rplus = r;
        else rplus = rmin;
        
        x = rplus/(particle.eq_dist*particle.radius);
//         xm6 = 1.0/ipow(x, 6);
        xm6 = 1.0/ipow(x, 3);
        xm6 = xm6*xm6;
    
        if (x <= 1.0) return(12.0*h1*xm6*(1.0 - xm6)/x);
        else if (x <= a)
        {
            y1 = ipow(a - 1.0, 3);
            return(6.0*h1*(x - 1.0)*(a - x)/y1);
        }
        else if (x <= phi)
        {
            y1 = ipow(phi - a, 3);
            return(6.0*h2*(x - a)*(x - phi)/y1);
        }
        else return(12.0*h2*phi6*(1.0 - phi6*xm6)*xm6/x);
    }
}

void dipole_lj_force(double r, double ca, double sa, double ca_rel, double sa_rel, double force[2], t_particle particle)
{
    int i;
    double rmin = 0.01, rplus, ratio = 1.0, a, aprime, f1, f2;
    
    if (r > REPEL_RADIUS*MU)
    {
        force[0] = 0.0;
        force[1] = 0.0;
    }
    else
    {
        if (r > rmin) rplus = r;
        else rplus = rmin;
    
//         ratio = ipow(particle.eq_dist*particle.radius/rplus, 6);
        ratio = ipow(particle.eq_dist*particle.radius/rplus, 3);
        ratio = ratio*ratio;
            
        a = 1.0 + 0.25*ca_rel;
        aprime = -0.25*sa_rel;
        
        f1 = ratio*(a - ratio)/rplus;
        f2 = ratio*aprime/rplus;
        
        force[0] = f1*ca - f2*sa;
        force[1] = f1*sa + f2*ca;
    }
}

void quadrupole_lj_force(double r, double ca, double sa, double ca_rel, double sa_rel, double force[2], t_particle particle)
{
    int i;
    double rmin = 0.01, rplus, ratio = 1.0, a, aprime, f1, f2, ca2, sa2, x, y, dplus, dminus;
    static int first = 1;
    static double a0, b0, aplus, aminus;

    if (first)
    {
        dplus = cos(0.2*PI)*cos(0.1*PI);
//         dminus = 0.8*dplus;
        dminus = QUADRUPOLE_RATIO*dplus;
        aplus = ipow(1.0/dplus, 6);
        aminus = ipow(1.0/dminus, 6);
//         aminus = ipow(cos(0.2*PI)*(0.25 + 0.5*sin(0.1*PI)), 6);
        a0 = 0.5*(aplus + aminus);
        b0 = 0.5*(aplus - aminus);
        first = 0;
    }
    
    if (r > REPEL_RADIUS*particle.radius)
    {
        force[0] = 0.0;
        force[1] = 0.0;
    }
    else
    {
        if (r > rmin) rplus = r;
        else rplus = rmin;
    
//         ratio = ipow(particle.eq_dist*particle.radius/rplus, 6);
        ratio = ipow(particle.eq_dist*particle.radius/rplus, 3);
        ratio = ratio*ratio;
        
        /* cos(2*phi) and sin(2*phi) */
        ca2 = ca_rel*ca_rel - sa_rel*sa_rel;
        sa2 = 2.0*ca_rel*sa_rel;
            
        a = a0 + b0*ca2;
//         if (a == 0.0) a = 1.0e-10; 
        aprime = -2.0*b0*sa2;
                
        f1 = ratio*(a - ratio)/rplus;
        f2 = ratio*aprime/rplus;
        
        force[0] = f1*ca - f2*sa;
        force[1] = f1*sa + f2*ca;
    }
}


void quadrupole_lj_force2(double r, double ca, double sa, double ca_rel, double sa_rel, double force[2], t_particle particle)
{
    int i;
    double rmin = 0.01, rplus, ratio = 1.0, a, aprime, f1, f2, ca2, sa2, x, y, eqdist;
    static int first = 1;
    static double aplus, aminus, a0, b0;

    if (first)
    {
        aplus = ipow(cos(0.2*PI)*cos(0.1*PI), 6);
        aminus = 0.1*aplus;
//         aminus = 0.0;
//         aminus = -2.0*ipow(cos(0.2*PI)*(0.5*sin(0.1*PI)), 6);
//         aminus = ipow(cos(0.2*PI)*(0.25 + 0.5*sin(0.1*PI)), 6);
        a0 = 0.5*(aplus + aminus);
        b0 = 0.5*(aplus - aminus);
        first = 0;
    }
    
    if (r > REPEL_RADIUS*particle.radius)
    {
        force[0] = 0.0;
        force[1] = 0.0;
    }
    else
    {        
        if (r > rmin) rplus = r;
        else rplus = rmin;
    
        /* correct distance */
//         ratio = ipow(particle.eq_dist*particle.radius/rplus, 6);
        ratio = ipow(particle.eq_dist*particle.radius/rplus, 3);
        ratio = ratio*ratio;
        
        /* cos(2*phi) and sin(2*phi) */
        
        ca2 = ca_rel*ca_rel - sa_rel*sa_rel;
        sa2 = 2.0*ca_rel*sa_rel;
            
        a = a0 + b0*ca2;
        if (a == 0.0) a = 1.0e-10; 
        aprime = -2.0*b0*sa2;
        
//         f1 = ratio*(a - ratio)/rplus;
//         f2 = ratio*aprime/rplus;
        
        f1 = ratio*(aplus - ratio)/(rplus);
        f2 = ratio*(aminus - ratio)/(rplus);
        
//         force[0] = f1*ca_rel - f2*sa_rel;
//         force[1] = f1*sa_rel + f2*ca_rel;

        force[0] = f1*ca - f2*sa;
        force[1] = f1*sa + f2*ca;
    }
}

double water_torque(double r, double ca, double sa, double ca_rel, double sa_rel, double ck_rel, double sk_rel)
/* compute torque of water molecule #k on water molecule #j (for interaction I_LJ_WATER) - OLD VERSION */
{
    double c1p, c1m, c2p, c2m, s2p, s2m, s21, s21p, s21m, c21, c21p, c21m, torque;
    double r2, rd, rd2, rr[3][3];
    static double cw = -0.5, sw = 0.866025404, delta = 1.5*MU, d2 = 2.25*MU*MU;
    int i, j;
    
    c1p = ck_rel*cw - sk_rel*sw;
    c1m = ck_rel*cw + sk_rel*sw;
    c2p = ca_rel*cw - sa_rel*sw;
    c2m = ca_rel*cw + sa_rel*sw;
    s2p = sa_rel*cw + ca_rel*sw;
    s2m = sa_rel*cw - ca_rel*sw;
    
    s21 = sa_rel*ck_rel - ca_rel*sk_rel;
    c21 = ca_rel*ck_rel + sa_rel*sk_rel;
    
    s21p = s21*cw - c21*sw;
    s21m = s21*cw + c21*sw;
    c21p = c21*cw + s21*sw;
    c21m = c21*cw - s21*sw;
    
    r2 = r*r;
    rd = 2.0*r*delta;
    rd2 = r2 + d2;
    
    rr[0][0] = r;
    rr[0][1] = sqrt(rd2 + rd*c2p);
    rr[0][2] = sqrt(rd2 + rd*c2m);
    rr[1][0] = sqrt(rd2 - rd*c1p);
    rr[2][0] = sqrt(rd2 - rd*c1m);
    
    rr[1][1] = sqrt(r2 + rd*(c2p - c1p) + 2.0*d2*(1.0 - c21));
    rr[1][2] = sqrt(r2 + rd*(c2m - c1p) + 2.0*d2*(1.0 - c21m));
    rr[2][1] = sqrt(r2 + rd*(c2p - c1m) + 2.0*d2*(1.0 - c21p));
    rr[2][2] = sqrt(r2 + rd*(c2m - c1m) + 2.0*d2*(1.0 - c21));
    
    for (i=0; i<3; i++) for (j=0; j<3; j++) 
    {   
        if (rr[i][j] < 1.0e-4) rr[i][j] = 1.0e-4;
        rr[i][j] = rr[i][j]*rr[i][j]*rr[i][j];
    }
    
    torque = rd*(s2p/rr[0][1] + s2m/rr[0][2]);
    torque += -0.5*rd*(s2p/rr[1][1] + s2p/rr[2][1] + s2m/rr[1][2] + s2m/rr[2][2]);
    torque += d2*(s21/rr[1][1] + s21/rr[2][2] + s21m/rr[1][2] + s21p/rr[2][1]);
    
    return(torque);
    
}

double water_force(double r, double ca, double sa, double ca_rel, double sa_rel, double ck_rel, double sk_rel, double f[2])
/* compute force and torque of water molecule #k on water molecule #j (for interaction I_LJ_WATER) */
{
    double x1[3], y1[3], x2[3], y2[3], rr[3][3], dx[3][3], dy[3][3], fx[3][3], fy[3][3], m[3][3], torque = 0.0;
    static double cw[3], sw[3], q[3], d[3], delta = 1.25*MU, dmin = 0.5*MU, fscale = 1.0;
    int i, j;
    static int first = 1;
    
    if (first)
    {
        cw[0] = 1.0; cw[1] = -0.5; cw[2] = -0.5;
        sw[0] = 0.0; sw[1] = 866025404; sw[2] = -866025404;     /* sines and cosines of angles */
        q[0] = -2.0; q[1] = 1.0; q[2] = 1.0;                    /* charges */
        d[0] = 0.5*delta; d[1] = delta; d[2] = delta;           /* distances to center */
        first = 0;
    }
    
    /* positions of O and H atoms */
    for (i=0; i<3; i++)
    {
        x1[i] = d[i]*(ca_rel*cw[i] - sa_rel*sw[i]);
        y1[i] = d[i]*(ca_rel*sw[i] + sa_rel*cw[i]);
        x2[i] = r + d[i]*(ck_rel*cw[i] - sk_rel*sw[i]);
        y2[i] = d[i]*(ck_rel*sw[i] + sk_rel*cw[i]);
    }
    
    /* relative positions */
    for (i=0; i<3; i++) for (j=0; j<3; j++) 
    {
        dx[i][j] = x2[j] - x1[i];
        dy[i][j] = y2[j] - y1[i];
        rr[i][j] = module2(dx[i][j], dy[i][j]);
        if (rr[i][j] < dmin) rr[i][j] = dmin;
        rr[i][j] = ipow(rr[i][j],3);
//         rr[i][j] = rr[i][j]*rr[i][j]*rr[i][j];
    }
    
    /* forces between particles */
    for (i=0; i<3; i++) for (j=0; j<3; j++) 
    {
        fx[i][j] = -q[i]*q[j]*dx[i][j]/rr[i][j];
        fy[i][j] = -q[i]*q[j]*dy[i][j]/rr[i][j];
    }
    
    /* torques between particles */
    for (i=0; i<3; i++) for (j=0; j<3; j++) 
    {
        m[i][j] = x1[i]*fy[i][j] - y1[i]*fx[i][j];
    }
    
    /* total force */
    f[0] = 0.0;
    f[1] = 0.0;
    for (i=0; i<3; i++) for (j=0; j<3; j++) 
    {
        f[0] += fscale*fx[i][j];
        f[1] += fscale*fy[i][j];
        torque += fscale*m[i][j];
    }
    
    return(torque);
}

int compute_particle_interaction(int i, int k, double force[2], double *torque, t_particle* particle, double distance, double krepel, double ca, double sa, double ca_rel, double sa_rel)
/* compute repelling force and torque of particle #k on particle #i */
/* returns 1 if distance between particles is smaller than NBH_DIST_FACTOR*MU */
{
    double x1, y1, x2, y2, r, f, angle, aniso, fx, fy, ff[2], dist_scaled, spin_f, ck, sk, ck_rel, sk_rel, alpha, amp;
    static double dxhalf = 0.5*(BCXMAX - BCXMIN), dyhalf = 0.5*(BCYMAX - BCYMIN);
    int wwrapx, wwrapy, twrapx, twrapy;
    
    if (BOUNDARY_COND == BC_GENUS_TWO)
    {
        dxhalf *= 0.75;
        dyhalf *= 0.75;
    }
    
    x1 = particle[i].xc;
    y1 = particle[i].yc;
    x2 = particle[k].xc;
    y2 = particle[k].yc;
        
    wwrapx = ((BOUNDARY_COND == BC_KLEIN)||(BOUNDARY_COND == BC_BOY)||(BOUNDARY_COND == BC_GENUS_TWO))&&(vabs(x2 - x1) > dxhalf);
    wwrapy = ((BOUNDARY_COND == BC_BOY)||(BOUNDARY_COND == BC_GENUS_TWO))&&(vabs(y2 - y1) > dyhalf);
    twrapx = ((BOUNDARY_COND == BC_KLEIN)||(BOUNDARY_COND == BC_BOY))&&(vabs(x2 - x1) > dxhalf);
    twrapy = (BOUNDARY_COND == BC_BOY)&&(vabs(y2 - y1) > dyhalf);
        
    switch (particle[k].interaction) {
        case (I_COULOMB):
        {
            f = -krepel/(1.0e-8 + distance*distance);
            force[0] = f*ca;
            force[1] = f*sa; 
            break;
        }
        case (I_LENNARD_JONES):
        {
            f = krepel*lennard_jones_force(distance, particle[k]);
            force[0] = f*ca;
            force[1] = f*sa;   
            break;
        }
        case (I_LJ_DIRECTIONAL):
        {
            aniso_lj_force(distance, ca, sa, ca_rel, sa_rel, ff, particle[k]);
            force[0] = krepel*ff[0];
            force[1] = krepel*ff[1];
            break;
        }
        case (I_LJ_PENTA):
        {
            penta_lj_force(distance, ca, sa, ca_rel, sa_rel, ff, particle[k]);
            force[0] = krepel*ff[0];
            force[1] = krepel*ff[1];
            break;
        }
        case (I_GOLDENRATIO):
        {
            f = krepel*golden_ratio_force(distance, particle[k]);
            force[0] = f*ca;
            force[1] = f*sa;   
            break;
        }
        case (I_LJ_DIPOLE):
        {
            dipole_lj_force(distance, ca, sa, ca_rel, sa_rel, ff, particle[k]);
            force[0] = krepel*ff[0];
            force[1] = krepel*ff[1];
            break;
        }
        case (I_LJ_QUADRUPOLE):
        {
            quadrupole_lj_force(distance, ca, sa, ca_rel, sa_rel, ff, particle[k]);
            force[0] = krepel*ff[0];
            force[1] = krepel*ff[1];
            break;
        }
        case (I_LJ_WATER):
        {
            f = krepel*lennard_jones_force(distance, particle[k]);
            force[0] = f*ca;
            force[1] = f*sa;   
            break;
        }
        case (I_VICSEK):
        {
            force[0] = 0.0;
            force[1] = 0.0;
            break;
        }
        case (I_VICSEK_REPULSIVE):
        {
            f = krepel*coulomb_force(distance, particle[k]);
//             f = krepel*harmonic_force(distance, particle[k]);
//             f = krepel*lennard_jones_force(distance, particle[k]);
            force[0] = f*ca;
            force[1] = f*sa;   
            break;
        }
        case (I_VICSEK_SPEED):
        {
            f = cos(0.5*(particle[k].angle - particle[i].angle));
            force[0] = f*KSPRING_VICSEK*(particle[k].vx - particle[i].vx);
            force[1] = f*KSPRING_VICSEK*(particle[k].vy - particle[i].vy);
            break;
        }
        case (I_VICSEK_SHARK):
        {
            if (particle[k].type == particle[i].type)
            {
                f = cos(0.5*(particle[k].angle - particle[i].angle));
                force[0] = f*KSPRING_VICSEK*(particle[k].vx - particle[i].vx);
                force[1] = f*KSPRING_VICSEK*(particle[k].vy - particle[i].vy);
            }
            else if (particle[i].type != 2)
            {
                f = krepel*coulomb_force(distance, particle[k]);
                force[0] = f*ca;
                force[1] = f*sa;
            }
            else 
            {
                if (VICSEK_REPULSION > 0.0)
                {
//                     f = VICSEK_REPULSION*harmonic_force(distance, particle[k]);
                    f = VICSEK_REPULSION*coulomb_force(distance, particle[k]);
                    force[0] = f*ca;
                    force[1] = f*sa;   
                }
                else
                {
                    force[0] = 0.0;
                    force[1] = 0.0;
                }
            }
            break;
        }
    }

    if (ROTATION) 
    {   
        dist_scaled = distance/(particle[i].spin_range*particle[i].radius);
        switch (particle[k].interaction) {
            case (I_LJ_WATER):
            {
                ck = cos(particle[k].angle);
                sk = sin(particle[k].angle);
                ck_rel = ca*ck + sa*sk;
                sk_rel = sa*ck - ca*sk;
//                 *torque = (-3.0*ca_rel*sk_rel + 2.0*sa_rel*ck_rel)/(1.0e-12 + dist_scaled*dist_scaled*dist_scaled);
//                 *torque = water_torque(distance, ca, sa, ca_rel, sa_rel, ck_rel, sk_rel);
//                 *torque = (0.5*sin(angle) + 0.5*sin(2.0*angle) - 0.45*sin(3.0*angle))/(1.0e-12 + dist_scaled*dist_scaled*dist_scaled);
                
                *torque = water_force(distance, ca, sa, ca_rel, sa_rel, ck_rel, sk_rel, ff);
                force[0] += ff[0];
                force[1] += ff[1];
//                 printf("force = (%.3lg, %.3lg)\n", ff[0], ff[1]);
                break;
            }
            case (I_VICSEK):
            {
                if (dist_scaled > 1.0) *torque = 0.0;
                else if (twrapx||twrapy) *torque = sin(-particle[k].angle - particle[i].angle);
                else *torque = sin(particle[k].angle - particle[i].angle);
                break;
            }
            case (I_VICSEK_REPULSIVE):
            {
                if (dist_scaled > 1.0) *torque = 0.0;
                else if (twrapx||twrapy) *torque = sin(-particle[k].angle - particle[i].angle);
                else *torque = sin(particle[k].angle - particle[i].angle);
                break;
            }
            case (I_VICSEK_SPEED):
            {
                if (dist_scaled > 1.0) *torque = 0.0;
                else if (twrapx||twrapy) *torque = sin(-particle[k].angle - particle[i].angle);
                else *torque = sin(particle[k].angle - particle[i].angle);
                break;
            }
            case (I_VICSEK_SHARK):
            {
                if (dist_scaled > 10.0) *torque = 0.0;
                else if (particle[k].type == particle[i].type)  /* fish adjusting direction */
                {
                    if (twrapx||twrapy) *torque = sin(-particle[k].angle - particle[i].angle);
                    else *torque = sin(particle[k].angle - particle[i].angle);
                }
                else if (particle[k].type == 2)     /* fish fleeing a shark */
                {
                    alpha = argument(ca,sa);
                    particle[i].angle = alpha + PI;
                    *torque = 0.0;
                }
                else    /* shark tracking fish */
                {
                    *torque = cos(particle[k].angle)*sa - sin(particle[k].angle)*ca;
                }
                    
                break;
            }
            default: 
            {
                spin_f = particle[i].spin_freq;
                if (twrapx||twrapy) *torque = sin(spin_f*(-particle[k].angle - particle[i].angle))/(1.0e-8 + dist_scaled*dist_scaled);
                else 
                *torque = sin(spin_f*(particle[k].angle - particle[i].angle))/(1.0e-8 + dist_scaled*dist_scaled);
            }
        }
        
        if (particle[i].type == particle[k].type) 
        {
            if (particle[i].type == 0) *torque *= KTORQUE;
            else *torque *= KTORQUE_B;
        }
        else *torque *= KTORQUE_DIFF;
    }
    else *torque = 0.0;

    if ((distance < NBH_DIST_FACTOR*particle[i].radius)&&(k != i)) return(1);
//     if ((distance < NBH_DIST_FACTOR*particle[i].radius)) return(1);
    else return(0);
}


int compute_repelling_force(int i, int j, double force[2], double *torque, t_particle* particle, double krepel)
/* compute repelling force of neighbour #j on particle #i */
/* returns 1 if distance between particles is smaller than NBH_DIST_FACTOR*MU */
{
    double distance, ca, sa, cj, sj, ca_rel, sa_rel, f[2], ff[2], torque1, ck, sk, ck_rel, sk_rel;
    static double distmin = 10.0*((XMAX - XMIN)/HASHX + (YMAX - YMIN)/HASHY);
    int interact, k;
        
    if (BOUNDARY_COND == BC_GENUS_TWO) distmin *= 2.0;
    
    k = particle[i].hashneighbour[j];
        
    distance = module2(particle[i].deltax[j], particle[i].deltay[j]);
    
    /* for monitoring purposes */
//     if (distance > distmin) 
//     {
//         printf("i = %i, hashcell %i, j = %i, hashcell %i\n", i, particle[i].hashcell, k, particle[k].hashcell);
//         printf("X = (%.3lg, %.3lg)\n", particle[i].xc, particle[i].yc);
//         printf("Y = (%.3lg, %.3lg) d = %.3lg\n", particle[k].xc, particle[k].yc, distance);
//     }
        
    if ((distance == 0.0)||(i == k))
    {
        force[0] = 0.0;
        force[1] = 0.0;
        *torque = 0.0;
        return(1);
    }
    else if (distance > REPEL_RADIUS*particle[i].radius) 
    {
        force[0] = 0.0;
        force[1] = 0.0;
        *torque = 0.0;
        return(0);
    }
    else
    {
        /* to avoid numerical problems, assign minimal value to distance */
        if (distance < 0.1*particle[i].radius) distance = 0.1*particle[i].radius;
            
        ca = (particle[i].deltax[j])/distance;
        sa = (particle[i].deltay[j])/distance;
        
        /* compute relative angle in case particles can rotate */
        if (ROTATION)
        {
            cj = cos(particle[j].angle);
            sj = sin(particle[j].angle);
            ca_rel = ca*cj + sa*sj;
            sa_rel = sa*cj - ca*sj;
        }
        else
        {
            ca_rel = ca;
            sa_rel = sa;
        }
    
        interact = compute_particle_interaction(i, k, f, torque, particle, distance, krepel, ca, sa, ca_rel, sa_rel);
        
        if (SYMMETRIZE_FORCE)
        {
            torque1 = *torque;
//             compute_particle_interaction(k, i, ff, torque, particle, distance, krepel, ca, sa, ca_rel, sa_rel);
            ck = cos(particle[j].angle);
            sk = sin(particle[j].angle);
            ck_rel = ca*ck + sa*sk;
            sk_rel = sa*ck - ca*sk;
            compute_particle_interaction(k, i, ff, torque, particle, distance, krepel, -ca, -sa, -ck_rel, -sk_rel);
            force[0] = 0.5*(f[0] - ff[0]);
            force[1] = 0.5*(f[1] - ff[1]);
            *torque = 0.5*(torque1 - *torque);
//             *torque = 0.5*(*torque + torque1);
        }
        else
        {
            force[0] = f[0];
            force[1] = f[1];
        }
        
//         printf("force = (%.3lg, %.3lg), torque = %.3lg\n", f[0], f[1], *torque);
        return(interact);
    }
}


int resample_particle(int n, int maxtrials, t_particle particle[NMAXCIRCLES])
/* resample y coordinate of particle n, returns 1 if no collision is created */
{
    double x, y, dist, dmin = 10.0;
    int i, j, closeby = 0, success = 0, trials = 0; 
    
    while ((!success)&&(trials < maxtrials)) 
    {
        success = 1;
        x = particle[n].xc - MU*(double)rand()/RAND_MAX;
        y = 0.95*(BCYMIN + (BCYMAX - BCYMIN)*(double)rand()/RAND_MAX);
        i = 0;
        while ((success)&&(i<ncircles)) 
        {
            if ((i!=n)&&(particle[i].active))
            {
//             dist = module2(x - particle[i].xc, y - particle[i].yc);
                for (j=-1; j<2; j++)
                {
                    dist = module2(x - particle[i].xc - (double)j*(BCXMAX - BCXMIN), y - particle[i].yc);
                    if (dist < dmin) dmin = dist;
                }
                if (dmin < SAFETY_FACTOR*MU) success = 0;
            }
            i++;
        }
        trials++;
//         printf("Trial no %i - (%.3lg, %.3lg)\t", trials, x, y);
    }
    
    if (success)
    {
        printf("\nTrial %i succesful\n", trials);
        printf("Moving particle %i from (%.3lg, %.3lg) to (%.3lg, %.3lg)\n\n", n, particle[n].xc, particle[n].yc, x, y);
        particle[n].xc = x;
        particle[n].yc = y;
        particle[n].vx = V_INITIAL*gaussian();
        particle[n].vy = V_INITIAL*gaussian();
//         particle[n].vy = V_INITIAL*(double)rand()/RAND_MAX;
        return(1);
    }
    else  
    {
        printf("\nCannot move particle %i\n\n", n);
        return(0);
    }
    
}

int add_particle(double x, double y, double vx, double vy, double mass, short int type, t_particle particle[NMAXCIRCLES])
{
    int i, closeby = 0;
    double dist;
    
    /* test distance to other particles */
    for (i=0; i<ncircles; i++)
    {
        dist = module2(x - particle[i].xc, y - particle[i].yc);
        if ((particle[i].active)&&(dist < SAFETY_FACTOR*MU)) closeby = 1;
    }
    
    if ((closeby)||(ncircles >= NMAXCIRCLES)) 
    {
        printf("Cannot add particle at (%.3lg, %.3lg)\n", x, y);
        return(0);
    }
    else
    {
        i = ncircles;
        
        particle[i].type = type;
        
        particle[i].xc = x;
        particle[i].yc = y;
        particle[i].radius = MU*sqrt(mass);
        particle[i].active = 1;
        particle[i].neighb = 0;
        particle[i].diff_neighb = 0;
        particle[i].thermostat = 1;

        particle[i].energy = 0.0;
        particle[i].emean = 0.0;
        particle[i].dirmean = 0.0;

        if (RANDOM_RADIUS) particle[i].radius = particle[i].radius*(0.75 + 0.5*((double)rand()/RAND_MAX));
        
        particle[i].mass_inv = 1.0/mass;
        if (particle[i].type == 0) particle[i].inertia_moment_inv = 1.0/PARTICLE_INERTIA_MOMENT;
        else particle[i].inertia_moment_inv = 1.0/PARTICLE_INERTIA_MOMENT;
                
        particle[i].vx = vx;
        particle[i].vy = vy;
        particle[i].energy = (particle[i].vx*particle[i].vx + particle[i].vy*particle[i].vy)*particle[i].mass_inv;
        
        particle[i].angle = DPI*(double)rand()/RAND_MAX;
        particle[i].omega = 0.0;
        
//         if (particle[i].type == 1)
//         {
//             particle[i].interaction = INTERACTION_B;
//             particle[i].eq_dist = EQUILIBRIUM_DIST_B;
//             particle[i].spin_range = SPIN_RANGE_B;
//             particle[i].spin_freq = SPIN_INTER_FREQUENCY_B;            
//         }
    
        if ((PLOT == P_NUMBER)||(PLOT_B == P_NUMBER))
            particle[i].color_hue = 360.0*(double)(i%N_PARTICLE_COLORS)/(double)N_PARTICLE_COLORS;
        
        ncircles++;
        
        printf("Added particle at (%.3lg, %.3lg)\n", x, y);
        printf("Number of particles: %i\n", ncircles);

        return(1);
    }
}

double neighbour_color(int nnbg)
{
    if (nnbg > 7) nnbg = 7;
    switch(nnbg){
        case (7): return(340.0);
        case (6): return(310.0);
        case (5): return(260.0);
        case (4): return(200.0);
        case (3): return(140.0);
        case (2): return(100.0);
        case (1): return(70.0);
        default:  return(30.0);
    }   
}


void compute_entropy(t_particle particle[NMAXCIRCLES], double entropy[2])
{
    int i, nleft1 = 0, nleft2 = 0;
    double p1, p2, x;
    static int first = 1, ntot1 = 0, ntot2 = 0;
    static double log2;
    
    if (first)
    {
        log2 = log(2.0);
        for (i=0; i<ncircles; i++) if (particle[i].type == 0) ntot1++;
        else ntot2++;
        first = 0;
    }
    
    for (i=0; i<ncircles; i++)
    {
        if (POSITION_Y_DEPENDENCE) x = particle[i].yc;
        else x = particle[i].xc;
        if (particle[i].type == 0) 
        {
            if (x < 0.0) nleft1++;
        }
        else
        {
            if (x < 0.0) nleft2++;
        }            
    }
    p1 = (double)nleft1/(double)ntot1;
    p2 = (double)nleft2/(double)ntot2;
    printf("Type 1: nleft = %i, ntot = %i, p = %.3lg\n", nleft1, ntot1, p1);
    printf("Type 2: nleft = %i, ntot = %i, p = %.3lg\n", nleft2, ntot2, p2);
    if ((p1==0.0)||(p1==1.0)) entropy[0] = 0.0;
    else entropy[0] = -(p1*log(p1) + (1.0-p1)*log(1.0-p1)/log2);
    if ((p2==0.0)||(p2==1.0)) entropy[1] = 0.0;
    else entropy[1] = -(p2*log(p2) + (1.0-p2)*log(1.0-p2)/log2);
}


void compute_particle_colors(t_particle particle, int plot, double rgb[3], double rgbx[3], double rgby[3], double *radius, int *width)
{
    double ej, angle, hue, huex, huey, lum;
    int i;
    
    switch (plot) {
        case (P_KINETIC): 
        {
            ej = particle.energy;
            if (ej > 0.0) 
            {
                hue = ENERGY_HUE_MIN + (ENERGY_HUE_MAX - ENERGY_HUE_MIN)*ej/PARTICLE_EMAX;
                if (hue > ENERGY_HUE_MIN) hue = ENERGY_HUE_MIN;
                if (hue < ENERGY_HUE_MAX) hue = ENERGY_HUE_MAX;
            }
            *radius = particle.radius;
            *width = BOUNDARY_WIDTH;
            break;
        }
        case (P_NEIGHBOURS): 
        {
            hue = neighbour_color(particle.neighb);
            *radius = particle.radius;
            *width = BOUNDARY_WIDTH;
            break;
        }
        case (P_BONDS):
        {
            hue = neighbour_color(particle.neighb);
            *radius = particle.radius;
            *width = 1;
            break;
        }
        case (P_ANGLE):
        {
            angle = particle.angle;
            hue = angle*particle.spin_freq/DPI;
            hue -= (double)((int)hue);
            huex = (DPI - angle)*particle.spin_freq/DPI;
            huex -= (double)((int)huex);
            angle = PI - angle;
            if (angle < 0.0) angle += DPI;
            huey = angle*particle.spin_freq/DPI;
            huey -= (double)((int)huey);
            hue = PARTICLE_HUE_MIN + (PARTICLE_HUE_MAX - PARTICLE_HUE_MIN)*(hue);
            huex = PARTICLE_HUE_MIN + (PARTICLE_HUE_MAX - PARTICLE_HUE_MIN)*(huex);
            huey = PARTICLE_HUE_MIN + (PARTICLE_HUE_MAX - PARTICLE_HUE_MIN)*(huey);
            *radius = particle.radius;
            *width = BOUNDARY_WIDTH;
            break;
        }
        case (P_TYPE):
        {
            hue = type_hue(particle.type);
            *radius = particle.radius;
            *width = BOUNDARY_WIDTH;
            break;
        }
        case (P_DIRECTION): 
        {
            hue = argument(particle.vx, particle.vy);
            if (hue > DPI) hue -= DPI;
            if (hue < 0.0) hue += DPI;
            hue = PARTICLE_HUE_MIN + (PARTICLE_HUE_MAX - PARTICLE_HUE_MIN)*(hue)/DPI;
            *radius = particle.radius;
            *width = BOUNDARY_WIDTH;
            break;
        }
        case (P_DIRECT_ENERGY): 
        {
            hue = argument(particle.vx, particle.vy);
            if (hue > DPI) hue -= DPI;
            if (hue < 0.0) hue += DPI;
            hue = PARTICLE_HUE_MIN + (PARTICLE_HUE_MAX - PARTICLE_HUE_MIN)*(hue)/DPI;
            if (particle.energy < 0.1*PARTICLE_EMAX) lum = 10.0*particle.energy/PARTICLE_EMAX;
            else lum = 1.0;
            *radius = particle.radius;
            *width = BOUNDARY_WIDTH;
            break;
        }
        case (P_ANGULAR_SPEED): 
        {
            hue = 160.0*(1.0 + tanh(SLOPE*particle.omega));
//                printf("omega = %.3lg, hue = %.3lg\n", particle[j].omega, hue);
            *radius = particle.radius;
            *width = BOUNDARY_WIDTH;
            break;
        }
        case (P_DIFF_NEIGHB): 
        {
            hue = (double)(particle.diff_neighb+1)/(double)(particle.neighb+1);
//                hue = PARTICLE_HUE_MIN + (PARTICLE_HUE_MAX - PARTICLE_HUE_MIN)*hue;
            hue = 180.0*(1.0 + hue);
            *radius = particle.radius;
            *width = BOUNDARY_WIDTH;
            break;
        }
        case (P_THERMOSTAT):
        {
            if (particle.thermostat) hue = 30.0;
            else hue = 270.0;
            *radius = particle.radius;
            *width = BOUNDARY_WIDTH;
            break;
        }
        case (P_INITIAL_POS):
        {
            hue = particle.color_hue;
            *radius = particle.radius;
            *width = BOUNDARY_WIDTH;
            break;
        }
        case (P_NUMBER):
        {
            hue = particle.color_hue;
            *radius = particle.radius;
            *width = BOUNDARY_WIDTH;
            break;
        }
        case (P_EMEAN): 
        {
            ej = particle.emean;
            if (ej > 0.0) 
            {
                hue = ENERGY_HUE_MIN + (ENERGY_HUE_MAX - ENERGY_HUE_MIN)*ej/PARTICLE_EMAX;
                if (hue > ENERGY_HUE_MIN) hue = ENERGY_HUE_MIN;
                if (hue < ENERGY_HUE_MAX) hue = ENERGY_HUE_MAX;
            }
            *radius = particle.radius;
            *width = BOUNDARY_WIDTH;
            break;
        }
        case (P_DIRECT_EMEAN): 
        {
            hue = particle.dirmean + COLOR_HUESHIFT*PI;
            if (hue > DPI) hue -= DPI;
            hue = PARTICLE_HUE_MIN + (PARTICLE_HUE_MAX - PARTICLE_HUE_MIN)*(hue)/DPI;
            ej = particle.emean;
            if (ej < 0.1*PARTICLE_EMAX) lum = 10.0*ej/PARTICLE_EMAX;
            else lum = 1.0;
            *radius = particle.radius;
            *width = BOUNDARY_WIDTH;
            break;
        }
        case (P_NOPARTICLE): 
        {
            hue = 0.0;
            lum = 1.0;
            *radius = particle.radius;
            *width = BOUNDARY_WIDTH;
            break;
        }
    }
    
    switch (plot) {
        case (P_KINETIC):  
        {
//             hsl_to_rgb_turbo(hue, 0.9, 0.5, rgb);
//             hsl_to_rgb_turbo(hue, 0.9, 0.5, rgbx);
//             hsl_to_rgb_turbo(hue, 0.9, 0.5, rgby);
            hsl_to_rgb_palette(hue, 0.9, 0.5, rgb, COLOR_PALETTE_EKIN);
            hsl_to_rgb_palette(hue, 0.9, 0.5, rgbx, COLOR_PALETTE_EKIN);
            hsl_to_rgb_palette(hue, 0.9, 0.5, rgby, COLOR_PALETTE_EKIN);
            break;
        }
        case (P_BONDS):  
        {
            hsl_to_rgb_turbo(hue, 0.9, 0.5, rgb);
            hsl_to_rgb_turbo(hue, 0.9, 0.5, rgbx);
            hsl_to_rgb_turbo(hue, 0.9, 0.5, rgby);
            break;
        }
        case (P_DIRECTION): 
        {
            hsl_to_rgb_palette(hue, 0.9, 0.5, rgb, COLOR_PALETTE_DIRECTION);
            hsl_to_rgb_palette(hue, 0.9, 0.5, rgbx, COLOR_PALETTE_DIRECTION);
            hsl_to_rgb_palette(hue, 0.9, 0.5, rgby, COLOR_PALETTE_DIRECTION);
            break;
        }
        case (P_ANGLE): 
        {
            hsl_to_rgb_palette(hue, 0.9, 0.5, rgb, COLOR_PALETTE_ANGLE);
            hsl_to_rgb_palette(hue, 0.9, 0.5, rgbx, COLOR_PALETTE_ANGLE);
            hsl_to_rgb_palette(hue, 0.9, 0.5, rgby, COLOR_PALETTE_ANGLE);
            break;
        }
        case (P_DIRECT_ENERGY): 
        {
            hsl_to_rgb_twilight(hue, 0.9, 0.5, rgb);
            hsl_to_rgb_twilight(hue, 0.9, 0.5, rgbx);
            hsl_to_rgb_twilight(hue, 0.9, 0.5, rgby);
            for (i=0; i<3; i++)
            {
                rgb[i] *= lum;
                rgbx[i] *= lum;
                rgby[i] *= lum;
            }
            break;
        }
        case (P_DIFF_NEIGHB):  
        {
            hsl_to_rgb_twilight(hue, 0.9, 0.5, rgb);
            hsl_to_rgb_twilight(hue, 0.9, 0.5, rgbx);
            hsl_to_rgb_twilight(hue, 0.9, 0.5, rgby);
            break;
        }
        case (P_INITIAL_POS): 
        {
            hsl_to_rgb_palette(hue, 0.9, 0.5, rgb, COLOR_PALETTE_INITIAL_POS);
            hsl_to_rgb_palette(hue, 0.9, 0.5, rgbx, COLOR_PALETTE_INITIAL_POS);
            hsl_to_rgb_palette(hue, 0.9, 0.5, rgby, COLOR_PALETTE_INITIAL_POS);
            break;
        }
        case (P_EMEAN):  
        {
            hsl_to_rgb_palette(hue, 0.9, 0.5, rgb, COLOR_PALETTE_EKIN);
            hsl_to_rgb_palette(hue, 0.9, 0.5, rgbx, COLOR_PALETTE_EKIN);
            hsl_to_rgb_palette(hue, 0.9, 0.5, rgby, COLOR_PALETTE_EKIN);
            break;
        }
        case (P_DIRECT_EMEAN): 
        {
            hsl_to_rgb_palette(hue, 0.9, 0.5, rgb, COLOR_PALETTE_DIRECTION);
            hsl_to_rgb_palette(hue, 0.9, 0.5, rgbx, COLOR_PALETTE_DIRECTION);
            hsl_to_rgb_palette(hue, 0.9, 0.5, rgby, COLOR_PALETTE_DIRECTION);
            for (i=0; i<3; i++)
            {
                rgb[i] *= lum;
                rgbx[i] *= lum;
                rgby[i] *= lum;
            }
            break;
        }
        default: 
        {
            hsl_to_rgb(hue, 0.9, 0.5, rgb);
            hsl_to_rgb(hue, 0.9, 0.5, rgbx);
            hsl_to_rgb(hue, 0.9, 0.5, rgby);
        }
    }
}


void set_segment_group_color(int group, double lum, double rgb[3])
{
    switch (group) {
        case (1):
        {
            hsl_to_rgb_palette(270.0, 0.9, 0.5, rgb, COLOR_PALETTE);
            break;
        }
        case (2):
        {
            hsl_to_rgb_palette(90.0, 0.9, 0.5, rgb, COLOR_PALETTE);
            break;
        }
        default:
        {
            rgb[0] = 1.0;
            rgb[1] = 1.0;
            rgb[2] = 1.0;
        }
    }
    
    glColor3f(lum*rgb[0], lum*rgb[1], lum*rgb[2]);
}

void draw_altitude_lines()
{
    int i, i1, i2;
    double x, y;
    
    glColor3f(0.5, 0.5, 0.5);
    glLineWidth(1);
    
    i1 = (int)(YMIN + ytrack) - 1.0;
    i2 = (int)(YMAX + ytrack) + 1.0;
    
    for (i = i1; i < i2; i++)
    {
        y = (double)i - ytrack;
        draw_line(BCXMIN, y, XMAX - 1.8, y);
    }
    
    i1 = (int)(XMIN + xtrack) - 1.0;
    i2 = (int)(XMAX - 1.8 + xtrack) + 1.0;
    
    for (i = i1; i < i2; i++)
    {
        x = (double)i - xtrack;
        if (x < XMAX - 1.8) draw_line(x, YMIN, x, YMAX);
    }
}

void draw_one_triangle(t_particle particle, int same_table[9*HASHMAX], int p, int q, int nsame)
{
    double x, y, dx, dy;
    int k;
    
    x = particle.xc + (double)p*(BCXMAX - BCXMIN);
    y = particle.yc + (double)q*(BCYMAX - BCYMIN);
    
    if (TRACK_SEGMENT_GROUPS) 
    {
        x -= xtrack;
        y -= ytrack;
    }
            
    glBegin(GL_TRIANGLE_FAN);
    glVertex2d(x, y);

    for (k=0; k<nsame; k++) 
    {
        dx = particle.deltax[same_table[k]];
        dy = particle.deltay[same_table[k]];
        if (module2(dx, dy) < particle.radius*NBH_DIST_FACTOR)
            glVertex2d(x + dx, y + dy);
    }
    glEnd();
}


void draw_triangles(t_particle particle[NMAXCIRCLES], int plot)
/* fill triangles between neighboring particles */
{
    int i, j, k, p, q, t0, tj, tmax, nsame = 0, same_table[9*HASHMAX], width;
    double rgb[3], hue, dx, dy, x, y, radius, rgbx[3], rgby[3];
    
//     printf("Number of nbs: ");
    for (i=0; i<ncircles; i++) if (particle[i].active)
    {
        nsame = 0;
        t0 = particle[i].type;
        
        /* find neighbours of same type */
        for (j=0; j<particle[i].hash_nneighb; j++)
        {
            k = particle[i].hashneighbour[j];
            if ((k!=i)&&((plot != P_TYPE)||(particle[k].type == t0))) 
            {
                same_table[nsame] = j;
                nsame++;
            }
        }
                
        /* draw the triangles */
        if (nsame > 1)
        {
            if (plot == P_TYPE)
            {
                hue = type_hue(t0);
                hsl_to_rgb(hue, 0.9, 0.3, rgb);
            }
            else
            {
                compute_particle_colors(particle[i], plot, rgb, rgbx, rgby, &radius, &width);
            }
                
            glColor3f(rgb[0], rgb[1], rgb[2]);
            if (PERIODIC_BC) for (p=-1; p<2; p++) for (q=-1; q<2; q++)
                draw_one_triangle(particle[i], same_table, p, q, nsame);
            else draw_one_triangle(particle[i], same_table, 0, 0, nsame);    
        }
    }
    printf("\n");
}


void draw_one_particle_links(t_particle particle)
/* draw links of one particle */
{
    int i, j, k;
    double x1, x2, y1, y2, length, linkcolor,periodx, periody, xt1, yt1, xt2, yt2;

    glLineWidth(LINK_WIDTH);
//     if (particle.active)
//     {
//             radius = particle[j].radius;
    for (k = 0; k < particle.hash_nneighb; k++)
    {
        x1 = particle.xc;
        if (CENTER_VIEW_ON_OBSTACLE) x1 -= xshift;
        y1 = particle.yc;
        x2 = x1 + particle.deltax[k];
        y2 = y1 + particle.deltay[k];
                
        length = module2(particle.deltax[k], particle.deltay[k])/particle.radius;
                
        if (COLOR_BONDS)
        {
            if (length < 1.5) linkcolor = 1.0;
            else linkcolor = 1.0 - 0.75*(length - 1.5)/(NBH_DIST_FACTOR - 1.5);
            glColor3f(linkcolor, linkcolor, linkcolor);
        }
                
        if (length < 1.0*NBH_DIST_FACTOR) draw_line(x1, y1, x2, y2);
    }
}


int draw_special_particle(t_particle particle, double xc1, double yc1, double radius, double angle, int nsides, double rgb[3], short int fill)
/* draw special particles shapes, for chemical reactions */
{
    double x1, y1, x2, y2, omega, r;
    int wsign, i, j;
    
    switch(RD_REACTION)
    {
        case (CHEM_AAB):
        {
            if (particle.type == 2) 
            {
                for (wsign = -1; wsign <= 1; wsign+=2)
                {
                    x1 = xc1 + (double)wsign*0.7*radius*cos(particle.angle);
                    y1 = yc1 + (double)wsign*0.7*radius*sin(particle.angle);
                    if (fill) draw_colored_polygon(x1, y1, 0.7*radius, nsides, angle + APOLY*PID, rgb);
                    else draw_polygon(x1, y1, 0.7*radius, nsides, angle + APOLY*PID);
                }
                return(0);
            }
            break;
        }
        case (CHEM_ABC):
        {
            if (particle.type == 3) 
            {
                for (wsign = -1; wsign <= 1; wsign+=2)
                {
                    x1 = xc1 + (double)wsign*0.7*radius*cos(particle.angle);
                    y1 = yc1 + (double)wsign*0.7*radius*sin(particle.angle);
                    if (wsign == 1) 
                    {
                        if (fill) draw_colored_polygon(x1, y1, 1.2*MU_B, nsides, angle + APOLY*PID, rgb);
                        else draw_polygon(x1, y1, 1.2*MU_B, nsides, angle + APOLY*PID);
                    }
                    else 
                    {
                        if (fill) draw_colored_polygon(x1, y1, 1.2*MU, nsides, angle + APOLY*PID, rgb);
                        else draw_polygon(x1, y1, 1.2*MU, nsides, angle + APOLY*PID);
                    }
                }
                return(0);
            }
            break;
        }
        case (CHEM_A2BC):
        {
            if (particle.type == 3) 
            {
                draw_colored_polygon(xc1, yc1, 1.2*MU_B, nsides, angle + APOLY*PID, rgb);
                for (wsign = -1; wsign <= 1; wsign+=2)
                {
                    x1 = xc1 + 1.5*radius*cos(particle.angle + 0.6*(double)wsign*PI);
                    y1 = yc1 + 1.5*radius*sin(particle.angle + 0.6*(double)wsign*PI);
                    if (fill) draw_colored_polygon(x1, y1, 1.2*MU, nsides, angle + APOLY*PID, rgb);
                    else draw_polygon(x1, y1, 1.2*MU, nsides, angle + APOLY*PID);
                }
                return(0);
            }
            break;
        }
        case (CHEM_CATALYSIS):
        {
            if (particle.type == 3) 
            {
                for (wsign = -1; wsign <= 1; wsign+=2)
                {
                    x1 = xc1 + (double)wsign*0.7*radius*cos(particle.angle);
                    y1 = yc1 + (double)wsign*0.7*radius*sin(particle.angle);
                    if (fill) draw_colored_polygon(x1, y1, 1.2*MU, nsides, angle + APOLY*PID, rgb);
                    else draw_polygon(x1, y1, 1.2*MU, nsides, angle + APOLY*PID);
                }
                return(0);
            }
            break;
        }
        case (CHEM_AUTOCATALYSIS):
        {
            if (particle.type == 2)
            {
                for (wsign = -1; wsign <= 1; wsign+=2)
                {
                    x1 = xc1 + (double)wsign*MU*0.7*cos(particle.angle);
                    y1 = yc1 + (double)wsign*MU*0.7*sin(particle.angle);
                    if (fill) draw_colored_polygon(x1, y1, MU, nsides, angle + APOLY*PID, rgb);
                    else draw_polygon(x1, y1, MU, nsides, angle + APOLY*PID);
                }
                return(0);
            }
            break;
        }
        case (CHEM_BAA):
        {
            if (particle.type == 2) 
            {
                for (wsign = -1; wsign <= 1; wsign+=2)
                {
                    x1 = xc1 + (double)wsign*1.2*radius*cos(particle.angle);
                    y1 = yc1 + (double)wsign*1.2*radius*sin(particle.angle);
                    if (fill) draw_colored_polygon(x1, y1, 0.9*radius, nsides, angle + APOLY*PID, rgb);
                    else draw_polygon(x1, y1, 0.9*radius, nsides, angle + APOLY*PID);
                }
                return(0);
            }
            break;
        }
        case (CHEM_AABAA):
        {
            if (particle.type == 2) 
            {
                for (wsign = -1; wsign <= 1; wsign+=2)
                {
                    x1 = xc1 + (double)wsign*0.7*radius*cos(particle.angle);
                    y1 = yc1 + (double)wsign*0.7*radius*sin(particle.angle);
                    if (fill) draw_colored_polygon(x1, y1, 0.9*radius, nsides, angle + APOLY*PID, rgb);
                    else draw_polygon(x1, y1, 0.9*radius, nsides, angle + APOLY*PID);
                }
                return(0);
            }
            break;
        }
        case (CHEM_POLYMER):
        {
            if (particle.type >= 3) 
            {
                omega = DPI/(double)(particle.type - 2);
                draw_colored_polygon(xc1, yc1, 1.2*MU_B, nsides, angle + APOLY*PID, rgb);
                for (i=0; i<particle.type-2; i++)
                {
                    x1 = xc1 + 1.5*MU_B*cos(particle.angle + (double)i*omega);
                    y1 = yc1 + 1.5*MU_B*sin(particle.angle + (double)i*omega);
                    if (fill) draw_colored_polygon(x1, y1, 1.2*MU, nsides, angle + APOLY*PID, rgb);
                    else draw_polygon(x1, y1, 1.2*MU, nsides, angle + APOLY*PID);
                }
                return(0);
            }
            break;
        }
        case (CHEM_POLYMER_DISS):
        {
            if (particle.type >= 3) 
            {
                omega = DPI/(double)(particle.type - 2);
                draw_colored_polygon(xc1, yc1, 1.2*MU_B, nsides, angle + APOLY*PID, rgb);
                for (i=0; i<particle.type-2; i++)
                {
                    x1 = xc1 + 1.5*MU_B*cos(particle.angle + (double)i*omega);
                    y1 = yc1 + 1.5*MU_B*sin(particle.angle + (double)i*omega);
                    if (fill) draw_colored_polygon(x1, y1, 1.2*MU, nsides, angle + APOLY*PID, rgb);
                    else draw_polygon(x1, y1, 1.2*MU, nsides, angle + APOLY*PID);
                }
                return(0);
            }
            break;
        }
        case (CHEM_POLYMER_STEP):
        {
            if (particle.type >= 2) 
            {
                omega = DPI/(double)(particle.type - 1);
                draw_colored_polygon(xc1, yc1, MU, nsides, angle + APOLY*PID, rgb);
                for (i=0; i<particle.type-1; i++)
                {
                    x1 = xc1 + 1.5*MU*cos(particle.angle + (double)i*omega);
                    y1 = yc1 + 1.5*MU*sin(particle.angle + (double)i*omega);
                    if (fill) draw_colored_polygon(x1, y1, MU, nsides, angle + APOLY*PID, rgb);
                    else draw_polygon(x1, y1, MU, nsides, angle + APOLY*PID);
                }
                return(0);
            }
            break;
        }
        case (CHEM_CATALYTIC_A2D):
        {
            if (particle.type == 3)
            {
                for (wsign = -1; wsign <= 1; wsign+=2)
                {
                    x1 = xc1 + (double)wsign*MU*0.7*cos(particle.angle);
                    y1 = yc1 + (double)wsign*MU*0.7*sin(particle.angle);
                    if (wsign == -1) r = MU; 
                    else r = MU_B;
                    if (fill) draw_colored_polygon(x1, y1, r, nsides, angle + APOLY*PID, rgb);
                    else draw_polygon(x1, y1, r, nsides, angle + APOLY*PID);
                }
                return(0);
            }
            else if (particle.type == 4)
            {
                for (wsign = -1; wsign <= 1; wsign+=2)
                {
                    x1 = xc1 + (double)wsign*MU*0.7*cos(particle.angle);
                    y1 = yc1 + (double)wsign*MU*0.7*sin(particle.angle);
                    if (fill) draw_colored_polygon(x1, y1, MU, nsides, angle + APOLY*PID, rgb);
                    else draw_polygon(x1, y1, MU, nsides, angle + APOLY*PID);
                }
                return(0);
            }
            break;
        }
        case (CHEM_ABCAB):
        {
            if (particle.type == 3) 
            {
                for (wsign = -1; wsign <= 1; wsign+=2)
                {
                    x1 = xc1 + (double)wsign*0.7*radius*cos(particle.angle);
                    y1 = yc1 + (double)wsign*0.7*radius*sin(particle.angle);
                    if (wsign == 1) 
                    {
                        if (fill) draw_colored_polygon(x1, y1, 1.2*MU_B, nsides, angle + APOLY*PID, rgb);
                        else draw_polygon(x1, y1, 1.2*MU_B, nsides, angle + APOLY*PID);
                    }
                    else 
                    {
                        if (fill) draw_colored_polygon(x1, y1, 1.2*MU, nsides, angle + APOLY*PID, rgb);
                        else draw_polygon(x1, y1, 1.2*MU, nsides, angle + APOLY*PID);
                    }
                }
                return(0);
            }
            break;
        }
        case (CHEM_ABCDABC):
        {
            if (particle.type == 3) 
            {
                for (wsign = -1; wsign <= 1; wsign+=2)
                {
                    x1 = xc1 + (double)wsign*0.7*radius*cos(particle.angle);
                    y1 = yc1 + (double)wsign*0.7*radius*sin(particle.angle);
                    if (wsign == 1) 
                    {
                        if (fill) draw_colored_polygon(x1, y1, MU_B, nsides, angle + APOLY*PID, rgb);
                        else draw_polygon(x1, y1, MU_B, nsides, angle + APOLY*PID);
                    }
                    else 
                    {
                        if (fill) draw_colored_polygon(x1, y1, MU, nsides, angle + APOLY*PID, rgb);
                        else draw_polygon(x1, y1, MU, nsides, angle + APOLY*PID);
                    }
                }
                return(0);
            }
            else if (particle.type == 4)
            {
                draw_colored_polygon(xc1, yc1, 1.2*MU_B, nsides, angle + APOLY*PID, rgb);
                for (wsign = -1; wsign <= 1; wsign+=2)
                {
                    x1 = xc1 + 0.7*radius*cos(particle.angle + 0.6*(double)wsign*PI);
                    y1 = yc1 + 0.7*radius*sin(particle.angle + 0.6*(double)wsign*PI);
                    if (fill) draw_colored_polygon(x1, y1, MU, nsides, angle + APOLY*PID, rgb);
                    else draw_polygon(x1, y1, MU, nsides, angle + APOLY*PID);
                }
                return(0);
            }
            break;
        }
        case (CHEM_BZ):
        {
            if ((particle.type >= 2)&&(particle.type <= 4)) 
            {
                omega = DPI/(double)(particle.type - 1);
                if (fill) draw_colored_polygon(xc1, yc1, MU, nsides, angle + APOLY*PID, rgb);
                else draw_polygon(xc1, yc1, MU, nsides, angle + APOLY*PID);
                for (i=0; i<particle.type-1; i++)
                {
                    x1 = xc1 + 1.5*MU*cos(particle.angle + (double)i*omega);
                    y1 = yc1 + 1.5*MU*sin(particle.angle + (double)i*omega);
                    if (fill) draw_colored_polygon(x1, y1, MU_B, nsides, angle + APOLY*PID, rgb);
                    else draw_polygon(x1, y1, MU_B, nsides, angle + APOLY*PID);
                }
                return(0);
            }
            else if (particle.type == 5)
            {
                if (fill) draw_colored_polygon(xc1, yc1, MU_B, 4, angle + APOLY*PID, rgb);
                else draw_polygon(xc1, yc1, MU_B, 4, angle + APOLY*PID);
                return(0);
            }
            else if (particle.type == 6)
            {
                if (fill) draw_colored_polygon(xc1, yc1, 1.5*MU_B, 6, angle + APOLY*PID, rgb);
                else draw_polygon(xc1, yc1, 1.5*MU_B, 6, angle + APOLY*PID);
                return(0);
            }
            else if (particle.type == 7)
            {
                for (i=-1; i<2; i+=2)
                {
                    x1 = xc1 + (double)i*1.5*MU*cos(particle.angle);
                    y1 = yc1 + (double)i*1.5*MU*sin(particle.angle);
                    if (fill) draw_colored_polygon(x1, y1, MU, nsides, angle + APOLY*PID, rgb);
                    else draw_polygon(x1, y1, MU, nsides, angle + APOLY*PID);
                    for (j=-1; j<2; j++)
                    {
                        x2 = x1 + 1.5*MU*cos(particle.angle + (double)j*PID);
                        y2 = y1 + 1.5*MU*sin(particle.angle + (double)j*PID);
                        if (fill) draw_colored_polygon(x2, y2, MU_B, nsides, angle + APOLY*PID, rgb);
                        else draw_polygon(x2, y2, MU_B, nsides, angle + APOLY*PID);
                    }
                
                }
                
                return(0);
            }
            else if (particle.type == 8)
            {
                x1 = xc1 + 1.5*MU*cos(particle.angle);
                y1 = yc1 + 1.5*MU*sin(particle.angle);
                if (fill) draw_colored_polygon(x1, y1, MU_B, 4, angle + APOLY*PID, rgb);
                else draw_polygon(x1, y1, MU_B, 4, angle + APOLY*PID);
                x1 = xc1 - 1.5*MU*cos(particle.angle);
                y1 = yc1 - 1.5*MU*sin(particle.angle);
                if (fill) draw_colored_polygon(x1, y1, MU_B, 4, angle + APOLY*PID, rgb);
                else draw_polygon(x1, y1, MU_B, 4, angle + APOLY*PID);
                return(0);
            }
            break;
        }
        case (CHEM_BRUSSELATOR):
        {
            if (particle.type == 4)
            {
                if (fill) draw_colored_polygon(xc1, yc1, MU_B, nsides, angle + APOLY*PID, rgb);
                else draw_polygon(xc1, yc1, MU_B, nsides, angle + APOLY*PID);
                x1 = xc1 + 1.2*MU_B*cos(particle.angle);
                y1 = yc1 + 1.2*MU_B*sin(particle.angle);
                if (fill) draw_colored_polygon(x1, y1, MU, nsides, angle + APOLY*PID, rgb);
                else draw_polygon(x1, y1, MU, nsides, angle + APOLY*PID);
                return(0);
            }
            break;
        }
        case (CHEM_ABDACBE):
        {
            if (particle.type == 4) 
            {
                for (wsign = -1; wsign <= 1; wsign+=2)
                {
                    x1 = xc1 + (double)wsign*0.7*radius*cos(particle.angle);
                    y1 = yc1 + (double)wsign*0.7*radius*sin(particle.angle);
                    if (wsign == 1) 
                    {
                        if (fill) draw_colored_polygon(x1, y1, 1.2*MU_B, nsides, angle + APOLY*PID, rgb);
                        else draw_polygon(x1, y1, 1.2*MU_B, nsides, angle + APOLY*PID);
                    }
                    else 
                    {
                        if (fill) draw_colored_polygon(x1, y1, 1.2*MU, nsides, angle + APOLY*PID, rgb);
                        else draw_polygon(x1, y1, 1.2*MU, nsides, angle + APOLY*PID);
                    }
                }
                return(0);
            }
            break;
        }
    }
    return(1);
}

void draw_one_particle(t_particle particle, double xc, double yc, double radius, double angle, int nsides, double width, double rgb[3])
/* draw one of the particles */ 
{
    double ca, sa, x1, x2, y1, y2, xc1, yc1, wangle, newradius = radius;
    int wsign, cont = 1, draw = 1;
    
    if (CENTER_VIEW_ON_OBSTACLE) xc1 = xc - xshift;
    else xc1 = xc;
    if (TRACK_SEGMENT_GROUPS) 
    {
        xc1 -= xtrack;
        yc1 = yc - ytrack;
    }
    else yc1 = yc;
    glColor3f(rgb[0], rgb[1], rgb[2]);
    
    /* specific shapes for chemical reactions */
    if (REACTION_DIFFUSION) cont = draw_special_particle(particle, xc1, yc1, radius, angle, nsides, rgb, 1);    
   
    if ((particle.interaction == I_LJ_QUADRUPOLE)||(particle.interaction == I_LJ_DIPOLE)||(particle.interaction == I_VICSEK)||(particle.interaction == I_VICSEK_REPULSIVE)||(particle.interaction == I_VICSEK_SPEED)||(particle.interaction == I_VICSEK_SHARK)) 
        draw_colored_rhombus(xc1, yc1, radius, angle + APOLY*PID, rgb);
    else if (cont) 
    {
        if (nsides == NSEG) draw_colored_circle_precomp(xc1, yc1, radius, rgb);
        else draw_colored_polygon(xc1, yc1, radius, nsides, angle + APOLY*PID, rgb);
    }
    
    /* draw crosses on particles of second type */
    if ((TWO_TYPES)&&(DRAW_CROSS))
        if (particle.type == 1)
        {
            if (ROTATION) angle = angle + APOLY*PID;
            else angle = APOLY*PID;
            ca = cos(angle);
            sa = sin(angle);
            glLineWidth(3);
            glColor3f(0.0, 0.0, 0.0);
            x1 = xc1 - MU_B*ca;
            y1 = yc1 - MU_B*sa;
            x2 = xc1 + MU_B*ca;
            y2 = yc1 + MU_B*sa;
            draw_line(x1, y1, x2, y2);
            x1 = xc1 - MU_B*sa;
            y1 = yc1 + MU_B*ca;
            x2 = xc1 + MU_B*sa;
            y2 = yc1 - MU_B*ca;
            draw_line(x1, y1, x2, y2);
        }
        
    glLineWidth(width);
    glColor3f(1.0, 1.0, 1.0);
    if (REACTION_DIFFUSION) cont = draw_special_particle(particle, xc1, yc1, radius, angle, nsides, rgb, 0);

    if ((particle.interaction == I_LJ_QUADRUPOLE)||(particle.interaction == I_LJ_DIPOLE)||(particle.interaction == I_VICSEK)||(particle.interaction == I_VICSEK_REPULSIVE)||(particle.interaction == I_VICSEK_SPEED)||(particle.interaction == I_VICSEK_SHARK)) 
        draw_rhombus(xc1, yc1, radius, angle + APOLY*PID);
    else if (cont) 
    {
        if (nsides == NSEG) draw_circle_precomp(xc1, yc1, radius);
        else draw_polygon(xc1, yc1, radius, nsides, angle + APOLY*PID); 
    }
    
    if (particle.interaction == I_LJ_WATER) for (wsign = -1; wsign <= 1; wsign+=2)
    {
        wangle = particle.angle + (double)wsign*DPI/3.0;
        x1 = xc1 + particle.radius*cos(wangle);
        y1 = yc1 + particle.radius*sin(wangle);
        draw_colored_polygon(x1, y1, 0.5*radius, nsides, angle + APOLY*PID, rgb);
        glColor3f(1.0, 1.0, 1.0);
        draw_polygon(x1, y1, 0.5*radius, nsides, angle + APOLY*PID);
    }
}

void draw_collisions(t_collision *collisions, int ncollisions)
/* draw discs where collisions happen */
{
    int i, j;
    double rgb[3], lum, x, y, x1, y1;
    
    
    
    for (i=0; i<ncollisions; i++) if (collisions[i].time > 0)
    {
        lum = (double)collisions[i].time/(double)COLLISION_TIME;
        if (collisions[i].color == 0.0) for (j=0; j<3; j++) rgb[j] = lum;
        else hsl_to_rgb_palette(collisions[i].color, 0.9, lum, rgb, COLOR_PALETTE);
        x = collisions[i].x;
        y = collisions[i].y;
        
        if (CENTER_VIEW_ON_OBSTACLE) x1 = x - xshift;
        else x1 = x;
        if (TRACK_SEGMENT_GROUPS) 
        {
            x1 -= xtrack;
            y1 = y - ytrack;
        }
        else y1 = y;
    
        draw_colored_polygon(x1, y1, 5.0*MU, NSEG, 0.0, rgb);
        collisions[i].time--;
    }
    
}

void draw_trajectory(t_tracer trajectory[TRAJECTORY_LENGTH*N_TRACER_PARTICLES], int traj_position, int traj_length)
/* draw tracer particle trajectory */
{
    int i, j, time;
    double x1, x2, y1, y2, rgb[3], lum;
    
//     blank();
    glLineWidth(TRAJECTORY_WIDTH);
    printf("drawing trajectory\n");
    
    for (j=0; j<N_TRACER_PARTICLES; j++)
    {
//         if (j == 0) hsl_to_rgb(HUE_TYPE1, 0.9, 0.5, rgb);
//         else if (j == 1) hsl_to_rgb(HUE_TYPE2, 0.9, 0.5, rgb);
//         else hsl_to_rgb(HUE_TYPE3, 0.9, 0.5, rgb);
        set_type_color(j+2, 0.5, rgb);
        glColor3f(rgb[0], rgb[1], rgb[2]);
    
        if (traj_length < TRAJECTORY_LENGTH) 
            for (i=0; i < traj_length-1; i++)
            {
                x1 = trajectory[j*TRAJECTORY_LENGTH + i].xc;
                x2 = trajectory[j*TRAJECTORY_LENGTH + i+1].xc;
                y1 = trajectory[j*TRAJECTORY_LENGTH + i].yc;
                y2 = trajectory[j*TRAJECTORY_LENGTH + i+1].yc;
            
                time = traj_length - i;
                lum = 1.0 - (double)time/(double)TRAJECTORY_LENGTH;
                glColor3f(lum*rgb[0], lum*rgb[1], lum*rgb[2]);
        
                if (module2(x2 - x1, y2 - y1) < 0.25*(YMAX - YMIN)) draw_line(x1, y1, x2, y2);
            
//             printf("(x1, y1) = (%.3lg, %.3lg), (x2, y2) = (%.3lg, %.3lg)\n", x1, y1, x2, y2);
            }
        else 
        {
            for (i = traj_position + 1; i < traj_length-1; i++)
            {
                x1 = trajectory[j*TRAJECTORY_LENGTH + i].xc;
                x2 = trajectory[j*TRAJECTORY_LENGTH + i+1].xc;
                y1 = trajectory[j*TRAJECTORY_LENGTH + i].yc;
                y2 = trajectory[j*TRAJECTORY_LENGTH + i+1].yc;
        
                time = traj_position + traj_length - i;
                lum = 1.0 - (double)time/(double)TRAJECTORY_LENGTH;
                glColor3f(lum*rgb[0], lum*rgb[1], lum*rgb[2]);
        
                if (module2(x2 - x1, y2 - y1) < 0.1*(YMAX - YMIN)) draw_line(x1, y1, x2, y2);
            }
            for (i=0; i < traj_position-1; i++)
            {
                x1 = trajectory[j*TRAJECTORY_LENGTH + i].xc;
                x2 = trajectory[j*TRAJECTORY_LENGTH + i+1].xc;
                y1 = trajectory[j*TRAJECTORY_LENGTH + i].yc;
                y2 = trajectory[j*TRAJECTORY_LENGTH + i+1].yc;
        
                time = traj_position - i;
                lum = 1.0 - (double)time/(double)TRAJECTORY_LENGTH;
                glColor3f(lum*rgb[0], lum*rgb[1], lum*rgb[2]);
        
                if (module2(x2 - x1, y2 - y1) < 0.1*(YMAX - YMIN)) draw_line(x1, y1, x2, y2);
            }
        }
    }
}
    

void draw_particles(t_particle particle[NMAXCIRCLES], int plot, double beta, t_collision *collisions, int ncollisions)
{
    int i, j, k, m, width, nnbg, nsides;
    double ej, hue, huex, huey, rgb[3], rgbx[3], rgby[3], radius, x1, y1, x2, y2, angle, ca, sa, length, linkcolor, sign = 1.0, angle1, signy = 1.0, periodx, periody, x, y, lum, logratio;
    char message[100];
    
//     if (!TRACER_PARTICLE) blank();
    if (plot == P_NOPARTICLE) blank();
    glColor3f(1.0, 1.0, 1.0);
    
    /* show region of partial thermostat */
    if ((PARTIAL_THERMO_COUPLING)&&(PARTIAL_THERMO_REGION == TH_INBOX))
    {
        if (INCREASE_BETA)
        {
            logratio = log(beta/BETA)/log(0.5*BETA_FACTOR);
            if (logratio > 1.0) logratio = 1.0;
            else if (logratio < 0.0) logratio = 0.0;
            if (BETA_FACTOR > 1.0) hue = PARTICLE_HUE_MAX - (PARTICLE_HUE_MAX - PARTICLE_HUE_MIN)*logratio;
            else hue = PARTICLE_HUE_MIN - (PARTICLE_HUE_MIN - PARTICLE_HUE_MAX)*logratio;
        }
        else hue = 0.25*PARTICLE_HUE_MIN + 0.75*PARTICLE_HUE_MAX;
        erase_area_hsl_turbo(0.0, YMIN, 2.0*PARTIAL_THERMO_WIDTH, PARTIAL_THERMO_HEIGHT*(YMAX - YMIN),  hue, 0.9, 0.15);
        
    }
    
    /* draw "altitude lines" */
    if (ALTITUDE_LINES) draw_altitude_lines();
    
    /* draw the bonds first */
    if ((DRAW_BONDS)||(plot == P_BONDS))
    {
        glLineWidth(LINK_WIDTH);
        for (j=0; j<ncircles; j++) if (particle[j].active) draw_one_particle_links(particle[j]);
    }
                
    /* fill triangles between particles */
    if (FILL_TRIANGLES) draw_triangles(particle, plot);
    
    /* draw collision discs */
    if ((REACTION_DIFFUSION)&&(ncollisions > 0)) draw_collisions(collisions, ncollisions);
    
    /* determine particle color and size */
    for (j=0; j<ncircles; j++) if (particle[j].active)
    {
        compute_particle_colors(particle[j], plot, rgb, rgbx, rgby, &radius, &width);
       
        switch (particle[j].interaction) {
            case (I_LJ_DIRECTIONAL): 
            {   
                nsides = 4;
                break;
            }
            case (I_LJ_PENTA): 
            {
                nsides = 5;
                break;
            }
            case (I_LJ_QUADRUPOLE):
            {
                nsides = 4;
                break;
            }
            case (I_LJ_WATER):
            {
                nsides = NSEG;
                radius *= 0.75;
                break;
            }
            default: nsides = NSEG;
        }

        angle = particle[j].angle + APOLY*DPI;
        
        /* in case of periodic b.c., draw translates of particles */
        if ((PERIODIC_BC)&&(plot != P_NOPARTICLE))
        {
            x1 = particle[j].xc;
            y1 = particle[j].yc;
            
            for (i=-1; i<2; i++)
                for (k=-1; k<2; k++)
                    draw_one_particle(particle[j], x1 + (double)i*(BCXMAX - BCXMIN), y1 + (double)k*(BCYMAX - BCYMIN), radius, angle, nsides, width, rgb);
        }
        else if ((BOUNDARY_COND == BC_KLEIN)&&(plot != P_NOPARTICLE))
        {
            x1 = particle[j].xc;
            y1 = particle[j].yc;
            
            for (i=-2; i<3; i++)
            {
                if (vabs(i) == 1) sign = -1.0;
                else sign = 1.0;
                angle1 = angle*sign;
                for (k=-1; k<2; k++)
                    draw_one_particle(particle[j], x1 + (double)i*(BCXMAX - BCXMIN), sign*(y1 + (double)k*(BCYMAX - BCYMIN)), 
                                      radius, angle1, nsides, width, rgb);
            }
        }
        else if ((BOUNDARY_COND == BC_BOY)&&(plot != P_NOPARTICLE))
        {
            x1 = particle[j].xc;
            y1 = particle[j].yc;
            
            for (i=-1; i<2; i++) for (k=-1; k<2; k++)
            {
                if (vabs(i) == 1) sign = -1.0;
                else sign = 1.0;
                if (vabs(k) == 1) signy = -1.0;
                else signy = 1.0;
                if (signy == 1.0) angle1 = angle*sign;
                else angle1 = PI - angle;
                if (sign == -1.0) draw_one_particle(particle[j], signy*(x1 + (double)i*(BCXMAX - BCXMIN)), 
                    sign*(y1 + (double)k*(BCYMAX - BCYMIN)), radius, angle1, nsides, width, rgbx);
                else if (signy == -1.0) draw_one_particle(particle[j], signy*(x1 + (double)i*(BCXMAX - BCXMIN)), 
                    sign*(y1 + (double)k*(BCYMAX - BCYMIN)), radius, angle1, nsides, width, rgby);
                else draw_one_particle(particle[j], signy*(x1 + (double)i*(BCXMAX - BCXMIN)), 
                    sign*(y1 + (double)k*(BCYMAX - BCYMIN)), radius, angle1, nsides, width, rgb);
            }
        }
        else if ((BOUNDARY_COND == BC_GENUS_TWO)&&(plot != P_NOPARTICLE))
        {
            x1 = particle[j].xc;
            y1 = particle[j].yc;
            
            if (x1 < 0.0) periody = BCYMAX - BCYMIN;
            else periody = 0.5*(BCYMAX - BCYMIN);
            
            if (y1 < 0.0) periodx = BCXMAX - BCXMIN;
            else periodx = 0.5*(BCXMAX - BCXMIN);
            
            if ((x1 < 0.0)&&(y1 < 0.0))
                for (i=-1; i<2; i++)
                    for (k=-1; k<2; k++)
                    {
                        x = x1 + (double)i*periodx;
                        y = y1 + (double)k*periody;
                        draw_one_particle(particle[j], x, y, radius, angle, nsides, width, rgb);
                    }
            else if ((x1 < 0.0)&&(y1 >= 0.0))
                for (i=-1; i<2; i++)
                    for (k=-1; k<2; k++)
                    {
                        x = x1 + (double)i*periodx;
                        y = y1 + (double)k*periody;
                        if (x < 1.2*particle[j].radius)
                            draw_one_particle(particle[j], x, y, radius, angle, nsides, width, rgb);
                    }
            else if ((x1 >= 0.0)&&(y1 < 0.0))
                for (i=-1; i<2; i++)
                    for (k=-1; k<2; k++)
                    {
                        x = x1 + (double)i*periodx;
                        y = y1 + (double)k*periody;
                         if (y < 1.2*particle[j].radius)
                            draw_one_particle(particle[j], x, y, radius, angle, nsides, width, rgb);
                    }
        }
        else if (plot != P_NOPARTICLE)
            draw_one_particle(particle[j], particle[j].xc, particle[j].yc, radius, angle, nsides, width, rgb);
                
        
    }
    
//     /* draw spin vectors */
    if ((DRAW_SPIN)||(DRAW_SPIN_B))
    {
        glLineWidth(width);
        for (j=0; j<ncircles; j++) 
            if ((particle[j].active)&&(((DRAW_SPIN)&&(particle[j].type == 0))||((DRAW_SPIN_B)&&(particle[j].type == 1))))
        {
//             x1 = particle[j].xc - 2.0*MU*cos(particle[j].angle);
//             y1 = particle[j].yc - 2.0*MU*sin(particle[j].angle);
            x1 = particle[j].xc;
//             if (CENTER_VIEW_ON_OBSTACLE) x1 -= xshift;
            y1 = particle[j].yc;
            x2 = particle[j].xc + 2.0*MU*cos(particle[j].angle);
//             if (CENTER_VIEW_ON_OBSTACLE) x2 -= xshift;
            y2 = particle[j].yc + 2.0*MU*sin(particle[j].angle);
            draw_line(x1, y1, x2, y2);
        }
    }
}


void draw_container(double xmin, double xmax, t_obstacle obstacle[NMAXOBSTACLES], t_segment segment[NMAXSEGMENTS], int wall)
/* draw the container, for certain boundary conditions */
{
    int i, j;
    double rgb[3], x, phi, angle, dx, dy, ybin, x1, x2, h;
    char message[100];
    
    /* draw fixed obstacles */
    if (ADD_FIXED_OBSTACLES)
    {
        glLineWidth(CONTAINER_WIDTH);
        hsl_to_rgb(300.0, 0.1, 0.5, rgb);
        for (i = 0; i < nobstacles; i++)
            draw_colored_circle_precomp(obstacle[i].xc - xtrack, obstacle[i].yc - ytrack, obstacle[i].radius, rgb);
        glColor3f(1.0, 1.0, 1.0);
        for (i = 0; i < nobstacles; i++)
            draw_circle_precomp(obstacle[i].xc - xtrack, obstacle[i].yc - ytrack, obstacle[i].radius);
    }
    if (ADD_FIXED_SEGMENTS)
    {
        glLineWidth(CONTAINER_WIDTH);
        glColor3f(1.0, 1.0, 1.0);
        for (i = 0; i < nsegments; i++) if (segment[i].active)
        {
            if (COLOR_SEG_GROUPS) set_segment_group_color(segment[i].group, 1.0, rgb);
            draw_line(segment[i].x1 - xtrack, segment[i].y1 - ytrack, segment[i].x2 - xtrack, segment[i].y2 - ytrack);
        }
    }

    switch (BOUNDARY_COND) {
        case (BC_SCREEN): 
        {
            /* do nothing */
            break;
        }
        case (BC_RECTANGLE): 
        {
            glColor3f(1.0, 1.0, 1.0);
            glLineWidth(CONTAINER_WIDTH);
            
            draw_line(BCXMIN, BCYMIN, BCXMAX, BCYMIN);
            draw_line(BCXMIN, BCYMAX, BCXMAX, BCYMAX);
            
            if (!SYMMETRIC_DECREASE) draw_line(BCXMAX, BCYMIN,  BCXMAX, BCYMAX);

            draw_line(xmin, BCYMIN, xmin, BCYMAX);
//             draw_line(XMIN, 0.5*(BCYMIN + BCYMAX), xmin, 0.5*(BCYMIN + BCYMAX));
            
            if (SYMMETRIC_DECREASE)
            {
                draw_line(xmax, BCYMIN, xmax, BCYMAX);
                draw_line(XMAX, 0.5*(BCYMIN + BCYMAX), xmax, 0.5*(BCYMIN + BCYMAX));
            }
            
            break;
        }
        case (BC_CIRCLE):
        {
            glLineWidth(CONTAINER_WIDTH);
            hsl_to_rgb(300.0, 0.1, 0.5, rgb);
            for (i=-1; i<2; i++)
            {
                if (CENTER_VIEW_ON_OBSTACLE) x = 0.0;
                else x = xmin + (double)i*(OBSXMAX - OBSXMIN);
                
                draw_colored_circle_precomp(x, 0.0, OBSTACLE_RADIUS, rgb);
                glColor3f(1.0, 1.0, 1.0);
                draw_circle_precomp(x, 0.0, OBSTACLE_RADIUS);
                
                glColor3f(0.0, 0.0, 0.0);
                sprintf(message, "Mach %.3f", xspeed/20.0);
//                 sprintf(message, "Speed %.2f", xspeed);
                write_text(x-0.17, -0.025, message); 
            }
            break;
        }
        case (BC_PERIODIC_CIRCLE):
        {
            glLineWidth(CONTAINER_WIDTH);
            hsl_to_rgb(300.0, 0.1, 0.5, rgb);
            for (i=-1; i<2; i++)
            {
                if (CENTER_VIEW_ON_OBSTACLE) x = 0.0;
                else x = xmin + (double)i*(OBSXMAX - OBSXMIN);
                
                draw_colored_circle_precomp(x, 0.0, OBSTACLE_RADIUS, rgb);
                glColor3f(1.0, 1.0, 1.0);
                draw_circle_precomp(x, 0.0, OBSTACLE_RADIUS);
                
                glColor3f(0.0, 0.0, 0.0);
                sprintf(message, "Mach %.2f", xspeed/20.0);
//                 sprintf(message, "Speed %.2f", xspeed);
                write_text(x-0.17, -0.025, message); 
            }
            break;
        }
        case (BC_PERIODIC_TRIANGLE):
        {
            glLineWidth(CONTAINER_WIDTH);
            hsl_to_rgb(300.0, 0.1, 0.5, rgb);
            for (i=-1; i<2; i++)
            {
                if (CENTER_VIEW_ON_OBSTACLE) x = 0.0;
                else x = xmin + (double)i*(OBSXMAX - OBSXMIN);
                
                x1 = x + OBSTACLE_RADIUS;
                x2 = x - OBSTACLE_RADIUS;
                h = 2.0*OBSTACLE_RADIUS*tan(APOLY*PID);
                draw_colored_triangle(x1, 0.0, x2, h, x2, -h, rgb);
                glColor3f(1.0, 1.0, 1.0);
                draw_triangle(x1, 0.0, x2, h, x2, -h);
                
                glColor3f(0.0, 0.0, 0.0);
                sprintf(message, "Mach %.2f", xspeed*3.0/40.0);
                write_text(x-0.25, -0.025, message); 
            }
            break;
        }
        case (BC_PERIODIC_FUNNEL):
        {
            glLineWidth(CONTAINER_WIDTH);
            hsl_to_rgb(300.0, 0.1, 0.5, rgb);
            for (i=-1; i<2; i++)
            {
                if (CENTER_VIEW_ON_OBSTACLE) x = 0.0;
                else x = xmin + (double)i*(OBSXMAX - OBSXMIN);
                
                for (j=-1; j<2; j+=2)
                {
                    draw_colored_circle_precomp(x, (double)j*(FUNNEL_WIDTH + OBSTACLE_RADIUS), OBSTACLE_RADIUS, rgb);
                    glColor3f(1.0, 1.0, 1.0);
                    draw_circle_precomp(x, (double)j*(FUNNEL_WIDTH + OBSTACLE_RADIUS), OBSTACLE_RADIUS);
                }
                
                glColor3f(0.0, 0.0, 0.0);
                sprintf(message, "Mach %.2f", xspeed/20.0);
                write_text(x-0.17, 0.75, message); 
            }
            break;
        }
        case (BC_RECTANGLE_LID): 
        {
            glColor3f(1.0, 1.0, 1.0);
            glLineWidth(CONTAINER_WIDTH);
            
            draw_line(BCXMIN, BCYMIN, BCXMAX, BCYMIN);
            draw_line(BCXMIN, BCYMIN, BCXMIN, BCYMAX);
            draw_line(BCXMAX, BCYMIN, BCXMAX, BCYMAX);
            
            hsl_to_rgb(300.0, 0.1, 0.5, rgb);
            draw_colored_rectangle(BCXMIN + 0.05, ylid, BCXMAX - 0.05, ylid + LID_WIDTH, rgb);
            glColor3f(1.0, 1.0, 1.0);
            draw_rectangle(BCXMIN + 0.05, ylid, BCXMAX - 0.05, ylid + LID_WIDTH);
            
            break;
        }
        case (BC_RECTANGLE_WALL): 
        {
            glColor3f(1.0, 1.0, 1.0);
            glLineWidth(CONTAINER_WIDTH);
            
            draw_rectangle(BCXMIN, BCYMIN, BCXMAX, BCYMAX);
            
            draw_line(0.5*(BCXMIN+BCXMAX), BCYMAX, 0.5*(BCXMIN+BCXMAX), BCYMAX + 0.5*WALL_WIDTH);
            draw_line(0.5*(BCXMIN+BCXMAX), BCYMIN, 0.5*(BCXMIN+BCXMAX), BCYMIN - 0.5*WALL_WIDTH);
            
            if (wall)
            {
                hsl_to_rgb(300.0, 0.1, 0.5, rgb);
                draw_colored_rectangle(xwall - 0.5*WALL_WIDTH, BCYMIN + 0.025, xwall + 0.5*WALL_WIDTH, BCYMAX - 0.025, rgb);
                glColor3f(1.0, 1.0, 1.0);
                draw_rectangle(xwall - 0.5*WALL_WIDTH, BCYMIN + 0.025, xwall + 0.5*WALL_WIDTH, BCYMAX - 0.025);
            }
            
            break;
        }
        case (BC_EHRENFEST):
        {
            glLineWidth(CONTAINER_WIDTH);
            glColor3f(1.0, 1.0, 1.0);
            phi = asin(EHRENFEST_WIDTH/EHRENFEST_RADIUS);
            glBegin(GL_LINE_LOOP);
            for (i=0; i<=NSEG; i++)
            {
                angle = -PI + phi + (double)i*2.0*(PI - phi)/(double)NSEG;
                glVertex2d(1.0 + EHRENFEST_RADIUS*cos(angle), EHRENFEST_RADIUS*sin(angle));
            }
            for (i=0; i<=NSEG; i++)
            {
                angle = phi + (double)i*2.0*(PI - phi)/(double)NSEG;
                glVertex2d(-1.0 + EHRENFEST_RADIUS*cos(angle), EHRENFEST_RADIUS*sin(angle));
            }
            glEnd();
            break;
        }
        case (BC_SCREEN_BINS): 
        {
            glLineWidth(CONTAINER_WIDTH);
            glColor3f(1.0, 1.0, 1.0);
            dy = (YMAX - YMIN)/((double)NGRIDX + 3);
            dx = dy/cos(PI/6.0);
            ybin = 2.75*dy;

            for (i=-1; i<=NGRIDX; i++) 
            {
                x = ((double)i - 0.5*(double)NGRIDX + 0.5)*dx;
                draw_line(x, YMIN, x, YMIN + ybin);
            }
            break;
        }
        case (BC_GENUS_TWO):
        {
            hsl_to_rgb(300.0, 0.1, 0.5, rgb);
            draw_colored_rectangle(0.0, 0.0, BCXMAX, BCYMAX, rgb);
            break;
        }
        default: 
        {
            /* do nothing */
        }
    }
    
}

void print_omega(double angle, double angular_speed, double fx, double fy)
{
    char message[100];
    double rgb[3], y1, frac, absa;
    static double xleftbox, xlefttext, xrightbox, xrighttext, y = YMAX - 0.1, ymin = YMIN + 0.05;
    static int first = 1;
    
    if (first)
    {
        xrightbox = XMAX - 0.54;
        xrighttext = xrightbox - 0.48;
        first = 0;
    }
    
    y1 = y;
    
    if (PRINT_ANGLE)
    {
        erase_area_hsl(xrightbox, y + 0.025, 0.42, 0.05, 0.0, 0.9, 0.0);
        glColor3f(1.0, 1.0, 1.0);
        angle = angle*360.0/DPI;
        absa = vabs(angle);
        frac = absa - (double)((int)absa);
        sprintf(message, "Angle = %3d.%d degrees", (int)absa, (int)(10.0*frac));
        write_text(xrighttext + 0.1, y, message);
        y1 -= 0.1;
    }
    if (PRINT_OMEGA)
    {
        erase_area_hsl(xrightbox, y1 + 0.025, 0.42, 0.05, 0.0, 0.9, 0.0);
        glColor3f(1.0, 1.0, 1.0);
        sprintf(message, "Angular speed = %.4f", angular_speed);
        write_text(xrighttext + 0.1, y1, message);
        y1 -= 0.1;
    }
    if (PRINT_SEGMENTS_FORCE)
    {
        erase_area_hsl(xrightbox, y1 + 0.025, 0.42, 0.05, 0.0, 0.9, 0.0);
        glColor3f(1.0, 1.0, 1.0);
        sprintf(message, "Fx = %.4f", fx);
        write_text(xrighttext + 0.1, y1, message);
        
        y1 -= 0.1;
        erase_area_hsl(xrightbox, y1 + 0.025, 0.42, 0.05, 0.0, 0.9, 0.0);
        glColor3f(1.0, 1.0, 1.0);
        sprintf(message, "Fy = %.4f", fy);
        write_text(xrighttext + 0.1, y1, message);
    }
}

void compute_segments_force(t_lj_parameters *params, t_segment segment[NMAXSEGMENTS])
{
    int i;
    double fx = 0.0, fy = 0.0;
    
    for (i=0; i<nsegments; i++) if (segment[i].active)
    {
        fx += segment[i].fx;
        fy += segment[i].fy;
    }
    params->bdry_fx = fx;
    params->bdry_fy = fy;
        
}

void print_parameters(t_lj_parameters params, short int left, double pressure[N_PRESSURES], short int refresh)
{
    char message[100];
    int i, j, k;
    double density, hue, rgb[3], logratio, x, y, meanpress[N_PRESSURES], phi, sphi, dphi, pprint, mean_temp, lengthcontainer, boundary_force, fx, fy, r1, r2;
    static double xbox, xtext, xmid, xmidtext, xxbox, xxtext, pressures[N_P_AVERAGE], meanpressure = 0.0, maxpressure = 0.0, mean_fx, mean_fy;
    static double press[N_PRESSURES][N_P_AVERAGE], temp[N_T_AVERAGE], scale;
    static int first = 1, i_pressure, i_temp;
    
    if (first)
    {
//         scale = (XMAX - XMIN)/4.0;
        scale = (YMAX - YMIN)/2.5;
        if (left)
        {
            xbox = XMIN + 0.45*scale;
            xtext = XMIN + 0.12*scale;
            xxbox = XMAX - 0.39*scale;
            xxtext = XMAX - 0.73*scale;
       }
        else
        {
            xbox = XMAX - 0.41*scale;
            xtext = XMAX - 0.73*scale;
            xxbox = XMIN + 0.4*scale;
            xxtext = XMIN + 0.08*scale;
        }
        xmid = 0.5*(XMIN + XMAX) - 0.1*scale;
        xmidtext = xmid - 0.3*scale;
        for (i=0; i<N_P_AVERAGE; i++) pressures[i] = 0.0;
        if (RECORD_PRESSURES) for (j=0; j<N_PRESSURES; j++) 
        {
            meanpress[j] = 0.0;
            for (i=0; i<N_P_AVERAGE; i++) press[j][i] = 0.0;
        }
        i_pressure = 0;
        i_temp = 0;
        for (i=0; i<N_T_AVERAGE; i++) temp[i] = 0.0;
        mean_fx = 0.0;
        mean_fy = 0.0;
        r1 = 0.005;
        r2 = 1.0 - r1;
        
        first = 0;
    }
    
    lengthcontainer = params.xmaxcontainer - params.xmincontainer;
    boundary_force = params.fboundary/(double)(ncircles*NVID);
    
    /* table of pressures */
    pressures[i_pressure] = boundary_force/(lengthcontainer + INITYMAX - INITYMIN);
    if (RECORD_PRESSURES) 
    {
        for (j=0; j<N_PRESSURES; j++) press[j][i_pressure] = pressure[j];
    }
    i_pressure++;
    if (i_pressure == N_P_AVERAGE) i_pressure = 0;
    
    for (i=0; i<N_P_AVERAGE; i++) meanpressure += pressures[i];
    meanpressure = meanpressure/(double)N_P_AVERAGE;
    
    if (RECORD_PRESSURES) for (j=0; j<N_PRESSURES; j++) 
    {
        meanpress[j] = 0.0;
        for (i=0; i<N_P_AVERAGE; i++) meanpress[j] += press[j][i];
        meanpress[j] = meanpress[j]/(double)N_P_AVERAGE;
    }
    
    
//     if (RECORD_PRESSURES) 
//         for (j=0; j<N_PRESSURES; j++) meanpress[j] = 
//     
//     for (j=0; j<N_PRESSURES; j++) printf("Mean pressure[%i] = %.5lg\n", j, meanpress[j]);

    if (RECORD_PRESSURES) 
    {
        for (j=0; j<N_PRESSURES; j++) if (meanpress[j] > maxpressure) maxpressure = meanpress[j];
        printf("Max pressure = %.5lg\n\n", maxpressure);
    }
    
    y = YMAX - 0.1*scale;
    if ((INCREASE_BETA)||(PRINT_TEMPERATURE))  /* print temperature */
    {
        logratio = log(params.beta/BETA)/log(0.5*BETA_FACTOR);
        if (logratio > 1.0) logratio = 1.0;
        else if (logratio < 0.0) logratio = 0.0;
        if (BETA_FACTOR > 1.0) hue = PARTICLE_HUE_MAX - (PARTICLE_HUE_MAX - PARTICLE_HUE_MIN)*logratio;
        else hue = PARTICLE_HUE_MIN - (PARTICLE_HUE_MIN - PARTICLE_HUE_MAX)*logratio;
        if (PRINT_LEFT) erase_area_hsl_turbo(xbox, y + 0.025*scale, 0.5*scale, 0.05*scale, hue, 0.9, 0.5);
        else erase_area_hsl_turbo(xmid + 0.1, y + 0.025*scale, 0.45*scale, 0.05*scale, hue, 0.9, 0.5);
        if ((hue < 90)||(hue > 270)) glColor3f(1.0, 1.0, 1.0);
        else glColor3f(0.0, 0.0, 0.0);
        sprintf(message, "Temperature %.2f", 1.0/params.beta);
        if (PRINT_LEFT) write_text(xtext, y, message);
        else write_text(xmidtext, y, message);
//         y -= 0.1;
        
//         erase_area_hsl(xxbox, y + 0.025, 0.37, 0.05, 0.0, 0.9, 0.0);
//         glColor3f(1.0, 1.0, 1.0);
//         sprintf(message, "Pressure %.3f", meanpressure);
//         write_text(xxtext, y, message);
    }
    if (DECREASE_CONTAINER_SIZE)  /* print density */
    {
        density = (double)ncircles/((lengthcontainer)*(INITYMAX - INITYMIN));
        erase_area_hsl(xbox, y + 0.025*scale, 0.37*scale, 0.05*scale, 0.0, 0.9, 0.0);
        glColor3f(1.0, 1.0, 1.0);
        sprintf(message, "Density %.3f", density);
        write_text(xtext, y, message);
        
        erase_area_hsl(xmid, y + 0.025*scale, 0.37*scale, 0.05*scale, 0.0, 0.9, 0.0);
        glColor3f(1.0, 1.0, 1.0);
        sprintf(message, "Temperature %.2f", params.mean_energy);
        write_text(xmidtext, y, message);

        erase_area_hsl(xxbox, y + 0.025*scale, 0.37*scale, 0.05*scale, 0.0, 0.9, 0.0);
        glColor3f(1.0, 1.0, 1.0);
        sprintf(message, "Pressure %.3f", meanpressure);
        write_text(xxtext, y, message);

    }   
    else if (INCREASE_KREPEL)  /* print force constant */
    {
        erase_area_hsl(xbox, y + 0.025*scale, 0.22*scale, 0.05*scale, 0.0, 0.9, 0.0);
        glColor3f(1.0, 1.0, 1.0);
        sprintf(message, "Force %.0f", params.krepel);
        write_text(xtext + 0.28, y, message);
    }   
    
    if (RECORD_PRESSURES) 
    {
        y = FUNNEL_WIDTH + OBSTACLE_RADIUS;
        for (i=0; i<N_PRESSURES; i++)
        {
            phi = DPI*(double)i/(double)N_PRESSURES;
            hue = PARTICLE_HUE_MIN + (PARTICLE_HUE_MAX - PARTICLE_HUE_MIN)*meanpress[i]/MAX_PRESSURE;
            if (hue > PARTICLE_HUE_MIN) hue = PARTICLE_HUE_MIN;
            if (hue < PARTICLE_HUE_MAX) hue = PARTICLE_HUE_MAX;
            hsl_to_rgb_turbo(hue, 0.9, 0.5, rgb);
        
            dphi = DPI/(double)N_PRESSURES;
//         x = 0.95*OBSTACLE_RADIUS*cos(phi);
            sphi = sin(phi);
            if (sphi < 0.0) draw_colored_sector(0.0, y, 0.95*OBSTACLE_RADIUS, OBSTACLE_RADIUS, phi, phi + dphi, rgb, 10);
            else draw_colored_sector(0.0, -y, 0.95*OBSTACLE_RADIUS, OBSTACLE_RADIUS, phi, phi + dphi, rgb, 10);
        }
    
        glColor3f(1.0, 1.0, 1.0);
        for (i=-1; i<2; i++) 
        {
            k = N_PRESSURES/4 + i*N_PRESSURES/9;
            phi = DPI*(double)k/(double)N_PRESSURES;
            pprint = 0.0;
            for (j=-2; j<3; j++) pprint += meanpress[k + j];
            sprintf(message, "p = %.0f", pprint*200.0/MAX_PRESSURE);
            write_text(0.85*OBSTACLE_RADIUS*cos(phi) - 0.1, -y + 0.85*OBSTACLE_RADIUS*sin(phi), message);
        }
    }
    
    if ((PARTIAL_THERMO_COUPLING)&&(!INCREASE_BETA)&&(!EXOTHERMIC))
    {
        printf("Temperature %i in average: %.3lg\n", i_temp, params.mean_energy);
        temp[i_temp] = params.mean_energy;
        i_temp++;
        if (i_temp >= N_T_AVERAGE) i_temp = 0;
        
        mean_temp = 0.0;
        for (i=0; i<N_T_AVERAGE; i++) mean_temp += temp[i];
        mean_temp = mean_temp/N_T_AVERAGE;
        
        hue = PARTICLE_HUE_MIN + 0.5*(PARTICLE_HUE_MAX - PARTICLE_HUE_MIN)*mean_temp/PARTICLE_EMAX;
        if (hue < PARTICLE_HUE_MAX) hue = PARTICLE_HUE_MAX;
        erase_area_hsl_turbo(xbox, y + 0.025*scale, 0.37*scale, 0.05*scale, hue, 0.9, 0.5);
        if ((hue < 90)||(hue > 270)) glColor3f(1.0, 1.0, 1.0);
        else glColor3f(0.0, 0.0, 0.0);
        sprintf(message, "Temperature %.2f", mean_temp);
        write_text(xtext, y, message);
    }
    
    if (INCREASE_GRAVITY)
    {
        erase_area_hsl(xmid, y + 0.025*scale, 0.22*scale, 0.05*scale, 0.0, 0.9, 0.0);
        glColor3f(1.0, 1.0, 1.0);
        sprintf(message, "Gravity %.2f", params.gravity/GRAVITY);
        write_text(xmidtext + 0.1, y, message);
    } 
    
    if (CHANGE_RADIUS)
    {
        erase_area_hsl(xmid, y + 0.025*scale, 0.3*scale, 0.05*scale, 0.0, 0.9, 0.0);
        glColor3f(1.0, 1.0, 1.0);
        sprintf(message, "Radius %.4f", params.radius);
        write_text(xmidtext + 0.05, y, message);
    }  
    
    if (PRINT_SEGMENTS_FORCE)
    {
        glColor3f(1.0, 1.0, 1.0);
        if (refresh)
        {
            fx = 0.01*params.bdry_fx/(double)params.nactive;
            fy = 0.01*params.bdry_fy/(double)params.nactive;
            /* average boundary force */
            mean_fx = r2*mean_fx + r1*fx;
            mean_fy = r2*mean_fy + r1*fy;
        }
        
        if ((PRINT_ANGLE)||(PRINT_OMEGA)) 
            draw_arrow(0.0, 0.0, FORCE_FACTOR*mean_fx, FORCE_FACTOR*mean_fy, 15.0, 0.1);
        else
        {
            erase_area_hsl(xmid, 0.0, 0.2*scale, 0.05*scale, 0.0, 0.9, 0.0);
            glColor3f(1.0, 1.0, 1.0);
            sprintf(message, "Fy = %.2f", mean_fy);
            write_text(xmidtext + 0.15*scale, -0.02*scale, message);
            if (mean_fx*mean_fx + mean_fy*mean_fy > 5.0*FORCE_FACTOR) 
                draw_arrow(0.0, 0.0, FORCE_FACTOR*mean_fx, FORCE_FACTOR*mean_fy, 15.0, 0.1);
        }
    }
    if ((PRINT_ANGLE)||(PRINT_OMEGA)) print_omega(params.angle, params.omega, mean_fx, mean_fy); 
    
}


void print_ehrenfest_parameters(t_particle particle[NMAXCIRCLES], double pleft, double pright)
{
    char message[100];
    int i, j, nleft1 = 0, nleft2 = 0, nright1 = 0, nright2 = 0;
    double density, hue, rgb[3], logratio, y, shiftx = 0.3, xmidplus, xmidminus;
    static double xleftbox, xlefttext, xmidbox, xmidtext, xrightbox, xrighttext, pressures[500][2], meanpressure[2];
    static int first = 1, i_pressure, naverage = 500, n_pressure;
    
    if (first)
    {
        xleftbox = -0.85;
        xlefttext = xleftbox - 0.5;
        xrightbox = 1.0;
        xrighttext = xrightbox - 0.45;
//         xmid = 0.5*(XMIN + XMAX) - 0.1;
//         xmidtext = xmid - 0.24;
        
        meanpressure[0] = 0.0;
        meanpressure[1] = 0.0;
        for (i=0; i<naverage; i++) 
        {
            pressures[i][0] = 0.0;
            pressures[i][1] = 0.0;
        }
        i_pressure = 0;
        n_pressure = 0;
        first = 0;
    }
    
    if (BOUNDARY_COND == BC_EHRENFEST)
    {
        xmidplus = 1.0 - EHRENFEST_RADIUS;
        xmidminus = -1.0 + EHRENFEST_RADIUS;
    }
    else 
    {
        xmidplus = xwall;
        xmidminus = xwall;
    }
    
     /* table of pressures */
    pressures[i_pressure][0] = pleft;
    pressures[i_pressure][1] = pright;
    i_pressure++;
    if (i_pressure == naverage) i_pressure = 0;
    if (n_pressure < naverage - 1) n_pressure++;
    
    for (i=0; i<n_pressure; i++) 
        for (j=0; j<2; j++)
            meanpressure[j] += pressures[i][j];
    for (j=0; j<2; j++) meanpressure[j] = meanpressure[j]/(double)n_pressure;
    
    
    for (i = 0; i < ncircles; i++) if (particle[i].active)
    {
        if (particle[i].xc < xmidminus)
        {
            if (particle[i].type == 0) nleft1++;
            else nleft2++;
        }
        else if (particle[i].xc > xmidplus) 
        {
            if (particle[i].type == 0) nright1++;
            else nright2++;
        }
    }
    
    y = YMIN + 0.05;
    
    erase_area_hsl(xleftbox - shiftx, y + 0.025, 0.22, 0.05, 0.0, 0.9, 0.0);
    hsl_to_rgb(HUE_TYPE0, 0.9, 0.5, rgb);
    glColor3f(rgb[0], rgb[1], rgb[2]);
    sprintf(message, "%i particles", nleft1);
    write_text(xlefttext + 0.28 - shiftx, y, message);
    
    erase_area_hsl(xleftbox + shiftx, y + 0.025, 0.22, 0.05, 0.0, 0.9, 0.0);
    hsl_to_rgb(HUE_TYPE1, 0.9, 0.5, rgb);
    glColor3f(rgb[0], rgb[1], rgb[2]);
    sprintf(message, "%i particles", nleft2);
    write_text(xlefttext + 0.28 + shiftx, y, message);
    
    erase_area_hsl(xrightbox - shiftx, y + 0.025, 0.22, 0.05, 0.0, 0.9, 0.0);
    hsl_to_rgb(HUE_TYPE0, 0.9, 0.5, rgb);
    glColor3f(rgb[0], rgb[1], rgb[2]);
    sprintf(message, "%i particles", nright1);
    write_text(xrighttext + 0.28 - shiftx, y, message);
    
    erase_area_hsl(xrightbox + shiftx, y + 0.025, 0.22, 0.05, 0.0, 0.9, 0.0);
    hsl_to_rgb(HUE_TYPE1, 0.9, 0.5, rgb);
    glColor3f(rgb[0], rgb[1], rgb[2]);
    sprintf(message, "%i particles", nright2);
    write_text(xrighttext + 0.28 + shiftx, y, message);
    y = YMAX - 0.1;
    
    erase_area_hsl(xleftbox - 0.1, y + 0.025, 0.22, 0.05, 0.0, 0.9, 0.0);
    hsl_to_rgb_turbo(HUE_TYPE1, 0.9, 0.5, rgb);
    glColor3f(rgb[0], rgb[1], rgb[2]);
    sprintf(message, "Pressure %.2f", 0.001*meanpressure[0]/(double)ncircles);
    write_text(xlefttext + 0.25, y, message);
    
    erase_area_hsl(xrightbox - 0.1, y + 0.025, 0.22, 0.05, 0.0, 0.9, 0.0);
    hsl_to_rgb_turbo(HUE_TYPE0, 0.9, 0.5, rgb);
    glColor3f(rgb[0], rgb[1], rgb[2]);
    sprintf(message, "Pressure %.2f", 0.001*meanpressure[1]/(double)ncircles);
    write_text(xrighttext + 0.2, y, message);

}

void count_particle_number(t_particle *particle, int *particle_numbers, int time)
{
    int type, i, n;
    
    n = time*(RD_TYPES+1);
    
    for (type = 0; type < RD_TYPES+1; type++)
        particle_numbers[n + type] = 0;
    
    for (i=0; i<ncircles; i++)
        if (particle[i].active) 
            particle_numbers[n + particle[i].type]++;
}

void print_particle_number(int npart)
{
    char message[100];
    double y = YMAX - 0.1;
    static double xleftbox, xlefttext;
    static int first = 1;
    
    if (first)
    {
        xleftbox = XMIN + 0.5;
        xlefttext = xleftbox - 0.5;
        first = 0;
    }
    
    erase_area_hsl(xleftbox, y + 0.025, 0.22, 0.05, 0.0, 0.9, 0.0);
    glColor3f(1.0, 1.0, 1.0);
    if (npart > 1) sprintf(message, "%i particles", npart);
    else sprintf(message, "%i particle", npart);
    write_text(xlefttext + 0.28, y, message);
}

void print_particle_types_number(t_particle *particle, int ntypes)
{
    int i, ntype[10];
    char message[100];
    double y = YMAX - 0.1, rgb[3];
    static double xleftbox, xlefttext;
    static int first = 1;
    
    if (first)
    {
        xleftbox = XMAX - 0.5;
        xlefttext = xleftbox - 0.45;
        first = 0;
    }
    
    for (i=0; i<10; i++) ntype[i] = 0;
    
    for (i=0; i<ncircles; i++)
        if (particle[i].active) ntype[particle[i].type]++;
    
    for (i=1; i<ntypes+1; i++)
    {
        erase_area_hsl(xleftbox, y + 0.025, 0.22, 0.05, 0.0, 0.9, 0.0);
        set_type_color(i, 0.5, rgb);
        glColor3f(rgb[0], rgb[1], rgb[2]);
        sprintf(message, "%i particles", ntype[i]);
        write_text(xlefttext + 0.28, y, message);
        y -= 0.12;
    }
}

void print_entropy(double entropy[2])
{
    char message[100];
    double rgb[3];
    static double xleftbox, xlefttext, xrightbox, xrighttext, y = YMAX - 0.1, ymin = YMIN + 0.05;
    static int first = 1;
    
    if (first)
    {
        xleftbox = XMIN + 0.4;
        xlefttext = xleftbox - 0.55;
        xrightbox = XMAX - 0.39;
        xrighttext = xrightbox - 0.55;
       first = 0;
    }
    
    if (POSITION_Y_DEPENDENCE)
    {
        erase_area_hsl(xrightbox, ymin + 0.025, 0.35, 0.05, 0.0, 0.9, 0.0);
        hsl_to_rgb_turbo(HUE_TYPE1, 0.9, 0.5, rgb);
        glColor3f(rgb[0], rgb[1], rgb[2]);
        sprintf(message, "Entropy = %.4f", entropy[1]);
        write_text(xrighttext + 0.28, ymin, message);
    }
    else
    {
        erase_area_hsl(xleftbox, y + 0.025, 0.35, 0.05, 0.0, 0.9, 0.0);
        hsl_to_rgb_turbo(HUE_TYPE1, 0.9, 0.5, rgb);
        glColor3f(rgb[0], rgb[1], rgb[2]);
        sprintf(message, "Entropy = %.4f", entropy[1]);
        write_text(xlefttext + 0.28, y, message);
    }

    erase_area_hsl(xrightbox, y + 0.025, 0.35, 0.05, 0.0, 0.9, 0.0);
    hsl_to_rgb_turbo(HUE_TYPE0, 0.9, 0.5, rgb);
    glColor3f(rgb[0], rgb[1], rgb[2]);
    sprintf(message, "Entropy = %.4f", entropy[0]);
    write_text(xrighttext + 0.28, y, message);

}

void print_segments_speeds(double vx[2], double vy[2])
{
    char message[100];
    double rgb[3], y;
    static double xleftbox, xlefttext, xrightbox, xrighttext, ymin = YMIN + 0.05, scale;
    static int first = 1;
    
    if (first)
    {
        scale = (XMAX - XMIN)/4.0;
        xleftbox = XMIN + 0.3*scale;
        xlefttext = xleftbox - 0.22*scale;
        xrightbox = XMAX - 0.39*scale;
        xrighttext = xrightbox - 0.3*scale;
        first = 0;
    }
    
    y = YMAX - 0.1*scale;
    erase_area_hsl(xrightbox, y + 0.025*scale, 0.25*scale, 0.05*scale, 0.0, 0.9, 0.0);
    glColor3f(1.0, 1.0, 1.0);
    sprintf(message, "Vx = %.2f", vx[0]);
    write_text(xrighttext + 0.1, y, message);

    y -= 0.1*scale;
    erase_area_hsl(xrightbox, y + 0.025*scale, 0.25*scale, 0.05*scale, 0.0, 0.9, 0.0);
    glColor3f(1.0, 1.0, 1.0);
    sprintf(message, "Vy = %.2f", vy[0]);
    write_text(xrighttext + 0.1, y, message);
    
    if (TWO_OBSTACLES)
    {
        y = YMAX - 0.1*scale;
        erase_area_hsl(xleftbox, y + 0.025*scale, 0.2*scale, 0.05*scale, 0.0, 0.9, 0.0);
        glColor3f(1.0, 1.0, 1.0);
        sprintf(message, "Vx = %.2f", vx[1]);
        write_text(xlefttext + 0.1, y, message);

        y -= 0.1*scale;
        erase_area_hsl(xleftbox, y + 0.025*scale, 0.2*scale, 0.05*scale, 0.0, 0.9, 0.0);
        glColor3f(1.0, 1.0, 1.0);
        sprintf(message, "Vy = %.2f", vy[1]);
        write_text(xlefttext + 0.1, y, message);
    }
}

void print_segment_group_speeds(t_group_segments *segment_group)
{
    char message[100];
    double rgb[3], y, av_vx = 0.0, av_vy = 0.0, av_omega = 0.0, xbox, xtext;
    int i, group, groupshift;
    static double xleftbox, xlefttext, xrightbox, xrighttext, ymin = YMIN + 0.05, scale;
    static double vx[100*NMAXGROUPS], vy[100*NMAXGROUPS], omega[100*NMAXGROUPS], inv_t_window, speed_ratio;
    static int first = 1, position[NMAXGROUPS], t_window = 50;
    
    if (first)
    {
        scale = (XMAX - XMIN)/4.0;
        xleftbox = XMIN + 0.3*scale;
        xlefttext = xleftbox - 0.22*scale;
        xrightbox = XMAX - 0.39*scale;
        xrighttext = xrightbox - 0.3*scale;
        speed_ratio = (double)(25*NVID)*DT_PARTICLE;
        inv_t_window = 1.0/(double)t_window;
        for (group=1; group<NMAXGROUPS; group++) position[group] = 0;
        
        first = 0;
    }
    
    /* compute time averages */
    for (group = 1; group < ngroups; group++)
    {
        groupshift = (group-1)*100;
        vx[groupshift + position[group]] = segment_group[group].vx*speed_ratio;
        vy[groupshift + position[group]] = segment_group[group].vy*speed_ratio;
        omega[groupshift + position[group]] = segment_group[group].omega*speed_ratio;
        position[group]++;
        if (position[group] >=t_window) position[group] = 0;
        for (i=0; i<t_window; i++)
        {
            av_vx += vx[groupshift + i];
            av_vy += vy[groupshift + i];
            av_omega += omega[groupshift + i];
        }
        av_vx *= inv_t_window;
        av_vy *= inv_t_window;
        av_omega *= inv_t_window;
        
        if (ngroups > 2)
        {
            xbox = xleftbox + (xrightbox - xleftbox)*(group-1)/(ngroups-2);
            xtext = xlefttext + (xrighttext - xlefttext)*(group-1)/(ngroups-2);
        }
        else
        {
            xbox = xrightbox - 0.2;
            xtext = xrighttext - 0.2;            
        }
        
//         printf("xbox = %.2f, xtext = %.2f, av_vy = %.2f\n", xbox, xtext, av_vy);
        
        y = YMAX - 0.1*scale;
        erase_area_hsl(xbox, y + 0.025*scale, 0.25*scale, 0.05*scale, 0.0, 0.9, 0.0);
        set_segment_group_color(group, 1.0, rgb);
        sprintf(message, "Vx = %.4f", av_vx);
        write_text(xtext + 0.1, y, message);

        y -= 0.1*scale;
        erase_area_hsl(xbox, y + 0.025*scale, 0.25*scale, 0.05*scale, 0.0, 0.9, 0.0);
        set_segment_group_color(group, 1.0, rgb);
        sprintf(message, "Vy = %.4f", av_vy);
        write_text(xtext + 0.1, y, message);

        y -= 0.1*scale;
        erase_area_hsl(xbox, y + 0.025*scale, 0.3*scale, 0.05*scale, 0.0, 0.9, 0.0);
        set_segment_group_color(group, 1.0, rgb);
        sprintf(message, "V = %.4f", module2(av_vx, av_vy));
        write_text(xtext + 0.1, y, message);
        
        y -= 0.1*scale;
        erase_area_hsl(xbox, y + 0.025*scale, 0.3*scale, 0.05*scale, 0.0, 0.9, 0.0);
        set_segment_group_color(group, 1.0, rgb);
        sprintf(message, "Omega = %.4f", av_omega);
        write_text(xtext + 0.1, y, message);
    }
}

void print_particles_speeds(t_particle particle[NMAXCIRCLES])
{
    char message[100];
    double y = YMAX - 0.1, vx = 0.0, vy = 0.0;
    int i, nactive = 0;
    static double xleftbox, xlefttext, xrightbox, xrighttext, ymin = YMIN + 0.05;
    static int first = 1;
    
    if (first)
    {
        xrightbox = XMAX - 0.39;
        xrighttext = xrightbox - 0.55;
        first = 0;
    }
    
    for (i=0; i<NMAXCIRCLES; i++) if (particle[i].active)
    {
        nactive++;
        vx += particle[i].vx;
        vy += particle[i].vy;
    }
    
    vx = vx/(double)nactive;
    vy = vy/(double)nactive;

    erase_area_hsl(xrightbox, y + 0.025, 0.35, 0.05, 0.0, 0.9, 0.0);
    glColor3f(1.0, 1.0, 1.0);
    sprintf(message, "Average vx = %.4f", vx);
    write_text(xrighttext + 0.1, y, message);

    y -= 0.1;
    
    erase_area_hsl(xrightbox, y + 0.025, 0.35, 0.05, 0.0, 0.9, 0.0);
    glColor3f(1.0, 1.0, 1.0);
    sprintf(message, "Average vy = %.4f", vy);
    write_text(xrighttext + 0.1, y, message);
}

double compute_boundary_force(int j, t_particle particle[NMAXCIRCLES], t_obstacle obstacle[NMAXOBSTACLES], 
                              t_segment segment[NMAXSEGMENTS], 
                              double xleft, double xright, double *pleft, double *pright, double pressure[N_PRESSURES], int wall)
{
    int i, k;
    double xmin, xmax, ymin, ymax, padding, r, rp, r2, cphi, sphi, f, fperp = 0.0, x, y, xtube, distance, dx, dy, width, ybin, angle, x1, x2, h, ytop, norm, dleft, dplus, dminus, tmp_pleft = 0.0, tmp_pright = 0.0, proj, pscal, pvect, pvmin;
    
    /* compute force from fixed circular obstacles */
    if (ADD_FIXED_OBSTACLES) for (i=0; i<nobstacles; i++)
    {
        x = particle[j].xc - obstacle[i].xc;
        y = particle[j].yc - obstacle[i].yc;
        distance = module2(x, y);
        if (distance < 1.0e-7) distance = 1.0e-7;
        cphi = x/distance;
        sphi = y/distance;
        r2 = obstacle[i].radius + particle[j].radius;
        
        if (distance < r2)
        {
            f = KSPRING_OBSTACLE*(r2 - distance);
            particle[j].fx += f*cphi;
            particle[j].fy += f*sphi;
        }
    }
    /* compute force from fixed linear obstacles */
    particle[j].close_to_boundary = 0;
    if (ADD_FIXED_SEGMENTS) for (i=0; i<nsegments; i++) if (segment[i].active)
    {
        x = particle[j].xc;
        y = particle[j].yc;
        proj = (segment[i].ny*(x - segment[i].x1) - segment[i].nx*(y - segment[i].y1))/segment[i].length;
        if ((proj > 0.0)&&(proj < 1.0))
        {
//             distance = segment[i].nx*x + segment[i].ny*y - segment[i].c;
            distance = segment[i].nx*x + segment[i].ny*y - segment[i].c;
            r = 1.5*particle[j].radius;
            if (vabs(distance) < r)
            {
                particle[j].close_to_boundary = 1;
                
                f = KSPRING_OBSTACLE*(r - distance);
                particle[j].fx += f*segment[i].nx;
                particle[j].fy += f*segment[i].ny;
                
                
            
                if ((MOVE_BOUNDARY)||(MOVE_SEGMENT_GROUPS)||(PRINT_SEGMENTS_FORCE))
                {
                    segment[i].fx -= f*segment[i].nx;
                    segment[i].fy -= f*segment[i].ny;
                    segment[i].torque -= (x - segment[i].xc)*f*segment[i].ny - (y - segment[i].yc)*f*segment[i].nx;
//                     printf("Segment %i: f = (%.3lg, %.3lg)\n", i, segment[i].fx, segment[i].fy);
                }
            }
            if ((VICSEK_INT)&&(vabs(distance) < 1.5*r))
            {
                pvmin = 2.0;
                pvect = cos(particle[j].angle)*segment[i].ny - sin(particle[j].angle)*segment[i].nx;
                if ((pvect > 0.0)&&(pvect < pvmin)) pvect = pvmin;
                else if ((pvect < 0.0)&&(pvect > -pvmin)) pvect = -pvmin;
//                 particle[j].torque += KTORQUE_BOUNDARY*pvect;
                particle[j].torque += KTORQUE_BOUNDARY*pvect*(1.5*r - vabs(distance));
            }
        }
        
        /* compute force from concave corners */
        if (segment[i].concave)
        {
            distance = module2(x - segment[i].x1, y - segment[i].y1);
            angle = argument(x - segment[i].x1, y - segment[i].y1);
            if (angle < segment[i].angle1) angle += DPI;
            
            /* added 24/9/22 */
            else if (angle > segment[i].angle2) angle -= DPI;
            
            r = 1.5*particle[j].radius;
            
            if ((distance < r)&&(angle > segment[i].angle1)&&(angle < segment[i].angle2))
            {
                f = KSPRING_OBSTACLE*(r - distance);
                particle[j].fx += f*cos(angle);
                particle[j].fy += f*sin(angle);
                if ((MOVE_BOUNDARY)||(MOVE_SEGMENT_GROUPS))
                {
                    segment[i].fx -= f*cos(angle);
                    segment[i].fy -= f*sin(angle);
                    segment[i].torque -= (x - segment[i].xc)*f*sin(angle) - (y - segment[i].yc)*f*cos(angle);
                }
            }
        }
    }

    switch(BOUNDARY_COND){
        case (BC_SCREEN):
        {
            /* add harmonic force outside screen */
            if (particle[j].xc > XMAX) particle[j].fx -= KSPRING_BOUNDARY*(particle[j].xc - XMAX);
            else if (particle[j].xc < XMIN) particle[j].fx += KSPRING_BOUNDARY*(XMIN - particle[j].xc);
            if (particle[j].yc > YMAX) particle[j].fy -= KSPRING_BOUNDARY*(particle[j].yc - YMAX);
            else if (particle[j].yc < YMIN) particle[j].fy += KSPRING_BOUNDARY*(YMIN - particle[j].yc);
//             if (particle[j].xc > BCXMAX) particle[j].fx -= KSPRING_BOUNDARY*(particle[j].xc - BCXMAX);
//             else if (particle[j].xc < BCXMIN) particle[j].fx += KSPRING_BOUNDARY*(BCXMIN - particle[j].xc);
//             if (particle[j].yc > BCYMAX) particle[j].fy -= KSPRING_BOUNDARY*(particle[j].yc - BCYMAX);
//             else if (particle[j].yc < BCYMIN) particle[j].fy += KSPRING_BOUNDARY*(BCYMIN - particle[j].yc);
            return(fperp);
        }
        case (BC_RECTANGLE):
        {
            /* add harmonic force outside rectangular box */
            padding = MU + 0.01;
            xmin = xleft + padding;
            xmax = xright - padding;
            ymin = BCYMIN + padding;
            ymax = BCYMAX - padding;
            
            if (particle[j].xc > xmax) 
            {
                fperp = KSPRING_BOUNDARY*(particle[j].xc - xmax);
                particle[j].fx -= fperp;
            }
            else if (particle[j].xc < xmin) 
            {
                fperp = KSPRING_BOUNDARY*(xmin - particle[j].xc);
                particle[j].fx += fperp;
            }
            if (particle[j].yc > ymax) 
            {
                fperp = KSPRING_BOUNDARY*(particle[j].yc - ymax);
                particle[j].fy -= fperp;
            }
            else if (particle[j].yc < ymin) 
            {
                fperp = KSPRING_BOUNDARY*(ymin - particle[j].yc);
                particle[j].fy += fperp;
            }
//             if (particle[j].xc > xmax) particle[j].fx -= KSPRING_BOUNDARY*(particle[j].xc - xmax);
//             else if (particle[j].xc < xmin) particle[j].fx += KSPRING_BOUNDARY*(xmin - particle[j].xc);
//             if (particle[j].yc > ymax) particle[j].fy -= KSPRING_BOUNDARY*(particle[j].yc - ymax);
//             else if (particle[j].yc < ymin) particle[j].fy += KSPRING_BOUNDARY*(ymin - particle[j].yc);
            
            return(fperp);
        }
        case (BC_CIRCLE):
        {
            /* add harmonic force outside screen */
            if (particle[j].xc > BCXMAX) particle[j].fx -= KSPRING_BOUNDARY*(particle[j].xc - BCXMAX);
            else if (particle[j].xc < BCXMIN) particle[j].fx += KSPRING_BOUNDARY*(BCXMIN - particle[j].xc);
            if (particle[j].yc > BCYMAX) particle[j].fy -= KSPRING_BOUNDARY*(particle[j].yc - BCYMAX);
            else if (particle[j].yc < BCYMIN) particle[j].fy += KSPRING_BOUNDARY*(BCYMIN - particle[j].yc);
            
            /* add harmonic force from obstacle */
            for (i=-2; i<2; i++) 
            {
                x = xleft + (double)i*(OBSXMAX - OBSXMIN);
                if (vabs(particle[j].xc - x) < 1.1*OBSTACLE_RADIUS)
                {
                    r = module2(particle[j].xc - x, particle[j].yc);
                    if (r < 1.0e-5) r = 1.0e-05;
                    cphi = (particle[j].xc - x)/r;
                    sphi = particle[j].yc/r;
                    padding = MU + 0.03;
                    if (r < OBSTACLE_RADIUS + padding)
                    {
                        f = KSPRING_OBSTACLE*(OBSTACLE_RADIUS + padding - r);
                        particle[j].fx += f*cphi;
                        particle[j].fy += f*sphi;
                    }
                }
            }
            return(fperp);
        }
        case (BC_PERIODIC_CIRCLE):
        {
            x = xleft;
            if (vabs(particle[j].xc - x) < 1.1*OBSTACLE_RADIUS)
            {
                r = module2(particle[j].xc - x, particle[j].yc);
                if (r < 1.0e-5) r = 1.0e-05;
                cphi = (particle[j].xc - x)/r;
                sphi = particle[j].yc/r;
                padding = MU + 0.03;
                if (r < OBSTACLE_RADIUS + padding)
                {
                    f = KSPRING_OBSTACLE*(OBSTACLE_RADIUS + padding - r);
                    particle[j].fx += f*cphi;
                    particle[j].fy += f*sphi;
                }
            }
            return(f);
        }
        case (BC_PERIODIC_TRIANGLE):
        {
            x = xleft;
            x1 = x + OBSTACLE_RADIUS;
            x2 = x - OBSTACLE_RADIUS;
            h = 2.0*OBSTACLE_RADIUS*tanh(APOLY*PID);
            padding = MU + 0.03;
//             ytop = 0.5*h*(1.0 - (particle[j].xc - x)/OBSTACLE_RADIUS);
            if ((vabs(particle[j].xc - x) < OBSTACLE_RADIUS + padding)&&(vabs(particle[j].yc < h + padding)))
            {
                /* signed distances to side of triangle */
                dleft = x2 - particle[j].xc;
                norm = module2(h, 2.0*OBSTACLE_RADIUS);
                
                if (particle[j].yc >= 0.0)
                {
                    dplus = (h*particle[j].xc + 2.0*OBSTACLE_RADIUS*particle[j].yc - h*(x+OBSTACLE_RADIUS))/norm;
                    if ((dleft < padding)&&(dleft > dplus))     /* left side is closer */
                    {
                        f = KSPRING_OBSTACLE*(padding - dleft);
                        particle[j].fx -= f;
                    }
                    else if (dplus < padding)   /* top side is closer */
                    {
                        f = KSPRING_OBSTACLE*(padding - dplus);
                        particle[j].fx += f*h/norm;
                        particle[j].fy += 2.0*f*OBSTACLE_RADIUS/norm;
                    }
                }
                else   
                {
                    dminus = (h*particle[j].xc - 2.0*OBSTACLE_RADIUS*particle[j].yc - h*(x+OBSTACLE_RADIUS))/norm;
                    if ((dleft < padding)&&(dleft > dminus))     /* left side is closer */
                    {
                        f = KSPRING_OBSTACLE*(padding - dleft);
                        particle[j].fx -= f;
                    }
                    else if (dminus < padding)   /* bottom side is closer */
                    {
                        f = KSPRING_OBSTACLE*(padding - dminus);
                        particle[j].fx += f*h/norm;
                        particle[j].fy += -2.0*f*OBSTACLE_RADIUS/norm;
                    }
                }
                /* force from tip of triangle */
                r = module2(particle[j].xc - x1, particle[j].yc);
                if (r < 0.5*padding)
                {   
                    if (r < 1.0e-5) r = 1.0e-05;
                    cphi = (particle[j].xc - x1)/r;
                    sphi = particle[j].yc/r;
                    f = KSPRING_OBSTACLE*(0.5*padding - r);
                    particle[j].fx += f*cphi;
                    particle[j].fy += f*sphi;
                }
            }
            return(f);
        }
        case (BC_PERIODIC_FUNNEL):
        {
            x = xleft;
            padding = MU + 0.02;
            if (vabs(particle[j].yc) > FUNNEL_WIDTH - padding) for (i=-1; i<2; i+=2)
            {
                r = module2(particle[j].xc - x, particle[j].yc - (double)i*(FUNNEL_WIDTH + OBSTACLE_RADIUS));
                if (r < 1.0e-5) r = 1.0e-05;
                cphi = (particle[j].xc - x)/r;
                sphi = (particle[j].yc - (double)i*(FUNNEL_WIDTH + OBSTACLE_RADIUS))/r;
                if (r < OBSTACLE_RADIUS + padding)
                {
                    f = KSPRING_OBSTACLE*(OBSTACLE_RADIUS + padding - r);
                    particle[j].fx += f*cphi;
                    particle[j].fy += f*sphi;
                    if (RECORD_PRESSURES) 
                    {
                        angle = argument(cphi, sphi);
                        if (angle < 0.0) angle += DPI;
                        k = (int)((double)N_PRESSURES*angle/DPI);
                        if (k >= N_PRESSURES) k = N_PRESSURES - 1;
                        pressure[k] += f;
                    }
                }
            }
            return(f);
        }
        case (BC_RECTANGLE_LID):
        {
            r = particle[j].radius;
            
            if (particle[j].yc < BCYMIN + r) particle[j].fy += KSPRING_BOUNDARY*(BCYMIN + r - particle[j].yc);
            else if (particle[j].yc > ylid - r) 
            {
                fperp = KSPRING_BOUNDARY*(particle[j].yc - ylid + r);
                particle[j].fy -= fperp;
            }
            if (particle[j].yc < BCYMAX + r)
            {
                if (particle[j].xc > BCXMAX - r) particle[j].fx -= KSPRING_BOUNDARY*(particle[j].xc - BCXMAX + r);
                else if (particle[j].xc < BCXMIN + r) particle[j].fx += KSPRING_BOUNDARY*(BCXMIN + r - particle[j].xc);
            }
            return(fperp);
        }
        case (BC_RECTANGLE_WALL):
        {
            padding = particle[j].radius + 0.01;
            xmin = BCXMIN + padding;
            xmax = BCXMAX - padding;
            ymin = BCYMIN + padding;
            ymax = BCYMAX - padding;
            
            if (particle[j].xc > xmax) 
            {
                fperp = KSPRING_BOUNDARY*(particle[j].xc - xmax);
                particle[j].fx -= fperp;
                tmp_pright += fperp;
            }
            else if (particle[j].xc < xmin) 
            {
                fperp = KSPRING_BOUNDARY*(xmin - particle[j].xc);
                particle[j].fx += fperp;
                tmp_pleft += fperp;
            }
            if (particle[j].yc > ymax) 
            {
                fperp = KSPRING_BOUNDARY*(particle[j].yc - ymax);
                particle[j].fy -= fperp;
                if (particle[j].xc > xwall) tmp_pright += fperp;
                else tmp_pleft += fperp;
            }
            else if (particle[j].yc < ymin) 
            {
                fperp = KSPRING_BOUNDARY*(ymin - particle[j].yc);
                particle[j].fy += fperp;
                if (particle[j].xc > xwall) tmp_pright += fperp;
                else tmp_pleft += fperp;
            }
            
            if (wall)
            {
                *pleft += tmp_pleft/(2.0*(BCYMAX - BCYMIN) + 2.0*(xwall - BCXMIN));
                *pright += tmp_pright/(2.0*(BCYMAX - BCYMIN) + 2.0*(BCXMAX - xwall));
            }
            else 
            {
                *pleft += tmp_pleft/(2.0*(BCYMAX - BCYMIN + BCXMAX - BCXMIN));
                *pright += tmp_pright/(2.0*(BCYMAX - BCYMIN + BCXMAX - BCXMIN));
            }
            
            if ((wall)&&(vabs(particle[j].xc - xwall) < 0.5*WALL_WIDTH + padding))
            {
                if (particle[j].xc > xwall)
                {
                    fperp = -KSPRING_BOUNDARY*(xwall + 0.5*WALL_WIDTH + padding - particle[j].xc);
                    particle[j].fx -= fperp;
                    *pright -= fperp/(BCYMAX - BCYMIN);
                }
                else
                {
                    fperp = KSPRING_BOUNDARY*(particle[j].xc - xwall + 0.5*WALL_WIDTH + padding);
                    particle[j].fx -= fperp;
                    *pleft += fperp/(BCYMAX - BCYMIN);
                }
                return(fperp);
            }
            
            return(0.0);
        }
        case (BC_EHRENFEST):
        {
            rp = particle[j].radius;
            xtube = 1.0 - sqrt(EHRENFEST_RADIUS*EHRENFEST_RADIUS - EHRENFEST_WIDTH*EHRENFEST_WIDTH);
            distance = 0.0;
            /* middle tube */
            if (vabs(particle[j].xc) <= xtube) 
            {
                if (particle[j].yc > EHRENFEST_WIDTH - rp) 
                {
                    distance = particle[j].yc - EHRENFEST_WIDTH;
                    particle[j].fy -= KSPRING_BOUNDARY*(distance + rp);
                }
                else if (particle[j].yc < -EHRENFEST_WIDTH + rp) 
                {
                    distance = - EHRENFEST_WIDTH - particle[j].yc;
                    particle[j].fy += KSPRING_BOUNDARY*(distance + rp);
                }
            }
            /* right container */
            else if (particle[j].xc > 0.0)
            {
                r = module2(particle[j].xc - 1.0, particle[j].yc);
                if ((r > EHRENFEST_RADIUS - rp)&&((particle[j].xc > 1.0)||(vabs(particle[j].yc) > EHRENFEST_WIDTH)))
                {
                    cphi = (particle[j].xc - 1.0)/r;
                    sphi = particle[j].yc/r;
                    f = KSPRING_BOUNDARY*(EHRENFEST_RADIUS - r - rp);
                    particle[j].fx += f*cphi;
                    particle[j].fy += f*sphi;
                    *pright -= f;
                }
            }
            /* left container */
            else 
            {
                r = module2(particle[j].xc + 1.0, particle[j].yc);
                if ((r > EHRENFEST_RADIUS - rp)&&((particle[j].xc < -1.0)||(vabs(particle[j].yc) > EHRENFEST_WIDTH)))
                {
                    cphi = (particle[j].xc + 1.0)/r;
                    sphi = particle[j].yc/r;
                    f = KSPRING_BOUNDARY*(EHRENFEST_RADIUS - r - rp);
                    particle[j].fx += f*cphi;
                    particle[j].fy += f*sphi;
                    *pleft -= f;
                }
            }
            
            /* add force from "corners" */
            if ((vabs(particle[j].xc) - xtube < rp)&&(vabs(particle[j].yc) - EHRENFEST_WIDTH < rp))
            {
                for (i=-1; i<=1; i+=2)
                    for (k=-1; k<=1; k+=2)
                    {
                        distance = module2(particle[j].xc - (double)i*xtube, particle[j].yc - (double)k*EHRENFEST_WIDTH);
                        if (distance < rp)
                        {
                            cphi = (particle[j].xc - (double)i*xtube)/distance;
                            sphi = (particle[j].yc - (double)k*EHRENFEST_WIDTH)/distance;
                            f = KSPRING_BOUNDARY*(rp - distance);
                            particle[j].fx += f*cphi;
                            particle[j].fy += f*sphi;
                        }
                    }
            }
            return(fperp);
        }
        case (BC_SCREEN_BINS):
        {
            /* add harmonic force outside screen */
            if (particle[j].xc > XMAX) particle[j].fx -= KSPRING_BOUNDARY*(particle[j].xc - XMAX);
            else if (particle[j].xc < XMIN) particle[j].fx += KSPRING_BOUNDARY*(XMIN - particle[j].xc);
            if (particle[j].yc > YMAX + 10.0*MU) particle[j].fy -= KSPRING_BOUNDARY*(particle[j].yc - YMAX - 10.0*MU);
            else if (particle[j].yc < YMIN) particle[j].fy += KSPRING_BOUNDARY*(YMIN - particle[j].yc);
                        
            /* force from the bins */
            dy = (YMAX - YMIN)/((double)NGRIDX + 3);
            dx = dy/cos(PI/6.0);
            rp = particle[j].radius;
            width = rp + 0.05*dx;
            ybin = 2.75*dy;
            
            if (particle[j].yc < YMIN + ybin) for (i=-1; i<=NGRIDX; i++) 
            {
                x = ((double)i - 0.5*(double)NGRIDX + 0.5)*dx;
                distance = vabs(particle[j].xc - x);
                if (distance < width)
                {
                    if (particle[j].xc > x) particle[j].fx += KSPRING_BOUNDARY*(width - distance);
                    else particle[j].fx -= KSPRING_BOUNDARY*(width - distance);
                }
            }
            else if (particle[j].yc < YMIN + ybin + particle[j].radius) for (i=-1; i<=NGRIDX; i++) 
            {
                x = ((double)i - 0.5*(double)NGRIDX + 0.5)*dx;
                distance = module2(particle[j].xc - x, particle[j].yc - YMIN - ybin);
                if (distance < rp)
                {
                    if (distance < 1.0e-8) distance = 1.0e-8;
                    cphi = (particle[j].xc - x)/distance;
                    sphi = (particle[j].yc - YMIN - ybin)/distance;
                    f = KSPRING_BOUNDARY*(rp - distance);
                    particle[j].fx += f*cphi;
                    particle[j].fy += f*sphi;
                }
            }
            return(fperp);
        }
        case (BC_ABSORBING):
        {
            /* add harmonic force outside screen */
            padding = 0.1;
            x = particle[j].xc;
            y = particle[j].yc;
            
            if ((x > BCXMAX + padding)||(x < BCXMIN - padding)||(y > BCYMAX + padding)||(y < BCYMIN - padding)) 
            {
                particle[j].active = 0;
                particle[j].vx = 0.0;
                particle[j].vy = 0.0;
                particle[j].xc = BCXMAX + 2.0*padding;
                particle[j].yc = BCYMAX + 2.0*padding;
            }
            
            return(fperp);
        }
        case (BC_REFLECT_ABS):
        {
            /* add harmonic force outside screen */
            padding = 0.1;
            x = particle[j].xc;
            y = particle[j].yc;
            
            if ((x > BCXMAX + padding)||(y > BCYMAX + padding)||(y < BCYMIN - padding)) 
            {
                particle[j].active = 0;
                particle[j].vx = 0.0;
                particle[j].vy = 0.0;
                particle[j].xc = BCXMAX + 2.0*padding;
                particle[j].yc = BCYMAX + 2.0*padding;
            }
             else if (particle[j].yc < BCYMIN) particle[j].fy += KSPRING_BOUNDARY*(BCYMIN - particle[j].yc);
            
            return(fperp);
        }
        case (BC_REFLECT_ABS_BOTTOM):
        {
            /* add harmonic force outside screen */
            padding = MU;
            x = particle[j].xc;
            y = particle[j].yc;
            
            if (y < BCYMIN)
            {
                particle[j].active = 0;
                particle[j].vx = 0.0;
                particle[j].vy = 0.0;
                particle[j].xc = BCXMAX + 2.0*padding;
                particle[j].yc = BCYMIN - 2.0*padding;
            }
            else 
            {
                if (y > BCYMAX + padding) particle[j].fy -= KSPRING_BOUNDARY*(y - BCYMAX - padding);
                if (x > BCXMAX - padding) particle[j].fx -= KSPRING_BOUNDARY*(x - BCXMAX + padding);
                else if (x < BCXMIN + padding) particle[j].fx += KSPRING_BOUNDARY*(BCXMIN - x + padding);
            }
                        
            return(0.0);
        }
    }
}


void compute_particle_force(int j, double krepel, t_particle particle[NMAXCIRCLES], t_hashgrid hashgrid[HASHX*HASHY])
/* compute force from other particles on particle j */
{
    int i0, j0, m0, k, m, q, close;
    double fx = 0.0, fy = 0.0, force[2], torque = 0.0, torque_ij, x, y;
    
    particle[j].neighb = 0;
    if (REACTION_DIFFUSION) particle[j].diff_neighb = 0;

    /* NEW */
    for (k=0; k<particle[j].hash_nneighb; k++) 
        if (particle[particle[j].hashneighbour[k]].active)
    {
        close = compute_repelling_force(j, k, force, &torque_ij, particle, krepel);
        fx += force[0];
        fy += force[1];
        torque += torque_ij;
        if (close) 
        {
            particle[j].neighb++;
            if (REACTION_DIFFUSION&&(particle[j].type != particle[particle[j].hashneighbour[k]].type)) 
                particle[j].diff_neighb++;
        }
    }
                
    particle[j].fx += fx;
    particle[j].fy += fy;
    particle[j].torque += torque;
}


int reorder_particles(t_particle particle[NMAXCIRCLES], double py[NMAXCIRCLES], double pangle[NMAXCIRCLES])
/* keep only active particles, beta */
{
    int i, k, new = 0, nactive = 0;
    
    for (i=0; i<ncircles; i++) if (particle[i].active)
    {
        particle[new].xc = particle[i].xc;
        particle[new].yc = particle[i].yc;
        particle[new].radius = particle[i].radius;
        particle[new].angle = particle[i].angle;
        particle[new].active = particle[i].active;
        particle[new].energy = particle[i].energy;
        particle[new].vx = particle[i].vx;
        particle[new].vy = particle[i].vy;
        particle[new].omega = particle[i].omega;
        particle[new].mass_inv = particle[i].mass_inv;
        particle[new].inertia_moment_inv = particle[i].inertia_moment_inv;
        particle[new].fx = particle[i].fx;
        particle[new].fy = particle[i].fy;
        particle[new].thermostat = particle[i].thermostat;
        particle[new].hashcell = particle[i].hashcell;
        particle[new].neighb = particle[i].neighb;
        particle[new].diff_neighb = particle[i].diff_neighb;
        particle[new].hash_nneighb = particle[i].hash_nneighb;
        particle[new].type = particle[i].type;
        particle[new].interaction = particle[i].interaction;
        particle[new].eq_dist = particle[i].eq_dist;
        particle[new].spin_range = particle[i].spin_range;
        particle[new].spin_freq = particle[i].spin_freq;
        
        py[new] = py[i];
        pangle[new] = pangle[i];
        
        for (k=0; k<9*HASHMAX; k++)
        {
            particle[new].hashneighbour[k] = particle[i].hashneighbour[k];
            particle[new].deltax[k] = particle[i].deltax[k];
            particle[new].deltay[k] = particle[i].deltay[k];
        }
        
        new++;
        nactive++;
    }
    
    return(nactive);
}


int initialize_configuration(t_particle particle[NMAXCIRCLES], t_hashgrid hashgrid[HASHX*HASHY], t_obstacle obstacle[NMAXOBSTACLES], double px[NMAXCIRCLES], double py[NMAXCIRCLES], double pangle[NMAXCIRCLES], int tracer_n[N_TRACER_PARTICLES],
t_segment segment[NMAXSEGMENTS])
/* initialize all particles, obstacles, and the hashgrid */
{
    int i, j, k, n, nactive = 0;
    double x, y, h, xx, yy, rnd;
    
    for (i=0; i < ncircles; i++) 
    {
        /* set particle type */
        particle[i].type = 0;
        if ((TWO_TYPES)&&((double)rand()/RAND_MAX > TYPE_PROPORTION)) 
        {
            particle[i].type = 2;
            particle[i].radius = MU_B;
        }
        if ((INTERACTION == I_VICSEK_SHARK)&&(i==1)) 
        {
            particle[i].type = 2;
            particle[i].radius = MU_B;
        }
        
        particle[i].neighb = 0;
        particle[i].diff_neighb = 0;
        particle[i].thermostat = 1;
        particle[i].close_to_boundary = 0;
        particle[i].emean = 0.0;
        particle[i].dirmean = 0.0;

//         particle[i].energy = 0.0;
//         y = particle[i].yc;
//         if (y >= YMAX) y -= particle[i].radius;
//         if (y <= YMIN) y += particle[i].radius;

        if (RANDOM_RADIUS) particle[i].radius = particle[i].radius*(0.75 + 0.5*((double)rand()/RAND_MAX));
        
        if (particle[i].type == 0)
        {
            particle[i].interaction = INTERACTION;
            particle[i].eq_dist = EQUILIBRIUM_DIST;
            particle[i].spin_range = SPIN_RANGE;
            particle[i].spin_freq = SPIN_INTER_FREQUENCY;
            particle[i].mass_inv = 1.0/PARTICLE_MASS;
            particle[i].inertia_moment_inv = 1.0/PARTICLE_INERTIA_MOMENT;
        }
        else 
        {
            particle[i].interaction = INTERACTION_B;
            particle[i].eq_dist = EQUILIBRIUM_DIST_B;
            particle[i].spin_range = SPIN_RANGE_B;
            particle[i].spin_freq = SPIN_INTER_FREQUENCY_B;
            particle[i].mass_inv = 1.0/PARTICLE_MASS_B;
            particle[i].inertia_moment_inv = 1.0/PARTICLE_INERTIA_MOMENT_B;
        }
                
        particle[i].vx = V_INITIAL*gaussian();
        particle[i].vy = V_INITIAL*gaussian();
        
        if ((INTERACTION == I_VICSEK_SHARK)&&(i==1)) 
        {
            particle[i].vx *= 1000.0;
            particle[i].vy *= 1000.0;
        }
        
        particle[i].energy = (particle[i].vx*particle[i].vx + particle[i].vy*particle[i].vy)*particle[i].mass_inv;
        particle[i].emean = particle[i].energy;
        particle[i].dirmean = 0.0;
        
        px[i] = particle[i].vx;
        py[i] = particle[i].vy;
        
        if (ROTATION) 
        {
            particle[i].angle = DPI*(double)rand()/RAND_MAX;
            particle[i].omega = OMEGA_INITIAL*gaussian();
            if (COUPLE_ANGLE_TO_THERMOSTAT) 
                particle[i].energy += particle[i].omega*particle[i].omega*particle[i].inertia_moment_inv;
        }
        else 
        {
            particle[i].angle = 0.0;
            particle[i].omega = 0.0;
        }
        pangle[i] = particle[i].omega;
        
        if ((PLOT == P_INITIAL_POS)||(PLOT_B == P_INITIAL_POS))
        {
            switch (INITIAL_POS_TYPE) {
                case (IP_X):
                {
                    particle[i].color_hue = 360.0*(particle[i].xc - INITXMIN)/(INITXMAX - INITXMIN);
                    break;
                }
                 case (IP_Y):
                {
                    particle[i].color_hue = 360.0*(particle[i].yc - INITYMIN)/(INITYMAX - INITYMIN);
                    break;
                }
            }
            
        }
        else if ((PLOT == P_NUMBER)||(PLOT_B == P_NUMBER))
            particle[i].color_hue = 360.0*(double)(i%N_PARTICLE_COLORS)/(double)N_PARTICLE_COLORS;
    }
    /* initialize dummy values in case particles are added */
    for (i=ncircles; i < NMAXCIRCLES; i++) 
    {
        particle[i].type = 0;
        particle[i].active = 0;
        particle[i].neighb = 0;
        particle[i].thermostat = 0;
        particle[i].energy = 0.0;
        particle[i].emean = 0.0;
        particle[i].dirmean = 0.0;
        particle[i].mass_inv = 1.0/PARTICLE_MASS;
        particle[i].inertia_moment_inv = 1.0/PARTICLE_INERTIA_MOMENT;
        particle[i].vx = 0.0;
        particle[i].vy = 0.0;
        px[i] = 0.0;
        py[i] = 0.0;        
        particle[i].angle = DPI*(double)rand()/RAND_MAX;
        particle[i].omega = 0.0;
        pangle[i] = 0.0;
        particle[i].interaction = INTERACTION;
        particle[i].eq_dist = EQUILIBRIUM_DIST;
        particle[i].spin_range = SPIN_RANGE;
        particle[i].spin_freq = SPIN_INTER_FREQUENCY;
        particle[i].close_to_boundary = 0;
   }
    
    /* add particles at the bottom as seed */
    if (PART_AT_BOTTOM) for (i=0; i<=NPART_BOTTOM; i++)
    {   
        x = XMIN + (double)i*(XMAX - XMIN)/(double)NPART_BOTTOM;
        y = YMIN + 2.0*MU;
        add_particle(x, y, 0.0, 0.0, MASS_PART_BOTTOM, 0, particle);
    }
    if (PART_AT_BOTTOM) for (i=0; i<=NPART_BOTTOM; i++)
    {   
        x = XMIN + (double)i*(XMAX - XMIN)/(double)NPART_BOTTOM;
        y = YMIN + 4.0*MU;
        add_particle(x, y, 0.0, 0.0, MASS_PART_BOTTOM, 0, particle);
    }
    
    /* add larger copies of particles (for Ehrenfest model)*/
    if (EHRENFEST_COPY)
    {
        for (i=0; i < ncircles; i++) 
        {
            n = ncircles + i;
            particle[n].xc = -particle[i].xc;
            particle[n].yc = particle[i].yc;
            particle[n].vx = -0.5*particle[i].vx;
            particle[n].vy = 0.5*particle[i].vy;
            px[n] = -0.5*px[i];
            py[n] = 0.5*py[i];        
            particle[n].energy = particle[i].energy;
            particle[n].radius = 2.0*particle[i].radius;
            particle[n].type = 2;
            particle[n].mass_inv = 1.25*particle[i].mass_inv;
            particle[n].thermostat = 1;
            particle[n].interaction = particle[i].interaction;
            particle[n].eq_dist = 0.45*particle[i].eq_dist;
            
            if ((double)rand()/RAND_MAX > 0.6) particle[n].active = 1;
        }
        ncircles *= 2;
    }
    
    /* change type of tracer particle */
    if (TRACER_PARTICLE) for (j=0; j<N_TRACER_PARTICLES; j++)
    {
        i = 0;
        if (j%2==0) xx = 1.0;
        else xx = -1.0;
        
        if (j/2 == 0) yy = -0.5;
        else yy = 0.5;
        
//         if (j%2 == 1) yy = -yy;
//         while ((!particle[i].active)||(module2(particle[i].xc, particle[i].yc) > 0.5)) i++;
        while ((!particle[i].active)||(module2(particle[i].xc + xx, particle[i].yc - yy) > 0.4)) i++;
        tracer_n[j] = i;
        particle[i].type = 2 + j;
        particle[i].radius *= 1.5;
        particle[i].mass_inv *= 1.0/TRACER_PARTICLE_MASS;
        particle[i].vx *= 0.1;
        particle[i].vy *= 0.1;
        particle[i].thermostat = 0;
        px[i] *= 0.1;
        py[i] *= 0.1;
    }
    
    /* position-dependent particle type */
    if (POSITION_DEPENDENT_TYPE) for (i=0; i<ncircles; i++)
        if (((!POSITION_Y_DEPENDENCE)&&((particle[i].xc - POSITION_DEP_X)*POSITION_DEP_SIGN < 0.0))||((POSITION_Y_DEPENDENCE)&&(particle[i].yc*POSITION_DEP_SIGN < 0.0))) 
        {
            particle[i].type = 2;
            particle[i].mass_inv = 1.0/PARTICLE_MASS_B;
            particle[i].radius = MU_B;
        }
    
    /* inactivate particles in obstacle */
    printf("Inactivating particles inside obstacles\n");
    if ((BOUNDARY_COND == BC_CIRCLE)||(BOUNDARY_COND == BC_PERIODIC_CIRCLE))
    {
        for (i=0; i< ncircles; i++)
            if ((module2(particle[i].xc - OBSTACLE_XMIN, particle[i].yc) < 1.2*OBSTACLE_RADIUS)) 
                particle[i].active = 0;
    }
    else if (BOUNDARY_COND == BC_PERIODIC_FUNNEL)
    {
        for (i=0; i< ncircles; i++)
            for (k=-1; k<2; k+=2)
                if ((module2(particle[i].xc, particle[i].yc - (double)k*(FUNNEL_WIDTH + OBSTACLE_RADIUS)) < OBSTACLE_RADIUS + 2.0*MU)) 
                {
                    printf("Inactivating particle at (%.3lg, %.3lg)\n", particle[i].xc, particle[i].yc);
                    particle[i].active = 0;
                }
    }
    else if (BOUNDARY_COND == BC_PERIODIC_TRIANGLE)
    {
        h = 2.0*OBSTACLE_RADIUS*tan(APOLY*PID);
        for (i=0; i< ncircles; i++)
            if ((vabs(particle[i].xc) < 1.1*OBSTACLE_RADIUS + 2.0*MU)
                &&(2.0*OBSTACLE_RADIUS*vabs(particle[i].yc) < h*(OBSTACLE_RADIUS + 2.0*MU - particle[i].xc)))
                {
                    printf("Inactivating particle at (%.3lg, %.3lg)\n", particle[i].xc, particle[i].yc);
                    particle[i].active = 0;
                }            
    }
    else if (BOUNDARY_COND == BC_EHRENFEST)
    {
        for (i=0; i< ncircles; i++)
            if (module2(vabs(particle[i].xc) -1.0, particle[i].yc) > EHRENFEST_RADIUS) 
                particle[i].active = 0;
    }
    else if (BOUNDARY_COND == BC_RECTANGLE_WALL)
    {
        for (i=0; i< ncircles; i++)
            if (vabs(particle[i].xc - xwall) < WALL_WIDTH) 
                particle[i].active = 0;
    }
    else if (BOUNDARY_COND == BC_GENUS_TWO)
    {
        for (i=0; i< ncircles; i++)
            if ((particle[i].xc > 0.0)&&(particle[i].yc > 0.0)) 
                particle[i].active = 0;
    }
    
    
    if (ADD_FIXED_OBSTACLES)
    {
        for (i=0; i< ncircles; i++) for (j=0; j < nobstacles; j++)
            if (module2(particle[i].xc - obstacle[j].xc, particle[i].yc - obstacle[j].yc) < OBSTACLE_RADIUS + particle[i].radius)
                particle[i].active = 0;
    }
    
    /* case of segment obstacles */
    if (ADD_FIXED_SEGMENTS) for (i=0; i< ncircles; i++)
            if (!in_segment_region(particle[i].xc, particle[i].yc, segment)) 
                particle[i].active = 0;
            
    /* case of reaction-diffusion equation/chemical reactions */
    if (REACTION_DIFFUSION) for (i=0; i< ncircles; i++)
    {
        switch (RD_INITIAL_COND) {
            case (IC_UNIFORM):
            {
                particle[i].type = 1;
                break;
            }
            case (IC_UNIFORM2):
            {
                particle[i].type = 2;
                particle[i].radius = MU_B;
                particle[i].omega = OMEGA_INITIAL*gaussian();
                break;
            }
            case (IC_RANDOM_UNIF):
            {
                particle[i].type = 1 + (int)(RD_TYPES*(double)rand()/(double)RAND_MAX);
                break;
            }
            case (IC_RANDOM_TWO):
            {
                if ((double)rand()/(double)RAND_MAX < TYPE_PROPORTION) particle[i].type = 1;
                else 
                {
                    particle[i].type = 2;
                    particle[i].radius = MU_B;
                    particle[i].mass_inv = 1.0/PARTICLE_MASS_B;
                }
                break;
            }
            case (IC_CIRCLE):
            {
                if (module2(particle[i].xc,particle[i].yc) < LAMBDA) particle[i].type = 1;
                else 
                {
                    particle[i].type = 2;
                    particle[i].radius = MU_B;
                    particle[i].mass_inv = 1.0/PARTICLE_MASS_B;
                }
                break;
            }
            case (IC_CATALYSIS):
            {
                if ((particle[i].xc > 0.0)||((double)rand()/(double)RAND_MAX < TYPE_PROPORTION))
                {
                    particle[i].type = 1;
                    particle[i].radius = MU;
                    particle[i].mass_inv = 1.0/PARTICLE_MASS;
                }
                else 
                {
                    particle[i].type = 2;
                    particle[i].radius = MU_B;
                    particle[i].mass_inv = 1.0/PARTICLE_MASS_B;
                }
                break;
            }
            case (IC_LAYERS):
            {
                if (particle[i].yc > 0.0)
                {
                    particle[i].type = 1;
                    particle[i].radius = MU;
                    particle[i].eq_dist = EQUILIBRIUM_DIST;
                    particle[i].mass_inv = 1.0/PARTICLE_MASS;
                }
                else if (particle[i].yc < -LAMBDA)
                {
                    particle[i].type = 2;
                    particle[i].radius = MU_B;
                    particle[i].eq_dist = EQUILIBRIUM_DIST_B;
                    particle[i].mass_inv = 1.0/PARTICLE_MASS_B;
                }
                else particle[i].active = 0;
                break;
            }
            case (IC_BZ):
            {
                rnd = (double)rand()/(double)RAND_MAX;
                if (rnd < 60.0/180.0) 
                {
                    particle[i].type = 4;
                    particle[i].radius = MU;
                    particle[i].eq_dist = EQUILIBRIUM_DIST;
                    particle[i].mass_inv = 1.0/(PARTICLE_MASS + 3.0*PARTICLE_MASS_B);
                }
                else if (rnd < 100.0/180.0) 
                {
                    particle[i].type = 5;
                    particle[i].radius = MU_B;
                    particle[i].eq_dist = EQUILIBRIUM_DIST;
                    particle[i].mass_inv = 1.0/(3.5*PARTICLE_MASS);
                }
                else 
                {
                    particle[i].type = 6;
                    particle[i].radius = 1.5*MU_B;
                    particle[i].eq_dist = EQUILIBRIUM_DIST;
                    particle[i].mass_inv = 0.2/(PARTICLE_MASS);
                }
                break;
            }
            case (IC_SIGNX):
            {
                if (particle[i].xc < 0.0) particle[i].type = 1;
                else 
                {
                    particle[i].type = 2;
                    particle[i].radius = MU_B;
                    particle[i].mass_inv = 1.0/PARTICLE_MASS_B;
                }
                break;
            }
            case (IC_TWOROCKETS):
            {
                if (vabs(particle[i].xc) < SEGMENTS_X0) particle[i].type = 1;
                else 
                {
                    particle[i].type = 2;
                    particle[i].radius = MU_B;
                    particle[i].mass_inv = 1.0/PARTICLE_MASS_B;
                }
                break;
            }
            case (IC_TWOROCKETS_TWOFUELS):
            {
                if (vabs(particle[i].xc) < SEGMENTS_X0) particle[i].type = 1;
                else 
                {
                    if (particle[i].xc < 0) particle[i].type = 2;
                    else particle[i].type = 3;
                    particle[i].radius = MU_B;
                    particle[i].mass_inv = 1.0/PARTICLE_MASS_B;
                }
                break;
            }
        }
    }   
    
    /* keep only active particles */
//     ncircles = reorder_particles(particle, py, pangle);
            
    /* count number of active particles */
    for (i=0; i< ncircles; i++) nactive += particle[i].active;
    printf("%i active particles\n", nactive);
    
//     for (i=0; i<ncircles; i++) printf("particle %i of type %i\n", i, particle[i].type);
    
    return(nactive);
}


int add_particles(t_particle particle[NMAXCIRCLES], double px[NMAXCIRCLES], double py[NMAXCIRCLES], int nadd_particle)
/* add several particles to the system */
{
    static int i = 0;
    double x, y; 
    
//             add_particle(XMIN + 0.1, 0.0, 50.0, 0.0, 3.0, 0, particle);
//             px[ncircles - 1] = particle[ncircles - 1].vx;
//             py[ncircles - 1] = particle[ncircles - 1].vy;
//             particle[ncircles - 1].radius = 1.5*MU;
//             j = 0;
//             while (module2(particle[j].xc,particle[j].yc) > 0.7) j = rand()%ncircles;
//             x =  particle[j].xc + 2.5*MU;
//             y =  particle[j].yc;
            
//             x =  XMIN + (XMAX - XMIN)*rand()/RAND_MAX;
//             y =  YMAX + 0.01*rand()/RAND_MAX;
//             add_particle(x, y, 0.0, 0.0, 1.0, 0, particle);

//             x =  XMIN + 0.25*(XMAX - XMIN);
//             y =  YMAX + 0.01;
//             prop = 1.0 - (double)nadd_particle/5.0;
//             vx = 100.0*prop;
//             add_particle(x, y, vx, -10.0, 5.0*prop, 0, particle);
//             particle[ncircles - 1].radius = 10.0*MU*prop;
//             particle[ncircles - 1].eq_dist = 2.0;
//             particle[ncircles - 1].thermostat = 0;
//             px[ncircles - 1] = particle[ncircles - 1].vx;
//             py[ncircles - 1] = particle[ncircles - 1].vy;
//     add_particle(MU*(2.0*rand()/RAND_MAX - 1.0), YMAX + 2.0*MU, 0.0, 0.0, PARTICLE_MASS, 0, particle);
     
//     add_particle(XMIN - 0.5*MU, 0.0, 50.0 + 5.0*(double)i, 0.0, 2.0*PARTICLE_MASS, 0, particle);
    
//     x = INITXMIN + (INITXMAX - INITXMIN)*(double)rand()/(double)RAND_MAX;
//     y = INITYMIN + (INITYMAX - INITYMIN)*(double)rand()/(double)RAND_MAX;
//     x = BCXMIN + (BCXMAX - BCXMIN)*(double)rand()/(double)RAND_MAX;
//     y = YMAX + 0.5*(BCYMAX - YMAX)*(double)rand()/(double)RAND_MAX;
    
    printf("Adding a particle\n\n");
    
    x = ADDXMIN + (ADDXMAX - ADDXMIN)*(double)rand()/(double)RAND_MAX;
    y = ADDYMIN + 0.5*(ADDYMAX - ADDYMIN)*(double)rand()/(double)RAND_MAX;
    add_particle(x, y, 0.0, 0.0, PARTICLE_MASS, 0, particle);
//     add_particle(BCXMIN + 0.1, 0.5*(BCYMIN + BCYMAX), 200.0, 0.0, PARTICLE_MASS, 0, particle);
//     i++;
    
    particle[ncircles - 1].radius = MU;
    particle[ncircles - 1].eq_dist = EQUILIBRIUM_DIST;
    particle[ncircles - 1].thermostat = 0;
    px[ncircles - 1] = particle[ncircles - 1].vx;
    py[ncircles - 1] = particle[ncircles - 1].vy;
    return (nadd_particle + 1);
}


void center_momentum(double p[NMAXCIRCLES])
{
    int i;
    double ptot = 0.0, pmean;
    
    for (i=0; i<ncircles; i++) ptot += p[i];
    pmean = ptot/(double)ncircles;
    for (i=0; i<ncircles; i++) p[i] -= pmean;
}


int floor_momentum(double p[NMAXCIRCLES])
{
    int i, floor = 0;
    double ptot = 0.0, pmean;
    
    for (i=0; i<ncircles; i++) 
    {
        if (p[i] > PMAX) 
        {
            p[i] = PMAX;
            floor = 1;
        }
        else if (p[i] < -PMAX) 
        {
            p[i] = -PMAX;
            floor = 1;
        }
    }
    if (floor) printf("Flooring momentum\n");
    return (floor); 
}

int partial_thermostat_coupling(t_particle particle[NMAXCIRCLES], double xmin, t_segment segment[NMAXSEGMENTS])
/* only couple particles satisfying condition PARTIAL_THERMO_REGION to thermostat */
{
    int condition, i, nthermo = 0;
    double x, y, height;
    static double maxheight;
    static int first = 1;
    
    if (first)
    {
        maxheight = YMIN + PARTIAL_THERMO_HEIGHT*(YMAX - YMIN);
        first = 0;
    }
    
    if (PARTIAL_THERMO_REGION == TH_LAYER_TYPE2)
    {
        height = YMIN;
        for (i=0; i<ncircles; i++)
        {
            y = particle[i].yc;
            if ((particle[i].active)&&(particle[i].type == 2)&&(y > height)&&(y <= maxheight)) 
                height = y;
        }
        if (height > maxheight) height = maxheight;
        printf("Thermostat region y > %.3lg, max height = %.3lg\n", height, maxheight);
    }
    
    for (i=0; i<ncircles; i++) 
    {
        switch (PARTIAL_THERMO_REGION) {
            case (TH_VERTICAL):
            {
                condition = (particle[i].xc > xmin);
                break;
            }
            case (TH_INSEGMENT):
            {
                condition = (in_segment_region(particle[i].xc, particle[i].yc, segment));
//                 condition = (in_segment_region(particle[i].xc - xsegments[0], particle[i].yc - ysegments[0]));
// //                 condition = (in_segment_region(particle[i].xc - xsegments[0], particle[i].yc - ysegments[0]));
//                 if (TWO_OBSTACLES) 
//                     condition = condition||(in_segment_region(particle[i].xc - xsegments[1], particle[i].yc - ysegments[1]));
                break;
            }
            case (TH_INBOX):
            {
                x = particle[i].xc;
                y = particle[i].yc;
                condition = ((y < YMIN + PARTIAL_THERMO_HEIGHT*(YMAX - YMIN))&&(vabs(x) < PARTIAL_THERMO_WIDTH*XMAX));
                break;
            }
            case (TH_LAYER):
            {
                y = particle[i].yc;
                condition = (y > PARTIAL_THERMO_HEIGHT);
                break;
            }
            case (TH_LAYER_TYPE2):
            {
                y = particle[i].yc;
                condition = (y > height);
                break;
            }
            default: condition = 1;
        }
        if (condition) 
        {
            particle[i].thermostat = 1;
            nthermo++;
        }
        else particle[i].thermostat = 0;
    }
    return(nthermo);
}

double compute_mean_energy(t_particle particle[NMAXCIRCLES])
{
    int i, nactive = 0;
    double total_energy = 0.0;
    
    for (i=0; i<ncircles; i++) if (particle[i].active)
    {
        total_energy += particle[i].energy;
        nactive++;
    }    
    return(total_energy/(double)nactive);
}


void compute_inverse_masses(double inv_masses[RD_TYPES+1])
/* compute inverse masses of molecules in chemical reactions */
{
    int type; 
    double mass;
    
    /* default values that apply in most cases */
    inv_masses[1] = 1.0/PARTICLE_MASS;
    inv_masses[2] = 1.0/PARTICLE_MASS_B;
    
    switch (RD_REACTION) {
        case (CHEM_AAB): 
        {
            inv_masses[2] = 0.5/PARTICLE_MASS;
            break;
        }
        case (CHEM_ABC): 
        {
            inv_masses[3] = 1.0/(PARTICLE_MASS + PARTICLE_MASS_B);
            break;
        }
        case (CHEM_A2BC): 
        {
            inv_masses[3] = 1.0/(2.0*PARTICLE_MASS + PARTICLE_MASS_B);
            break;
        }
        case (CHEM_CATALYSIS): 
        {
            inv_masses[3] = 0.5/PARTICLE_MASS;
            break;
        }
        case (CHEM_BAA): 
        {
            inv_masses[1] = 2.0/PARTICLE_MASS_B;
            break;
        }
        case (CHEM_AABAA): 
        {
            inv_masses[2] = 0.5/PARTICLE_MASS;
            break;
        }
        case (CHEM_POLYMER):
        {
            for (type = 3; type < RD_TYPES+1; type++)
            {
                mass = PARTICLE_MASS_B + (double)(type-2)*PARTICLE_MASS;
                inv_masses[type] = 1.0/(PARTICLE_MASS + mass);
            }
            break;
        }
        case (CHEM_POLYMER_DISS):
        {
            for (type = 3; type < RD_TYPES+1; type++)
            {
                mass = PARTICLE_MASS_B + (double)(type-2)*PARTICLE_MASS;
                inv_masses[type] = 1.0/(PARTICLE_MASS + mass);
            }
            break;
        }
        case (CHEM_POLYMER_STEP):
        {
            for (type = 2; type < RD_TYPES+1; type++)
                inv_masses[type] = 1.0/((double)type*PARTICLE_MASS);
            break;
        }
        case (CHEM_CATALYTIC_A2D): 
        {
            inv_masses[3] = 1.0/(1.0/PARTICLE_MASS + 1.0/PARTICLE_MASS_B);
            inv_masses[4] = 0.5/PARTICLE_MASS;
            break;
        }
        case (CHEM_ABCAB): 
        {
            inv_masses[3] = 1.0/(PARTICLE_MASS + PARTICLE_MASS_B);
            break;
        }
        case (CHEM_ABCDABC): 
        {
            inv_masses[3] = 1.0/(PARTICLE_MASS + PARTICLE_MASS_B);
            inv_masses[4] = 1.0/(2.0*PARTICLE_MASS + PARTICLE_MASS_B);
            break;
        }
        case (CHEM_BZ):
        {
            for (type = 2; type <= 4; type++)
            {
                mass = PARTICLE_MASS + (double)(type-2)*PARTICLE_MASS_B;
                inv_masses[type] = 1.0/(mass);
            }
//             inv_masses[5] = 1.0/(3.5*PARTICLE_MASS);
//             inv_masses[6] = 0.2/PARTICLE_MASS;
            inv_masses[5] = 0.5/(PARTICLE_MASS);
            inv_masses[6] = 0.5/PARTICLE_MASS;
            inv_masses[7] = 0.5*inv_masses[3];
            inv_masses[8] = 0.5*inv_masses[5];
            break;
        }
        case (CHEM_BRUSSELATOR):
        {
            inv_masses[3] = inv_masses[1];
            inv_masses[4] = 1.0/(PARTICLE_MASS + PARTICLE_MASS_B);
            inv_masses[5] = inv_masses[1];
            inv_masses[6] = inv_masses[1];
            break;
        }
        case (CHEM_ABDACBE): 
        {
            inv_masses[3] = inv_masses[1];
            inv_masses[4] = 1.0/(PARTICLE_MASS + PARTICLE_MASS_B);
            inv_masses[5] = inv_masses[1];
            break;
        }
    }
}

void compute_radii(double radii[RD_TYPES+1])
/* compute radii of molecules in chemical reactions */
{
    int type; 
    
    /* default values that apply in most cases */
    radii[1] = MU;
    radii[2] = MU_B;
    
    switch (RD_REACTION) {
        case (CHEM_ABC): 
        {
            radii[3] = 1.5*MU;
            break;
        }
        case (CHEM_A2BC): 
        {
            radii[3] = 1.5*MU;
            break;
        }
        case (CHEM_CATALYSIS): 
        {
            radii[3] = MU_B;
            break;
        }
        case (CHEM_POLYMER):
        {
            for (type = 3; type < RD_TYPES+1; type++) radii[type] = 1.5*MU;
            break;
        }
        case (CHEM_POLYMER_DISS):
        {
            for (type = 3; type < RD_TYPES+1; type++) radii[type] = 1.5*MU;
            break;
        }
        case (CHEM_POLYMER_STEP):
        {
            for (type = 2; type < RD_TYPES+1; type++) radii[type] = 1.2*MU;
            break;
        }
        case (CHEM_CATALYTIC_A2D): 
        {
            radii[3] = MU_B;
            radii[4] = MU*1.2;
            break;
        }
        case (CHEM_ABCAB): 
        {
            radii[3] = 1.5*MU;
            break;
        }
        case (CHEM_ABCDABC): 
        {
            radii[3] = 1.5*MU;
            radii[4] = 1.5*MU;
            break;
        }
        case (CHEM_BZ): 
        {
//             for (type = 2; type <= 8; type++) radii[type] = MU;
            for (type = 2; type <= 4; type++) radii[type] = MU;
            radii[5] = MU_B;
            radii[6] = 1.5*MU_B;
            radii[7] = 1.2*radii[3];
            radii[8] = 1.2*radii[5];
            break;
        }
        case (CHEM_BRUSSELATOR): 
        {
            radii[3] = MU;
            radii[4] = 1.2*MU_B;
            radii[5] = MU;
            radii[6] = MU;
            break;
        }
        case (CHEM_ABDACBE): 
        {
            radii[3] = MU;
            radii[4] = 1.2*MU_B;
            radii[5] = MU;
            break;
        }
    }
}

void adapt_speed_exothermic(t_particle *particle, double delta_ekin)
/* change the particle speed in case of exothermic reaction */
{
    double old_energy, e_ratio;
    
    old_energy = (particle->vx*(particle->vx) + particle->vy*(particle->vy))*(particle->mass_inv);
    if (old_energy + delta_ekin > 0.0) e_ratio = sqrt((old_energy + delta_ekin)/old_energy);
    else e_ratio = 0.0;
    particle->vx *= e_ratio;
    particle->vy *= e_ratio;
}

int chem_merge_AAB(int i, int newtype, t_particle particle[NMAXCIRCLES], t_collision *collisions, int ncollisions, double inv_masses[RD_TYPES+1], double radii[RD_TYPES+1])
/* merging particle i with a particle of same type into a particle of type newtype */
/* particular case of chem_merge */
{
    int j, k, type1;
    short int search = 1;
    double distance, rx, ry, mr, newmass_inv;
    
    type1 = particle[i].type;
    newmass_inv = inv_masses[newtype];
    mr = inv_masses[newtype]/inv_masses[type1];
    
    for (j=0; j<particle[i].hash_nneighb; j++) 
    {
        search = 1;
        k = particle[i].hashneighbour[j];
        if ((search)&&(particle[k].active)&&(particle[k].type == type1))
        {
            distance  = module2(particle[i].deltax[j], particle[i].deltay[j]);
            if ((distance < REACTION_DIST*MU)&&((double)rand()/RAND_MAX < REACTION_PROB))
            {
                particle[i].type = newtype;
                particle[i].radius = radii[newtype];
                rx = particle[i].xc - particle[k].xc;
                ry = particle[i].yc - particle[k].yc;
                particle[i].angle = argument(rx, ry);
                particle[i].omega = rx*particle[i].vy - ry*particle[i].vx;
                particle[i].omega -= rx*particle[k].vy - ry*particle[k].vx;
                if (CENTER_COLLIDED_PARTICLES)
                {
                    particle[i].xc = 0.5*(particle[i].xc + particle[k].xc);
                    particle[i].yc = 0.5*(particle[i].yc + particle[k].yc);
                }
                particle[i].vx = mr*(particle[i].vx + particle[k].vx);
                particle[i].vy = mr*(particle[i].vy + particle[k].vy);
                particle[i].mass_inv = newmass_inv;
                
                particle[k].active = 0;
                
                collisions[ncollisions].x = particle[i].xc;
                collisions[ncollisions].y = particle[i].yc;
                collisions[ncollisions].time = COLLISION_TIME;
                collisions[ncollisions].color = 0.0;
                
                if (ncollisions < NMAXCOLLISIONS - 1) ncollisions++;
                else printf("Too many collisions\n");
                
                search = 0;
            }
        }
    }
    return(ncollisions); 
}

int chem_merge(int i, int type2, int newtype, t_particle particle[NMAXCIRCLES], t_collision *collisions, int ncollisions, double inv_masses[RD_TYPES+1], double radii[RD_TYPES+1])
/* merging of particle i and a particle of type type2 into a particle of type newtype */
{
    int j, k, type1;
    short int search = 1;
    double distance, rx, ry, m2, mr1, mr2, newmass_inv, old_energy, e_ratio;
        
    type1 = particle[i].type;
    newmass_inv = inv_masses[newtype];
    mr1 = newmass_inv/inv_masses[type1];
    mr2 = newmass_inv/inv_masses[type2];
    
    for (j=0; j<particle[i].hash_nneighb; j++) 
    {
        k = particle[i].hashneighbour[j];
        if ((particle[k].active)&&(particle[k].type == type2))
        {
            distance  = module2(particle[i].deltax[j], particle[i].deltay[j]);
            if ((distance < REACTION_DIST*MU)&&((double)rand()/RAND_MAX < REACTION_PROB))
            {
                particle[i].type = newtype;
                particle[i].radius = radii[newtype];
                rx = particle[i].xc - particle[k].xc;
                ry = particle[i].yc - particle[k].yc;
                particle[i].angle = argument(rx, ry);
                particle[i].omega = rx*particle[i].vy - ry*particle[i].vx;
                particle[i].omega -= rx*particle[k].vy - ry*particle[k].vx;
                if (CENTER_COLLIDED_PARTICLES)
                {
                    particle[i].xc = 0.5*(particle[i].xc + particle[k].xc);
                    particle[i].yc = 0.5*(particle[i].yc + particle[k].yc);
                }
                particle[i].vx = mr1*particle[i].vx + mr2*particle[k].vx;
                particle[i].vy = mr1*particle[i].vy + mr2*particle[k].vy;
                particle[i].mass_inv = newmass_inv;
                
                if (EXOTHERMIC) adapt_speed_exothermic(&particle[i], DELTA_EKIN);
                
                particle[k].active = 0;
                
                collisions[ncollisions].x = particle[i].xc;
                collisions[ncollisions].y = particle[i].yc;
                collisions[ncollisions].time = COLLISION_TIME;
                collisions[ncollisions].color = 0.0;
                
                if (ncollisions < 2*NMAXCOLLISIONS - 1) ncollisions++;
                else printf("Too many collisions\n");
            }
        }
    }
    return(ncollisions); 
}


int chem_multi_merge(int i, int ni, int type2, int newtype, t_particle particle[NMAXCIRCLES], t_collision *collisions, int ncollisions, double inv_masses[RD_TYPES+1], double radii[RD_TYPES+1])
/* merging of ni particles of the same type as particle i (including particle i) */
/* and one particle of type type2 into a particle of type newtype */
{
    int j, k, type1, n1, n2, k1[10], k2;
    short int search = 1;
    double distance, xg, yg, rx, ry, rx1[10], ry1[10], rx2, ry2, mr1, mr2, newmass_inv;
        
    if (ni > 10)
    {
        printf("Error: need to increase size of k1 table in chem_multi_merge()\n");
        exit(1);
    }
    
    type1 = particle[i].type;
    newmass_inv = inv_masses[newtype];
    mr1 = newmass_inv/inv_masses[type1];
    mr2 = newmass_inv/inv_masses[type2];
    
    n1 = 0;
    n2 = 0;
    for (j=0; j<particle[i].hash_nneighb; j++) 
    {
        k = particle[i].hashneighbour[j];
        if ((particle[k].active)&&(particle[k].type == type1))
        {
            distance  = module2(particle[i].deltax[j], particle[i].deltay[j]);
            if ((distance < REACTION_DIST*MU)&&((double)rand()/RAND_MAX < REACTION_PROB))
            {
                k1[n1] = k;
                n1++;
            }
        }
        else if ((particle[k].active)&&(particle[k].type == type2))
        {
            distance  = module2(particle[i].deltax[j], particle[i].deltay[j]);
            if ((distance < REACTION_DIST*MU)&&((double)rand()/RAND_MAX < REACTION_PROB))
            {
                k2 = k;
                n2 = 1;
            }
        }
    }
                
    if ((n1 == ni-1)&&(n2 == 1))
    {
        particle[i].type = radii[newtype];
        particle[i].radius *= 1.5;
        xg = particle[i].xc;
        yg = particle[i].yc;
        for (n1 = 0; n1 < ni-1; n1++)
        {
            xg += particle[k1[n1]].xc;
            yg += particle[k1[n1]].yc;
        }
        xg = xg*mr1 + particle[k2].xc*mr2;
        yg = yg*mr1 + particle[k2].yc*mr2;
        
        rx = particle[i].xc - xg;
        ry = particle[i].yc - yg;
        for (n1 = 0; n1 < ni-1; n1++)
        {
            rx1[n1] = particle[k1[n1]].xc - xg;
            ry1[n1] = particle[k1[n1]].yc - yg;
        }
        rx2 = particle[k2].xc - xg;
        ry2 = particle[k2].yc - yg;
        particle[i].angle = argument(rx2, ry2);
        for (n1 = 0; n1 < ni-1; n1++) particle[i].angle += argument(rx1[n1], ry1[n1]);
        if (CENTER_COLLIDED_PARTICLES)
        {
            particle[i].xc = xg;
            particle[i].yc = yg;
        }
        particle[i].vx = mr1*particle[i].vx + mr2*particle[k2].vx;
        particle[i].vy = mr1*particle[i].vy + mr2*particle[k2].vy;
        for (n1 = 0; n1 < ni-1; n1++)
        {
            particle[i].vx += mr1*particle[k1[n1]].vx; 
            particle[i].vy += mr1*particle[k1[n1]].vx;
        }
        
        particle[i].omega = mr1*(rx*particle[i].vy - ry*particle[i].vx);
        for (n1 = 0; n1 < ni-1; n1++) 
            particle[i].omega += mr1*(rx1[n1]*particle[k1[n1]].vy - ry1[n1]*particle[k1[n1]].vx);
        particle[i].omega += mr2*(rx2*particle[k2].vy - ry2*particle[k2].vx);
        particle[i].mass_inv = newmass_inv;
                
        for (k==0; k<ni; k++) particle[k1[k]].active = 0;
                
        collisions[ncollisions].x = particle[i].xc;
        collisions[ncollisions].y = particle[i].yc;
        collisions[ncollisions].time = COLLISION_TIME;
        collisions[ncollisions].color = 0.0;
                            
        if (ncollisions < 2*NMAXCOLLISIONS - 1) ncollisions++;
        else printf("Too many collisions\n");
    }
            
    return(ncollisions); 
}


int chem_catalytic_merge(int i, int type_catalyst, int newtype, t_particle particle[NMAXCIRCLES], t_collision *collisions, int ncollisions, double inv_masses[RD_TYPES+1], double radii[RD_TYPES+1])
/* merging of 2 particles of the same type as particle i (including particle i) */
/* into a particle of type newtype, provided a particle of type type_catalyst is present */
{
    int j, k, type1, n1, n2, k1, k2;
    short int search = 1;
    double distance, xg, yg, rx, rx1, ry, ry1, mr, newmass_inv;
        
    type1 = particle[i].type;
    newmass_inv = inv_masses[newtype];
    mr = inv_masses[newtype]/inv_masses[type1];
    
    n1 = 0;
    n2 = 0;
    for (j=0; j<particle[i].hash_nneighb; j++) 
    {
        k = particle[i].hashneighbour[j];
        if ((particle[k].active)&&(particle[k].type == type1))
        {
            distance  = module2(particle[i].deltax[j], particle[i].deltay[j]);
            if ((distance < REACTION_DIST*MU)&&((double)rand()/RAND_MAX < REACTION_PROB))
            {
                k1 = k;
                n1 = 1;
            }
        }
        else if ((particle[k].active)&&(particle[k].type == type_catalyst))
        {
            distance  = module2(particle[i].deltax[j], particle[i].deltay[j]);
            if ((distance < REACTION_DIST*MU)&&((double)rand()/RAND_MAX < REACTION_PROB))
            {
                k2 = k;
                n2 = 1;
            }
        }
    }
                
    if ((n1 == 1)&&(n2 == 1))
    {
        particle[i].type = newtype;
        particle[i].radius  = radii[newtype];
        xg = 0.5*(particle[i].xc + particle[k1].xc + particle[k2].xc);
        yg = 0.5*(particle[i].yc + particle[k1].yc + particle[k2].yc);
        rx = particle[i].xc - xg;
        ry = particle[i].yc - yg;
        rx1 = particle[k1].xc - xg;
        ry1 = particle[k1].yc - yg;
        particle[i].angle = argument(rx1, ry1);
        if (CENTER_COLLIDED_PARTICLES)
        {
            particle[i].xc = xg;
            particle[i].yc = yg;
        }
        particle[i].vx = mr*(particle[i].vx + particle[k1].vx);
        particle[i].vy = mr*(particle[i].vy + particle[k1].vy);
        particle[i].omega = mr*(rx*particle[i].vy - ry*particle[i].vx);
        particle[i].omega += mr*(rx1*particle[k1].vy - ry1*particle[k1].vx);
        particle[i].mass_inv = newmass_inv;
                
        particle[k1].active = 0;
                
        collisions[ncollisions].x = particle[i].xc;
        collisions[ncollisions].y = particle[i].yc;
        collisions[ncollisions].time = COLLISION_TIME;
        collisions[ncollisions].color = 0.0;
                            
        if (ncollisions < 2*NMAXCOLLISIONS - 1) ncollisions++;
        else printf("Too many collisions\n");
    }
            
    return(ncollisions); 
}


int chem_split(int i, int newtype1, int newtype2, t_particle particle[NMAXCIRCLES], t_collision *collisions, int ncollisions, double inv_masses[RD_TYPES+1], double radii[RD_TYPES+1])
/* split particle i into particles of type newtype1 and newtype2 */
{
    int j, k, oldtype;
    short int success;
    double distance, rx, ry, xg, yg, normv;
    
    normv = module2(particle[i].vx, particle[i].vy); 
    rx = -MU*particle[i].vy/normv;
    ry = MU*particle[i].vx/normv;
                    
    xg = particle[i].xc;
    yg = particle[i].yc;
    particle[i].xc += rx;
    particle[i].yc += ry;
    oldtype = particle[i].type;
       
    /* test whether there is room to put two particles */
    success = add_particle(xg - rx, yg - ry, particle[i].vx, particle[i].vy, 1.0/inv_masses[newtype1], newtype1, particle);
    if (success)
    {
        particle[i].type = newtype2;
        particle[i].mass_inv = inv_masses[newtype2];
        particle[i].radius = radii[newtype2];
        
        particle[ncircles-1].type = newtype1;
        particle[ncircles-1].mass_inv = inv_masses[newtype1];
        particle[ncircles-1].radius = radii[newtype1];
        
        if (EXOTHERMIC) 
        {
            adapt_speed_exothermic(&particle[i], -0.5*DELTA_EKIN);
            adapt_speed_exothermic(&particle[ncircles-1], -0.5*DELTA_EKIN);
        }
            
        collisions[ncollisions].x = particle[i].xc;
        collisions[ncollisions].y = particle[i].yc;
        collisions[ncollisions].time = COLLISION_TIME;
        collisions[ncollisions].color = type_hue(oldtype);
            
        if (ncollisions < NMAXCOLLISIONS - 1) ncollisions++;
        else printf("Too many collisions\n");
    }
    else
    {
        particle[i].xc -= rx;
        particle[i].yc -= ry;
    }
    
    return(ncollisions);
}

int chem_transfer(int i, int type2, int newtype1, int newtype2, t_particle particle[NMAXCIRCLES], t_collision *collisions, int ncollisions, double inv_masses[RD_TYPES+1], double radii[RD_TYPES+1], double reac_dist)
/* reaction of particle i and a particle of type type2 into a particles of types newtype1 and newtype2 */
{
    int j, k, type1;
    short int search = 1;
    double distance, rx, ry, mass_inv1, mass_inv2, newmass_inv1, newmass_inv2, old_energy, e_ratio, omega;
        
    type1 = particle[i].type;
    mass_inv1 = inv_masses[type1];
    mass_inv2 = inv_masses[type2];
    newmass_inv1 = inv_masses[newtype1];
    newmass_inv2 = inv_masses[newtype2];
    
    for (j=0; j<particle[i].hash_nneighb; j++) 
    {
        k = particle[i].hashneighbour[j];
        if ((particle[k].active)&&(particle[k].type == type2))
        {
            distance  = module2(particle[i].deltax[j], particle[i].deltay[j]);
            if ((distance < reac_dist*MU)&&((double)rand()/RAND_MAX < REACTION_PROB))
            {
                particle[i].type = newtype1;
                particle[i].radius = radii[newtype1];
                particle[k].type = newtype2;
                particle[k].radius = radii[newtype2];
                
                /* momentum conservation does not determine outgoing velocities completely */
                /* so make some arbitrary/random choices here */
                rx = particle[i].xc - particle[k].xc;
                ry = particle[i].yc - particle[k].yc;
                particle[i].angle = argument(rx, ry);
                particle[k].angle = argument(rx, ry) + PI;
                
                omega = gaussian()*OMEGA_INITIAL;
                particle[i].omega = omega;
                particle[k].omega = -omega;
                
                particle[i].vx *= newmass_inv1/mass_inv1;
                particle[i].vy *= newmass_inv1/mass_inv1;
                particle[i].mass_inv = newmass_inv1;
                particle[k].vx *= newmass_inv2/mass_inv2;
                particle[k].vy *= newmass_inv2/mass_inv2;
                particle[k].mass_inv = newmass_inv2;
                
                if (EXOTHERMIC) 
                {
                    adapt_speed_exothermic(&particle[i], DELTA_EKIN);
                    adapt_speed_exothermic(&particle[k], DELTA_EKIN);
                }
                
                collisions[ncollisions].x = 0.5*(particle[i].xc + particle[k].xc);
                collisions[ncollisions].y = 0.5*(particle[i].yc + particle[k].yc);
                collisions[ncollisions].time = COLLISION_TIME;
                collisions[ncollisions].color = 0.0;
                
                if (ncollisions < 2*NMAXCOLLISIONS - 1) ncollisions++;
                else printf("Too many collisions\n");
            }
        }
    }
    return(ncollisions); 
}

int chem_convert(int i, int newtype, t_particle particle[NMAXCIRCLES], t_collision *collisions, int ncollisions, double inv_masses[RD_TYPES+1], double radii[RD_TYPES+1], double reac_prob)
/* convert particle i into new type with probability reac_prob */
{
    int oldtype;
    
    oldtype = particle[i].type;
//     if (oldtype == 5) printf("Converting type %i to type %i\n", oldtype, newtype); 
    
    if ((double)rand()/RAND_MAX < reac_prob)
    {
        particle[i].type = newtype;
        particle[i].radius = radii[newtype];
        particle[i].mass_inv = inv_masses[newtype];
        
        collisions[ncollisions].x = particle[i].xc;
        collisions[ncollisions].y = particle[i].yc;
        collisions[ncollisions].time = COLLISION_TIME;
        collisions[ncollisions].color = type_hue(oldtype);
                
        if (ncollisions < 2*NMAXCOLLISIONS - 1) ncollisions++;
        else printf("Too many collisions\n");
    }
    return(ncollisions); 
}

int chem_catalytic_convert(int i, int type2, int newtype, t_particle particle[NMAXCIRCLES], t_collision *collisions, int ncollisions, double inv_masses[RD_TYPES+1], double radii[RD_TYPES+1], double reaction_dist, double reaction_prob)
/* if 2 particles of the same type as particle i (including particle i) are close  */
/* to a particle of type type2, convert particle type from type2 to newtype */
{
    int j, k, type1, n1, n2, k1, k2, oldtype;
    short int search = 1;
    double distance;
        
    type1 = particle[i].type;
    
    oldtype = particle[i].type;
    
    n1 = 0;
    n2 = 0;
    for (j=0; j<particle[i].hash_nneighb; j++) 
    {
        k = particle[i].hashneighbour[j];
        if ((particle[k].active)&&(particle[k].type == type1))
        {
            distance  = module2(particle[i].deltax[j], particle[i].deltay[j]);
            if ((distance < reaction_dist*MU)&&((double)rand()/RAND_MAX < reaction_prob))
            {
                k1 = k;
                n1 = 1;
            }
        }
        else if ((particle[k].active)&&(particle[k].type == type2))
        {
            distance  = module2(particle[i].deltax[j], particle[i].deltay[j]);
            if ((distance < reaction_dist*MU)&&((double)rand()/RAND_MAX < reaction_prob))
            {
                k2 = k;
                n2 = 1;
            }
        }
    }
                
    if ((n1 == 1)&&(n2 == 1))
    {
        particle[k2].type = newtype;
        particle[k2].radius  = radii[newtype];
        particle[k2].mass_inv = inv_masses[newtype];
                
        collisions[ncollisions].x = particle[i].xc;
        collisions[ncollisions].y = particle[i].yc;
        collisions[ncollisions].time = COLLISION_TIME;
        collisions[ncollisions].color = type_hue(oldtype);
                            
        if (ncollisions < 2*NMAXCOLLISIONS - 1) ncollisions++;
        else printf("Too many collisions\n");
    }
            
    return(ncollisions); 
}

int update_types(t_particle particle[NMAXCIRCLES], t_collision *collisions, int ncollisions, int *particle_numbers, int time, double *delta_e)
/* update the types in case of reaction-diffusion equation */
{
    int i, j, k, n, n3, n4, type, oldncollisions, delta_n;
    double distance, rnd, p1, p2;
    static double inv_masses[RD_TYPES+1], radii[RD_TYPES+1];
    static int first = 1;
    
    if (first)  /* compute total mass and mass ratios */
    {
        compute_inverse_masses(inv_masses);
        compute_radii(radii);
        first = 0;
    }
    
    if (EXOTHERMIC) oldncollisions = ncollisions;
    
    switch (RD_REACTION) {
        case (CHEM_RPS): 
        {
            for (i=0; i<ncircles; i++)
                for (j=0; j<particle[i].hash_nneighb; j++)
                {
                    k = particle[i].hashneighbour[j];
                    if ((particle[k].type == particle[i].type + 1)||((particle[i].type == RD_TYPES)&&(particle[k].type == 1)))
                    {
                        distance  = module2(particle[i].deltax[j], particle[i].deltay[j]);
                        if ((distance < EQUILIBRIUM_DIST)&&((double)rand()/RAND_MAX < REACTION_PROB))
                            particle[k].type = particle[i].type;
//                         printf("Changed particle type to %i\n", particle[k].type);
                    }
                }
            return(0);
        }
        case (CHEM_AAB):
        {
            for (i=0; i<ncircles; i++) 
                if ((particle[i].active)&&(particle[i].type == 1))
                    ncollisions = chem_merge_AAB(i, 2, particle, collisions, ncollisions, inv_masses, radii);
            return(ncollisions);
        }
        case (CHEM_ABC):
        {
            for (i=0; i<ncircles; i++) if ((particle[i].active)&&(particle[i].type == 1))
                ncollisions = chem_merge(i, 2, 3, particle, collisions, ncollisions, inv_masses, radii);
            printf("%i collisions\n", ncollisions);
            delta_n = ncollisions - oldncollisions; 
            printf("delta_n = %i\n", delta_n);
            if (EXOTHERMIC) *delta_e = (double)(delta_n)*DELTA_EKIN;
            return(ncollisions);
        }
        case (CHEM_A2BC):
        {
            for (i=0; i<ncircles; i++) if ((particle[i].active)&&(particle[i].type == 1))
                ncollisions = chem_multi_merge(i, 2, 2, 3, particle, collisions, ncollisions, inv_masses, radii);
            printf("%i collisions\n", ncollisions);
            return(ncollisions);
        }
        case (CHEM_CATALYSIS):
        {
            for (i=0; i<ncircles; i++) if ((particle[i].active)&&(particle[i].type == 1))
                ncollisions = chem_catalytic_merge(i, 2, 3, particle, collisions, ncollisions, inv_masses, radii);
            printf("%i collisions\n", ncollisions);
            return(ncollisions);
        }
        case (CHEM_AUTOCATALYSIS):
        {
            for (i=0; i<ncircles; i++) if ((particle[i].active)&&(particle[i].type == 1))
                ncollisions = chem_catalytic_merge(i, 2, 2, particle, collisions, ncollisions, inv_masses, radii);
            printf("%i collisions\n", ncollisions);
            return(ncollisions);
        }
        case (CHEM_BAA):
        {
            for (i=0; i<ncircles; i++) 
                if ((particle[i].active)&&(particle[i].type == 2)&&((double)rand()/RAND_MAX < DISSOCIATION_PROB))
                    ncollisions = chem_split(i, 1, 1, particle, collisions, ncollisions, inv_masses, radii);
            printf("%i collisions\n", ncollisions);
            return(ncollisions);
        }
        case (CHEM_AABAA):
        {
            *delta_e = 0.0;
            for (i=0; i<ncircles; i++)
            {
                oldncollisions = ncollisions;
                if ((particle[i].active)&&(particle[i].type == 1))
                {
                    ncollisions = chem_merge_AAB(i, 2, particle, collisions, ncollisions, inv_masses, radii);
                    if (EXOTHERMIC) *delta_e += (double)(ncollisions - oldncollisions)*DELTA_EKIN;
                }
                else if ((particle[i].active)&&(particle[i].type == 2)&&((double)rand()/RAND_MAX < DISSOCIATION_PROB))
                {
                    ncollisions = chem_split(i, 1, 1, particle, collisions, ncollisions, inv_masses, radii);
                    if (EXOTHERMIC) *delta_e -= (double)(ncollisions - oldncollisions)*DELTA_EKIN;
                }
            }

            printf("Delta_E = %.3lg\n", *delta_e);
            printf("%i collisions\n", ncollisions);
            return(ncollisions);
        }
        case (CHEM_POLYMER):
        {
            for (i=0; i<ncircles; i++)  if ((particle[i].active)&&(particle[i].type == 1))
                for (k=2; k<RD_TYPES; k++)
                    ncollisions = chem_merge(i, k, k+1, particle, collisions, ncollisions, inv_masses, radii);   

            printf("%i collisions\n", ncollisions);
            return(ncollisions);            
        }
        case (CHEM_POLYMER_DISS):
        {
            for (i=0; i<ncircles; i++)  if (particle[i].active)
            {
                if (particle[i].type == 1) for (k=2; k<RD_TYPES; k++)
                    ncollisions = chem_merge(i, k, k+1, particle, collisions, ncollisions, inv_masses, radii);
                else if ((particle[i].type > 2)&&((double)rand()/RAND_MAX < DISSOCIATION_PROB)) 
                    ncollisions = chem_split(i, 1, particle[i].type - 1, particle, collisions, ncollisions, inv_masses, radii);
            }
            
            printf("%i collisions\n", ncollisions);
            return(ncollisions);            
        }
        case (CHEM_POLYMER_STEP):
        {
            for (i=0; i<ncircles; i++)  if (particle[i].active)
            {
                type = particle[i].type;
                for (k=1; k<RD_TYPES+1-type; k++)
                    ncollisions = chem_merge(i, k, k+type, particle, collisions, ncollisions, inv_masses, radii);
                if ((type > 1)&&((double)rand()/RAND_MAX < DISSOCIATION_PROB)) 
                {
                    k = rand()%(type-1) + 1;
//                     printf("Splitting type %i into type %i and type %i\n", type, k, type - k);
                    ncollisions = chem_split(i, k, type - k, particle, collisions, ncollisions, inv_masses, radii);
                }
            }
            
            printf("%i collisions\n", ncollisions);
            return(ncollisions);            
        }
        case (CHEM_CATALYTIC_A2D):
        {
            for (i=0; i<ncircles; i++) if ((particle[i].active)&&(particle[i].type == 1))
            {
                ncollisions = chem_merge(i, 2, 3, particle, collisions, ncollisions, inv_masses, radii);
                ncollisions = chem_transfer(i, 3, 4, 2, particle, collisions, ncollisions, inv_masses, radii, 1.05*REACTION_DIST);
            }
            printf("%i collisions\n", ncollisions);
            return(ncollisions);
        }
        case (CHEM_ABCAB):
        {
            *delta_e = 0.0;
            for (i=0; i<ncircles; i++)
            {
                oldncollisions = ncollisions;
                if ((particle[i].active)&&(particle[i].type == 1))
                {
                    ncollisions = chem_merge(i, 2, 3, particle, collisions, ncollisions, inv_masses, radii);
                    if (EXOTHERMIC) *delta_e += (double)(ncollisions - oldncollisions)*DELTA_EKIN;
                }
                else if ((particle[i].active)&&(particle[i].type == 3)&&((double)rand()/RAND_MAX < DISSOCIATION_PROB))
                {
                    ncollisions = chem_split(i, 1, 2, particle, collisions, ncollisions, inv_masses, radii);
                    if (EXOTHERMIC) *delta_e -= (double)(ncollisions - oldncollisions)*DELTA_EKIN;
                }
            }
            printf("Delta_E = %.3lg\n", *delta_e);
            printf("%i collisions\n", ncollisions);
            return(ncollisions);
        }
        case (CHEM_ABCDABC):
        {
            *delta_e = 0.0;
            for (i=0; i<ncircles; i++)
            {
                oldncollisions = ncollisions;
                if ((particle[i].active)&&(particle[i].type == 1))
                {
                    ncollisions = chem_merge(i, 2, 3, particle, collisions, ncollisions, inv_masses, radii);
                    if (EXOTHERMIC) *delta_e += (double)(ncollisions - oldncollisions)*DELTA_EKIN;
                    
                    ncollisions = chem_merge(i, 3, 4, particle, collisions, ncollisions, inv_masses, radii);
                    if (EXOTHERMIC) *delta_e += (double)(ncollisions - oldncollisions)*DELTA_EKIN;
                }
                else if ((particle[i].active)&&(particle[i].type == 3)&&((double)rand()/RAND_MAX < DISSOCIATION_PROB))
                {
                    ncollisions = chem_split(i, 1, 2, particle, collisions, ncollisions, inv_masses, radii);
                    if (EXOTHERMIC) *delta_e -= (double)(ncollisions - oldncollisions)*DELTA_EKIN;
                }
                else if ((particle[i].active)&&(particle[i].type == 4)&&((double)rand()/RAND_MAX < DISSOCIATION_PROB))
                {
                    ncollisions = chem_split(i, 1, 3, particle, collisions, ncollisions, inv_masses, radii);
                    if (EXOTHERMIC) *delta_e -= (double)(ncollisions - oldncollisions)*DELTA_EKIN;
                }
            }
            printf("Delta_E = %.3lg\n", *delta_e);
            printf("%i collisions\n", ncollisions);
            return(ncollisions);
        }
        case (CHEM_BZ):
        {
            for (i=0; i<ncircles; i++) if (particle[i].active)
            {
                oldncollisions = ncollisions;
                switch (particle[i].type) {
                    case (1): 
                    {
                        ncollisions = chem_transfer(i, 4, 2, 3, particle, collisions, ncollisions, inv_masses, radii, REACTION_DIST);
                        ncollisions = chem_transfer(i, 3, 2, 2, particle, collisions, ncollisions, inv_masses, radii, REACTION_DIST);
                        break;
                    }
                    case (2):
                    {
//                         particle[i].active = 0;
                        ncollisions = chem_transfer(i, 2, 1, 3, particle, collisions, ncollisions, inv_masses, radii, 0.95*REACTION_DIST);
                        break;
                    }
                    case (3):
                    {
                        ncollisions = chem_transfer(i, 3, 2, 4, particle, collisions, ncollisions, inv_masses, radii, REACTION_DIST);
                        ncollisions = chem_transfer(i, 4, 7, 8, particle, collisions, ncollisions, inv_masses, radii, REACTION_DIST);
                        break;
                    }
                    case (5):
                    {
                        rnd = (double)rand()/RAND_MAX;
                        if (rnd < 0.2)
                            ncollisions = chem_merge(i, 6, 1, particle, collisions, ncollisions, inv_masses, radii);
                        else /*if (rnd < 0.9)*/
                            ncollisions = chem_merge(i, 6, 6, particle, collisions, ncollisions, inv_masses, radii);
//                         else if (rnd < 0.6)
//                             ncollisions = chem_transfer(i, 6, 1, 1, particle, collisions, ncollisions, inv_masses, radii, REACTION_DIST);
                        break;
                    }
                    case (7):
                    {
                        ncollisions = chem_split(i, 3, 3, particle, collisions, ncollisions, inv_masses, radii);
                        break;
                    }
                    case (8):
                    {
                        ncollisions = chem_split(i, 5, 5, particle, collisions, ncollisions, inv_masses, radii);
                        break;
                    }
                }
            }
            return(ncollisions);
        }
        case (CHEM_BRUSSELATOR):
        {
            for (i=0; i<ncircles; i++) if (particle[i].active)
            {
                oldncollisions = ncollisions;
                switch (particle[i].type) {
                    case (1): /* A -> X */
                    {
                        ncollisions = chem_convert(i, 3, particle, collisions, ncollisions, inv_masses, radii, DISSOCIATION_PROB);
                        break;
                    }
                    case (2): /* B + X -> Y + D */
                    {
                        ncollisions = chem_transfer(i, 3, 4, 5, particle, collisions, ncollisions, inv_masses, radii, REACTION_DIST);
                        break;
                    }
                    case (3): 
                    {
                        /* X -> B, modified from Brusselator */
                        ncollisions = chem_convert(i, 2, particle, collisions, ncollisions, inv_masses, radii, 0.5*DISSOCIATION_PROB);
                        /* 2X + Y -> 3X */
                        ncollisions = chem_catalytic_convert(i, 4, 3, particle, collisions, ncollisions, inv_masses, radii, 3.0*REACTION_DIST, 1.0);
                        break;
                    }
                    case (5): /* D -> A or B if concentration of X or Y is small */ 
                    {
                        n = time*(RD_TYPES+1);
                        n3 = particle_numbers[n+3];
                        n4 = particle_numbers[n+4];
                        p1 = 1.0/((double)(n3*n3/10+1));
                        p2 = 5.0*DISSOCIATION_PROB + 1.0/((double)(n4*n4/200+1));
                        rnd = (double)rand()/RAND_MAX;
                        if (rnd < p1/(p1+p2))
                            ncollisions = chem_convert(i, 1, particle, collisions, ncollisions, inv_masses, radii, p1);
                        else /*if (n4 < 1000)*/
                            ncollisions = chem_convert(i, 2, particle, collisions, ncollisions, inv_masses, radii, p2);
                    }
                }
            }
            return(ncollisions);
        }
        case (CHEM_ABDACBE):
        {
            for (i=0; i<ncircles; i++) 
            {
//                 oldncollisions = ncollisions;
                if ((particle[i].active)&&(particle[i].type == 1))
                {
                    ncollisions = chem_merge(i, 2, 4, particle, collisions, ncollisions, inv_masses, radii);
//                     if (EXOTHERMIC) *delta_e += (double)(ncollisions - oldncollisions)*DELTA_EKIN;
                    
                    ncollisions = chem_transfer(i, 3, 2, 5, particle, collisions, ncollisions, inv_masses, radii, REACTION_DIST);
//                     if (EXOTHERMIC) *delta_e += (double)(ncollisions - oldncollisions)*DELTA_EKIN;
                }
            }
            printf("%i collisions\n", ncollisions);
            delta_n = ncollisions - oldncollisions; 
            printf("delta_n = %i\n", delta_n);
            if (EXOTHERMIC) *delta_e = (double)(delta_n)*DELTA_EKIN;
            return(ncollisions);
        }
    }
}
    
double plot_coord(double x, double xmin, double xmax)
{
    return(xmin + x*(xmax - xmin));
}


void draw_speed_plot(t_group_data *group_speeds, int i)
/* draw plot of obstacle speeds as a function of time */
{
    int j, group;
    char message[100];
    static double xmin, xmax, ymin, ymax, xmid, ymid, dx, dy, plotxmin, plotxmax, plotymin, plotymax;
    double pos[2], x1, y1, x2, y2, rgb[3];
    static int first = 1, gshift = INITIAL_TIME + NSTEPS;
    
    if (first)
    {
//         xmin = XMAX - 1.35;
//         xmax = XMAX - 0.05;
//         ymin = YMAX - 1.8;
//         ymax = YMAX - 0.5;
                
        xmin = XMAX - 1.8;
        xmax = XMAX - 0.066;
        ymin = YMAX - 2.5;
        ymax = YMAX - 0.76;

        xmid = 0.5*(xmin + xmax);
        ymid = 0.5*(ymin + ymax);
        
        dx = 0.5*(xmax - xmin);
        dy = 0.5*(ymax - ymin);
        
        plotxmin = xmin + 0.05;
        plotxmax = xmax - 0.1;
        plotymin = ymin + 0.07;
        plotymax = ymax - 0.15;
        
        first = 0;
    }
    
//     rgb[0] = 1.0; rgb[1] = 1.0; rgb[2] = 1.0;
    
    glLineWidth(2);
    
    /* plot angular speed */
    for (group=1; group<ngroups; group++)
    {
        set_segment_group_color(group, 0.5, rgb);
        x1 = plotxmin;
        y1 = plotymin;
        for (j=0; j<i; j++)
        {
            x2 = plot_coord((double)j/(double)NSTEPS, plotxmin, plotxmax);
            y2 = plot_coord(group_speeds[(group-1)*gshift + j].omega/VMAX_PLOT_SPEEDS, plotymin, plotymax);
        
            draw_line(x1, y1, x2, y2);
            x1 = x2;
            y1 = y2;
        }
    
        sprintf(message, "omega");
        write_text_fixedwidth(plotxmin - 0.22 + (double)(group-1)*0.3, plotymax + 0.16, message);
    }

    /* plot speed of obstacles */
    for (group=1; group<ngroups; group++)
    {
        set_segment_group_color(group, 1.0, rgb);
        x1 = plotxmin;
        y1 = plotymin;
        for (j=0; j<i; j++)
        {
            x2 = plot_coord((double)j/(double)NSTEPS, plotxmin, plotxmax);
            y2 = plot_coord(group_speeds[(group-1)*gshift + j].vy/VMAX_PLOT_SPEEDS, plotymin, plotymax);
        
            draw_line(x1, y1, x2, y2);
            x1 = x2;
            y1 = y2;
        }
        
        sprintf(message, "vy");
        write_text_fixedwidth(plotxmin - 0.22 + (double)(group-1)*0.3, plotymax + 0.25, message);
    }
    
    glColor3f(1.0, 1.0, 1.0);
    
    /* axes and labels */
    draw_line(plotxmin, plotymin, plotxmax + 0.05, plotymin);
    draw_line(plotxmin, plotymin, plotxmin, plotymax + 0.1);
    
    for (j=1; j<=(int)(10.0*VMAX_PLOT_SPEEDS); j++)
    {
        y1 = plot_coord((double)j/(10.0*VMAX_PLOT_SPEEDS), plotymin, plotymax);
        draw_line(plotxmin - 0.02, y1, plotxmin + 0.02, y1);
    }
    
    sprintf(message, "%.1f", VMAX_PLOT_SPEEDS);
    write_text_fixedwidth(plotxmin - 0.28, y1 - 0.025, message);
//     write_text_fixedwidth(plotxmin - 0.22, y1 - 0.025, message);
        
    sprintf(message, "time");
    write_text_fixedwidth(plotxmax - 0.13, plotymin - 0.12, message);
//     write_text_fixedwidth(plotxmax - 0.1, plotymin - 0.08, message);
}


void draw_trajectory_plot(t_group_data *group_speeds, int i)
/* draw plot of obstacle speeds as a function of time */
{
    int j, group;
    char message[100];
    static double xmin, xmax, ymin, ymax, xmid, ymid, dx, dy, plotxmin, plotxmax, plotymin, plotymax, scalex, scaley, yinitial[NMAXGROUPS];
    double pos[2], x0, y0, x1, y1, x2, y2, rgb[3];
    static int first = 1, gshift = INITIAL_TIME + NSTEPS;
    
    if (first)
    {
        xmin = XMAX - 1.8;
        xmax = XMAX - 0.066;
        ymin = YMIN + 0.1;
        ymax = YMIN + 1.9;

        xmid = 0.5*(xmin + xmax);
        ymid = 0.5*(ymin + ymax);
        
        dx = 0.5*(xmax - xmin);
        dy = 0.5*(ymax - ymin);
        
        plotxmin = xmin + 0.05;
        plotxmax = xmax - 0.1;
        plotymin = ymin + 0.07;
        plotymax = ymax - 0.15;
        
        scalex = 6.0;
        scaley = 3.0;
        
        for (group = 1; group < ngroups; group++)
            yinitial[group] = group_speeds[(group-1)*gshift].yc;
        
        first = 0;
    }
    
    if (ALTITUDE_LINES)
    {
        glLineWidth(1);
        glColor3f(0.5, 0.5, 0.5);
        for (j=1; j<(int)scaley + 1; j++)
        {
            y1 = plot_coord((double)j/scaley, plotymin, plotymax);
            if (y1 < plotymax) draw_line(plotxmin, y1, plotxmax + 0.05, y1);
        }
    }
    
    glLineWidth(2);
    
    printf("ngroups = %i\n", ngroups);
    
    /* plot trajectories */
    for (group=1; group<ngroups; group++)
    {
        set_segment_group_color(group, 0.75, rgb);
        x1 = group_speeds[(group-1)*gshift].xc;
        y1 = group_speeds[(group-1)*gshift].yc - yinitial[group];
            
        for (j=0; j<i-1; j++)
        {
            x0 = group_speeds[(group-1)*gshift + j].xc;
            y0 = group_speeds[(group-1)*gshift + j].yc - yinitial[group];
            
            if (y0 > scaley) scaley = y0;
            
            x2 = plot_coord(0.5 + x0/scalex, plotxmin, plotxmax);
            y2 = plot_coord(y0/scaley + 0.05, plotymin, plotymax);
            
//             printf("yinitial = %.3lg, (x0, y0) = (%.3lg, %.3lg), (x2, y2) = (%.3lg, %.3lg)\n", yinitial, x0, y0, x2, y2); 
        
            if ((j>0)&&(module2(x2-x1, y2-y1) < 0.1)) draw_line(x1, y1, x2, y2);
            x1 = x2;
            y1 = y2;
        }
        if (i>0) draw_colored_circle_precomp(x1, y1, 0.015, rgb);
    }

    glColor3f(1.0, 1.0, 1.0);
    
    /* axes and labels */
    draw_line(plotxmin, plotymin, plotxmax + 0.05, plotymin);
    draw_line(plotxmin, plotymin, plotxmin, plotymax + 0.1);
    
    for (j=1; j<(int)scaley; j++)
    {
        y1 = plot_coord((double)j/scaley, plotymin, plotymax);
        draw_line(plotxmin - 0.02, y1, plotxmin + 0.02, y1);
    }
    
    sprintf(message, "%i", (int)scaley - 1);
    write_text_fixedwidth(plotxmin - 0.15, y1 - 0.025, message);
    
//     sprintf(message, "%.1f", VMAX_PLOT_SPEEDS);
    sprintf(message, "y");
    write_text_fixedwidth(plotxmin - 0.1, plotymax - 0.005, message);
        
    sprintf(message, "x");
    write_text_fixedwidth(plotxmax - 0.1, plotymin - 0.1, message);
}


void write_size_text(double x, double y, char *message, int largewin)
{
    if (largewin) write_text(x, y, message);
    else write_text_fixedwidth(x, y, message);                     
}

void draw_particle_nb_plot(int *particle_numbers, int i)
/* draw plot of number of particles as a function of time */
{
    int j, type, power;
    char message[100];
    static double xmin, xmax, ymin, ymax, xmid, ymid, dx, dy, plotxmin, plotxmax, plotymin, plotymax;
    double pos[2], x1, y1, x2, y2, rgb[3];
    static int first = 1, gshift = INITIAL_TIME + NSTEPS, nmax, ygrad, largewin;
    
    if (first)
    {                
        xmin = XMAX - 1.05;
        xmax = XMAX - 0.05;
        ymin = YMAX - 1.0;
        ymax = YMAX - 0.05;

        xmid = 0.5*(xmin + xmax);
        ymid = 0.5*(ymin + ymax);
        
        dx = 0.5*(xmax - xmin);
        dy = 0.5*(ymax - ymin);
        
        plotxmin = xmin + 0.18;
        plotxmax = xmax - 0.1;
        plotymin = ymin + 0.07;
        plotymax = ymax - 0.15;
        
        nmax = 0;
        for (j=1; j<RD_TYPES+1; j++) 
            if (particle_numbers[RD_TYPES+1+j] > nmax) nmax = particle_numbers[RD_TYPES+1+j];
        
        nmax *= PARTICLE_NB_PLOT_FACTOR;
            
        power = (int)(log((double)nmax)/log(10.0)) - 1.0;
        ygrad = (int)ipow(10.0, power);
        
        largewin = (WINWIDTH > 1280);
        
        first = 0;
    }
    
    erase_area_hsl(xmid, ymid, dx, dy, 0.0, 0.9, 0.0);
    
    glLineWidth(2);
    
    /* plot particle number */
    for (type=1; type<RD_TYPES+1; type++)
    {
        set_type_color(type, 1.0, rgb);
        x1 = plotxmin;
        y1 = plotymin;
        glBegin(GL_LINE_STRIP);
        for (j=1; j<i; j++)
        {
            x2 = plot_coord((double)j/(double)NSTEPS, plotxmin, plotxmax);
            y2 = plot_coord(particle_numbers[(RD_TYPES+1)*j+type]/(double)nmax, plotymin, plotymax);
            glVertex2d(x2, y2);
//             draw_line(x1, y1, x2, y2);
//             x1 = x2;
//             y1 = y2;
        }
        glEnd(); 
        
        sprintf(message, "%5d molecules", particle_numbers[(RD_TYPES+1)*i+type]);
//         write_size_text(xmax - 0.5, ymax - 0.04 - 0.06*(double)type, message, largewin);
//         if (largewin) write_text(xmax - 0.5, ymax - 0.1 - 0.1*(double)type, message);
//         else 
        write_text_fixedwidth(xmax - 0.5, ymax - 0.06*(double)type, message);
    }
    
    glColor3f(1.0, 1.0, 1.0);
    
    /* axes and labels */
    draw_line(plotxmin, plotymin, plotxmax + 0.05, plotymin);
    draw_line(plotxmin, plotymin, plotxmin, plotymax + 0.1);
    
    glLineWidth(1);
    j = 1;
    while (j*ygrad <= nmax + 5)
    {
        y1 = plot_coord((double)(j*ygrad)/(double)nmax, plotymin, plotymax);
        draw_line(plotxmin - 0.015, y1, plotxmin + 0.015, y1);
        
        if (j%10 == 0)
        {
            draw_line(plotxmin - 0.025, y1, plotxmin + 0.025, y1);
            sprintf(message, "%i", j*ygrad);
//             write_size_text(plotxmin - 0.15, y1 - 0.015, message, largewin);
//             if (largewin) write_text(plotxmin - 0.17, y1 - 0.015, message);
//             else 
            write_text_fixedwidth(plotxmin - 0.15, y1 - 0.015, message);
        }
        else if ((j%5 == 0)&&(nmax<5*ygrad))
        {
            sprintf(message, "%i", j*ygrad);
//             if (j*ygrad >= 1000) write_size_text(plotxmin - 0.15, y1 - 0.015, message, largewin);
//             else write_size_text(plotxmin - 0.12, y1 - 0.015, message, largewin);
            if (j*ygrad >= 1000) write_text_fixedwidth(plotxmin - 0.15, y1 - 0.015, message);
            else write_text_fixedwidth(plotxmin - 0.12, y1 - 0.015, message);
        }
        j++;
    }
    
    sprintf(message, "time");
//     write_size_text(plotxmax - 0.1, plotymin - 0.05, message, largewin);
    write_text_fixedwidth(plotxmax - 0.1, plotymin - 0.05, message);
}

void init_segment_group(t_segment segment[NMAXSEGMENTS], int group, t_group_segments segment_group[NMAXGROUPS])
/* initialize center of mass and similar data of grouped segments */
{
    int i, nseg_group = 0;
    double xc = 0.0, yc = 0.0;
    
    for (i=0; i<nsegments; i++) if (segment[i].group == group) 
    {
        xc += 0.5*(segment[i].x1 + segment[i].x2);
        yc += 0.5*(segment[i].y1 + segment[i].y2);
        
//         printf("Segment group data %i: z1 = (%.3lg, %.3lg)\n z2 = (%.3lg, %.3lg)\n zc = (%.3lg, %.3lg)\n\n", group, 
//                segment[i].x1, segment[i].y1, segment[i].x2, segment[i].y2, xc, yc);
        
        nseg_group++;
    }
    
    if (nseg_group == 0) nseg_group = 1;
    
    segment_group[group].xc = xc/(double)nseg_group;
    segment_group[group].yc = yc/(double)nseg_group;
    segment_group[group].angle = 0.0;
    segment_group[group].vx = 0.0;
    segment_group[group].vy = 0.0;
    segment_group[group].omega = 0.0;
    segment_group[group].mass = SEGMENT_GROUP_MASS;
    segment_group[group].moment_inertia = SEGMENT_GROUP_I;
    
    printf("Segment group data %i: (%.3lg, %.3lg)\n", group, segment_group[group].xc, segment_group[group].yc);
    
}


void reset_energy(t_particle particle[NMAXCIRCLES], double px[NMAXCIRCLES], double py[NMAXCIRCLES], double totalenergy, double emean)
/* decrease energy in case of blow-up */
{
    int i;
    double vratio, emax;
    
    emax = 10.0*emean/(double)ncircles;
    
//     printf("Warning: blow-up, resetting energy from %.5lg to %.5lg\n\n", totalenergy, emean);
    printf("Warning: blow-up, resetting energy of some particles\n");
    
//     vratio = sqrt(emean/totalenergy);
    
    for (i=0; i < ncircles; i++) if (particle[i].energy > emax)
    {
        printf("Particle %i at (%.3lg, %.3lg) has energy %.5lg, resetting to 0\n", i, particle[i].xc, particle[i].yc, particle[i].energy);
        particle[i].vx = 0.0;
        particle[i].vy = 0.0;
        px[i] = 0.0;
        py[i] = 0.0;
        particle[i].energy = 0.0;
    }
    printf("\n\n"); 
    
}
