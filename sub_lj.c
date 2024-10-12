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
  TIFFSetField(file, TIFFTAG_IMAGEWIDTH, (uint32_t) width);
  TIFFSetField(file, TIFFTAG_IMAGELENGTH, (uint32_t) height);
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

void erase_rect_rgb(double xmin, double ymin, double xmax, double ymax, double rgb[3])
{
    double pos[2];
    
    glColor3f(rgb[0], rgb[1], rgb[2]);
    glBegin(GL_QUADS);
    glVertex2d(xmin, ymin);
    glVertex2d(xmax, ymin);
    glVertex2d(xmax, ymax);
    glVertex2d(xmin, ymax);
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

void erase_rectangle_hsl_turbo(double xmin, double ymin, double xmax, double ymax, double h, double s, double l)
{
    double pos[2], rgb[3];
    
    hsl_to_rgb_turbo(h, s, l, rgb);
    erase_rect_rgb(xmin, ymin, xmax, ymax, rgb);
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

void draw_colored_triangle_turbo(double x1, double y1, double x2, double y2, double x3, double y3, double h, double s, double l)
{    
    double rgb[3];
    
    hsl_to_rgb_turbo(h, s, l, rgb);
    draw_colored_triangle(x1, y1, x2, y2, x3, y3, rgb);   
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

void draw_rotated_rect(double x, double y, double l, double w, double angle)
{
    int i;
    double ca, sa, lca, lsa, wca, wsa;
    
    ca = cos(angle);
    sa = sin(angle);
    lca = l*ca;
    lsa = l*sa;
    wca = w*ca;
    wsa = w*sa;
        
    glBegin(GL_LINE_LOOP);
    glVertex2d(x + lca - wsa, y + lsa + wca);
    glVertex2d(x - lca - wsa, y - lsa + wca);
    glVertex2d(x - lca + wsa, y - lsa - wca);
    glVertex2d(x + lca + wsa, y + lsa - wca);
    glEnd();
}

void draw_colored_rotated_rect(double x, double y, double l, double w, double angle, double rgb[3])
{
    int i;
    double ca, sa, lca, lsa, wca, wsa;
    
    ca = cos(angle);
    sa = sin(angle);
    lca = l*ca;
    lsa = l*sa;
    wca = w*ca;
    wsa = w*sa;
    
    glColor3f(rgb[0], rgb[1], rgb[2]);
    glBegin(GL_TRIANGLE_FAN);
    glVertex2d(x, y);
    glVertex2d(x + lca - wsa, y + lsa + wca);
    glVertex2d(x - lca - wsa, y - lsa + wca);
    glVertex2d(x - lca + wsa, y - lsa - wca);
    glVertex2d(x + lca + wsa, y + lsa - wca);
    glVertex2d(x + lca - wsa, y + lsa + wca);
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


double partner_color(int np)
{
    switch(np){
        case (0): return(340.0);
        case (1): return(260.0);
        case (2): return(210.0);
        case (3): return(140.0);
        case (4): return(70.0);
        default:  return(20.0);
//         case (0): return(70.0);
//         case (1): return(200.0);
//         case (2): return(280.0);
//         case (3): return(140.0);
//         case (4): return(320.0);
//         default:  return(20.0);
    }   
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
        case (1): return(HUE_TYPE1);
        case (2): return(HUE_TYPE2);
        case (3): return(HUE_TYPE3);
        case (4): return(HUE_TYPE4);
        case (5): return(HUE_TYPE5);
        case (6): return(HUE_TYPE6);
        case (7):
        {
            if (RD_REACTION == CHEM_BZ) return(HUE_TYPE2);
            else return(HUE_TYPE7);
        }
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
    
    if (COUNT_PARTNER_TYPE) hue = partner_color(type-1);
    else hue = type_hue(type);
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
            dx = (INITXMAX - INITXMIN)/((double)NGRIDX);
            dy = (INITYMAX - INITYMIN)/((double)NGRIDY);
            for (i = 0; i < NGRIDX; i++)
                for (j = 0; j < NGRIDY; j++)
                {
                    n = NGRIDY*i + j;
                    particles[n].xc = INITXMIN + ((double)i - 0.5)*dx;   
                    particles[n].yc = INITYMIN + ((double)j - 0.5)*dy;
                    particles[n].radius = MU;
                    /* activate only circles that intersect the domain */
                    if ((particles[n].yc < INITYMAX + MU)&&(particles[n].yc > INITYMIN - MU)&&(particles[n].xc < INITXMAX + MU)&&(particles[n].xc > INITXMIN - MU)) particles[n].active = 1;
                    else particles[n].active = 0;
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
                    if (particles[n].yc > INITYMAX) particles[n].yc += INITYMIN - INITYMAX;
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
    
    for (i=0; i<ncircles; i++) 
    {
        particles[i].coulomb = 1;
        particles[i].reactive = 1;
    }
}

void add_particle_config(t_particle particles[NMAXCIRCLES], double xmin, double xmax, double ymin, double ymax, int nx, int ny, short int reactive, double radius)
/* add particles to configuration */
{
    int i, j, k, n, ncirc0, n_p_active, ncandidates = PDISC_CANDIDATES, naccepted, newcircles; 
    double dx, dy, p, phi, r, r0, ra[5], sa[5], height, x, y = 0.0, gamma, dpoisson, xx[4], yy[4];
    short int active_poisson[NMAXCIRCLES], far;
    
    dpoisson = PDISC_DISTANCE*radius;
    
    switch (CIRCLE_PATTERN_B) {
        case (C_SQUARE):
        {
            n = ncircles;
            dx = (xmax - xmin)/((double)nx);
            dy = (ymax - ymin)/((double)ny);
            for (i = 0; i < nx; i++)
                for (j = 0; j < NGRIDY; j++)
                {
                    particles[n].xc = xmin + ((double)i + 0.5)*dx;
                    particles[n].yc = ymin + ((double)j + 0.5)*dy;
                    particles[n].radius = MU;
                    particles[n].active = 1;
                    particles[n].added = 1;
                    particles[n].coulomb = 1;
                    particles[n].reactive = reactive;
                    n++;
                    ncircles++;
                }
            break;
        }
        case (C_HEX):   /* TODO */
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
                        particles[ncircles].added = 1;
                        particles[ncircles].coulomb = 1;
                        particles[ncircles].reactive = reactive;
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
    int i, j, n, jmin, jmax, nx, ny, ntot; 
    double x, y, dx, dy, width, lpocket, xmid = 0.5*(BCXMIN + BCXMAX), radius, c;
    
    /* set default rotation to 0 */
    for (i=0; i<NMAXOBSTACLES; i++) 
    {
        obstacle[i].omega = 0.0;
        obstacle[i].angle = 0.0;
        obstacle[i].oscillate = 0;
    }
    
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
        case (O_HEX):
        {
            dx = (XMAX - XMIN)/(double)NOBSX;
            dy = (YMAX - YMIN)/(double)NOBSY;
            n = 0;
            
            if ((ADD_FIXED_SEGMENTS)&&(SEGMENT_PATTERN == S_CYLINDER))
            /* pseudo-cylindrical boundary conditions (e.g. for Hall effect) */
            {
                jmin = 1; 
                jmax = NOBSY-1;
            }
            else
            {
                jmin = 0;
                jmax = NOBSY;
            }
            
            for (i=0; i<NOBSX; i++)
                for (j=jmin; j<jmax; j++)
                {
                    obstacle[n].xc = XMIN + ((double)i + 0.25)*dx;
                    obstacle[n].yc = YMIN + ((double)j + 0.5)*dy;
                    if (j%2 == 1) obstacle[n].xc += 0.5*dx;
                    if (obstacle[n].xc > XMAX) obstacle[n].xc += (XMAX - XMIN);
                    obstacle[n].radius = OBSTACLE_RADIUS;
                    obstacle[n].active = 1;
                    n++;
                }
            nobstacles = n;
//             printf("Added %i obstacles\n", nobstacles);
//             for (n=0; n<nobstacles; n++) printf("Obstacle %i at (%.3f, %.3f)\n", n, obstacle[n].xc, obstacle[n].yc);
            break;
        }
        case (O_SIDES):
        {
            n = 0;
            for (i = 0; i < 9; i++)
                for (j = 0; j < 5; j++)
                    if ((i == 0)||(i == 8)||(j == 0)||(j == 4))
                {
                    obstacle[n].xc = BCXMIN + 0.125*((double)i)*(BCXMAX - BCXMIN);
                    obstacle[n].yc = BCYMIN + 0.25*((double)j)*(BCYMAX - BCYMIN);
                    obstacle[n].radius = OBSTACLE_RADIUS;
                    obstacle[n].active = 1;
                    n++;
                }
            nobstacles = n;
            break;
        }
        case (O_SIDES_B):
        {
            n = 0;
            nx = 16;
            ny = 8;
            dx = (BCXMAX - BCXMIN)/(double)nx;
            dy = (BCYMAX - BCYMIN)/(double)ny;
            for (i = 0; i < nx+1; i++)
                for (j = 0; j < ny+1; j++)
                    if ((i == 0)||(i == nx)||(j == 0)||(j == ny))
                    {
                        obstacle[n].xc = BCXMIN + (double)i*dx;
                        obstacle[n].yc = BCYMIN + (double)j*dy;
                        obstacle[n].radius = OBSTACLE_RADIUS;
                        obstacle[n].active = 1;
                        n++;
                    }
            nobstacles = n;
            break;
        }
        case (O_SIEVE):
        {
            n = 0;
            width = 1.0;
            
            dx = width/12.0;
            dy = dx*0.6;
            for (i = 0; i < 13; i++)
            {
                obstacle[n].xc = -1.0 + (double)i*dx;
                obstacle[n].yc = 0.7 - (double)i*dy;
                obstacle[n].radius = 0.95*OBSTACLE_RADIUS;
                obstacle[n].omega = 1.25*OBSTACLE_OMEGA;
                obstacle[n].active = 1;
                n++;
            }
            for (i = 0; i < 13; i++)
            {
                obstacle[n].xc = -1.0 + (double)i*dx;
                obstacle[n].yc = 0.2 - (double)i*dy;
                obstacle[n].radius = 1.2*OBSTACLE_RADIUS;
                obstacle[n].omega = OBSTACLE_OMEGA;
                obstacle[n].active = 1;
                n++;
            }
            for (i = 0; i < 13; i++)
            {
                obstacle[n].xc = -1.0 + (double)i*dx;
                obstacle[n].yc = -0.3 - (double)i*dy;
                obstacle[n].radius = 1.6*OBSTACLE_RADIUS;
                obstacle[n].omega = 0.75*OBSTACLE_OMEGA;
                obstacle[n].active = 1;
                n++;
            }
                        
            nobstacles = n;
            break;
        }
        case (O_SIEVE_B):
        {
            n = 0;
            width = 1.2;
            
            dx = width/14.0;
            dy = dx*0.6;
            for (i = 0; i < 15; i++)
            {
                obstacle[n].xc = -1.2 + (double)i*dx;
                obstacle[n].yc = 0.85 - (double)i*dy;
                obstacle[n].radius = 0.9*OBSTACLE_RADIUS;
                obstacle[n].omega = 1.25*OBSTACLE_OMEGA;
                obstacle[n].active = 1;
                n++;
            }
            for (i = 0; i < 15; i++)
            {
                obstacle[n].xc = -1.2 + (double)i*dx;
                obstacle[n].yc = 0.35 - (double)i*dy;
                obstacle[n].radius = 1.2*OBSTACLE_RADIUS;
                obstacle[n].omega = OBSTACLE_OMEGA;
                obstacle[n].active = 1;
                n++;
            }
            for (i = 0; i < 15; i++)
            {
                obstacle[n].xc = -1.2 + (double)i*dx;
                obstacle[n].yc = -0.15 - (double)i*dy;
                obstacle[n].radius = 1.6*OBSTACLE_RADIUS;
                obstacle[n].omega = 0.75*OBSTACLE_OMEGA;
                obstacle[n].active = 1;
                n++;
            }
            
            if (RATTLE_OBSTACLES) for (i = 0; i < n; i++)
            {
                obstacle[i].oscillate = 1;
                obstacle[i].period = 8;
                obstacle[i].amplitude = 0.0015;
                obstacle[i].phase = (double)i*DPI/10.0;
            }
                        
            nobstacles = n;
            break;
        }
        case (O_SIEVE_LONG):
        {
            n = 0;
            ntot = 36;
            width = 1.2;
            
            dx = (XMAX - XMIN)/(double)ntot;
            dy = dx*0.3;
            for (i = 0; i < ntot + 1; i++)
            {
                obstacle[n].xc = XMIN + (double)i*dx;
                obstacle[n].yc = 0.6 - (double)i*dy;
                obstacle[n].radius = OBSTACLE_RADIUS*(2.4 - 1.6*(double)i/(double)ntot);
                obstacle[n].omega = OBSTACLE_OMEGA;
                obstacle[n].active = 1;
                n++;
            }
            
            if (RATTLE_OBSTACLES) for (i = 0; i < n; i++)
            {
                obstacle[i].oscillate = 1;
                obstacle[i].period = 8;
                obstacle[i].amplitude = 0.0018;
                obstacle[i].phase = (double)i*DPI/10.0;
            }
                        
            nobstacles = n;
            break;
        }
        case (O_SIEVE_LONG_B):
        {
            n = 0;
            ntot = 45;
            width = 1.2;
            
            dx = (XMAX - XMIN)/(double)ntot;
            dy = dx*0.2;
            x = XMIN;
            y = 0.4;
            c = 0.019;
            
            while (x < XMAX)
            {
                obstacle[n].xc = x;
                obstacle[n].yc = y;
                obstacle[n].radius = OBSTACLE_RADIUS;
                obstacle[n].omega = OBSTACLE_OMEGA;
                obstacle[n].active = 1;
                x += dx;
                y -= dy;
                dx += c*dx;
                dy += c*dy;
                n++;
            }
            
            if (RATTLE_OBSTACLES) for (i = 0; i < n; i++)
            {
                obstacle[i].oscillate = 1;
                obstacle[i].period = 8;
                obstacle[i].amplitude = 0.0018;
                obstacle[i].phase = (double)i*DPI/10.0;
            }
                        
            nobstacles = n;
            break;
        }
        default: 
        {
            printf("Function init_obstacle_config not defined for this pattern \n");
        }
    }
    
    if (CHARGE_OBSTACLES)
        for (n=0; n<nobstacles; n++) obstacle[n].charge = OBSTACLE_CHARGE;
        
    for (n=0; n<nobstacles; n++) 
    {
        obstacle[n].omega0 = obstacle[n].omega;
        obstacle[n].xc0 = obstacle[n].xc;
        obstacle[n].yc0 = obstacle[n].yc;
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
        
        
void add_tree_to_segments(double x, double y, double r, double lfoot[2], double rfoot[2], t_segment segment[NMAXSEGMENTS], int group)
/* add segments forming a tree to linear obstacle configuration */
{
    int i, n = nsegments, nplus, nminus, nadd; 
    double angle, h, h1, h2, h3, x1, x2, xp1, xp2, dx1, dx2, xt, yt;
    
    angle = PI/3.0;
    h = r*tan(angle);
    h1 = 0.56*h;
    h2 = h1;
    dx1 = h1/tan(angle);
    dx2 = 0.5*dx1;
    x1 = r - dx1;
    xp1 = x1 + dx2;
    x2 = xp1 - dx1;
    xp2 = x2 + dx2;
    h3 = xp2*tan(angle);
    xt = 0.5*dx2;
    yt = dx1;
    
    nadd = 15;
    
    if (nsegments + nadd < NMAXSEGMENTS) 
    {
        /* bottom left */
        segment[n].x1 = x - r;
        segment[n].y1 = y;
        segment[n].angle1 = PID + angle;
        segment[n].angle2 = PI + PID;
        segment[n].concave = 1;
        n++;
        
        /* level 1 */
        segment[n].x1 = x - x1;
        segment[n].y1 = y + h1;
        segment[n].concave = 0;
        n++;
        
        segment[n].x1 = x - xp1;
        segment[n].y1 = y + h1;
        segment[n].angle1 = PID + angle;
        segment[n].angle2 = PI + PID;
        segment[n].concave = 1;
        n++;
        
        /* level 2 */
        segment[n].x1 = x - x2;
        segment[n].y1 = y + h1 + h2;
        segment[n].concave = 0;
        n++;
        
        segment[n].x1 = x - xp2;
        segment[n].y1 = y + h1 + h2;
        segment[n].angle1 = PID + angle;
        segment[n].angle2 = PI + PID;
        segment[n].concave = 1;
        n++;
        
        /* top */
        segment[n].x1 = x;
        segment[n].y1 = y + h1 + h2 + h3;
        segment[n].angle1 = PID - angle;
        segment[n].angle2 = PID + angle;
        segment[n].concave = 1;
        n++;
        
        /* level 2 */
        segment[n].x1 = x + xp2;
        segment[n].y1 = y + h1 + h2;
        segment[n].angle1 = -PID;
        segment[n].angle2 = angle;
        segment[n].concave = 1;
        n++;
        
        segment[n].x1 = x + x2;
        segment[n].y1 = y + h1 + h2;
        segment[n].concave = 0;
        n++;
        
        /* level 1 */
        segment[n].x1 = x + xp1;
        segment[n].y1 = y + h1;
        segment[n].angle1 = -PID;
        segment[n].angle2 = angle;
        segment[n].concave = 1;
        n++;
        
        segment[n].x1 = x + x1;
        segment[n].y1 = y + h1;
        segment[n].concave = 0;
        n++;
        
        /* bottom right */
        segment[n].x1 = x + r;
        segment[n].y1 = y;
        segment[n].angle1 = -PID;
        segment[n].angle2 = angle;
        segment[n].concave = 1;
        n++;
        
        /*stem */
        segment[n].x1 = x + xt;
        segment[n].y1 = y;
        segment[n].concave = 0;
        n++;
        
        segment[n].x1 = x + xt;
        segment[n].y1 = y - yt;
        segment[n].angle1 = -PID;
        segment[n].angle2 = 0.0;
        segment[n].concave = 1;
        n++;
        
        segment[n].x1 = x - xt;
        segment[n].y1 = y - yt;
        segment[n].angle1 = PI;
        segment[n].angle2 = 3.0*PID;
        segment[n].concave = 1;
        n++;
        
        segment[n].x1 = x - xt;
        segment[n].y1 = y;
        segment[n].concave = 0;
        n++;
        
        rfoot[0] = x + xt;
        rfoot[1] = y - yt;
        lfoot[0] = x - xt;
        lfoot[1] = y - yt;
        
        for (i=0; i<nadd-1; i++)
        {
            segment[nsegments+i].x2 = segment[nsegments+i+1].x1;
            segment[nsegments+i].y2 = segment[nsegments+i+1].y1;
        }
        segment[nsegments+nadd-1].x2 = segment[nsegments].x1;
        segment[nsegments+nadd-1].y2 = segment[nsegments].y1;
        
        for (i=0; i<nadd; i++)
        {
//             segment[nsegments+i].concave = 1;
            segment[nsegments+i].group = group;
        }
    
        nsegments += nadd;
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


int add_shovel_to_belt(double x, double y, double angle, double size, double width, t_segment segment[NMAXSEGMENTS], t_belt belt)
/* add showver to conveyor belt */
{
    int i, n, iminus, iplus;
    double sq2, angle1, angle2, cg, sg, t;
    
    if (nsegments + 6 > NMAXSEGMENTS) 
    {
        printf("NMAXSEGMENTS too small");
        return(0);
    }
    
    n = nsegments; 
//     sq2 = sqrt(2.0)/2.0;
    cg = cos(PI/6.0);
    sg = sin(PI/6.0);
    t = size + width*(cg - 1.0)/sg;
    
    segment[n].x1 = 0.0;
    segment[n].y1 = 0.0;
    
    segment[n+1].x1 = 0.0;
    segment[n+1].y1 = size;
    
    segment[n+2].x1 = size*sg;
    segment[n+2].y1 = size*(1.0+cg);
    
    segment[n+3].x1 = size*sg + width*cg;
    segment[n+3].y1 = size*(1.0+cg) - width*sg;
    
    segment[n+4].x1 = width;
    segment[n+4].y1 = size*(1.0 + cg) - width*sg - t*cg;
    
//     segment[n+2].x1 = size*sq2;
//     segment[n+2].y1 = size*(1.0+sq2);
//     
//     segment[n+3].x1 = size*sq2 + width*sq2;
//     segment[n+3].y1 = size*(1.0+sq2) - width*sq2;
//     
//     segment[n+4].x1 = width;
//     segment[n+4].y1 = size + width*(1.0 - sqrt(2.0));

    segment[n+5].x1 = width;
    segment[n+5].y1 = 0.0;    
    
    for (i=0; i<6; i++)
    {
        segment[n+i].x2 = segment[n+i+1].x1;
        segment[n+i].y2 = segment[n+i+1].y1;
    }
    segment[n+6].x2 = segment[n].x1;
    segment[n+6].y2 = segment[n].y1;
    
    for (i=n; i<n+6; i++) 
    {
        segment[i].inactivate = 0;
        segment[i].conveyor = 0;
        segment[i].concave = 1;
    }
    
    /* deal with concave corners */
    for (i=n; i<n+6; i++)
    {
        iminus = i-1;  
        iplus = i+1;  
        if (iminus < n) iminus = n + 5;
        if (iplus > n + 5) iplus = n;
        angle1 = argument(segment[iplus].x1 - segment[i].x1, segment[iplus].y1 - segment[i].y1) + PID;
        angle2 = argument(segment[i].x1 - segment[iminus].x1, segment[i].y1 - segment[iminus].y1) + PID;
        if (angle2 < angle1) angle2 += DPI;
        segment[i].angle1 = angle1;
        segment[i].angle2 = angle2;
    }
    
    for (i=n; i<n+6; i++)
    {
        translate_one_segment(segment, i, x, y);
        rotate_one_segment(segment, i, angle, x, y);
    }
    
    nsegments += 6;
    return(1);
}

int add_conveyor_belt(double x1, double y1, double x2, double y2, double width, double speed, int nshovels, t_segment segment[NMAXSEGMENTS], t_belt conveyor_belt[NMAXBELTS])
/* add segments forming a conveyor belt */
{
    double length, tx, ty, angle, angle1, angle2, beltlength, shovel_dist, shovel_pos, x, y, beta;
    int i, n, iplus, iminus, shovel;
    
    length = module2(x2 - x1, y2 - y1);
    angle = argument(x2 - x1, y2 - y1);
    if (angle < 0.0) angle += DPI;
    tx = (x2 - x1)/length;
    ty = (y2 - y1)/length;
    
    n = nsegments;
    
    if (nsegments + 14 > NMAXSEGMENTS) 
    {
        printf("NMAXSEGMENTS too small");
        return(0);
    }
    if (nbelts + 1 > NMAXBELTS) 
    {
        printf("NMAXBELTS too small");
        return(0);
    }

    segment[n].x1 = x1 - ty*width;
    segment[n].y1 = y1 + tx*width;
                        
    for (i=0; i<7; i++)
    {
        segment[n+1+i].x1 = x2 + width*sin(-angle + (double)i*PI/6.0);
        segment[n+1+i].y1 = y2 + width*cos(-angle + (double)i*PI/6.0);
    }
    for (i=7; i<14; i++)
    {
        segment[n+1+i].x1 = x1 + width*sin(-angle + (double)(i-1)*PI/6.0);
        segment[n+1+i].y1 = y1 + width*cos(-angle + (double)(i-1)*PI/6.0);
    }
            
    for (i=0; i<13; i++)
    {
        segment[n+i].x2 = segment[n+i+1].x1;
        segment[n+i].y2 = segment[n+i+1].y1;
    }
    segment[n+13].x2 = segment[n].x1;
    segment[n+13].y2 = segment[n].y1;
                        
    for (i=n; i<n+14; i++) 
    {
        segment[i].concave = 1;
        segment[i].inactivate = 0;
        segment[i].conveyor = 1;
        segment[i].conveyor_speed = speed;
        segment[i].align_torque = 1;
    }
    
    /* deal with concave corners */
    for (i=n; i<n+14; i++)
    {
        iminus = i-1;  
        iplus = i+1;  
        if (iminus < n) iminus = n + 13;
        if (iplus > n + 13) iplus = n;
        angle1 = argument(segment[iplus].x1 - segment[i].x1, segment[iplus].y1 - segment[i].y1) + PID;
        angle2 = argument(segment[i].x1 - segment[iminus].x1, segment[i].y1 - segment[iminus].y1) + PID;
        if (angle2 < angle1) angle2 += DPI;
        segment[i].angle1 = angle1;
        segment[i].angle2 = angle2;
    }
    
    nsegments += 14;

    conveyor_belt[nbelts].x1 = x1;
    conveyor_belt[nbelts].x2 = x2;
    conveyor_belt[nbelts].y1 = y1;
    conveyor_belt[nbelts].y2 = y2;
    conveyor_belt[nbelts].speed = speed;
    conveyor_belt[nbelts].width = width;
    conveyor_belt[nbelts].position = 0.0;
    conveyor_belt[nbelts].length = length;
    conveyor_belt[nbelts].angle = angle;
    conveyor_belt[nbelts].tx = tx;
    conveyor_belt[nbelts].ty = ty;
    conveyor_belt[nbelts].nshovels = nshovels; 
    nbelts++; 
    
    if (nshovels > NMAXSHOVELS)
    {
        printf("NMAXSHOVELS too small");
        conveyor_belt[nbelts].nshovels = 0; 
        return(0);
    }
    
    if (nshovels > 0)
    {
        beltlength = 2.0*length + DPI*width;
        shovel_dist = beltlength/(double)nshovels; 
        shovel_pos = 0.0;
        shovel = 0;
        
        while (shovel_pos < length)
        {
            x = x1 - width*ty + shovel_pos*tx;
            y = y1 + width*tx + shovel_pos*ty;
            conveyor_belt[nbelts-1].shovel_segment[shovel] = nsegments;
            add_shovel_to_belt(x, y, angle, width, 0.3*width, segment, conveyor_belt[nbelts-1]);
            conveyor_belt[nbelts-1].shovel_pos[shovel] = shovel_pos;
            shovel_pos += shovel_dist;
            shovel++;
        }
        while (shovel_pos < length + width*PI)
        {
            beta = (shovel_pos - length)/width;
            x = x2 + width*cos(angle + PID - beta);
            y = y2 + width*sin(angle + PID - beta);
            conveyor_belt[nbelts-1].shovel_segment[shovel] = nsegments;
            add_shovel_to_belt(x, y, angle - beta, width, 0.3*width, segment, conveyor_belt[nbelts-1]);
            conveyor_belt[nbelts-1].shovel_pos[shovel] = shovel_pos;
            shovel_pos += shovel_dist;
            shovel++;
        }
        while (shovel_pos < 2.0*length + width*PI)
        {
            x = x2 + width*ty - (shovel_pos - length - width*PI)*tx;
            y = y2 - width*tx - (shovel_pos - length - width*PI)*ty;
            conveyor_belt[nbelts-1].shovel_segment[shovel] = nsegments;
            add_shovel_to_belt(x, y, angle - PI, width, 0.3*width, segment, conveyor_belt[nbelts-1]);
            conveyor_belt[nbelts-1].shovel_pos[shovel] = shovel_pos;
            shovel_pos += shovel_dist;
            shovel++;
        }
        while (shovel_pos < beltlength)
        {
            beta = (shovel_pos - 2.0*length - width*PI)/width;
            x = x1 + width*cos(angle - PID - beta);
            y = y1 + width*sin(angle - PID - beta);
            conveyor_belt[nbelts-1].shovel_segment[shovel] = nsegments;
            add_shovel_to_belt(x, y, angle - PI - beta, width, 0.3*width, segment, conveyor_belt[nbelts-1]);
            conveyor_belt[nbelts-1].shovel_pos[shovel] = shovel_pos;
            shovel_pos += shovel_dist;
            shovel++;
        }
    }
    
    return(1);
}

void init_segment_config(t_segment segment[NMAXSEGMENTS], t_belt conveyor_belt[NMAXBELTS])
/* initialise linear obstacle configuration */
{
    int i, j, cycle = 0, iminus, iplus, nsides, n, concave = 1; 
    double angle, angle2, dangle, dx, width, height, a, b, length, xmid = 0.5*(BCXMIN + BCXMAX), lpocket, r, x, x1, y1, x2, y2, nozx, nozy, y, yy, dy, ca, sa, padding;
    double lfoot[NTREES][2], rfoot[NTREES][2];
    
    /* set default to no conveyor, no torque */
    for (i=0; i<NMAXSEGMENTS; i++) 
    {
        segment[i].conveyor = 0;
        segment[i].align_torque = 0;
    }
    
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
        case (S_HLINE_HOLE_SLOPED):
        {
            x = 0.15;
            x1 = XMAX + 1.0;
            width = 0.05;
            angle = APOLY*DPI;
            y1 = (x1 - x)*tan(angle);
                        
            segment[0].x1 = x1;
            segment[0].y1 = y1;
            
            segment[1].x1 = x1 + width*sin(angle);
            segment[1].y1 = y1 - width*cos(angle);
            
            segment[2].x1 = x + width*sin(angle);
            segment[2].y1 = -width*cos(angle);
            
            segment[3].x1 = x;
            segment[3].y1 = 0.0;
            
            for (i=0; i<3; i++)
            {
                segment[i].x2 = segment[i+1].x1;
                segment[i].y2 = segment[i+1].y1;
            }
            segment[3].x2 = segment[0].x1;
            segment[3].y2 = segment[0].y1;
            
            segment[4].x1 = -x1;
            segment[4].y1 = y1;
            
            segment[5].x1 = -x;
            segment[5].y1 = 0.0;

            segment[6].x1 = -x - width*sin(angle);
            segment[6].y1 = -width*cos(angle);
            
            segment[7].x1 = -x1 - width*sin(angle);
            segment[7].y1 = y1 - width*cos(angle);
            
            for (i=4; i<7; i++)
            {
                segment[i].x2 = segment[i+1].x1;
                segment[i].y2 = segment[i+1].y1;
            }
            segment[7].x2 = segment[4].x1;
            segment[7].y2 = segment[4].y1;
            
            nsegments = 8;
            
            for (i=0; i<nsegments; i++) segment[i].concave = 1;
            
            /* closing segment */
            segment[nsegments].x1 = -x;
            segment[nsegments].y1 = 0.0;
            segment[nsegments].x2 = x;
            segment[nsegments].y2 = 0.0;
            nsegments++;
        
            cycle = 0;
            concave = 1;
            
            for (i=0; i<nsegments; i++) 
            {
                segment[i].inactivate = 0;
            }
            segment[nsegments-1].inactivate = 1;
            
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
        case (S_BIN_LARGE):
        {
            add_rectangle_to_segments(LAMBDA, -MU, -LAMBDA, 0.0, segment, 0);
            add_rectangle_to_segments(LAMBDA, -MU, LAMBDA - MU, 0.25, segment, 0);
            add_rectangle_to_segments(-LAMBDA + MU, -MU, -LAMBDA, 0.25, segment, 0);
                        
            cycle = 0;
            concave = 0;    /* add_rectangle_to_segments already deals with concave corners */
            
            for (i=0; i<nsegments; i++) 
            {
                segment[i].group = 0;
                segment[i].inactivate = 0;
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
        case (S_COANDA):
        {
            y1 = SEGMENTS_Y0;
            width = 0.05;
            padding = 0.01;
            
            height = 0.25;
            dx = (XMAX - XMIN + 2.0*padding)/(double)(NPOLY-1);
            
            for (i=0; i<NPOLY; i++)
            {
                x = XMIN - padding + (double)i*dx;
                y = y1 + height*cos(PI*x/XMAX);
                
                segment[i].x1 = x;
                segment[i].y1 = y + width;
                segment[i].x2 = x + dx;
                segment[i].x2 = y1 + height*cos(PI*(x+dx)/XMAX) + width;
                
                segment[i].concave = 1;
            }
            
            for (i=0; i<NPOLY; i++)
            {
                x = XMAX + padding - (double)i*dx;
                y = y1 + height*cos(PI*x/XMAX);
                
                segment[NPOLY+i].x1 = x;
                segment[NPOLY+i].y1 = y - width;
                segment[NPOLY+i].x2 = x - dx;
                segment[NPOLY+i].x2 = y1 + height*cos(PI*(x-dx)/XMAX) - width;
                
                segment[NPOLY+i].concave = 1;
            }
            
            nsegments = 2*NPOLY;
            
            cycle = 1;
            concave = 1;    
            
            for (i=0; i<nsegments; i++) 
            {
                segment[i].inactivate = 0;
            }
            
            break;
        }
        case (S_COANDA_SHORT):
        {
            y1 = SEGMENTS_Y0;
            width = 0.05;
            padding = 0.1;
            x1 = XMIN + padding;
            
            height = 0.2;
            length = 0.85*(XMAX - XMIN - 2.0*padding);
            dx = length/(double)(NPOLY-1);
            
            for (i=0; i<NPOLY; i++)
            {
                x = x1 + (double)i*dx;
                y = y1 - height*cos(DPI*(x-x1)/length);
                yy = y1 - height*cos(DPI*(x+dx-x1)/length);
                add_rotated_angle_to_segments(x, y, x+dx, yy, width, 1, segment, 0);
            }
            
            cycle = 0;
            concave = 1;    
            
            for (i=0; i<nsegments; i++) 
            {
                segment[i].inactivate = 0;
            }
            
            break;
        }
        case (S_CYLINDER):
        {
            add_rectangle_to_segments(XMAX + 0.1, YMAX - 0.05, XMIN - 0.1, YMAX + 1.0, segment, 0);
            add_rectangle_to_segments(XMAX + 0.1, YMIN - 1.0, XMIN - 0.1, YMIN + 0.05, segment, 0);
            
            cycle = 0;
            concave = 1;    
            
            for (i=0; i<nsegments; i++) 
            {
                segment[i].inactivate = 0;
            }
            
            break;
        }
        case (S_TREE):
        {
            for (i=0; i<NTREES; i++)
            {
                x = XMIN + ((double)i - 0.5 + 0.1*(double)rand()/(double)RAND_MAX)*(XMAX-XMIN)*1.1/(double)NTREES; 
                y = YMIN + 0.1 + 0.014*(double)i*(YMAX - YMIN);
                add_tree_to_segments(x, y, LAMBDA*(1.0 +  0.4*(double)rand()/(double)RAND_MAX), lfoot[i], rfoot[i], segment, 0);
            }
            /* add ground */
            for (i=0; i<NTREES-1; i++)
            {
                segment[nsegments].x1 = rfoot[i][0];
                segment[nsegments].y1 = rfoot[i][1];
                segment[nsegments].x2 = lfoot[i+1][0];
                segment[nsegments].y2 = lfoot[i+1][1];
                segment[nsegments].concave = 0;
                nsegments++;
            }
            cycle = 0;
            concave = 0;
            ngroups = 1;
            
            for (i=0; i<nsegments; i++) 
            {
                segment[i].concave = 1;
                segment[i].inactivate = 0;
            }
            
            break;
        }
        case (S_CONE):
        {
            add_rotated_angle_to_segments(LAMBDA, LAMBDA, 0.05, -0.2, WALL_WIDTH, 0, segment, 0);
            add_rotated_angle_to_segments(-0.05, -0.2, -LAMBDA, LAMBDA, WALL_WIDTH, 0, segment, 0);
                        
            cycle = 0;
            concave = 1;    
            
            for (i=0; i<nsegments; i++) 
            {
                segment[i].group = 0;
                segment[i].inactivate = 0;
            }
            break;
        }
        case (S_CONVEYOR_BELT):
        {
            add_conveyor_belt(XMIN + 0.05, 0.0, 0.5, 0.0, 0.05, BELT_SPEED1, 0, segment, conveyor_belt);
                                
            cycle = 0;
            concave = 0;
            break;
        }
        case (S_TWO_CONVEYOR_BELTS):
        {
            add_conveyor_belt(XMIN + 0.05, 0.2, 0.8, 0.5, 0.05, BELT_SPEED1, 0, segment, conveyor_belt);
            add_conveyor_belt(-0.5, -0.3, XMAX - 0.05, -0.8, 0.05, BELT_SPEED2, 0, segment, conveyor_belt);
                                
            cycle = 0;
            concave = 0;
            break;
        }
        case (S_PERIODIC_CONVEYORS):
        {
            add_conveyor_belt(0.05, -0.3, XMAX - XMIN - 0.125, 0.4, 0.05, BELT_SPEED1, 0, segment, conveyor_belt);
            add_conveyor_belt(XMIN - XMAX + 0.05, -0.3, -0.125, 0.4, 0.05, BELT_SPEED1, 0, segment, conveyor_belt);
            add_conveyor_belt(-1.0, -0.6, 0.2, -0.6, 0.05, BELT_SPEED2, 0, segment, conveyor_belt);
                                
            cycle = 0;
            concave = 0;
            break;            
        }
        case (S_CONVEYOR_SHOVELS):
        {
            add_conveyor_belt(-1.0, -0.97, 1.0, 0.8, 0.05, BELT_SPEED1, 25, segment, conveyor_belt);
            add_rectangle_to_segments(WALL_WIDTH, YMIN - 0.05, -WALL_WIDTH, -0.5, segment, 0);
            
            cycle = 0;
            concave = 0;
            break;
        }
        case (S_CONVEYOR_MIXED):
        {
            add_conveyor_belt(XMIN - 0.32, 0.65, -0.2, 0.65, 0.05, BELT_SPEED2, 0, segment, conveyor_belt);
            add_conveyor_belt(XMAX - 0.32, 0.65, XMAX - XMIN - 0.2, 0.65, 0.05, BELT_SPEED2, 0, segment, conveyor_belt);
            add_conveyor_belt(XMIN + 0.3, 0.3, 0.0, 0.3, 0.05, -BELT_SPEED1, 0, segment, conveyor_belt);
            add_conveyor_belt(XMIN + 0.1, -0.05, -0.2, -0.05, 0.05, BELT_SPEED2, 0, segment, conveyor_belt);
            add_conveyor_belt(XMIN - 0.5, -0.4, 0.0, -0.4, 0.05, -BELT_SPEED1, 0, segment, conveyor_belt);
            add_conveyor_belt(XMAX - 0.5, -0.4, XMAX - XMIN, -0.4, 0.05, -BELT_SPEED1, 0, segment, conveyor_belt);
            add_conveyor_belt(XMIN - 0.3, -0.75, -0.25, -0.75, 0.05, BELT_SPEED2, 0, segment, conveyor_belt);
            add_conveyor_belt(XMAX - 0.7, -0.75, XMAX - XMIN - 0.25, -0.75, 0.05, BELT_SPEED2, 0, segment, conveyor_belt);
            
            add_conveyor_belt(-0.13, -0.97, 1.67, 0.95, 0.05, BELT_SPEED3, 24, segment, conveyor_belt);
            
            cycle = 0;
            concave = 0;
            break;
        }
        case (S_CONVEYOR_SIEVE):
        {
            add_conveyor_belt(XMIN, 0.8, -0.4, 0.8, 0.05, BELT_SPEED1, 0, segment, conveyor_belt);
            
            add_conveyor_belt(0.05, 0.0, XMAX - 0.3, 0.0, 0.05, BELT_SPEED1, 0, segment, conveyor_belt);
            add_conveyor_belt(0.05, -0.5, XMAX - 1.0, -0.5, 0.05, BELT_SPEED1, 0, segment, conveyor_belt);
            
            add_rectangle_to_segments(1.3, YMIN - 0.05, 1.3 - WALL_WIDTH, -0.3, segment, 0);
            add_rectangle_to_segments(0.6, YMIN - 0.05, 0.6 - WALL_WIDTH, -0.75, segment, 0);
            add_rectangle_to_segments(0.0, YMIN - 0.05, - WALL_WIDTH, -0.95, segment, 0);
            
            cycle = 0;
            concave = 0;
            break;
        }
        case (S_CONVEYOR_SIEVE_B):
        {
            add_conveyor_belt(-0.7, 0.8, XMAX + 0.1, 0.8, 0.05, -BELT_SPEED1, 0, segment, conveyor_belt);
            
            add_conveyor_belt(0.05, 0.0, XMAX - 0.3, 0.0, 0.05, BELT_SPEED2, 0, segment, conveyor_belt);
            add_conveyor_belt(0.05, -0.5, XMAX - 1.0, -0.5, 0.05, BELT_SPEED2, 0, segment, conveyor_belt);
            
            add_rectangle_to_segments(1.5, YMIN - 0.05, 1.5 - WALL_WIDTH, -0.15, segment, 0);
            add_rectangle_to_segments(0.75, YMIN - 0.05, 0.75 - WALL_WIDTH, -0.75, segment, 0);
            add_rectangle_to_segments(0.0, YMIN - 0.05, 0.0 - WALL_WIDTH, -0.95, segment, 0);
            add_rectangle_to_segments(-1.0, YMIN - 0.05, -1.0 - WALL_WIDTH, -0.95, segment, 0);
            
            cycle = 0;
            concave = 0;
            break;
        }
        case (S_CONVEYOR_SIEVE_LONG):
        {
            add_conveyor_belt(-1.4, 0.8, XMAX + 0.1, 0.8, 0.05, -BELT_SPEED1, 0, segment, conveyor_belt);
                       
            for (i=0; i<7; i++)
            {
                x = XMIN + (double)i*(XMAX-XMIN)/7.0; 
                add_rectangle_to_segments(x, YMIN - 0.05, x - WALL_WIDTH, -0.75, segment, 0);
            }
            
            cycle = 0;
            concave = 0;
            break;
        }
        case (S_MASS_SPECTROMETER):
        {
            for (i=0; i<7; i++)
            {                
                y = YMIN + (double)i*(YMAX-YMIN)/7.0;
                add_rectangle_to_segments(XMAX + 0.2, y, XMAX - 0.6, y + WALL_WIDTH, segment, 0);
            }
            
            cycle = 0;
            concave = 0;
            break;
        }
        case (S_WIND_FORCE):
        {
            for (i=0; i<10; i++)
            {                
                x = XMIN + (double)i*(XMAX-XMIN)/10.0; 
                add_rectangle_to_segments(x, YMIN - 0.05, x - WALL_WIDTH, WIND_YMIN - 0.1, segment, 0);
            }
            
            cycle = 0;
            concave = 0;
            break;
        }
        case (S_TEST_CONVEYORS):
        {
            add_conveyor_belt(-1.0, -0.7, 1.0, 0.7, 0.05, BELT_SPEED1, 25, segment, conveyor_belt);
//             add_conveyor_belt(-0.5, -0.3, 0.5, 0.3, 0.05, BELT_SPEED1, 100, segment, conveyor_belt);
            
                                
            cycle = 0;
            concave = 0;
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
        segment[i].nangle = argument(segment[i].nx, segment[i].ny);
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
    
    /* initialize averaged pressure */
    if (SHOW_SEGMENTS_PRESSURE) for (i=0; i<nsegments; i++) 
        segment[i].avrg_pressure = 0.0;
    
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
    double angle, dx, height, width, length, theta, lx, ly, x1, y1, x2, y2, padding, ca, sa, r;
    
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
        case (S_COANDA):
        {
            return(vabs(y - SEGMENTS_Y0 - 0.25*cos(PI*x/XMAX)) > 0.1);
        }
        case (S_COANDA_SHORT):
        {
            x1 = XMIN + 0.1;
            length = 0.85*(XMAX - XMIN - 0.2);
            if (x > x1 + length + 0.1) return(1);
            return(vabs(y - SEGMENTS_Y0 + 0.2*cos(DPI*(x-x1)/length)) > 0.05);
        }
        case (S_CYLINDER):
        {
            if (y > YMAX - 0.05 - MU) return(0);
            if (y < YMIN + 0.05 + MU) return(0);
            return(1);
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

void segment_interaction(t_particle part_i, t_particle part_k, double ca, double sa, double distance, double ffm[3], int charge)
/* compute interaction for I_SEGMENT interaction */
/* computes force and torque of particle k on particle i */
{
    int s;
    double x0, y0, x, y, fx, fy, cb, sb, cc, sc, cab, sab, ccb, scb, r;
    double width, torque, mu_i, mu_k, r3, endfact, cfact;
    
    fx = 0.0;
    fy = 0.0;
    torque = 0.0;
    mu_i = part_i.radius; 
    mu_k = part_k.radius;
    width = 0.2*mu_i;
    
    cb = cos(part_i.angle);
    sb = sin(part_i.angle);
    cc = cos(part_k.angle);
    sc = sin(part_k.angle);
    
    cab = ca*cb + sa*sb;
    sab = sa*cb - ca*sb;
    ccb = cc*cb + sc*sb;
    scb = sc*cb - cc*sb;
    
    /* relative coordinates of particle k in basis where particle i is horizontal */
    x0 = distance*cab;
    y0 = distance*sab;
    
    for (s=-1; s<=1; s+=2)
    {
        x = x0 + (double)s*mu_k*ccb;
        y = y0 + (double)s*mu_k*scb;
        
        if (vabs(y) < width)
        {
            if (vabs(x) < mu_i)
            {
//                 if (((double)s*scb < 0.0)&&(y > 0.0)) fy += y - width;
                if (y > 0.0) fy += y - width;
                else fy += y + width;
            }
            else if ((x >= mu_i)&&(x < mu_i+width))
            {
                r = module2(x-mu_i, y) + 1.0e-10;
                if (r < width)
                {
                    fx -= (width - r)*(x-mu_i)/r;
                    fy -= (width - r)*y/r;
                }
            }
            else if ((x <= -mu_i)&&(x > -mu_i-width))
            {
                r = module2(x+mu_i, y) + 1.0e-10;
                if (r < width)
                {
                    fx -= (width - r)*(x+mu_i)/r;
                    fy -= (width - r)*y/r;
                }
            }
            
            torque += x*fy - y*fx;
        }
        
        if (charge)     /* add Coulomb force between ends */
        {
            cfact = CHARGE*KREPEL/KSPRING_OBSTACLE;
            
            r = module2(x-mu_i,y) + 1.0e-2;
            r3 = r*r*r;
            fx += cfact*(double)s*(x+mu_i)/r3;
            fy += cfact*(double)s*y/r3;
            
            torque += mu_i*fx;
            
            r = module2(x+mu_i,y) + 1.0e-2;
            r3 = r*r*r;
            fx -= cfact*(double)s*(x+mu_i)/r3;
            fy -= cfact*(double)s*y/r3;
            
            torque -= mu_i*fy;
        }
    }
    
    ffm[0] = cb*fx - sb*fy;
    ffm[1] = sb*fx + cb*fy;
    ffm[2] = torque;
}

void polygon_interaction(t_particle part_i, t_particle part_k, double ca, double sa, double distance, double ffm[3], int charge)
/* compute interaction for I_POLYGON interaction */
/* computes force and torque of particle k on particle i */
{
    int s, k, s0, s1;
    double x0, y0, x, y, f, fx, fy, cb, sb, cab, sab, r, r0, z2, phi, ckt, skt, d, fx1, fy1, dmax2;
    double torque, mu_i, mu_k, r3, width, gamma_beta, rot, cr, sr, xx[NPOLY], yy[NPOLY], dd[NPOLY], z, z0, z02, d1, dmin;
    static double theta, twotheta, ctheta, stheta, cornerx[NPOLY+1], cornery[NPOLY+1], nx[NPOLY], ny[NPOLY];
    static int first = 1;
    
    if (first)
    {
        theta = PI/(double)NPOLY;
        twotheta = 2.0*theta;
        ctheta = cos(theta);
        stheta = sin(theta);
        
        for (s=0; s<NPOLY+1; s++)
        {
            cornerx[s] = cos((double)s*twotheta);
            cornery[s] = sin((double)s*twotheta);
        }
        
        for (s=0; s<NPOLY; s++)
        {
            nx[s] = cos((double)(2*s+1)*theta);
            ny[s] = sin((double)(2*s+1)*theta);
        }
        
        first = 0;
    }
    
    fx = 0.0;
    fy = 0.0;
    torque = 0.0;
    mu_i = part_i.radius; 
    mu_k = part_k.radius;
    dmax2 = 4.0*mu_i*mu_i;
    width = 0.1*mu_i;
    mu_i -= width;
    gamma_beta = part_k.angle - part_i.angle;
    r0 = mu_i*ctheta;
    z0 = mu_i*stheta;
    z02 = z0*z0;
    
    cb = cos(part_i.angle);
    sb = sin(part_i.angle);
    cab = ca*cb + sa*sb;
    sab = sa*cb - ca*sb;
    
    /* relative coordinates of particle k in basis where particle i is horizontal */
    x0 = distance*cab;
    y0 = distance*sab;
    
    /* find vertex distances */
    for (s = 0; s<NPOLY; s++)
    {
        rot = gamma_beta + (double)s*twotheta;
        cr = cos(rot);
        sr = sin(rot); 
        
        xx[s] = x0 + cr*mu_k;
        yy[s] = y0 + sr*mu_k;
        
        dd[s] = xx[s]*xx[s] + yy[s]*yy[s];
    }
        
    /* compute force from vertices */
    for (s = 0; s<NPOLY; s++) if (dd[s] < dmax2)
    {
        x = xx[s];
        y = yy[s];
            
        phi = argument(x, y);
        if (phi < 0.0) phi += DPI;
        
        k = (int)(phi/twotheta);
        
        d = x*nx[k] + y*ny[k];
        
        if (d < r0 + width)
        {                
            z = -x*ny[k] + y*nx[k];
            if (z*z < z02) 
            {
                f = r0 + width - d;
                if (f > width) f = width;
                fx1 = f*nx[k];
                fy1 = f*ny[k];
                fx -= fx1;
                fy -= fy1;
                torque += x*fy1 - y*fx1;
//                 fx -= f*nx[k];
//                 fy -= f*ny[k];
//                 torque += x*fy - y*fx;
            }
            else for (s1=0; s1<NPOLY; s1++)    /* corners */
            {
                d1 = module2(x - mu_i*cornerx[s1], y - mu_i*cornery[s1]) + 1.0e-8;
                fx1 = 0.0;
                fy1 = 0.0;
                if (d1 < width)
                {
                    f = width - d1;
                    if (f > width) f = width;
                    fx1 -= f*(x - mu_i*cornerx[s1])/d1;
                    fy1 -= f*(y - mu_i*cornery[s1])/d1;
                }
                fx += fx1;
                fy += fy1;
                torque += x*fy1 - y*fx1;
            }
        }
    }
    
    ffm[0] = cb*fx - sb*fy;
    ffm[1] = sb*fx + cb*fy;
    ffm[2] = torque;
}

int dna_charge(int type1, int type2, int coulomb)
/* returns -1 for matching nucleotide ends, 1 for non-matching, and 0 else */
{
    int min, max;
    
//     if (coulomb == 0) return(0);
    
    if (type1 == 0) return(0);
    if (type2 == 0) return(0);
    
    if ((type1 == 7)&&(type2 == 7)) return(2);
    
    if (type1 == type2) return(1);
    
    if (type1 > type2)
    {
        max = type1;
        min = type2;
    }
    else
    {
        max = type2;
        min = type1;
    }
    
    /* enzymes */
    if ((min >= 3)&&(max == 7)) return(-10);
    
    if (max - min >= 2) return(1);
    if (min == 1) return(-4*coulomb);
    if (min%2 == 1) return(-10*coulomb);
    return(1);
}

int compute_particle_interaction(int i, int k, double force[2], double *torque, t_particle* particle, double distance, double krepel, double ca, double sa, double ca_rel, double sa_rel)
/* compute repelling force and torque of particle #k on particle #i */
/* returns 1 if distance between particles is smaller than NBH_DIST_FACTOR*MU */
{
    double x1, y1, x2, y2, r, f, angle, aniso, fx, fy, ff[2], dist_scaled, spin_f, ck, sk, ck_rel, sk_rel, alpha, amp, charge, f1, ffm[3];
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
        case (I_COULOMB_LJ):
        {
            charge = particle[i].charge*particle[k].charge;
            f = -100.0*krepel*charge/(1.0e-12 + distance*distance);
            if (charge <= 0.0) 
                f += COULOMB_LJ_FACTOR*krepel*lennard_jones_force(distance, particle[k]);
            force[0] = f*ca;
            force[1] = f*sa;   
            break;
        }
        case (I_COULOMB_PENTA):
        {
            charge = particle[i].charge*particle[k].charge;
            if (charge < 0.0)
            {
                penta_lj_force(distance, ca, sa, ca_rel, sa_rel, ff, particle[k]);
                force[0] = -charge*krepel*ff[0];
                force[1] = -charge*krepel*ff[1];
            }
            else
            {
                f = krepel*lennard_jones_force(distance, particle[k]);
                force[0] = f*ca;
                force[1] = f*sa;
            }
            break;
        }
        case (I_COULOMB_IMAGINARY):
        {
            charge = particle[i].charge*particle[k].charge;
            f = -100.0*krepel*charge/(1.0e-12 + distance*distance);
            if (charge <= 0.0) 
                f1 = 0.01*krepel*lennard_jones_force(distance, particle[k]);
            else f1 = 0.0;
            force[0] = f1*ca - f*sa;
            force[1] = f1*sa + f*ca;   
            break;
        }
        case (I_DNA_CHARGED):
        {
            f = 0.2*krepel*lennard_jones_force(distance, particle[k]);
            charge = CHARGE*(double)dna_charge(particle[i].type, particle[k].type, particle[i].coulomb);
            if ((RD_REACTION == CHEM_DNA_ENZYME)||(RD_REACTION == CHEM_DNA_ENZYME_REPAIR))
            {
                /* make some interactions repulsive */
                if ((particle[i].reactive == 0)&&(particle[k].reactive == 0)&&(particle[i].added == particle[k].added))
                    charge = vabs(charge);
                
                /* TEST */
                /* make interactions repulsive between base-paired and other molecules */
//                 else if (particle[i].paired != particle[k].paired) 
//                     charge = vabs(charge);
            }
            f -= krepel*charge/(1.0e-12 + distance*distance);
            force[0] = f*ca;
            force[1] = f*sa;
            break;
        }
        case (I_DNA_CHARGED_B):
        {
            f = 0.2*krepel*lennard_jones_force(distance, particle[k]);
            charge = CHARGE*(double)dna_charge(particle[i].type, particle[k].type, particle[i].coulomb);
            if (particle[i].added != particle[k].added) charge *= 1.5;
            if ((RD_REACTION == CHEM_DNA_ENZYME)||(RD_REACTION == CHEM_DNA_ENZYME_REPAIR))
            {
                /* make some interactions repulsive */
                if ((particle[i].reactive == 0)&&(particle[k].reactive == 0)&&(particle[i].added == particle[k].added))
                    charge = vabs(charge);
                
                /* TEST */
                /* make interactions repulsive between base-paired and other molecules */
//                 else if (particle[i].paired != particle[k].paired) 
//                     charge = vabs(charge);
            }
            f -= krepel*charge/(1.0e-12 + distance*distance);
            force[0] = f*ca;
            force[1] = f*sa;
            break;
        }
        case (I_SEGMENT):
        {
            segment_interaction(particle[i], particle[k], ca, sa, distance, ffm, 0);
            force[0] = KSPRING_OBSTACLE*ffm[0];
            force[1] = KSPRING_OBSTACLE*ffm[1];
            *torque = KTORQUE*ffm[2];

            f = 0.1*krepel*lennard_jones_force(distance, particle[k]);
            force[0] += f*ca;
            force[1] += f*sa; 
            
            break;
        }
        case (I_SEGMENT_CHARGED):
        {
            segment_interaction(particle[i], particle[k], ca, sa, distance, ffm, 1);
            force[0] = KSPRING_OBSTACLE*ffm[0];
            force[1] = KSPRING_OBSTACLE*ffm[1];
            *torque = KTORQUE*ffm[2];
            
//             printf("Force between particles %i and %i: (%.3lg, %.3lg)\n", i, k, force[0], force[1]); 
            
//             f = 0.1*krepel*lennard_jones_force(distance, particle[k]);
//             force[0] += f*ca;
//             force[1] += f*sa; 
            
            break;
        }
        case (I_POLYGON):
        {
            polygon_interaction(particle[i], particle[k], ca, sa, distance, ffm, 0);
            force[0] = KSPRING_OBSTACLE*ffm[0];
            force[1] = KSPRING_OBSTACLE*ffm[1];
            *torque = KTORQUE*ffm[2];

            f = 0.1*krepel*lennard_jones_force(distance, particle[k]);
            force[0] += f*ca;
            force[1] += f*sa; 
            
            break;
        }
        case (I_POLYGON_ALIGN):
        {
            polygon_interaction(particle[i], particle[k], ca, sa, distance, ffm, 0);
            force[0] = KSPRING_OBSTACLE*ffm[0];
            force[1] = KSPRING_OBSTACLE*ffm[1];
            *torque = KTORQUE*ffm[2];
            
            f = 0.1*krepel*lennard_jones_force(distance, particle[k]);
            force[0] += f*ca;
            force[1] += f*sa; 
            
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
            case (I_SEGMENT):
            {
                /* Do nothing, torque already computed */
                break;
            }
            case (I_SEGMENT_CHARGED):
            {
                /* Do nothing, torque already computed */
                break;
            }
            case (I_POLYGON):
            {
                /* Do nothing, torque already computed */
                break;
            }
            case (I_POLYGON_ALIGN):
            {
                spin_f = particle[i].spin_freq;
                if (twrapx||twrapy) 
                    *torque += sin(spin_f*(-particle[k].angle - particle[i].angle))/(1.0e-8 + dist_scaled*dist_scaled);
                else 
                {
                    if (NPOLY%2 == 0) 
                        *torque += sin(spin_f*(particle[k].angle - particle[i].angle))/(1.0e-8 + dist_scaled*dist_scaled);
                    else 
                        *torque -= sin(spin_f*(particle[k].angle - particle[i].angle))/(1.0e-8 + dist_scaled*dist_scaled);

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


void print_partners(int i, t_particle particle[NMAXCIRCLES])
/* print partner list and properties, for debugging */
{
    int p;
    
    if (particle[i].active)
    {
//         printf("Particle %i, type %i, charge %.3lg\n", i, particle[i].type, particle[i].charge);
        fprintf(lj_log, "Particle %i, type %i, charge %.3lg\n", i, particle[i].type, particle[i].charge);
    
//         printf("(x, y) = (%.3lg, %.3lg), radius = %.3lg\n", particle[i].xc, particle[i].yc, particle[i].radius);
    
        printf("%i partners: ", particle[i].npartners);
        fprintf(lj_log, "%i partners: ", particle[i].npartners);
    
        for (p=0; p<particle[i].npartners; p++)
        {
            printf("%i  ", particle[i].partner[p]);
            fprintf(lj_log, "%i  ", particle[i].partner[p]);
        }
        
        printf("\n\n");
        fprintf(lj_log, "\n\n");
    }
}



void add_particle_inpair(int i, int nlist[NMAXPARTNERS], double angle0, int pairing_type, int npartners, int newtype, int rtype, double newr, double newmass, double newcharge, int narms, t_particle particle[NMAXCIRCLES])
/* add new particles in list nlist to old particle i, for paired particles */
{
    int n, k, l, p, p1, q, counter, oldtype, closeby = 0, different = 1, armlen, arm, kplus, kminus, kpplus, kmminus;
    short int no_change = 0;
    double angle, dist, oldr, oldmass_inv, oldcharge, x, y, alpha, beta, beta2, r;
        
    oldtype = particle[i].type;
    oldr = particle[i].radius;
    oldmass_inv = particle[i].mass_inv;
    oldcharge = particle[i].charge;
    
    if (pairing_type == POLY_SOAP_NMIX) 
    {
        if ((double)rand()/RAND_MAX < THIRD_TYPE_PROPORTION) no_change = 1;
    }
    
    armlen = npartners/narms;
    if (rtype == 0) rtype = 3 + rand()%4;

    particle[i].npartners = npartners;
    particle[i].p0 = i;
    particle[i].p1 = nlist[0];
    
    /* shift center of polymer for long chains */
    if ((pairing_type == POLY_SOAP_B)||(pairing_type == POLY_SOAP_N)||(pairing_type == POLY_SOAP_NMIX)||(pairing_type == POLY_PLUSMINUS))
    {
        dist = 0.5*PAIR_DRATIO*oldr*(double)(npartners);
        particle[i].xc -= dist*cos(angle);
        particle[i].yc -= dist*cos(angle);
    }
    /* set some common values for POLY_SEG_POLYGON pairing */
    else if (pairing_type == POLY_SEG_POLYGON)
    {
//         alpha = PI/(double)(npartners+1);
        alpha = PI/(double)npartners;
        r = MU/tan(alpha);
        particle[i].angle = 0.0;
//         particle[i].radius *= 0.75;
    }
    
    if (pairing_type == POLY_PLUSMINUS)
    {
        particle[i].type = newtype;
        particle[i].radius = newr;
    }
    
    if (pairing_type == POLY_DNA_ALT) particle[i].angle = angle0;
    
    /* randomize sign of charge */
    if ((pairing_type == POLY_STAR_CHARGED)||(pairing_type == POLY_POLYGON))
        if (rand()%2 == 1) newcharge *= -1.0;
    
    for (k=0; k<npartners; k++)
    {
        n = nlist[k];
        
        particle[n].type = newtype;
        
        particle[n].thermostat = 1;
        particle[n].charge = newcharge;
        particle[n].active = 1;
        particle[n].cluster = particle[i].cluster;
        
        particle[n].p0 = i;
        particle[n].p1 = nlist[0];
        
        /* set first partner */
        /* TODO; adapt mol_angle to other pairing types? */
//         particle[i].mol_angle = 1;
        if ((pairing_type == POLY_SOAP_B)||(pairing_type == POLY_SOAP_N)||(pairing_type == POLY_SOAP_NMIX)||(pairing_type == POLY_PLUSMINUS))
        {
            particle[i].npartners = 1;
            particle[i].partner[0] = nlist[0];
        }
        else if ((pairing_type == POLY_HYDRA)||(pairing_type == POLY_HYDRA_RIGID))
        {
            particle[i].radius = newr;
            particle[i].npartners = narms;
            for (l=0; l<narms; l++) particle[i].partner[l] = nlist[l*armlen];
//             particle[i].mol_angle = narms;
        }
        else particle[i].partner[k] = n;
        
//         particle[k].mol_angle = particle[i].mol_angle;
//         printf("Particle %i mol_angle = %i\n", i, particle[k].mol_angle);
            
        /* compute positions of new particles */
        switch (pairing_type) {
            case (POLY_STAR_CHARGED): 
            {
                particle[i].radius = 2.0*MU;
                angle = angle0 + DPI*(double)k/(double)npartners;
                dist = (oldr + newr)*PAIR_DRATIO;
                particle[n].radius = newr;
                particle[n].mass_inv = 1.0/newmass;
                if (k%2 == 1) particle[n].charge *= -1.0;
                break;
            }
            case (POLY_POLYGON): 
            {
                particle[i].radius = 2.0*MU;
                particle[n].mass_inv = 1.0/newmass;
                
                if (k < npartners/2)
                {
                    angle = angle0 + DPI*(double)(2*k+1)/(double)npartners;
                    dist = (oldr + newr)*PAIR_DRATIO*cos(DPI/(double)npartners)*1.2;
                    particle[n].radius = newr*0.75;
                    particle[n].charge = 0.0;
                }
                else
                {
                    angle = angle0 + DPI*(double)(2*k)/(double)npartners;
                    dist = (oldr + newr)*PAIR_DRATIO;
                    particle[n].radius = newr;
                    if (k%2 == 1) particle[n].charge *= -1.0;
                }
                break;
            }
            case (POLY_KITE): 
            {
                particle[i].radius = 1.25*MU*sin(PI/5.0);
                particle[i].charge = 0.0;
                particle[n].mass_inv = 1.0/newmass;
                
                switch (k){
                    case (0):
                    {
                        angle = PI;
                        dist = 0.5*MU;
                        particle[n].radius = 0.4*MU;
                        particle[n].charge = 0.0;
                        break;
                    }
                    case (1):
                    {
                        angle = 0.0;
                        dist = 2.0*MU*cos(DPI/5.0);
                        particle[n].radius = 0.3*MU;
                        break;
                    }
                    case (2):
                    {
                        angle = DPI/5.0;
                        dist = MU;
                        particle[n].radius = 0.3*MU;
                        particle[n].charge *= -1.0;
                        break;
                    }
                    case (3):
                    {
                        angle = PI;
                        dist = MU;
                        particle[n].radius = 0.3*MU;
                        break;
                    }
                    case (4):
                    {
                        angle = -DPI/5.0;
                        dist = MU;
                        particle[n].radius = 0.3*MU;
                        particle[n].charge *= -1.0;
                        break;
                    }
                }
                break;
            }
            case (POLY_DART): 
            {
                particle[i].radius = 0.3*MU;
                particle[i].charge = -newcharge;
                particle[n].mass_inv = 1.0/newmass;
                
                switch (k){
                    case (0):
                    {
                        angle = 3.0*PI/10.0;
                        dist = MU*sin(PI/5.0);
                        particle[n].radius = 0.3*MU;
                        particle[n].charge = 0.0;
                        break;
                    }
                    case (1):
                    {
                        angle = -3.0*PI/10.0;
                        dist = MU*sin(PI/5.0);
                        particle[n].radius = 0.3*MU;
                        particle[n].charge = 0.0;
                        break;
                    }                    
                    case (2):
                    {
                        angle = 0.0;
                        dist = MU;
                        particle[n].radius = 0.3*MU;
                        particle[n].charge *= -1.0;
                        break;
                    }
                    case (3):
                    {
                        angle = 3.0*PI/5.0;
                        dist = MU;
                        particle[n].radius = 0.3*MU;
                        break;
                    }
                    case (4):
                    {
                        angle = -3.0*PI/5.0;
                        dist = MU;
                        particle[n].radius = 0.3*MU;
                        break;
                    }
                }
                break;
            }
            case (POLY_SEG_POLYGON):
            {
                dist = r;
                angle = (double)(2*k)*alpha;
                particle[n].angle = angle + PID;
                particle[n].radius = MU;
                particle[n].charge = 0.0;
                fprintf(lj_log, "k = %i, angle = %.3lg, orientation = %.3lg\n", k, angle, particle[n].angle);
                break;
            }
            case (POLY_WATER): 
            {
                angle = angle0 + (double)(k+1)*PARTNER_ANGLE*PI/180.0;
                dist = oldr + newr;
                particle[n].radius = newr;
                particle[n].mass_inv = 1.0/newmass;
                break;
            }
            case (POLY_SOAP):
            {
                if (k%2 == 0) angle = angle0;
                else angle = angle0 + PI;
                dist = 2.0*PAIR_DRATIO*oldr*(double)(k/2+1);
                particle[i].charge = 0.0;
                if (k==npartners-1) 
                {
                    particle[n].charge = oldcharge;
                    particle[n].type = newtype;
                }
                else 
                {
                    particle[n].charge = 0.0;
                    particle[n].type = oldtype;
                }
                if (k == npartners-1) particle[n].radius = newr;
                else particle[n].radius = oldr;
                particle[n].mass_inv = oldmass_inv;
                break;
            }
            case (POLY_SOAP_B):
            {
                angle = angle0;
                dist = 2.0*PAIR_DRATIO*oldr*(double)(k+1);
                particle[i].charge = 0.0;
                if (k==npartners-1) 
                {
                    particle[n].charge = oldcharge;
                    particle[n].type = newtype;
                }
                else 
                {
                    particle[n].charge = 0.0;
                    particle[n].type = oldtype;
                }
                if (k == npartners-1) particle[n].radius = newr;
                else particle[n].radius = oldr;
                particle[n].mass_inv = oldmass_inv;
                break;
            }
            case (POLY_SOAP_N):
            {
                angle = angle0;
                dist = 2.0*PAIR_DRATIO*oldr*(double)(k+1);
                particle[i].charge = 0.0;
                if (k==npartners-1) 
                {
                    particle[n].charge = oldcharge;
                    particle[n].type = newtype;
                    dist += 2.0*PAIR_DRATIO*oldr;
                }
                else if (k==npartners-2) 
                {
                    particle[n].charge = -oldcharge;
                    particle[n].type = newtype;
                    dist += PAIR_DRATIO*oldr;
                }
                else 
                {
                    particle[n].charge = 0.0;
                    particle[n].type = oldtype;
                }
                if ((k == npartners-1)||(k == npartners-2)) particle[n].radius = newr;
                else particle[n].radius = oldr;
                particle[n].mass_inv = oldmass_inv;
                break;
            }
            case (POLY_SOAP_NMIX):
            {
                angle = angle0;
                dist = 2.0*PAIR_DRATIO*oldr*(double)(k+1);
                particle[i].charge = 0.0;
                if (no_change)  /* only modify ends with probability 1/2 */
                {
                    particle[n].charge = 0.0;
                    particle[n].type = oldtype + 3;
                    particle[n].radius = oldr;
                    particle[i].charge = 0.0;
                    particle[i].type = oldtype + 3;
                    particle[i].radius = oldr;
                }
                else
                {
                    if (k==npartners-1) 
                    {
                        particle[n].charge = oldcharge;
                        particle[n].type = newtype;
                        dist += 2.0*PAIR_DRATIO*oldr;
                    }
                    else if (k==npartners-2) 
                    {
                        particle[n].charge = -oldcharge;
                        particle[n].type = newtype;
                        dist += PAIR_DRATIO*oldr;
                    }
                    else 
                    {
                        particle[n].charge = 0.0;
                        particle[n].type = oldtype;
                    }
                    if ((k == npartners-1)||(k == npartners-2)) particle[n].radius = newr;
                    else particle[n].radius = oldr;
                }
                particle[n].mass_inv = oldmass_inv;
                break;
            }
            case (POLY_PLUSMINUS):
            {
                angle = angle0;
                dist = 2.0*PAIR_DRATIO*oldr*(double)(k+1);
                if (k==npartners-1) 
                {
                    particle[n].charge = oldcharge;
                    particle[n].type = newtype;
                }
                else 
                {
                    particle[n].charge = 0.0;
                    particle[n].type = oldtype;
                }
                if (k == npartners-1) particle[n].radius = newr;
                else particle[n].radius = oldr;
                particle[n].mass_inv = oldmass_inv;
                break;
            }
            case (POLY_HYDRA):
            {
                arm = k/armlen;
                p = k%armlen;
                angle = angle0 + DPI*(double)arm/(double)narms;
                dist = PAIR_DRATIO*(newr + oldr + 2.0*oldr*(double)p);
                if (p == armlen-1) 
                {
                    particle[n].charge = newcharge;
                    if ((ALTERNATE_POLY_CHARGE)||(arm%2 == 1)) particle[n].charge *= -1.0;
                    particle[n].type = newtype;
                    particle[n].radius = newr;
                    particle[n].mass_inv = 1.0/newmass;
                }
                else 
                {
                    particle[n].charge = 0.0;
                    particle[n].type = oldtype;
                    particle[n].radius = oldr;
                    particle[n].mass_inv = oldmass_inv;
                }
                break;
            }
            case (POLY_HYDRA_RIGID):    /* same as POLY_HYDRA */
            {
                arm = k/armlen;
                p = k%armlen;
                angle = angle0 + DPI*(double)arm/(double)narms;
                dist = PAIR_DRATIO*(newr + oldr + 2.0*oldr*(double)p);
                if (p == armlen-1) 
                {
                    particle[n].charge = newcharge;
                    if ((ALTERNATE_POLY_CHARGE)&&(arm%2 == 0)) particle[n].charge *= -1.0;
                    particle[n].type = newtype;
                    particle[n].radius = newr;
                    particle[n].mass_inv = 1.0/newmass;
                }
                else 
                {
                    particle[n].charge = 0.0;
                    particle[n].type = oldtype;
                    particle[n].radius = oldr;
                    particle[n].mass_inv = oldmass_inv;
                }
                break;
            }
            case (POLY_DNA):
            {
                angle = angle0;
                dist = PAIR_DRATIO*(newr + oldr);
                particle[n].radius = newr;
                particle[n].mass_inv = 1.0/newmass;
                particle[n].type = 0;
                
                switch (k) {
                    case (0):
                    {
                        angle -= PID;
                        particle[n].type = 2;
                        break;
                    }
                    case (1):
                    {
                        angle += PID;
                        particle[n].type = 2;
                        break;
                    }
                    default: 
                    {
                        dist += (double)(k-2)*2.0*PAIR_DRATIO*newr;
                        if (k == npartners - 1) particle[n].type = 3 + rand()%4;
                    }
                }
                break;
            }
            case (POLY_DNA_ALT):
            {
                angle = angle0;
                dist = PAIR_DRATIO*(newr + oldr);
                particle[n].radius = newr;
                particle[n].mass_inv = 1.0/newmass;
                particle[n].type = 0;
                particle[n].angle = angle0;
                
                switch (k) {
                    case (0):
                    {
                        angle -= PID;
                        break;
                    }
                    case (1):
                    {
                        angle -= PID;
                        dist += 2.0*PAIR_DRATIO*newr;
                        particle[n].type = 1;
                        break;
                    }
                    case (2):
                    {
                        angle += PID;
                        break;
                    }
                    case (3):
                    {
                        angle += PID;
                        dist += 2.0*PAIR_DRATIO*newr;
                        particle[n].type = 2;
                        break;
                    }
                    default: 
                    {
                        dist += (double)(k-4)*2.0*PAIR_DRATIO*newr;
                        if (k == npartners - 1) particle[n].type = 3 + rand()%4;
                    }
                }
                break;
            }
            case (POLY_DNA_DOUBLE):
            {
                angle = angle0;
                dist = PAIR_DRATIO*(newr + oldr);
                particle[n].radius = newr;
                particle[n].mass_inv = 1.0/newmass;
                particle[n].type = 0;
                particle[n].angle = angle0;
                beta = atan(newr/(dist + 2.0*PAIR_DRATIO*newr));
                beta2 = atan(1.5/((double)(NPARTNERS_DNA-8)*2.0*PAIR_DRATIO));
                
                switch (k) {
                    case (0):
                    {
                        angle -= PID;
                        break;
                    }
                    case (1):
                    {
                        angle -= (PID + DNA_RIGIDITY*beta);
                        dist += 2.0*PAIR_DRATIO*newr;
                        particle[n].type = 1;
//                         particle[n].type = 2;
                        break;
                    }
                    case (2):
                    {
                        angle -= (PID - DNA_RIGIDITY*beta);
                        dist += 2.0*PAIR_DRATIO*newr;
                        particle[n].type = 1;
//                         particle[n].type = 2;
                        break;
                    }
                    case (3):
                    {
                        angle += PID;
                        break;
                    }
                    case (4):
                    {
                        angle += PID + DNA_RIGIDITY*beta;
                        dist += 2.0*PAIR_DRATIO*newr;
                        particle[n].type = 2;
                        break;
                    }
                    case (5):
                    {
                        angle += PID - DNA_RIGIDITY*beta;
                        dist += 2.0*PAIR_DRATIO*newr;
                        particle[n].type = 2;
                        break;
                    }
                    case (NPARTNERS_DNA-2):
                    {
                        dist += (double)(NPARTNERS_DNA-8)*2.0*PAIR_DRATIO*newr;
                        angle += 0.5*beta2;
                        particle[n].type = rtype;
                        break;
                    }
                    case (NPARTNERS_DNA-1):
                    {
                        dist += (double)(NPARTNERS_DNA-8)*2.0*PAIR_DRATIO*newr;
                        angle -= 0.5*beta2;
                        particle[n].type = rtype;
                        break;
                    }
                    default: 
                    {
                        dist += (double)(k-6)*2.0*PAIR_DRATIO*newr;
                    }
                }
                break;
            }
            case (POLY_DNA_FLEX):
            {
                angle = angle0;
                dist = PAIR_DRATIO*(newr + oldr);
                particle[n].radius = newr;
                particle[n].mass_inv = 1.0/newmass;
                particle[n].type = 0;
                particle[n].angle = angle0;
                beta = atan(newr/(dist + 2.0*PAIR_DRATIO*newr));
                beta2 = atan(1.5/((double)(NPARTNERS_DNA-6)*2.0*PAIR_DRATIO));
                
                switch (k) {
                    case (0):
                    {
                        angle -= PID;
                        break;
                    }
                    case (1):
                    {
                        angle -= PID;
                        dist += 2.0*PAIR_DRATIO*newr;
                        particle[n].type = 1;
                        break;
                    }
                    case (2):
                    {
                        angle += PID;
                        break;
                    }
                    case (3):
                    {
                        angle += PID;
                        dist += 2.0*PAIR_DRATIO*newr;
                        particle[n].type = 2;
                        break;
                    }
                    case (NPARTNERS_DNA-2):
                    {
                        dist += (double)(NPARTNERS_DNA-6)*2.0*PAIR_DRATIO*newr;
                        angle += 0.5*beta2;
                        particle[n].type = rtype;
                        break;
                    }
                    case (NPARTNERS_DNA-1):
                    {
                        dist += (double)(NPARTNERS_DNA-6)*2.0*PAIR_DRATIO*newr;
                        angle -= 0.5*beta2;
                        particle[n].type = rtype;
                        break;
                    }
                    default: 
                    {
                        dist += (double)(k-4)*2.0*PAIR_DRATIO*newr;
                    }
                }
                break;
            }
            default: 
            {
                angle = angle0 + DPI*(double)k/(double)npartners;
                dist = oldr + newr;
                particle[n].radius = newr;
                particle[n].mass_inv = 1.0/newmass;
//                 particle[n].partner_eqd[0] = (MU + MU_B)*PAIR_DRATIO;
//                 particle[i].partner_eqd[k] = (MU + MU_B)*PAIR_DRATIO;
            }

        }
        particle[n].xc = particle[i].xc + dist*cos(angle);
        particle[n].yc = particle[i].yc + dist*sin(angle);
        
//         printf("Distance %.3lg, (xi,yi) = (%.3lg, %.3lg), (xn,yn) = (%.3lg, %.3lg)\n", 
//                dist, particle[i].xc, particle[i].yc, particle[n].xc, particle[n].yc);
    
        /* adjust list of partners */
        switch (pairing_type) {
            case (POLY_STAR):
            {
                particle[n].npartners = 1;
                particle[n].partner[0] = i;
                break;
            }
            case (POLY_STAR_CHARGED):
            {
                particle[n].npartners = 3;
                particle[n].partner[0] = i;
                if (k == 0) particle[n].partner[1] = nlist[npartners-1];
                else particle[n].partner[1] = nlist[k-1];
                if (k == npartners-1) particle[n].partner[2] = nlist[0];
                else particle[n].partner[2] = nlist[k+1];
                break;
            }
            case (POLY_POLYGON):
            {
                particle[n].npartners = npartners;
                particle[n].partner[0] = i;
                p = 0;
                p1 = 1;
                while (p < npartners)
                {
                    q = nlist[p];
                    if (q != n)
                    {
                        particle[n].partner[p1] = q;
                        p1++;
                    }
                    p++;
                }
                particle[n].npartners = p1;
                break;
            }
            case (POLY_KITE):
            {
                particle[n].npartners = npartners;
                particle[n].partner[0] = i;
                p = 0;
                p1 = 1;
                while (p < npartners)
                {
                    q = nlist[p];
                    if (q != n)
                    {
                        particle[n].partner[p1] = q;
                        p1++;
                    }
                    p++;
                }
                particle[n].npartners = p1;
                break;
            }
            case (POLY_SEG_POLYGON):
            {
                particle[n].npartners = 3;
                particle[n].partner[0] = i;
                if (k == 0) particle[n].partner[1] = nlist[npartners-1];
                else particle[n].partner[1] = nlist[k-1];
                if (k == npartners-1) particle[n].partner[2] = nlist[0];
                else particle[n].partner[2] = nlist[k+1];
                break;
            }
            case (POLY_SOAP_B):
            {
                if (k==npartners-1)
                {
                    particle[n].npartners = 1;
                    particle[n].partner[0] = nlist[npartners-2];
                }
                else if (k==0)
                {
                    particle[n].npartners = 2;
                    particle[n].partner[0] = i;
                    particle[n].partner[1] = nlist[1];
                }
                else
                {
                    particle[n].npartners = 2;
                    particle[n].partner[0] = nlist[k-1];
                    particle[n].partner[1] = nlist[k+1];
                }
                break;
            }
            case (POLY_SOAP_N):
            {
                if (k==npartners-1)
                {
                    particle[n].npartners = 1;
                    particle[n].partner[0] = nlist[npartners-2];
                }
                else if (k==0)
                {
                    particle[n].npartners = 2;
                    particle[n].partner[0] = i;
                    particle[n].partner[1] = nlist[1];
                }
                else
                {
                    particle[n].npartners = 2;
                    particle[n].partner[0] = nlist[k-1];
                    particle[n].partner[1] = nlist[k+1];
                }
                break;
            }
            case (POLY_SOAP_NMIX):
            {
                if (k==npartners-1)
                {
                    particle[n].npartners = 1;
                    particle[n].partner[0] = nlist[npartners-2];
                }
                else if (k==0)
                {
                    particle[n].npartners = 2;
                    particle[n].partner[0] = i;
                    particle[n].partner[1] = nlist[1];
                }
                else
                {
                    particle[n].npartners = 2;
                    particle[n].partner[0] = nlist[k-1];
                    particle[n].partner[1] = nlist[k+1];
                }
                break;
            }
            case (POLY_PLUSMINUS):
            {
                if (k==npartners-1)
                {
                    particle[n].npartners = 1;
                    particle[n].partner[0] = nlist[npartners-2];
                }
                else if (k==0)
                {
                    particle[n].npartners = 2;
                    particle[n].partner[0] = i;
                    particle[n].partner[1] = nlist[1];
                }
                else
                {
                    particle[n].npartners = 2;
                    particle[n].partner[0] = nlist[k-1];
                    particle[n].partner[1] = nlist[k+1];
                }
                break;
            }
            case (POLY_HYDRA):
            {
                p = k%armlen;
                if (p == armlen-1)
                {
                    particle[n].npartners = 1;
                    particle[n].partner[0] = nlist[k-1];
                }
                else if (p == 0)
                {
                    particle[n].npartners = 2;
                    particle[n].partner[0] = i;
                    particle[n].partner[1] = nlist[k+1];
                }
                else
                {
                    particle[n].npartners = 2;
                    particle[n].partner[0] = nlist[k-1];
                    particle[n].partner[1] = nlist[k+1];
                }
                break;
            }
            case (POLY_HYDRA_RIGID):
            {
                p = k%armlen;
                if (p == armlen-1)
                {
                    particle[n].npartners = 4;
                    particle[n].partner[0] = nlist[k-1];
                    particle[n].partner[1] = i;
                    particle[n].partner[2] = nlist[(k+armlen)%NPARTNERS];
                    q = k-armlen;
                    if (q < 0) q += NPARTNERS;
                    particle[n].partner[3] = nlist[q];
                }
                else if (p == 0)
                {
                    particle[n].npartners = 4;
                    particle[n].partner[0] = i;
                    particle[n].partner[1] = nlist[k+1];
                    particle[n].partner[2] = nlist[(k+armlen)%NPARTNERS];
                    q = k-armlen;
                    if (q < 0) q += NPARTNERS;
                    particle[n].partner[3] = nlist[q];
                }
                else
                {
                    particle[n].npartners = 2;
                    particle[n].partner[0] = nlist[k-1];
                    particle[n].partner[1] = nlist[k+1];
                }
                break;
            }
            default:
            {
                particle[n].npartners = npartners;
                particle[n].partner[0] = i;
                counter = 1;
                for (l=0; l<npartners; l++) if (l != k)
                {
                    particle[n].partner[counter] = nlist[l];
                    counter++;
                }
            }
         }
    }
        
    /* adjust equilibrium distances */
    if ((pairing_type == POLY_ALL))
    {
        for (l=0; l<particle[i].npartners; l++) 
            particle[i].partner_eqd[l] = (oldr + newr)*PAIR_DRATIO;
        for (k=0; k<npartners; k++)
        {
            n = nlist[k];
            for (l=0; l<particle[n].npartners; l++)
                particle[n].partner_eqd[l] = (oldr + newr)*PAIR_DRATIO;
        }
    }
    else
    {
        for (l=0; l<particle[i].npartners; l++)
        {
            p = particle[i].partner[l];
            dist = module2(particle[i].xc - particle[p].xc, particle[i].yc - particle[p].yc);
            particle[i].partner_eqd[l] = dist;
            particle[i].partner_eqa[l] = particle[p].angle - particle[i].angle;
        }
        for (k=0; k<npartners; k++)
        {
            n = nlist[k];
            for (l=0; l<particle[n].npartners; l++)
            {
                p = particle[n].partner[l];
                dist = module2(particle[n].xc - particle[p].xc, particle[n].yc - particle[p].yc);
                particle[n].partner_eqd[l] = dist;
                particle[n].partner_eqa[l] = particle[p].angle - particle[n].angle;
            }
        }
    }
    
    /* set molecule numbers */
    particle[i].molecule = i;
    for (k=0; k<npartners; k++) 
    {
        n = particle[i].partner[k];
        particle[n].molecule = i;
        particle[n].added = particle[i].added;
        particle[n].coulomb = particle[i].coulomb;
        particle[n].reactive = particle[i].reactive;
//         printf("i = %i, added = %i, coulomb = %i\t nb = %i, added = %i, coulomb = %i\n", i, particle[i].added, particle[i].coulomb, n, particle[n].added, particle[n].coulomb);
    }
    
    /* test for presence of other molecules that are too close */
    if (DEACIVATE_CLOSE_PAIRS)
    {
        closeby = 0;
        for (k=0; k<NMAXCIRCLES; k++) if (particle[k].active)
        {
            different = (k!=i);
            for (l=0; l<npartners; l++) different *= (k!=nlist[l]); 
            if (different) for (l=0; l<npartners; l++) 
            {
                n = nlist[l];
                for (p=-1; p<2; p++) for (q=-1; q<2; q++)
                {
                    x = particle[k].xc + p*(BCXMAX - BCXMIN);
                    y = particle[k].yc + q*(BCYMAX - BCYMIN);
                    dist = module2(x - particle[n].xc, y - particle[n].yc);
                    closeby += (dist < PAIR_SAFETY_FACTOR*(particle[l].radius + particle[k].radius));
                }
            }
        }
        if (closeby > 0) 
        {
            printf("Deactivating particle cluster %i\n", i);
            particle[i].active = 0;
            for (l=0; l<npartners; l++) particle[nlist[l]].active = 0;
        }
    }
}

void update_single_molecule_data(int mol, t_particle particle[NMAXCIRCLES], t_molecule molecule[NMAXCIRCLES])
/* update connection data for a single molecule */
{
    int np, p, pp, nq, q, qq, molq, npartners = 0, new, connection;
    
//     molecule[mol].npartners = 0;
    
    printf("Resetting partner list of molecule %i\n", mol);
    fprintf(lj_log, "[update_single_molecule_data] Resetting partner list of molecule %i\n", mol);
    
    fprintf(lj_log, "[update_single_molecule_data] Before reset, molecule %i has %i partners: ", mol, molecule[mol].npartners);
    for (p=0; p<molecule[mol].npartners; p++)
    {
        printf("%i(%i) ", molecule[mol].partner[p], molecule[mol].connection_type[p]);
        fprintf(lj_log, "%i(%i) ", molecule[mol].partner[p], molecule[mol].connection_type[p]);
    }
    printf("\n");
    fprintf(lj_log, "\n");
    
    np = molecule[mol].nparticles;
    
    for (p=0; p<np; p++)
    {
        pp = molecule[mol].particle[p];
        nq = particle[pp].npartners;
        
        for (q=0; q<nq; q++)
        {
            qq = particle[pp].partner[q];
            molq = particle[qq].molecule;
            if (molq != mol)
            {
                new = 1;
                fprintf(lj_log, "Testing molecule %i\n", molq);
                /* test if connection has already been found */
                for (connection = 0; connection < npartners; connection++)
                    if (molecule[mol].partner[connection] == molq) new = 0;
                if (new)
                {
                    molecule[mol].partner[npartners] = molq;
                    molecule[mol].connection_type[npartners] = particle[pp].type;
                    npartners++;
                    fprintf(lj_log, "Added molecule %i(%i) as partner %i\n", molq, particle[pp].type, npartners);
                }
                else fprintf(lj_log, "Molecule %i not added\n", molq);
            }
        }
    }
    molecule[mol].npartners = npartners;
    
    printf("Molecule %i has %i partners: ", mol, npartners);
    fprintf(lj_log, "[update_single_molecule_data] After reset, molecule %i has %i partners: ", mol, npartners);
    for (p=0; p<npartners; p++)
    {
        printf("%i(%i) ", molecule[mol].partner[p], molecule[mol].connection_type[p]);
        fprintf(lj_log, "%i(%i) ", molecule[mol].partner[p], molecule[mol].connection_type[p]);
    }
    printf("\n");
    fprintf(lj_log, "\n");
}


void init_molecule_data(t_particle particle[NMAXCIRCLES], t_molecule molecule[NMAXCIRCLES])
/* initialize date in molecule structure */
{
    int i, m, np;
    
    printf("Initializing molecule structure\n");
    #pragma omp parallel for private(m)
    for (m=0; m<NMAXCIRCLES; m++)
    {
        molecule[m].nparticles = 0;
        molecule[m].npartners = 0;
        molecule[m].added = 0;
    }
    
//     printf("1\n");
    
    for (i=0; i<ncircles; i++)
    {
//         printf("i = %i\n", i);
        m = particle[i].molecule;
        if (m + 1 > nmolecules) nmolecules = m + 1;
        np = molecule[m].nparticles;
//         printf("np = %i\n", np);
        if (np < NPARTNERS+1)
        {
            molecule[m].particle[np] = i;
            molecule[m].nparticles++;
        }
        molecule[m].added = particle[i].added;
    }  
    printf("Found %i molecules\n", nmolecules);
    fprintf(lj_log, "Found %i molecules\n", nmolecules);
    
    /* for debugging */
    for (m=0; m<nmolecules; m++)
    {
        printf("Molecule %i has %i particles: \n", m, molecule[m].nparticles);
        fprintf(lj_log, "Molecule %i has %i particles: \n", m, molecule[m].nparticles);
        for (i=0; i<molecule[m].nparticles; i++)
        {
            printf(" %i-%i |" , molecule[m].particle[i], particle[molecule[m].particle[i]].molecule);
            fprintf(lj_log, " %i-%i |" , molecule[m].particle[i], particle[molecule[m].particle[i]].molecule);
        }
        printf("\n");
        fprintf(lj_log, "\n");
        printf("Molecule %i.added = %i\n", m, molecule[m].added);
    }
}

void init_particle_pairs(t_particle particle[NMAXCIRCLES], t_molecule molecule[NMAXCIRCLES])
/* initialize data structure for paired particles */
{
    int i, k, l, n, p, q, counter, nlist[NMAXPARTNERS], rtype = 0, newrtype;
    double angle, dist;
    
    if (ncircles*NMAXPARTNERS > NMAXCIRCLES)
    {
        printf("Error: NMAXCIRCLES is too small\n");
        exit(1); 
    }
    
    nmolecules = ncircles;
    
    for (i=0; i<ncircles; i++) 
    {
        if ((particle[i].active)&&(particle[i].type == 0))
        {
            if (RANDOMIZE_ANGLE) angle = DPI*(double)rand()/(double)RAND_MAX;
            else if ((RD_INITIAL_COND == IC_DNA_POLYMERASE)||(RD_INITIAL_COND == IC_DNA_POLYMERASE_REC))
            {
                if (i%NGRIDY == NGRIDY-2) 
                {
                    angle = PID;
                    rtype = 3 + rand()%4;
                }
                else if (i%NGRIDY == NGRIDY-1) 
                {
                    angle = -PID;
//                     printf("old type = %i\n", rtype);
                    switch (rtype) {
                        case (3):
                        {
                            newrtype = 4;
                            break;
                        }
                        case (4):
                        {
                            newrtype = 3;
                            break;
                        }
                        case (5):
                        {
                            newrtype = 6;
                            break;
                        }
                        case (6):
                        {
                            newrtype = 5;
                            break;
                        }
                    }
//                     printf("new type = %i\n", newrtype); 
                    rtype = newrtype; 
                }
                printf("i = %i, type = %i\n", i, rtype); 
            }
            else angle = 0.0;
            for (k=0; k<NPARTNERS; k++) nlist[k] = ncircles*(k+1) + i;
            if (PAIR_TYPEB_PARTICLES) add_particle_inpair(i, nlist, angle, PAIRING_TYPE, NPARTNERS, 2, rtype, MU_C, PARTICLE_MASS_C, CHARGE_C, NARMS, particle);
            else add_particle_inpair(i, nlist, angle, PAIRING_TYPE, NPARTNERS, 2, rtype, MU_B, PARTICLE_MASS_B, CHARGE_B, NARMS, particle);
        }
    }
    ncircles += ncircles*NPARTNERS;
    
    if (PAIR_TYPEB_PARTICLES) 
    {
        for (i=0; i<ncircles; i++) 
        {
            if ((particle[i].active)&&(particle[i].type == 1))
            {
                if (RANDOMIZE_ANGLE) angle = DPI*(double)rand()/(double)RAND_MAX;
                else angle = 0.0;
                for (k=0; k<NPARTNERS_B; k++) nlist[k] = ncircles*(k+1) + i;
                add_particle_inpair(i, nlist, angle, PAIRING_TYPE_B, NPARTNERS_B, 3, rtype, MU_D, PARTICLE_MASS_D, CHARGE_D, NARMS_B, particle);
            }
        }
        ncircles += ncircles*NPARTNERS;
    }
    
    init_molecule_data(particle, molecule);
    printf("Done\n");
    
    /* find first two partner numbers, for P_MOL_ANGLE color scheme */
//     for (i=0; i<ncircles; i++) 
//     {
//         p = i;
//         for (k=0; k<particle[i].npartners; k++)
//             if (particle[i].partner[k] < p) p = particle[i].partner[k];
//         particle[i].p0 = p;
//         particle[i].p1 = particle[p].partner[0];
//         for (k=0; k<particle[i].npartners; k++)
//         {
//             q = particle[i].partner[k];
//             if (q > p) 
//             {
//                 particle[q].p0 = p;
//                 particle[q].p1 = particle[p].partner[0];
//             }
//         }
//     }
    
    printf("ncircles = %i\n", ncircles);
    
    for (i=0; i<ncircles; i++) print_partners(i, particle);
}

int add_particle(double x, double y, double vx, double vy, double mass, short int type, t_particle particle[NMAXCIRCLES])
{
    int i, k, l, n, closeby = 0, nlist[NMAXPARTNERS];
    double dist, angle;
    
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
        particle[i].radius = MU;
//         particle[i].radius = MU*sqrt(mass);
        particle[i].active = 1;
        particle[i].neighb = 0;
        particle[i].diff_neighb = 0;
        particle[i].thermostat = 1;
        particle[i].charge = CHARGE;

        particle[i].energy = 0.0;
        particle[i].emean = 0.0;
        particle[i].dirmean = 0.0;
        
        if (RANDOM_RADIUS) particle[i].radius = particle[i].radius*(RANDOM_RADIUS_MIN + RANDOM_RADIUS_RANGE*((double)rand()/RAND_MAX));
        
        particle[i].mass_inv = 1.0/mass;
        if (particle[i].type == 0) particle[i].inertia_moment_inv = 1.0/PARTICLE_INERTIA_MOMENT;
        else particle[i].inertia_moment_inv = 1.0/PARTICLE_INERTIA_MOMENT;
        
        if ((RANDOM_RADIUS)&&(ADAPT_MASS_TO_RADIUS > 0.0)) 
            particle[i].mass_inv *= 1.0/(1.0 + pow(particle[i].radius/MU, ADAPT_MASS_TO_RADIUS));
        if ((RANDOM_RADIUS)&&(ADAPT_DAMPING_TO_RADIUS > 0.0))
            particle[i].damping = 1.0 + ADAPT_DAMPING_FACTOR*pow(particle[i].radius/MU, ADAPT_DAMPING_TO_RADIUS);
       
                
        particle[i].vx = vx;
        particle[i].vy = vy;
        particle[i].energy = (particle[i].vx*particle[i].vx + particle[i].vy*particle[i].vy)*particle[i].mass_inv;
        
        particle[i].angle = DPI*(double)rand()/RAND_MAX;
        particle[i].omega = OMEGA_INITIAL*gaussian();
        
//         printf("Particle[%i].omega = %.4lg\n", i, particle[i].omega); 
        
//         if (particle[i].type == 1)
//         {
//             particle[i].interaction = INTERACTION_B;
//             particle[i].eq_dist = EQUILIBRIUM_DIST_B;
//             particle[i].spin_range = SPIN_RANGE_B;
//             particle[i].spin_freq = SPIN_INTER_FREQUENCY_B;            
//         }
    
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
        
        if ((PAIR_PARTICLES)&&(type == 0))
        {
            angle = DPI*(double)rand()/(double)RAND_MAX;
            for (k=0; k<NPARTNERS; k++) nlist[k] = ncircles + k + 1;
            if (PAIR_TYPEB_PARTICLES) add_particle_inpair(i, nlist, angle, PAIRING_TYPE, NPARTNERS, 2, 0, MU_C, PARTICLE_MASS_C, CHARGE_C, NARMS, particle);
            else add_particle_inpair(ncircles, nlist, angle, PAIRING_TYPE, NPARTNERS, 2, 0, MU_B, PARTICLE_MASS_B, CHARGE_B, NARMS, particle);
            
            ncircles += NPARTNERS+1;
        }
        else ncircles++;

        printf("Added particle at (%.3lg, %.3lg)\n", x, y);
        printf("Number of particles: %i\n", ncircles);
        
        return(1);
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


void compute_particle_colors(t_particle particle, t_cluster cluster[NMAXCIRCLES], int plot, double rgb[3], double rgbx[3], double rgby[3], double *radius, int *width, t_particle other_particle[NMAXCIRCLES])
{
    double ej, angle, hue, huex, huey, lum, x, y, ccluster;
    int i, k, p, q, cl;
    
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
//             hue = 180.0*(1.0 + hue);
            hue = 20.0 + 320.0*hue;
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
        case (P_EMEAN_DENSITY): 
        {
            ej = particle.emean;
            cl = particle.cluster;
            ej *= PARTICLE_MASS/cluster[cl].mass;
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
        case (P_LOG_EMEAN): 
        {
            ej = particle.emean;
            if (ej > 0.0) 
            {
                ej = log(ej/PARTICLE_EMIN)/log(PARTICLE_EMAX/PARTICLE_EMIN);
                if (ej < 0.0) ej = 0.0;
                else if (ej > 1.0) ej = 1.0;
                hue = ENERGY_HUE_MIN + (ENERGY_HUE_MAX - ENERGY_HUE_MIN)*ej;
//                 if (hue > ENERGY_HUE_MIN) hue = ENERGY_HUE_MIN;
//                 if (hue < ENERGY_HUE_MAX) hue = ENERGY_HUE_MAX;
            }
            *radius = particle.radius;
            *width = BOUNDARY_WIDTH;
            break;
        }
        case (P_DIRECT_EMEAN): 
        {
            hue = particle.dirmean + COLOR_HUESHIFT*PI;
//             printf("dirmean = %.3lg\n", particle.dirmean);
            if (hue > DPI) hue -= DPI;
            if (hue < 0.0) hue += DPI;
            hue = PARTICLE_HUE_MIN + (PARTICLE_HUE_MAX - PARTICLE_HUE_MIN)*(hue)/DPI;
            ej = particle.emean;
            if (ej < 0.5*PARTICLE_EMAX) lum = 2.0*ej/PARTICLE_EMAX;
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
        case (P_NPARTNERS):
        {
            hue = partner_color(particle.npartners);
            *radius = particle.radius;
            *width = BOUNDARY_WIDTH;
            break;
        }
        case (P_CHARGE):
        {
            hue = (-tanh(0.5*particle.charge)+1.0)*180.0;
            *radius = particle.radius;
            *width = BOUNDARY_WIDTH;
            break;
        }
        case (P_MOL_ANGLE):
        {
            p = particle.p0;
            q = particle.p1;
            x = other_particle[q].xc - other_particle[p].xc;
            y = other_particle[q].yc - other_particle[p].yc;
            /* deal with periodic boundary conditions */
            if (x > 0.5*(XMAX - XMIN)) x -= (XMAX - XMIN);
            else if (x < 0.5*(XMIN - XMAX)) x += (XMAX - XMIN);
            if (y > 0.5*(YMAX - YMIN)) y -= (YMAX - YMIN);
            else if (y < 0.5*(YMIN - YMAX)) y += (YMAX - YMIN);
            angle = argument(x, y)*MOL_ANGLE_FACTOR;
//             printf("Particle p = %i, mol_angle = %i\n", p, particle.mol_angle);
//             angle = argument(x, y)*(double)particle.mol_angle;
//             angle = argument(x, y)*(double)(other_particle[particle.p0].npartners);
            while (angle > DPI) angle -= DPI;
            while (angle < 0.0) angle += DPI;
            hue = PARTICLE_HUE_MIN + (PARTICLE_HUE_MAX - PARTICLE_HUE_MIN)*(angle)/DPI;
            *radius = particle.radius;
            *width = BOUNDARY_WIDTH;
            break;
        }
        case (P_CLUSTER):
        {
//             cluster = (double)(particle.cluster)/(double)(ncircles);
            ccluster = (double)(particle.cluster_color)/(double)(ncircles);
            ccluster -= (double)((int)ccluster);
            hue = PARTICLE_HUE_MIN + (PARTICLE_HUE_MAX - PARTICLE_HUE_MIN)*ccluster;
            *radius = particle.radius;
            *width = BOUNDARY_WIDTH;
            break;
        }
        case (P_CLUSTER_SIZE):
        {
//             cluster = (double)(particle.cluster)/(double)(ncircles);
            ccluster = 1.0 - 5.0/((double)particle.cluster_size + 4.0);
            hue = PARTICLE_HUE_MIN + (PARTICLE_HUE_MAX - PARTICLE_HUE_MIN)*ccluster;
            *radius = particle.radius;
            *width = BOUNDARY_WIDTH;
            break;
        }
        case (P_CLUSTER_SELECTED):
        {
            cl = particle.cluster;
            if (cluster[cl].selected) hue = COLOR_HUE_CLUSTER_SELECTED;
            else hue = COLOR_HUE_CLUSTER_NOT_SELECTED;
            *radius = particle.radius;
            *width = BOUNDARY_WIDTH;
            break;
        }
        case (P_COLLISION):
        {
            hue = (double)particle.collision;
            if (hue > 0.0) hue = atan(0.25*(0.03*hue + 1.0))/PID;
//             {
//                 hue += 10.0;
//                 hue *= 1.0/(40.0 + hue);
//             }
            hue = PARTICLE_HUE_MIN + (PARTICLE_HUE_MAX - PARTICLE_HUE_MIN)*hue;
            *radius = particle.radius;
            *width = BOUNDARY_WIDTH;
            break;
        }
        case (P_RADIUS):
        {
//             hue = atan(5.0*(particle.radius/MU - 0.75))/PID;
            hue = (particle.radius/MU - RANDOM_RADIUS_MIN)/RANDOM_RADIUS_RANGE;
//             hue = 0.5*(hue + 1.0);
            hue = PARTICLE_HUE_MIN + (PARTICLE_HUE_MAX - PARTICLE_HUE_MIN)*hue;
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
            hsl_to_rgb_palette(hue, 0.9, 0.5, rgb, COLOR_PALETTE_DIFFNEIGH);
            hsl_to_rgb_palette(hue, 0.9, 0.5, rgbx, COLOR_PALETTE_DIFFNEIGH);
            hsl_to_rgb_palette(hue, 0.9, 0.5, rgby, COLOR_PALETTE_DIFFNEIGH);
//             hsl_to_rgb_twilight(hue, 0.9, 0.5, rgb);
//             hsl_to_rgb_twilight(hue, 0.9, 0.5, rgbx);
//             hsl_to_rgb_twilight(hue, 0.9, 0.5, rgby);
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
        case (P_LOG_EMEAN):  
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
        case (P_CHARGE):  
        {
            hsl_to_rgb_palette(hue, 0.9, 0.5, rgb, COLOR_PALETTE_CHARGE);
            hsl_to_rgb_palette(hue, 0.9, 0.5, rgbx, COLOR_PALETTE_CHARGE);
            hsl_to_rgb_palette(hue, 0.9, 0.5, rgby, COLOR_PALETTE_CHARGE);
            break;
        }
        case (P_MOL_ANGLE):  
        {
            hsl_to_rgb_palette(hue, 0.9, 0.5, rgb, COLOR_PALETTE_ANGLE);
            hsl_to_rgb_palette(hue, 0.9, 0.5, rgbx, COLOR_PALETTE_ANGLE);
            hsl_to_rgb_palette(hue, 0.9, 0.5, rgby, COLOR_PALETTE_ANGLE);
            break;
        }
        case (P_CLUSTER):  
        {
            hsl_to_rgb_palette(hue, 0.9, 0.5, rgb, COLOR_PALETTE_CLUSTER);
            hsl_to_rgb_palette(hue, 0.9, 0.5, rgbx, COLOR_PALETTE_CLUSTER);
            hsl_to_rgb_palette(hue, 0.9, 0.5, rgby, COLOR_PALETTE_CLUSTER);
            break;
        }
        case (P_CLUSTER_SIZE):  
        {
            hsl_to_rgb_palette(hue, 0.9, 0.5, rgb, COLOR_PALETTE_CLUSTER_SIZE);
            hsl_to_rgb_palette(hue, 0.9, 0.5, rgbx, COLOR_PALETTE_CLUSTER_SIZE);
            hsl_to_rgb_palette(hue, 0.9, 0.5, rgby, COLOR_PALETTE_CLUSTER_SIZE);
            break;
        }
        case (P_CLUSTER_SELECTED):  
        {
            hsl_to_rgb_palette(hue, 0.9, 0.5, rgb, COLOR_PALETTE_CLUSTER_SELECTED);
            hsl_to_rgb_palette(hue, 0.9, 0.5, rgbx, COLOR_PALETTE_CLUSTER_SELECTED);
            hsl_to_rgb_palette(hue, 0.9, 0.5, rgby, COLOR_PALETTE_CLUSTER_SELECTED);
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

void set_segment_pressure_color(double pressure, double lum, double rgb[3])
{
    double hue;
    
//     if (pressure < 0.0) pressure = 0.0;
    hue = PARTICLE_HUE_MIN + (PARTICLE_HUE_MAX - PARTICLE_HUE_MIN)*pressure/SEGMENT_PMAX;
    if (hue > PARTICLE_HUE_MIN) hue = PARTICLE_HUE_MIN;
    else if (hue < PARTICLE_HUE_MAX) hue = PARTICLE_HUE_MAX;
        
//     hsl_to_rgb_turbo(hue, 0.9, lum, rgb);
    
    hsl_to_rgb_palette(hue, 0.9, lum, rgb, COLOR_PALETTE_PRESSURE);
    
    glColor3f(rgb[0], rgb[1], rgb[2]);
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


void draw_triangles(t_particle particle[NMAXCIRCLES], t_cluster cluster[NMAXCIRCLES], int plot)
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
                compute_particle_colors(particle[i], cluster, plot, rgb, rgbx, rgby, &radius, &width, particle);
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
//         case (CHEM_DNA_ENZYME):
//         {
//             if (particle.type == 7)
//             {
//                 /* TODO */
//             }
//                 return(0);
//             break;
//         }
    }
    return(1);
}

void draw_one_particle(t_particle particle, double xc, double yc, double radius, double angle, int nsides, double width, double rgb[3])
/* draw one of the particles */ 
{
    double x1, x2, y1, y2, xc1, yc1, wangle, newradius = radius, pradius;
    int wsign, cont = 1, draw = 1;
    static double ca, sa;
    static int first = 1;
    
    if (first)
    {
        angle = APOLY*PID;
        ca = cos(angle);
        sa = sin(angle);
        first = 0;
    }
    
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
    else if (particle.interaction == I_SEGMENT)
        draw_colored_rotated_rect(xc1, yc1, particle.radius, 0.1*particle.radius, angle + APOLY*PID, rgb);
    else if (particle.interaction == I_SEGMENT_CHARGED)
    {
        draw_colored_rotated_rect(xc1, yc1, particle.radius, 0.1*particle.radius, angle + APOLY*PID, rgb);
        /* TODO: draw ends */
    }
    else if (cont) 
    {
        if (nsides == NSEG) draw_colored_circle_precomp(xc1, yc1, radius, rgb);
        else draw_colored_polygon(xc1, yc1, radius, nsides, angle + APOLY*PID, rgb);
    }
    
    /* draw crosses/bars on charged particles */
    if (TWO_TYPES)
    {
        if ((DRAW_CROSS)&&(particle.charge > 0.0))
        {
            pradius = particle.radius;
            if (ROTATION) 
            {
                angle = angle + APOLY*PID;
                ca = cos(angle);
                sa = sin(angle);
            }
            glLineWidth(2);
            glColor3f(1.0, 1.0, 1.0);
            x1 = xc1 - pradius*ca;
            y1 = yc1 - pradius*sa;
            x2 = xc1 + pradius*ca;
            y2 = yc1 + pradius*sa;
            draw_line(x1, y1, x2, y2);
            x1 = xc1 - pradius*sa;
            y1 = yc1 + pradius*ca;
            x2 = xc1 + pradius*sa;
            y2 = yc1 - pradius*ca;
            draw_line(x1, y1, x2, y2);
        }
        if ((DRAW_MINUS)&&(particle.charge < 0.0))
        {
            pradius = particle.radius;
            if (ROTATION) 
            {
                angle = angle + APOLY*PID;
                ca = cos(angle);
                sa = sin(angle);
            }
            glLineWidth(2);
            glColor3f(1.0, 1.0, 1.0);
            x1 = xc1 - pradius*ca;
            y1 = yc1 - pradius*sa;
            x2 = xc1 + pradius*ca;
            y2 = yc1 + pradius*sa;
            draw_line(x1, y1, x2, y2);
        }
    }
        
    glLineWidth(width);
    glColor3f(1.0, 1.0, 1.0);
    if (REACTION_DIFFUSION) cont = draw_special_particle(particle, xc1, yc1, radius, angle, nsides, rgb, 0);

    if ((particle.interaction == I_LJ_QUADRUPOLE)||(particle.interaction == I_LJ_DIPOLE)||(particle.interaction == I_VICSEK)||(particle.interaction == I_VICSEK_REPULSIVE)||(particle.interaction == I_VICSEK_SPEED)||(particle.interaction == I_VICSEK_SHARK)) 
        draw_rhombus(xc1, yc1, radius, angle + APOLY*PID);
    else if (particle.interaction == I_SEGMENT)
        draw_rotated_rect(xc1, yc1, particle.radius, 0.1*particle.radius, angle + APOLY*PID);
    else if (particle.interaction == I_SEGMENT_CHARGED)
    {
        draw_rotated_rect(xc1, yc1, particle.radius, 0.1*particle.radius, angle + APOLY*PID);
        /* TODO: draw ends */
    }
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
    
//     printf("Particle radius %.3f, pradius %.3f, mass %.3f\n", radius, particle.radius, 1.0/particle.mass_inv);
//     sleep(1);
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
    
        draw_colored_polygon(x1, y1, COLLISION_RADIUS*MU, NSEG, 0.0, rgb);
        collisions[i].time--;
    }
    
}

void draw_trajectory(t_tracer trajectory[TRAJECTORY_LENGTH*N_TRACER_PARTICLES], int traj_position, int traj_length, t_particle *particle, t_cluster *cluster, int *tracer_n, int plot)
/* draw tracer particle trajectory */
{
    int i, j, time, p, width;
    double x1, x2, y1, y2, rgb[3], rgbx[3], rgby[3], radius, lum;
    
//     blank();
    glLineWidth(TRAJECTORY_WIDTH);
    printf("drawing trajectory\n");
    
    for (j=0; j<n_tracers; j++)
    {
//         if (j == 0) hsl_to_rgb(HUE_TYPE1, 0.9, 0.5, rgb);
//         else if (j == 1) hsl_to_rgb(HUE_TYPE2, 0.9, 0.5, rgb);
//         else hsl_to_rgb(HUE_TYPE3, 0.9, 0.5, rgb);
//         set_type_color(j+2, 0.5, rgb);
        
        compute_particle_colors(particle[tracer_n[j]], cluster, plot, rgb, rgbx, rgby, &radius, &width, particle);
        
        glColor3f(rgb[0], rgb[1], rgb[2]);
            
        if (traj_length < TRAJECTORY_LENGTH) 
            for (i=0; i < traj_length-1; i++)
            {
                x1 = trajectory[j*TRAJECTORY_LENGTH + i].xc;
                x2 = trajectory[j*TRAJECTORY_LENGTH + i+1].xc;
                y1 = trajectory[j*TRAJECTORY_LENGTH + i].yc;
                y2 = trajectory[j*TRAJECTORY_LENGTH + i+1].yc;
            
                time = traj_length - i;
                lum = 0.8 - TRACER_LUM_FACTOR*(double)time/(double)TRAJECTORY_LENGTH;
                if (lum < 0.0) lum = 0.0;
                glColor3f(lum*rgb[0], lum*rgb[1], lum*rgb[2]);
        
                if (module2(x2 - x1, y2 - y1) < 0.25*(YMAX - YMIN)) draw_line(x1, y1, x2, y2);
            
//                 printf("tracer = %i, (x1, y1) = (%.3lg, %.3lg), (x2, y2) = (%.3lg, %.3lg)\n", j, x1, y1, x2, y2);
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
        
                if (module2(x2 - x1, y2 - y1) < 0.25*(YMAX - YMIN)) draw_line(x1, y1, x2, y2);
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
        
                if (module2(x2 - x1, y2 - y1) < 0.25*(YMAX - YMIN)) draw_line(x1, y1, x2, y2);
            }
        }
    }
}


void color_background(t_particle particle[NMAXCIRCLES], int bg_color, t_hashgrid hashgrid[HASHX*HASHY])
/* color background according to particle properties */
{
    int i, j, k, n, p, q, m, nnb, number, avrg_fact;
    double rgb[3], hue, value, p1, p2, pp1, pp2, oldhue;
    static int first = 1, counter = 0;
    
    p1 = 0.75;
    p2 = 1.0 - p1;
    pp1 = 0.95;
    pp2 = 1.0 - pp1;
   
    glBegin(GL_QUADS);
    for (i=0; i<HASHX; i++)
        for (j=0; j<HASHY; j++)
        {
            n = mhash(i, j);
            if (first) 
            {
                hashgrid[n].hue1 = 180.0;
                hashgrid[n].hue2 = 180.0;
            }
            /* set two old values for option DOUBLE_MOVIE */
            if (DOUBLE_MOVIE)
            {
                if (counter) oldhue = hashgrid[n].hue1;
                else oldhue = hashgrid[n].hue2;
            }
            else oldhue = hashgrid[n].hue1;
            
            switch (bg_color) {
                case (BG_DENSITY):
                {
                    nnb = hashgrid[n].nneighb;
                    number = 0;
                    for (q=0; q<nnb; q++)
                    {
                        m = hashgrid[n].neighbour[q];
                        number += hashgrid[m].number;
                    }
                    number += hashgrid[n].number;
                    
//                     hue = 50.0*(double)hashgrid[n].number;
                    hue = 75.0*(double)number/(double)(nnb + 1);
                    hue = p1*oldhue + p2*hue;
                    rgb[0] = hue/360.0;
                    rgb[1] = hue/360.0;
                    rgb[2] = hue/360.0;
                    break;
                }
                case (BG_CHARGE):
                {
                    avrg_fact = 3;      /* weight of central cell in hashgrid average */
                    value = 0.0;
                    nnb = hashgrid[n].nneighb;
                    for (q=0; q<nnb; q++)
                    {
                        m = hashgrid[n].neighbour[q];
                        for (k=0; k<hashgrid[m].number; k++)
                        {
                            p = hashgrid[m].particles[k];
                            value += particle[p].charge;
                        }
                    }
                    /* hashcell n counts double */
                    for (k=0; k<hashgrid[n].number; k++)
                    {
                        p = hashgrid[n].particles[k];
                        value += (double)(avrg_fact-1)*particle[p].charge;
                    }
                    value *= 1.0/(double)(nnb + avrg_fact);
                    if (CHARGE_OBSTACLES) value += hashgrid[n].charge;
                    hue = (-tanh(0.5*value)+1.0)*180.0;
                    hue = p1*oldhue + p2*hue;
                    hsl_to_rgb_twilight(hue, 0.9, 0.5, rgb);
                    break;
                }
                case (BG_EKIN):
                {
                    value = 0.0;
                    nnb = hashgrid[n].nneighb;
                    for (q=0; q<nnb; q++)
                    {
                        m = hashgrid[n].neighbour[q];
                        for (k=0; k<hashgrid[m].number; k++)
                        {
                            p = hashgrid[m].particles[k];
                            value += particle[p].energy;
                        }
                    }
                    /* hashcell n counts double */
                    for (k=0; k<hashgrid[n].number; k++)
                    {
                        p = hashgrid[n].particles[k];
                        value += particle[p].energy;
                    }
                    value *= 1.0/(double)(nnb + 1);
                    hue = ENERGY_HUE_MIN + (ENERGY_HUE_MAX - ENERGY_HUE_MIN)*value/PARTICLE_EMAX;
                    if (hue > ENERGY_HUE_MIN) hue = ENERGY_HUE_MIN;
                    if (hue < ENERGY_HUE_MAX) hue = ENERGY_HUE_MAX;
                    hue = p1*oldhue + p2*hue;
                    hsl_to_rgb_palette(hue, 0.9, 0.5, rgb, COLOR_PALETTE_EKIN);
                    break;
                }
                case (BG_FORCE):
                {
                    value = 0.0;
                    nnb = hashgrid[n].nneighb;
                    for (q=0; q<nnb; q++)
                    {
                        m = hashgrid[n].neighbour[q];
                        for (k=0; k<hashgrid[m].number; k++)
                        {
                            p = hashgrid[m].particles[k];
                            value += module2(particle[p].fx, particle[p].fy);
                        }
                    }
                    /* hashcell n counts double */
                    for (k=0; k<hashgrid[n].number; k++)
                    {
                        p = hashgrid[n].particles[k];
                        value += module2(particle[p].fx, particle[p].fy);
                    }
                    value *= BG_FORCE_SLOPE/(double)(nnb + 1);
                    hue = (1.0 - tanh(value))*360.0;
                    hue = pp1*oldhue + pp2*hue;
                    hsl_to_rgb_turbo(hue, 0.9, 0.5, rgb);
                    break;
                }
            }
            if (DOUBLE_MOVIE)
            {
                if (counter) hashgrid[n].hue1 = hue;
                else hashgrid[n].hue2 = hue;
            }
            else hashgrid[n].hue1 = hue;
            
            glColor3f(rgb[0], rgb[1], rgb[2]);
            glVertex2d(hashgrid[n].x1, hashgrid[n].y1);
            glVertex2d(hashgrid[n].x2, hashgrid[n].y1);
            glVertex2d(hashgrid[n].x2, hashgrid[n].y2);
            glVertex2d(hashgrid[n].x1, hashgrid[n].y2);
        }
    glEnd();
    first = 0;
    counter = 1 - counter;
}
    

void draw_particles(t_particle particle[NMAXCIRCLES], t_cluster cluster[NMAXCIRCLES], int plot, double beta, t_collision *collisions, int ncollisions, int bg_color, t_hashgrid hashgrid[HASHX*HASHY], t_lj_parameters params)
{
    int i, j, k, m, width, nnbg, nsides, cl, p, q;
    double ej, hue, huex, huey, rgb[3], rgbx[3], rgby[3], radius, x1, y1, x2, y2, angle, ca, sa, length, linkcolor, sign = 1.0, angle1, signy = 1.0, periodx, periody, x, y, lum, logratio;
    char message[100];
    
//     if (!TRACER_PARTICLE) blank();
    if (plot == P_NOPARTICLE) blank();
    
    if ((COLOR_BACKGROUND)&&(bg_color > 0)) color_background(particle, bg_color, hashgrid);
    else if (!TRACER_PARTICLE) blank();
    
    glColor3f(1.0, 1.0, 1.0);
    
    /* show region of partial thermostat */
    if (PARTIAL_THERMO_COUPLING)
    {
        switch (PARTIAL_THERMO_REGION) {
            case (TH_INBOX):
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
                break;
            }
            case(TH_LAYER): 
            {
                hue = 0.75*PARTICLE_HUE_MIN + 0.25*PARTICLE_HUE_MAX;
                erase_area_hsl_turbo(0.0, YMIN, XMAX, PARTIAL_THERMO_HEIGHT - YMIN,  hue, 0.9, 0.15);
                break;
            }
            case (TH_RING):
            {
                hue = 0.75*PARTICLE_HUE_MIN + 0.25*PARTICLE_HUE_MAX;
                hsl_to_rgb_turbo(hue, 0.9, 0.15, rgb);
                draw_colored_circle(0.0, 0.0, PARTIAL_THERMO_WIDTH, 180, rgb);
                break;
            }
            case (TH_RING_EXPAND):
            {
                hue = 0.75*PARTICLE_HUE_MIN + 0.25*PARTICLE_HUE_MAX;
                hsl_to_rgb_turbo(hue, 0.9, 0.15, rgb);
                draw_colored_circle(0.0, 0.0, params.thermo_radius, 180, rgb);
                break;
            }
            case(TH_INIT): 
            {
                hue = 0.75*PARTICLE_HUE_MIN + 0.25*PARTICLE_HUE_MAX;
                erase_rectangle_hsl_turbo(INITXMIN, INITYMIN, INITXMAX, INITYMAX,  hue, 0.9, 0.15);
                break;
            }
            case(TH_THERMO): 
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
                erase_rectangle_hsl_turbo(THERMOXMIN, THERMOYMIN, THERMOXMAX, THERMOYMAX,  hue, 0.9, 0.15);
                break;
            }
            case (TH_CONE):
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
                
                draw_colored_triangle_turbo(0.05, -0.2, LAMBDA, LAMBDA, -LAMBDA, LAMBDA, hue, 0.9, 0.15);
                draw_colored_triangle_turbo(0.05, -0.2, -LAMBDA, LAMBDA, -0.04, -0.2, hue, 0.9, 0.15);
                break;
            }
        }
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
    if (FILL_TRIANGLES) draw_triangles(particle, cluster, plot);
    
    /* draw collision discs */
    if ((REACTION_DIFFUSION)&&(ncollisions > 0)) draw_collisions(collisions, ncollisions);
    
    /* determine particle color and size */
    for (j=0; j<ncircles; j++) if (particle[j].active)
    {
        compute_particle_colors(particle[j], cluster, plot, rgb, rgbx, rgby, &radius, &width, particle);
       
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
            case (I_COULOMB_PENTA): 
            {
                nsides = 5;
                break;
            }
            case (I_DNA_CHARGED):
            {
                if (particle[j].type == 7) nsides = 3;
                else nsides = NSEG;
                break;
            }
            case (I_SEGMENT): 
            {
                nsides = 2;
                break;
            }
            case (I_SEGMENT_CHARGED): 
            {
                nsides = 2;
                break;
            }
            case (I_POLYGON): 
            {
                nsides = NPOLY;
                break;
            }
            case (I_POLYGON_ALIGN): 
            {
                nsides = NPOLY;
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
    
    /* draw links between particles in cluster */
    if (DRAW_CLUSTER_LINKS)
    {
        glLineWidth(LINK_WIDTH);
        for (j=0; j<ncircles; j++) if (particle[j].active) 
        {
            cl = particle[j].cluster;
            x1 = particle[j].xc;
            y1 = particle[j].yc;
            for (p=0; p<cluster[cl].nparticles; p++)
            {
                q = cluster[cl].particle[p];
                x2 = particle[q].xc;
                y2 = particle[q].yc;
                draw_line(x1, y1, x2, y2);
            }
        }
    }
}


void draw_container(double xmin, double xmax, t_obstacle obstacle[NMAXOBSTACLES], t_segment segment[NMAXSEGMENTS], t_belt *conveyor_belt, int wall)
/* draw the container, for certain boundary conditions */
{
    int i, j;
    double rgb[3], x, y, r, phi, angle, dx, dy, ybin, x1, x2, y1, y2, h, bangle, blength, pos0, pos, bsegment, tx, ty, bwidth, ca, sa;
    char message[100];
    
    /* draw fixed obstacles */
    if (ADD_FIXED_OBSTACLES)
    {
        glLineWidth(CONTAINER_WIDTH);
        if (CHARGE_OBSTACLES) hsl_to_rgb(30.0, 0.1, 0.5, rgb);
        else hsl_to_rgb(300.0, 0.1, 0.5, rgb);
        for (i = 0; i < nobstacles; i++)
            draw_colored_circle_precomp(obstacle[i].xc - xtrack, obstacle[i].yc - ytrack, obstacle[i].radius, rgb);
        glColor3f(1.0, 1.0, 1.0);
        for (i = 0; i < nobstacles; i++)
        {
            x = obstacle[i].xc - xtrack;
            y = obstacle[i].yc - ytrack;
            r = obstacle[i].radius;
            draw_circle_precomp(x, y, r);
            if (CHARGE_OBSTACLES) 
            {
                draw_line(x - r, y, x + r, y);
                draw_line(x, y - r, x, y + r);
            }
            if (ROTATE_OBSTACLES)
            {
                ca = cos(obstacle[i].angle);
                sa = sin(obstacle[i].angle);
//                 printf("Obstacle %i angle = %.5lg\n", i, obstacle[i].angle);
                draw_line(x - r*ca, y - r*sa, x + r*ca, y + r*sa);
                draw_line(x + r*sa, y - r*ca, x - r*sa, y + r*ca);
            }
        }
    }
    if (ADD_FIXED_SEGMENTS)
    {
        glLineWidth(CONTAINER_WIDTH);
        glColor3f(1.0, 1.0, 1.0);
        for (i = 0; i < nsegments; i++) if (segment[i].active)
        {
            if (COLOR_SEG_GROUPS) set_segment_group_color(segment[i].group, 1.0, rgb);
            else if (SHOW_SEGMENTS_PRESSURE) set_segment_pressure_color(segment[i].avrg_pressure, 1.0, rgb);
            draw_line(segment[i].x1 - xtrack, segment[i].y1 - ytrack, segment[i].x2 - xtrack, segment[i].y2 - ytrack);
        }
        
        /* draw conveyor belt */
        if (ADD_CONVEYOR_FORCE) for (i=0; i<nbelts;i++)
        {
            bwidth = conveyor_belt[i].width;
            r = 0.7*bwidth;
            tx = conveyor_belt[i].tx;
            ty = conveyor_belt[i].ty;
            bangle = -conveyor_belt[i].position/(bwidth);
            blength = 2.0*conveyor_belt[i].length + DPI*conveyor_belt[i].width;
            dx = r*cos(bangle);
            dy = r*sin(bangle);
            
            x1 = conveyor_belt[i].x1;
            y1 = conveyor_belt[i].y1;
            draw_line(x1-dx, y1-dy, x1+dx, y1+dy);
            draw_line(x1+dy, y1-dx, x1-dy, y1+dx);
            draw_circle_precomp(x1, y1, r);

            x2 = conveyor_belt[i].x2;
            y2 = conveyor_belt[i].y2;
            draw_line(x2-dx, y2-dy, x2+dx, y2+dy);
            draw_line(x2+dy, y2-dx, x2-dy, y2+dx);
            draw_circle_precomp(x2, y2, r);
            
            draw_line(x1-r*ty, y1+r*tx, x2-r*ty, y2+r*tx);
            draw_line(x1+r*ty, y1-r*tx, x2+r*ty, y2-r*tx);
            
            bsegment = 10.0*blength;
            bsegment = blength/(double)((int)bsegment);
            
            if (conveyor_belt[i].position > 0.0)
                pos0 = conveyor_belt[i].position - bsegment*(double)((int)(conveyor_belt[i].position/bsegment));
            else 
            {
                pos = -conveyor_belt[i].position;
                pos0 = pos - bsegment*(double)((int)(pos/bsegment));
                pos0 = bsegment - pos0;
            }
            
            while (pos0 < conveyor_belt[i].length)
            {
                x = x1 + tx*pos0;
                y = y1 + ty*pos0;
                draw_line(x - r*ty, y + r*tx, x - bwidth*ty, y + bwidth*tx);
                pos0 += bsegment;
            }
            while (pos0 < conveyor_belt[i].length + PI*bwidth)
            {
                pos = (pos0 - conveyor_belt[i].length)/bwidth;
                angle = PID + conveyor_belt[i].angle - pos;
                ca = cos(angle);
                sa = sin(angle);
                draw_line(x2 + r*ca, y2 + r*sa, x2 + bwidth*ca, y2 + bwidth*sa);
                pos0 += bsegment;
            }
            while (pos0 < 2.0*conveyor_belt[i].length + PI*bwidth)
            {
                pos = pos0 - conveyor_belt[i].length - PI*bwidth;
                x = x2 - tx*pos;
                y = y2 - ty*pos;
                draw_line(x + r*ty, y - r*tx, x + bwidth*ty, y - bwidth*tx);
                pos0 += bsegment;
            }
            while (pos0 < 2.0*conveyor_belt[i].length + DPI*bwidth)
            {
                pos = (pos0 - 2.0*conveyor_belt[i].length - PI*bwidth)/bwidth;
                angle = -PID + conveyor_belt[i].angle - pos;
                ca = cos(angle);
                sa = sin(angle);
                draw_line(x1 + r*ca, y1 + r*sa, x1 + bwidth*ca, y1 + bwidth*sa);
                pos0 += bsegment;
            }
            
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
        
        y1 -= 0.1;
        erase_area_hsl(xrightbox, y1 + 0.025, 0.42, 0.05, 0.0, 0.9, 0.0);
        glColor3f(1.0, 1.0, 1.0);
        sprintf(message, "Fy/Fx = %.4f", fy/fx);
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
//         xmid = 0.5*(XMIN + XMAX) + 0.05*scale;
        xmidtext = xmid - 0.2*scale;
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
        write_text(xtext + 0.1, y, message);
        
        erase_area_hsl(xmid, y + 0.025*scale, 0.37*scale, 0.05*scale, 0.0, 0.9, 0.0);
        glColor3f(1.0, 1.0, 1.0);
        sprintf(message, "Temperature %.2f", params.mean_energy);
        write_text(xmidtext - 0.1, y, message);

        erase_area_hsl(xxbox, y + 0.025*scale, 0.37*scale, 0.05*scale, 0.0, 0.9, 0.0);
        glColor3f(1.0, 1.0, 1.0);
        sprintf(message, "Pressure %.3f", meanpressure);
        write_text(xxtext, y, message);

    }   
    else if (INCREASE_KREPEL)  /* print force constant */
    {
        erase_area_hsl(xbox, y + 0.025*scale, 0.35*scale, 0.05*scale, 0.0, 0.9, 0.0);
        glColor3f(1.0, 1.0, 1.0);
        sprintf(message, "Force constant %.2f", params.krepel);
        write_text(xtext + 0.03, y, message);
    }   
    else if (INCREASE_E)  /* print electric field */
    {
        erase_area_hsl(xbox, y + 0.025*scale, 0.27*scale, 0.05*scale, 0.0, 0.9, 0.0);
        glColor3f(1.0, 1.0, 1.0);
        sprintf(message, "E field = %.2f", 25.0*NVID*DT_PARTICLE*params.efield);
        write_text(xtext + 0.08, y, message);
    }   
    else if (INCREASE_B)  /* print magnetic field */
    {
        erase_area_hsl(xbox, y + 0.025*scale, 0.27*scale, 0.05*scale, 0.0, 0.9, 0.0);
        glColor3f(1.0, 1.0, 1.0);
        sprintf(message, "B field = %.2f", 25.0*NVID*DT_PARTICLE*params.bfield);
        write_text(xtext + 0.08, y, message);
    }   
    
    if (PRINT_NPARTICLES)   /* print number of particles */
    {
        erase_area_hsl(xbox, y + 0.025*scale, 0.27*scale, 0.05*scale, 0.0, 0.9, 0.0);
        glColor3f(1.0, 1.0, 1.0);
        sprintf(message, "%i electrons", params.nactive);
        write_text(xtext + 0.08, y, message);
    } 
    if (PRINT_TYPE_PROP)   /* print proportion of types */
    {
        erase_area_hsl(xbox, y + 0.025*scale, 0.27*scale, 0.05*scale, 0.0, 0.9, 0.0);
        glColor3f(1.0, 1.0, 1.0);
        sprintf(message, "%.1f %% anions", 100.0*params.prop);
        write_text(xtext + 0.08, y, message);
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
//         write_text(xmidtext + 0.1, y, message);
        write_text(xmidtext, y, message);
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
    
    if (COUNT_PARTNER_TYPE)
    {
        for (i=0; i<ncircles; i++) if (particle[i].active)
        {
            if (particle[i].type == 0)
                particle_numbers[n + particle[i].npartners + 1]++;
            else if ((particle[i].type == 2)&&(particle[i].npartners == 0))
                particle_numbers[n + 1]++;
        }
//         for (type = 0; type < RD_TYPES+1; type++) printf("type = %i, %i particles\n", type, particle_numbers[n + type]);
    }
    else for (i=0; i<ncircles; i++)
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

double compute_boundary_force(int j, t_particle particle[NMAXCIRCLES], t_obstacle obstacle[NMAXOBSTACLES], t_segment segment[NMAXSEGMENTS], double xleft, double xright, double *pleft, double *pright, double pressure[N_PRESSURES], int wall, double krepel)
{
    int i, k, corner;
    double xmin, xmax, ymin, ymax, padding, r, rp, r2, cphi, sphi, f, fperp = 0.0, x, y, xtube, distance, dx, dy, width, ybin, angle, x1, y1, x2, h, ytop, norm, dleft, dplus, dminus, tmp_pleft = 0.0, tmp_pright = 0.0, proj, pscal, pvect, pvmin, charge, d2, speed, pangle, distance1, sangle;
    static int first = 1;
    static double seqangle;
    
    if (first)
    {
        seqangle = PI/(double)NPOLY;
        first = 0;
    }
    
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
            
            if (ROTATE_OBSTACLES)
            {
                speed = particle[j].vx*sphi - particle[j].vy*cphi;
                f = KSPRING_BELT*(obstacle[i].omega*obstacle[i].radius/particle[j].mass_inv - speed);
                particle[j].fx += f*sphi;
                particle[j].fy -= f*cphi;
            }
        }
        
        if (CHARGE_OBSTACLES)
        {
            charge = obstacle[i].charge*particle[k].charge;
            d2 = distance*distance;
            f = KCOULOMB_OBSTACLE*charge/(1.0e-12 + d2);
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
            /* case of interacting polygons */
//             if (POLYGON_INTERACTION) 
//             {
// //                 distance = segment[i].nx*x + segment[i].ny*y - segment[i].c;
//                 pangle = DPI/(double)NPOLY;
//                 dx = 0.0;
//                 dy = 0.0;
//                 for (corner=0; corner<NPOLY; corner++)
//                 {
//                     angle = particle[j].angle + pangle*(double)corner;
//                     x1 = x + particle[j].radius*cos(angle);
//                     y1 = y + particle[j].radius*sin(angle);
//                     distance1 = segment[i].nx*x1 + segment[i].ny*y1 - segment[i].c;
//                     if (distance1 < distance)
//                     {
//                         x = x1;
//                         y = y1;
//                         distance = distance1;
//                         dx = x1 - particle[j].xc;
//                         dy = y1 - particle[j].yc;
//                     }
//                 }
// //                 printf("(dx, dy) = (%.3lg, %.3lg)\n", dx, dy);
//             }
            r = SEGMENT_FORCE_EQR*particle[j].radius;
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
                }
                
                if (SHOW_SEGMENTS_PRESSURE)
                {
                    segment[i].pressure = f/segment[i].length;
                    segment[i].avrg_pressure *= (1.0 - P_AVRG_FACTOR);
                    segment[i].avrg_pressure += P_AVRG_FACTOR*segment[i].pressure;
//                      printf("Segment %i: pressure = %.3lg\n", i, segment[i].avrg_pressure);
                }
                
                if (ADD_CONVEYOR_FORCE)
                {
                    speed = particle[j].vx*segment[i].ny - particle[j].vy*segment[i].nx;
                    f = KSPRING_BELT*(segment[i].conveyor_speed/particle[j].mass_inv - speed);
                    particle[j].fx += f*segment[i].ny;
                    particle[j].fy -= f*segment[i].nx;
//                     printf("Belt force (%.5lg, %.5lg)\n", f*segment[i].ny, -f*segment[i].nx); 
                }
                
                if ((POLYGON_INTERACTION)&&(segment[i].align_torque))
                {
                    sangle = segment[i].nangle;
                    particle[j].torque -=        KTORQUE_BOUNDARY*sin(particle[j].spin_freq*(particle[j].angle - sangle) - seqangle);
//                     particle[j].torque -= KTORQUE_BOUNDARY*(dx*f*segment[i].ny - dy*f*segment[i].nx);
//                     printf("Torque on particle %i = %.3lg\n", j, particle[j].torque);
                }
                
            }
            if ((VICSEK_INT)&&(vabs(distance) < SEGMENT_FORCE_EQR*r))
            {
                pvmin = 2.0;
                pvect = cos(particle[j].angle)*segment[i].ny - sin(particle[j].angle)*segment[i].nx;
                if ((pvect > 0.0)&&(pvect < pvmin)) pvect = pvmin;
                else if ((pvect < 0.0)&&(pvect > -pvmin)) pvect = -pvmin;
//                 particle[j].torque += KTORQUE_BOUNDARY*pvect;
                particle[j].torque += KTORQUE_BOUNDARY*pvect*(SEGMENT_FORCE_EQR*r - vabs(distance));
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
            
            r = SEGMENT_FORCE_EQR*particle[j].radius;
            
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
            padding = particle[j].radius + 0.01;
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

void dissociate_particles_findp(int i, int j, t_particle particle[NMAXCIRCLES])
/* dissociate partnered particles i and j */
{
    int np, nq, r, p, q; 
    
    np = particle[i].npartners;
    nq = particle[j].npartners;
    
    /* find which partner of j is particle i */
    p = 0;
    while ((particle[i].partner[p] != j)&&(p < np)) p++;
   
    /* remove partner p from partner list */
    for (q=p; q<np-1; q++)
        particle[i].partner[q] = particle[i].partner[q+1];
    particle[i].npartners--;
    
    /* find which partner of j is particle i */
    q = 0;
    while ((particle[j].partner[q] != i)&&(q < nq)) q++;
    
    /* remove partner q from partner list */
    for (r=q; r<nq-1; r++)
        particle[j].partner[r] = particle[j].partner[r+1];
    particle[j].npartners--;
}


void dissociate_particles(int i, int j, int p, t_particle particle[NMAXCIRCLES])
/* dissociate partnered particles i and j */
/* assuming j is pth partner of i */
{
    int np, nq, r, q; 
    
    np = particle[i].npartners;
    nq = particle[j].npartners;
   
    /* remove partner p from partner list */
    for (q=p; q<np-1; q++)
        particle[i].partner[q] = particle[i].partner[q+1];
    particle[i].npartners--;
    
    /* find which partner of j is particle i */
    q = 0;
    while ((particle[j].partner[q] != i)&&(q < nq)) q++;
    
    /* remove partner q from partner list */
    for (r=q; r<nq-1; r++)
        particle[j].partner[r] = particle[j].partner[r+1];
    particle[j].npartners--;
}

void dissociate_molecule(int i, int j, t_particle particle[NMAXCIRCLES])
/* dissociate all particles in molecule containing particles i and j */
{
    int np, nq, p, q, n;
    
    np = particle[i].npartners;
    nq = particle[j].npartners;
    
    particle[i].npartners = 0;
    particle[j].npartners = 0;
    
    particle[i].vx = 0.0;
    particle[i].vy = 0.0;
    particle[i].thermostat = 1;
    particle[j].vx = 0.0;
    particle[j].vy = 0.0;
    particle[j].thermostat = 1;
    
    for (p=0; p<np; p++) 
    {
        n = particle[i].partner[p];
        particle[n].npartners = 0;
        particle[n].vx = 0.0;
        particle[n].vy = 0.0;
        particle[n].thermostat = 1;
    }
    for (q=0; q<nq; q++) 
    {
        n = particle[j].partner[q];
        particle[n].npartners = 0;
        particle[n].vx = 0.0;
        particle[n].vy = 0.0;
        particle[n].thermostat = 1;
    }
}


void compute_partner_force(int j, int n, double eq_distance, double eq_angle, double f[2], double *torque, t_particle particle[NMAXCIRCLES])
/* compute force of partner particle n on particle j */
{
    double dx, dy, r, ca, sa, force, sangle, torque2;
    double x1, x2, y1, y2, rmax, alpha, cosa, rj, rn, phi;
    int p, q;
    
    dx = particle[n].xc - particle[j].xc;
    dy = particle[n].yc - particle[j].yc;
    
    if (dx > 0.5*(BCXMAX - BCXMIN)) dx -= (BCXMAX-BCXMIN);
    else if (dx < -0.5*(BCXMAX - BCXMIN)) dx += (BCXMAX-BCXMIN);
    if (dy > 0.5*(BCYMAX - BCYMIN)) dy -= (BCYMAX-BCYMIN);
    else if (dy < -0.5*(BCYMAX - BCYMIN)) dy += (BCYMAX-BCYMIN);
    
    
    
    if ((PAIR_PARTICLES)&&(PAIRING_TYPE == POLY_SEG_POLYGON))
    {
        rmax = 1.0*particle[j].radius;
        for (p=-1; p<=1; p+=2)
            for (q=-1; q<=1; q+=2)
            {
                x1 = particle[n].xc + (double)p*particle[n].radius*cos(particle[n].angle);
                y1 = particle[n].yc + (double)p*particle[n].radius*sin(particle[n].angle);
                
                x2 = particle[j].xc + (double)q*particle[n].radius*cos(particle[j].angle);
                y2 = particle[j].yc + (double)q*particle[n].radius*sin(particle[j].angle);
                
                r = module2(x2-x1, y2-y1);
                if ((r>0.0)&&(r < rmax))
                {
//                     f[0] += KSPRING_PAIRS*(x1-x2)/r;
//                     f[1] += KSPRING_PAIRS*(y1-y2)/r;
                    f[0] = 1.0e8*(x1-x2)/r;
                    f[1] = 1.0e8*(y1-y2)/r;
                }
            }
    }
//     else
    {
        r = module2(dx, dy);
        if (r < 1.0e-10) r = 1.0e-10;
        if (r > 2.0*eq_distance) r = 2.0*eq_distance;
//     if (r > 1.5*eq_distance) r = 1.5*eq_distance;
        ca = dx/r;
        sa = dy/r;
    
        /* TODO: adjust max distance */
        if (r < 1.5*eq_distance) force = KSPRING_PAIRS*(r - eq_distance);
        else 
        {
//         printf("Dissociating partners %i and %i because max distance exceeded\n", j, n);
            force = 0.0;
//         dissociate_particles_findp(j, n, particle);
//         dissociate_molecule(j, n, particle);
        }
    
        f[0] = force*ca;
        f[1] = force*sa;
        
        /* TEST */
        f[0] -= particle[j].vx*DAMPING;
        f[1] -= particle[j].vy*DAMPING;
    }
    
    if (ROTATION)
    {
        *torque = 0.0;
        if ((REACTION_DIFFUSION)&&(RD_REACTION == CHEM_POLYGON_AGGREGATION))
        {
//             r = module2(dx, dy);
//             if (r < 2.0*MU)
//             {
//             if (NPOLY%2 == 0) 
//                 sangle = sin((double)NPOLY*(particle[n].angle - particle[j].angle));
//             else sangle = sin((double)NPOLY*(particle[n].angle - particle[j].angle) - PI);
//                 sangle = sin((double)NPOLY*(particle[n].angle - particle[j].angle - eq_angle));
//                 *torque += KTORQUE_PAIRS*sangle;
            
                /* stabilize relative angle */
                q = particle[j].partner[0];
                if ((particle[j].npartners == 1)&&(particle[q].npartners == 1))
                {
                    sangle = sin((double)NPOLY*(particle[n].angle - particle[j].angle - eq_angle));
                    *torque += KTORQUE_PAIRS*sangle;
                    phi = argument(dx, dy);
                    sangle = sin((double)NPOLY*(phi - particle[j].angle) - PI);
                    *torque += KTORQUE_PAIRS*sangle;
                }
                
                /* TEST: add some damping to rotation */
                *torque -= particle[j].omega*DAMPING_ROT;
            
//             sangle = sin((double)NPOLY*(phi - particle[j].angle) - PI);
//             *torque += KTORQUE_PAIRS*sangle;
//             }
        }
        else if ((PAIR_PARTICLES)&&(PAIRING_TYPE == POLY_SEG_POLYGON))
        {
            /* TODO */
            sangle = sin(particle[n].angle - particle[j].angle - eq_angle);
//             if (sangle > 0.0) *torque = KTORQUE_PAIRS*(1.0 + sangle);
//             else *torque = KTORQUE_PAIRS*(-1.0 + sangle);
            *torque = KTORQUE_PAIRS*sangle;
            
            
//             *torque = 0.0;
            
//             if (COMPUTE_PAIR_TORQUE)
//             {
//                 torque2 = KTORQUE_PAIR_ANGLE*sangle;
//                 f[0] -= torque2*sa;
//                 f[1] += torque2*ca;
//             }
        }
        else
        {
            sangle = sin(SPIN_INTER_FREQUENCY*(particle[n].angle - particle[j].angle));
            if (sangle > 0.0) *torque = KTORQUE_PAIRS*(1.0 + sangle);
            else *torque = KTORQUE_PAIRS*(-1.0 + sangle);
        
            if (COMPUTE_PAIR_TORQUE)
            {
                torque2 = KTORQUE_PAIR_ANGLE*sangle;
                f[0] -= torque2*sa;
                f[1] += torque2*ca;
            }
        }
    }
}


void compute_particle_force(int j, double krepel, t_particle particle[NMAXCIRCLES], t_hashgrid hashgrid[HASHX*HASHY])
/* compute force from other particles on particle j */
{
    int i0, j0, m0, k, l, m, q, close, n, test, different_cluster = 1;
    double fx = 0.0, fy = 0.0, force[2], torque = 0.0, torque_ij, x, y;
    
    particle[j].neighb = 0;
    if (REACTION_DIFFUSION) particle[j].diff_neighb = 0;

    for (k=0; k<particle[j].hash_nneighb; k++) 
    {
        n = particle[j].hashneighbour[k];
        if (CLUSTER_PARTICLES)  /* test whether n is not in the same cluster */
                different_cluster = (particle[n].cluster != particle[j].cluster);
        else different_cluster = 1;
        
        if ((particle[n].active)&&(different_cluster))
        {
            if (PAIR_FORCE)
            {
                /* test whether n is not a partner particle */
                test = 1;
                for (l=0; l<particle[n].npartners; l++)
                    if (particle[n].partner[l] == j) test = 0;
                
                if (test)
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
            }
            else
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
        }
    }
    
    if (PAIR_FORCE) for (l=0; l<particle[j].npartners; l++)
    {
        n = particle[j].partner[l];
        if ((!CLUSTER_PARTICLES)||(particle[n].cluster != particle[j].cluster))
        {
            compute_partner_force(j, n, particle[j].partner_eqd[l], particle[j].partner_eqa[l], force, &torque_ij, particle);
            fx += force[0];
            fy += force[1];
            torque += torque_ij;
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


int twotype_config(int i, t_particle particle[NMAXCIRCLES])
/* assign different particle types */
{
    switch (TWOTYPE_CONFIG) {
        case (TTC_RANDOM): return((double)rand()/RAND_MAX > TYPE_PROPORTION);
        case (TTC_CHESSBOARD): 
        {
            switch (CIRCLE_PATTERN) {
                case (C_SQUARE): return((i%NGRIDY + i/NGRIDY)%2 == 0);
                case (C_HEX): return(i%2 == 0);
                default: return(i%2 == 0);
            }
        }
        case (TTC_COANDA): 
        {
            if (vabs(particle[i].yc) < LAMBDA) return(1);
            if (vabs(particle[i].yc) < LAMBDA + MU) return(-1);
            return(0);
        }
        default: return(0);
    }
}


int initialize_configuration(t_particle particle[NMAXCIRCLES], t_hashgrid hashgrid[HASHX*HASHY], t_obstacle obstacle[NMAXOBSTACLES], double px[NMAXCIRCLES], double py[NMAXCIRCLES], double pangle[NMAXCIRCLES], int tracer_n[N_TRACER_PARTICLES],
t_segment segment[NMAXSEGMENTS], t_molecule molecule[NMAXCIRCLES])
/* initialize all particles, obstacles, and the hashgrid */
{
    int i, j, k, n, type, nactive = 0, hashcell;
    double x, y, h, xx, yy, rnd, angle;
    
    for (i=0; i < ncircles; i++) 
    {
        /* set particle type */
        particle[i].type = 0;
//         if ((TWO_TYPES)&&((double)rand()/RAND_MAX > TYPE_PROPORTION)) 
//         if ((TWO_TYPES)&&(twotype_config(i, particle))) 
        if (TWO_TYPES)
        {
            type = twotype_config(i, particle);
            if (type == 1)
            {
                particle[i].type = 1;
//                 particle[i].type = 2;
                particle[i].radius = MU_B;
            }
            else if (type == -1) particle[i].active = 0;
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
        particle[i].damping = 1.0;
        particle[i].dirmean = 0.0;
        particle[i].paired = 0;
        particle[i].collision = 0;
//         particle[i].added = 0;
//         particle[i].coulomb = 1;

//         particle[i].energy = 0.0;
//         y = particle[i].yc;
//         if (y >= YMAX) y -= particle[i].radius;
//         if (y <= YMIN) y += particle[i].radius;

        if (RANDOM_RADIUS) particle[i].radius = particle[i].radius*(RANDOM_RADIUS_MIN + RANDOM_RADIUS_RANGE*((double)rand()/RAND_MAX));
        
        if (particle[i].type == 0)
        {
            particle[i].interaction = INTERACTION;
            particle[i].eq_dist = EQUILIBRIUM_DIST;
            particle[i].spin_range = SPIN_RANGE;
            particle[i].spin_freq = SPIN_INTER_FREQUENCY;
            particle[i].mass_inv = 1.0/PARTICLE_MASS;
            if ((RANDOM_RADIUS)&&(ADAPT_MASS_TO_RADIUS > 0.0)) 
                particle[i].mass_inv *= 1.0/(1.0 + pow(particle[i].radius/MU, ADAPT_MASS_TO_RADIUS));
            if ((RANDOM_RADIUS)&&(ADAPT_DAMPING_TO_RADIUS > 0.0))
                particle[i].damping = 1.0 + ADAPT_DAMPING_FACTOR*pow(particle[i].radius/MU, ADAPT_DAMPING_TO_RADIUS);
            particle[i].inertia_moment_inv = 1.0/PARTICLE_INERTIA_MOMENT;
            particle[i].charge = CHARGE;
        }
        else 
        {
            particle[i].interaction = INTERACTION_B;
            particle[i].eq_dist = EQUILIBRIUM_DIST_B;
            particle[i].spin_range = SPIN_RANGE_B;
            particle[i].spin_freq = SPIN_INTER_FREQUENCY_B;
            particle[i].mass_inv = 1.0/PARTICLE_MASS_B;
            particle[i].inertia_moment_inv = 1.0/PARTICLE_INERTIA_MOMENT_B;
            particle[i].charge = CHARGE_B;
        }
               
        switch (V_INITIAL_TYPE)
        {
            case (VI_RANDOM):
            {
                particle[i].vx = V_INITIAL*gaussian();
                particle[i].vy = V_INITIAL*gaussian();
                break;
            }
            case (VI_COANDA):
            {
                if (vabs(particle[i].yc) < LAMBDA) particle[i].vx = V_INITIAL;
                else particle[i].vx = 0.0;
                particle[i].vy = 0.0;
                break;
            }            
        }
        
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
        if ((PLOT == P_NUMBER)||(PLOT_B == P_NUMBER))
            particle[i].color_hue = 360.0*(double)(i%N_PARTICLE_COLORS)/(double)N_PARTICLE_COLORS;
//         if ((PLOT == P_CLUSTER)||(PLOT_B == P_CLUSTER))
        {
            if (PAIR_PARTICLES) particle[i].cluster = rand()%(ncircles*CLUSTER_COLOR_FACTOR);
            else particle[i].cluster = rand()%ncircles;
        }
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
        particle[i].damping = 1.0;
        particle[i].dirmean = 0.0;
        particle[i].mass_inv = 1.0/PARTICLE_MASS;
        particle[i].inertia_moment_inv = 1.0/PARTICLE_INERTIA_MOMENT;
        particle[i].charge = 0.0;
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
        particle[i].coulomb = 1;
        particle[i].added = 1;
        particle[i].reactive = 1;
        particle[i].paired = 0;
        particle[i].collision = 0;
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
    if (TRACER_PARTICLE) for (j=0; j<n_tracers; j++)
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
        
        n_tracers++;
    }
    
    /* position-dependent particle type */
    if (POSITION_DEPENDENT_TYPE) for (i=0; i<ncircles; i++)
        if (((!POSITION_Y_DEPENDENCE)&&((particle[i].xc - POSITION_DEP_X)*POSITION_DEP_SIGN < 0.0))||((POSITION_Y_DEPENDENCE)&&(particle[i].yc*POSITION_DEP_SIGN < 0.0))) 
        {
            particle[i].type = 2;
            particle[i].mass_inv = 1.0/PARTICLE_MASS_B;
            particle[i].radius = MU_B;
        }
        
    /* add copies in case of particle pairing */
    if (PAIR_PARTICLES) init_particle_pairs(particle, molecule);
    
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
        {
            if (module2(particle[i].xc - obstacle[j].xc, particle[i].yc - obstacle[j].yc) < OBSTACLE_RADIUS + particle[i].radius)
                particle[i].active = 0;
        
            hashcell = hash_cell(obstacle[j].xc, obstacle[j].yc);
            hashgrid[hashcell].charge = OBSTACLE_CHARGE;
        }
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
            case (IC_DNA_POLYMERASE):
            {
                /* do nothing? */
                break;
            }
            case (IC_DNA_POLYMERASE_REC):
            {
                /* do nothing? */
                break;
            }
            default:
            {
                /* do nothing */
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


int add_particles(t_particle particle[NMAXCIRCLES], double px[NMAXCIRCLES], double py[NMAXCIRCLES], int nadd_particle, int type, t_molecule molecule[NMAXCIRCLES], 
int tracer_n[N_TRACER_PARTICLES])
/* add several particles to the system */
{
    static int i = 0;
    double x, y, r, phi; 
    
    printf("Adding a particle\n\n");
    
    switch (ADD_REGION) {
        case (ADD_RECTANGLE) :
        {
            x = ADDXMIN + (ADDXMAX - ADDXMIN)*(double)rand()/(double)RAND_MAX;
            y = ADDYMIN + (ADDYMAX - ADDYMIN)*(double)rand()/(double)RAND_MAX;
            break;
        }
        case (ADD_RING) :
        {
            r = ADDRMIN + (ADDRMAX - ADDRMIN)*(double)rand()/(double)RAND_MAX;
            phi = DPI*(double)rand()/(double)RAND_MAX;
            x = r*cos(phi);
            y = r*sin(phi);
            break;
        }
    }
    
    add_particle(x, y, 0.0, V_INITIAL, PARTICLE_MASS, type, particle);
    
//     add_particle(x, y, V_INITIAL*(double)rand()/RAND_MAX, 0.0, PARTICLE_MASS, type, particle);
        
//     if (y > 0.0)
//         add_particle(x, y, 5.0*V_INITIAL*(double)rand()/RAND_MAX, -10.0*V_INITIAL, PARTICLE_MASS, type, particle);
//     else
//         add_particle(x, y, 5.0*V_INITIAL*(double)rand()/RAND_MAX, 10.0*V_INITIAL, PARTICLE_MASS, type, particle);        

    particle[ncircles - 1].eq_dist = EQUILIBRIUM_DIST;
    particle[ncircles - 1].thermostat = 1;
    px[ncircles - 1] = particle[ncircles - 1].vx;
    py[ncircles - 1] = particle[ncircles - 1].vy;
    
    if ((TRACER_PARTICLE)&&(n_tracers < N_TRACER_PARTICLES))
    {
        tracer_n[n_tracers] = ncircles - 1;
        n_tracers++;
        printf("%i tracer particles\n", n_tracers); 
    }
    
//     init_molecule_data(particle, molecule);

    return (ncircles);
    
    
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

int partial_thermostat_coupling(t_particle particle[NMAXCIRCLES], double xmin, t_segment segment[NMAXSEGMENTS], t_lj_parameters params)
/* only couple particles satisfying condition PARTIAL_THERMO_REGION to thermostat */
{
    int condition, i, nthermo = 0;
    double x, y, height, a, b;
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
            case (TH_RING):
            {
                x = particle[i].xc;
                y = particle[i].yc;
                condition = x*x + y*y > PARTIAL_THERMO_WIDTH*PARTIAL_THERMO_WIDTH;
                break;
            }
            case (TH_RING_EXPAND):
            {
                x = particle[i].xc;
                y = particle[i].yc;
                condition = x*x + y*y > params.thermo_radius*params.thermo_radius;
                break;
            }
            case (TH_INIT):
            {
                x = particle[i].xc;
                y = particle[i].yc;
                condition = ((x > INITXMIN)&&(x < INITXMAX)&&(y > INITYMIN)&&(y < INITYMAX));
                break;
            }
            case (TH_THERMO):
            {
                x = particle[i].xc;
                y = particle[i].yc;
                condition = ((x > THERMOXMIN)&&(x < THERMOXMAX)&&(y > THERMOYMIN)&&(y < THERMOYMAX));
                break;
            }
            case (TH_CONE):
            {
                x = particle[i].xc;
                y = particle[i].yc;
                a = (LAMBDA + 0.2)/(LAMBDA - 0.05);
                b = LAMBDA*(1.0 - a);
                condition = ((y < LAMBDA)&&(y > 0.0)&&(y > a*vabs(x) + b));
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


int partial_bfield(double x, double y)
/*  */
{
    switch (BFIELD_REGION) {
        case (BF_CONST): return(1);
        case (BF_SQUARE): return(vabs(x) < YMAX);
    }
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

// int chem_dissociate_water_old(int i, t_particle particle[NMAXCIRCLES], t_collision *collisions, int ncollisions, double reaction_prob)
// /* dissociation of water molecule i into a pair H+, OH- */
// {
//     int j, k, p, q, closeby = 0, reaction = 0;
//     double distance, move_factor, xnew, ynew, vxnew, vynew, mr1, mr2; 
//     
//     if ((particle[i].npartners == 2)&&((double)rand()/RAND_MAX < reaction_prob))
//     {
//         /* choose particle k to dissociate */
//         if (rand()%2 == 0)
//         {
//             j = particle[i].partner[0];
//             k = particle[i].partner[1];
//         }
//         else
//         {
//             j = particle[i].partner[1];
//             k = particle[i].partner[0];
//         }
//         
//         /* propose moving dissociated particles further away */
//         move_factor = 2.0;
//         xnew = particle[i].xc + move_factor*(particle[k].xc - particle[i].xc);
//         ynew = particle[i].yc + move_factor*(particle[k].yc - particle[i].yc);
//         
//         /* test distance to other particles */
//         closeby = 0;
//         for (p=0; p<particle[k].hash_nneighb; p++) 
//         {
//             q = particle[k].hashneighbour[p];
//             if ((q!=i)&&(q!=j)&&(q!=k)&&(particle[q].active))
//             {
//                 distance = module2(xnew - particle[q].xc, ynew - particle[q].yc);
//                 if (distance < SAFETY_FACTOR*MU_B) closeby = 1;
//             }
//         }
//     
//         if (closeby) printf("Cannot split molecule %i\n", i);
//         else
//         {
//             printf("Splitting molecule %i\n", i);
//                 
//             particle[i].npartners = 1;
//             particle[i].partner[0] = j;
//         
//             particle[j].npartners = 1;
//             particle[j].partner[0] = i;
//         
//             particle[k].npartners = 0;
//             particle[k].xc = xnew;
//             particle[k].yc = ynew;
//         
//             collisions[ncollisions].x = particle[i].xc;
//             collisions[ncollisions].y = particle[i].yc;
//             collisions[ncollisions].time = COLLISION_TIME;
//             collisions[ncollisions].color = type_hue(1);
//             
//             if (ncollisions < 2*NMAXCOLLISIONS - 1) ncollisions++;
//             else printf("Too many collisions\n");
//         }
//     }
//     else if (particle[i].npartners == 1)
//     {
//         mr1 = PARTICLE_MASS/(PARTICLE_MASS + 2.0*PARTICLE_MASS_B);
//         mr2 = 1.0 - mr1;
//         
// //         for (p=0; p<particle[i].hash_nneighb; p++) 
//         
//         p = 0;
//         while ((p<particle[i].hash_nneighb)&&(!reaction))
//         {
//             k = particle[i].hashneighbour[p];
//             if ((particle[k].active)&&(particle[k].type == 2)&&(particle[k].npartners == 0))
//             {
//                 distance  = module2(particle[i].deltax[p], particle[i].deltay[p]);
//                 if ((distance < REACTION_DIST*MU)&&((double)rand()/RAND_MAX < REACTION_PROB))
//                 {
//                     reaction = 1;
//                     printf("Merging molecule %i with particle %i\n", i, k);
//         
//                     j = particle[i].partner[0];
//                     particle[i].npartners = 2;
//                     particle[i].partner[1] = k;
//                     
//                     particle[j].npartners = 2;
//                     particle[j].partner[0] = i;     /* needed? */
//                     particle[j].partner[1] = k;
//                     
//                     particle[k].npartners = 2;
//                     particle[k].partner[0] = i;
//                     particle[k].partner[1] = j;
//                     
//                     distance = 2.0*sin(PARTNER_ANGLE*PI/360.0)*(MU + MU_B)*PAIR_DRATIO;
//                     particle[k].partner_eqd[1] = distance;
//                     particle[j].partner_eqd[1] = distance;
//                       
//                     /* equalize speeds */
//                     vxnew = mr1*particle[i].vx + mr2*particle[j].vx + mr2*particle[k].vx;
//                     vynew = mr1*particle[i].vy + mr2*particle[j].vy + mr2*particle[k].vy;
//                     
//                     particle[i].vx = vxnew;
//                     particle[i].vy = vynew;
//                     particle[j].vx = vxnew;
//                     particle[j].vy = vynew;
//                     particle[k].vx = vxnew;
//                     particle[k].vy = vynew;
//                     
//                     collisions[ncollisions].x = particle[i].xc;
//                     collisions[ncollisions].y = particle[i].yc;
//                     collisions[ncollisions].time = COLLISION_TIME;
//                     collisions[ncollisions].color = 0.0;
// //                     collisions[ncollisions].color = type_hue(2);
//                 
//                     if (ncollisions < 2*NMAXCOLLISIONS - 1) ncollisions++;
//                     else printf("Too many collisions\n");
//                 }
//             }
//             p++;
//         }
//     }
//     
//     return(ncollisions);
// }

int chem_dissociate_molecule(int i, t_particle particle[NMAXCIRCLES], t_collision *collisions, int ncollisions, double reaction_prob)
/* dissociation of water molecule i into a pair H+, OH- */
{
    int k, p, q, r, ndiss, closeby = 0, reaction = 0, diff, np, j[NMAXPARTNERS];
    double distance, move_factor, xnew, ynew, vxnew, vynew, mr1, mr2; 
    
    if ((double)rand()/RAND_MAX < reaction_prob)
    {
        np = particle[i].npartners;
    
        /* choose particle k to dissociate */
        ndiss = rand()%np;
        k = particle[i].partner[ndiss];
        
        /* new list of partners of i */
        for (p=0; p<ndiss; p++) j[p] = particle[i].partner[p];
        for (p=ndiss+1; p<np; p++) j[p-1] = particle[i].partner[p];
                
        /* propose moving dissociated particle further away */
        move_factor = 2.0 + (double)(np-2);
        xnew = particle[i].xc + move_factor*(particle[k].xc - particle[i].xc);
        ynew = particle[i].yc + move_factor*(particle[k].yc - particle[i].yc);
        
        /* test distance to other particles */
        closeby = 0;
        for (p=0; p<particle[k].hash_nneighb; p++) 
        {
            q = particle[k].hashneighbour[p];
            diff = ((q!=i)&&(q!=k)&&(particle[q].active));
            for (r=0; r<np-1; r++) if (q=j[r]) diff = 0;
            if (diff)
            {
                distance = module2(xnew - particle[q].xc, ynew - particle[q].yc);
                if (distance < SAFETY_FACTOR*MU_B) closeby = 1;
            }
        }
    
        if (closeby) printf("Cannot split molecule %i\n", i);
        else
        {
            printf("Splitting molecule %i\n", i);
                
            particle[i].npartners = np-1;
            for (p=0; p<np-1; p++) particle[i].partner[p] = j[p];
        
            for (p=0; p<np-1; p++) 
            {
                particle[j[p]].npartners = np-1;
                particle[j[p]].partner[0] = i;
                
                for (q=0; q<p; q++) particle[j[p]].partner[q+1] = j[q];
                for (q=p+1; q<np-1; q++) particle[j[p]].partner[q] = j[q];
            } 
            particle[k].npartners = 0;
            particle[k].xc = xnew;
            particle[k].yc = ynew;
        
            collisions[ncollisions].x = particle[i].xc;
            collisions[ncollisions].y = particle[i].yc;
            collisions[ncollisions].time = COLLISION_TIME;
            collisions[ncollisions].color = partner_color(np);
//             collisions[ncollisions].color = type_hue(1);
            
            if (ncollisions < 2*NMAXCOLLISIONS - 1) ncollisions++;
            else printf("Too many collisions\n");
        }
    }
    return(ncollisions);
}



int chem_merge_molecule(int i, int type2, int maxpartners, t_particle particle[NMAXCIRCLES], t_collision *collisions, int ncollisions, double reaction_prob)
/* merging of molecule i with another molecule of type type 2 */
/* having at most nmaxpartners partners already */
{
    int k, p, q, np, n, m, closeby = 0, reaction = 0, j[NMAXPARTNERS], r, kk;
    double distance, angle, move_factor, xnew, ynew, vxnew, vynew, mr1, mr2; 
    
    np = particle[i].npartners;
    mr1 = PARTICLE_MASS/(PARTICLE_MASS + (double)(np+1)*PARTICLE_MASS_B);
    mr2 = 1.0 - mr1;
        
    p = 0;
    while ((p<particle[i].hash_nneighb)&&(!reaction))
    {
        k = particle[i].hashneighbour[p];
        if ((particle[k].active)&&(particle[k].type == type2)&&(particle[k].npartners <= maxpartners))
        {
            distance  = module2(particle[i].deltax[p], particle[i].deltay[p]);
            if ((distance < REACTION_DIST*MU)&&((double)rand()/RAND_MAX < reaction_prob))
            {
                reaction = 1;
                printf("Merging molecule %i with particle %i\n", i, k);
    
                for (r=0; r<np; r++) j[r] = particle[i].partner[r];
                particle[i].npartners++;
                particle[i].partner[np] = k;
                
                for (r=0; r<np; r++)
                {
                    particle[j[r]].npartners = np+1;
                    particle[j[r]].partner[np] = k;
                }
                    
                particle[k].npartners = np+1;
                particle[k].partner[0] = i;
                particle[k].partner_eqd[0] = (MU + MU_B)*PAIR_DRATIO;
                for (kk=0; kk<np; kk++) particle[k].partner[kk+1] = j[kk];
                  
                if (np == 1)
                {
                    distance = 2.0*sin(PARTNER_ANGLE*PI/360.0)*(MU + MU_B)*PAIR_DRATIO;
                    particle[k].partner_eqd[1] = distance;
                    particle[j[0]].partner_eqd[1] = distance;                    
                }
                else    /* distribute particles uniformly around i */
                {
                    j[np] = k;
                    print_partners(i, particle);
                    for (q=0; q<np+1; q++)
                    {
                        angle = (double)q*DPI/(double)(np+1);
                        particle[j[q]].npartners = np+1;
                        particle[j[q]].partner[0] = i;
                        particle[j[q]].partner_eqd[0] = (MU + MU_B)*PAIR_DRATIO;
                        particle[i].partner_eqd[q] = (MU + MU_B)*PAIR_DRATIO;
                        particle[j[q]].xc = particle[i].xc + (MU + MU_B)*PAIR_DRATIO*cos(angle);
                        particle[j[q]].yc = particle[i].yc + (MU + MU_B)*PAIR_DRATIO*sin(angle);
                    }
                    for (q=0; q<np+1; q++)
                        for (r=0; r<np+1; r++) 
                        {
                            m = particle[j[q]].partner[r];
                            particle[j[q]].partner_eqd[r] = module2(particle[j[q]].xc - particle[m].xc, particle[j[q]].yc - particle[m].yc);
                        }                        
//                     for (q=0; q<np+1; q++) print_partners(j[q], particle);
                }
                    
                /* equalize speeds */
                vxnew = mr1*particle[i].vx + mr2*particle[k].vx;
                vynew = mr1*particle[i].vy + mr2*particle[k].vy;
                for (r=0; r<np; r++)
                {
                    vxnew += mr2*particle[j[r]].vx;
                    vynew += mr2*particle[j[r]].vy;
                }
                                
                particle[i].vx = vxnew;
                particle[i].vy = vynew;
                for (r=0; r<np; r++)
                {
                    particle[j[r]].vx = vxnew;
                    particle[j[r]].vy = vynew;
                }
                particle[k].vx = vxnew;
                particle[k].vy = vynew;
                    
                collisions[ncollisions].x = particle[i].xc;
                collisions[ncollisions].y = particle[i].yc;
                collisions[ncollisions].time = COLLISION_TIME;
                collisions[ncollisions].color = 0.0;
                
                if (ncollisions < 2*NMAXCOLLISIONS - 1) ncollisions++;
                else printf("Too many collisions\n");
            }
        }
        p++;
    }
    return(ncollisions);
}


int chem_multi_glue_molecule(int i, int type2, int maxpartners, int no_triangles, int no_cluster, int require_charge, int equalize_charge, t_particle particle[NMAXCIRCLES], t_molecule molecule[NMAXCIRCLES], t_collision *collisions, int ncollisions, double reaction_prob)
/* glue molecule containing particle i with another molecule of type type 2 */
/* having at most nmaxpartners partners already */
/* if no_triangles >= 1, particles that have a common partner are not paired */
/* if no_cluster = 1, particles are not paired if they belong to the same cluster */
{
    int k, p, q, p1, q1, p2, np, nq, n, m, closeby = 0, reaction = 0, jp[NMAXPARTNERS], jq[NMAXPARTNERS], r, kk, different, np2, moli, molk, mp, newpartner, pp, type1;
    double distance, angle, move_factor, xnew, ynew, vxnew, vynew, m1, m2, mr1, mr2, deltav, x, y; 
    short int charge_condition, triangle_condition; 
    
    type1 = particle[i].type;
    np = particle[i].npartners;
    moli = particle[i].molecule;
    if (np > maxpartners) return(ncollisions);
    
    m1 = 1.0/particle[i].mass_inv;
    for (p=0; p<np; p++) 
    {
        jp[p] = particle[i].partner[p];
        m1 += 1.0/particle[jp[p]].mass_inv;
    }
        
    p = 0;
    while ((p<particle[i].hash_nneighb)&&(!reaction))
    {
        k = particle[i].hashneighbour[p];
        nq = particle[k].npartners;
        molk = particle[k].molecule;
        charge_condition = (!require_charge)||(particle[i].charge*particle[k].charge < 0.0);
        triangle_condition = 1;
        
        /* exclude mergers between added particles */
        if (particle[i].added == particle[k].added)
        {
            triangle_condition = particle[i].reactive*particle[k].reactive;
            /* allow mergers if molecules are paired to same DNA strand */
            if ((particle[i].paired)&&(particle[k].paired)&&(particle[i].partner_molecule % NGRIDY == particle[k].partner_molecule % NGRIDY)) triangle_condition = 1;
        }
        
        if ((RD_REACTION == CHEM_DNA_ENZYME)||(RD_REACTION == CHEM_DNA_ENZYME_REPAIR))
        {
            /* exclude mergers of added molecules to original ones via backbone */
            if ((particle[i].added != particle[k].added)&&(type2 <= 2)&&(type1 <= 2))                 
                triangle_condition = 0;
            
            /* NEW TEST */
//             if (molk == moli) triangle_condition = 0;
            
            /* NEW TEST */
            /* avoid original backbones to form a loop */
            if ((!particle[i].added)&&(!particle[k].added)&&((molk-moli > NGRIDX/2)||(moli-molk > NGRIDX/2))) triangle_condition = 0;
        }
        
        if (no_triangles == 1)   /* do not pair i and k if they have a common partner */
        {
            for (q=0; q<nq; q++)
                for (p1=0; p1<np; p1++)
                    if (particle[k].partner[q] == particle[i].partner[p1]) triangle_condition = 0;
        }
        if (no_triangles == 2)  /* also exclude i or k being a partner of k or i */
        {
            for (q=0; q<nq; q++) if (particle[k].partner[q] == i) triangle_condition = 0;
            for (p1=0; p1<np; p1++) if (particle[i].partner[p1] == k) triangle_condition = 0;
        }
        if (no_triangles >= 3) /* exclude connections of different type between the same molecules */
        {
            if (molk != moli) for (p1=0; p1<np; p1++)
            {
                n = particle[i].partner[p1]; 
                for (q=0; q<particle[n].npartners; q++)
                {
                    m = particle[n].partner[q];
                    if (particle[m].molecule == molk)
                        if (particle[n].type != particle[i].type)
                            triangle_condition = 0;
                }
                
                if ((particle[n].molecule != moli)&&(particle[n].molecule != molk))
                    triangle_condition = 0;
                
                if (particle[n].molecule == moli) if (particle[n].type == particle[i].type)
                {
                    for (q=0; q<particle[n].npartners; q++)
                    {
                        m = particle[n].partner[q];
                        if ((particle[m].molecule != moli)&&(particle[m].molecule != molk))
                            triangle_condition = 0;
                    }
                }
            }
            
            if (molk != moli) for (q=0; q<nq; q++)
            {
                n = particle[k].partner[q];
                if ((particle[n].molecule != molk)&&(particle[n].molecule != moli))
                    triangle_condition = 0;
            }
        }
        if (no_triangles >= 4)   /* also exclude second-order partners */
        {
            for (q=0; q<nq; q++)
                for (p1=0; p1<np; p1++)
                {
                    n = particle[k].partner[q];
                    for (q1 = 0; q1 < particle[n].npartners; q1++) 
                        if (particle[n].partner[q1] == particle[i].partner[p1]) 
                            triangle_condition = 0;
                    
                    n = particle[i].partner[p1];
                    for (p2 = 0; p2 < particle[n].npartners; p2++) 
                        if (particle[n].partner[p2] == particle[k].partner[q]) 
                            triangle_condition = 0;
                }
        }
        if ((no_cluster)&&(particle[k].cluster = particle[i].cluster)) triangle_condition = 0;
        
        
        
        if ((particle[k].active)&&(charge_condition)&&(triangle_condition)&&(particle[k].type == type2)&&(nq <= maxpartners))
        {
            distance  = module2(particle[i].deltax[p], particle[i].deltay[p]);
            if ((distance < REACTION_DIST*MU)&&((double)rand()/RAND_MAX < reaction_prob))
            {
                m2 = 1.0/particle[k].mass_inv;
                for (q=0; q<nq; q++) 
                {
                    jq[q] = particle[k].partner[q];
                    m2 += 1.0/particle[jq[q]].mass_inv;
                }
                
                mr1 = m1/(m1 + m2);
                mr2 = 1.0 - mr1;   
                
                if ((np == 0)&&(nq == 0))
                {
                    reaction = 1;
                    printf("Merging molecule %i with particle %i\n", i, k);
                
                    distance = (particle[i].radius + particle[k].radius)*PAIR_DRATIO;

                    particle[i].npartners++;
                    particle[i].partner[0] = k;
                    particle[i].partner_eqd[0] = distance;    
                    
                    particle[k].npartners++;
                    particle[k].partner[0] = i;
                    particle[k].partner_eqd[0] = distance;
                    
                    /* equalize speeds */
                    vxnew = mr1*particle[i].vx + mr2*particle[k].vx;
                    vynew = mr1*particle[i].vy + mr2*particle[k].vy;
                    
                    particle[i].vx = vxnew;
                    particle[i].vy = vynew;
                    particle[k].vx = vxnew;
                    particle[k].vy = vynew;
                
                    /* equalize charges */
                    if (equalize_charge == 1) particle[k].charge = particle[i].charge;
                    else if (equalize_charge == 2)
                    {
                        particle[i].coulomb = 0;
                        particle[k].coulomb = 0;
                    }
                    
                    /* update molecule data */
                    mp = molecule[moli].npartners;
                    if (mp < NMAXPARTNERMOLECULES)
                    {
                        newpartner = 1;
                        for (pp=0; pp<mp; pp++)
                            if (molecule[moli].partner[pp] == molk) newpartner = 0;
                        if (newpartner)
                        {
                            printf("adding molecule %i as %ith partner of molecule %i\n", molk, mp, moli);
                            molecule[moli].partner[mp] = molk;
                            molecule[moli].connection_type[mp] = particle[i].type;
                            molecule[moli].npartners++;
                        }
                    }
                    mp = molecule[molk].npartners;
                    if (mp < NMAXPARTNERMOLECULES)
                    {
                        newpartner = 1;
                        for (pp=0; pp<mp; pp++)
                            if (molecule[molk].partner[pp] == moli) newpartner = 0;
                        if (newpartner)
                        {
                            printf("adding molecule %i as %ith partner of molecule %i\n", moli, mp, molk);
                            molecule[molk].partner[mp] = moli;
                            molecule[molk].connection_type[mp] = type2;
                            molecule[molk].npartners++;
                        }
                    }
                }
                else
                {
                    /* check if i and k belong to different molecules */
//                     different = 1;
//                     for (p1=0; p1<np; p1++) 
//                         if (k == jp[p1]) different = 0;
//                     for (q1=0; q1<nq; q1++) 
//                         if (i == jq[q1]) different = 0;
                        /* TODO: replace by moli != molk ? */
                    different = (moli != molk);
                    
                    if (different)
                    {
                        deltav = module2(particle[i].vx - particle[k].vx, particle[i].vy - particle[k].vy);
//                         printf("Delta v is %.3lg\n", deltav);
                    
                        if (deltav < DELTAVMAX) 
                        {
                            reaction = 1;
                            printf("Pairing molecules %i and %i containing particles %i and %i\n", moli, molk, i, k);
                        
                            distance = (particle[i].radius + particle[k].radius)*PAIR_DRATIO;

                            particle[i].npartners++;
                            particle[i].partner[np] = k;
                            particle[i].partner_eqd[np] = distance;    
                    
                            particle[k].npartners++;
                            particle[k].partner[nq] = i;
                            particle[k].partner_eqd[nq] = distance;
                    
                            /* equalize speeds */
                            vxnew = mr1*particle[i].vx + mr2*particle[k].vx;
                            vynew = mr1*particle[i].vy + mr2*particle[k].vy;
                    
                            particle[i].vx = vxnew;
                            particle[i].vy = vynew;
                            particle[k].vx = vxnew;
                            particle[k].vy = vynew;
                            
                            if (equalize_charge) particle[k].charge = particle[i].charge;
                        
                            for (p1=0; p1<np; p1++) for (q1=0; q1<nq; q1++)
                            {
                                particle[jp[p1]].vx = vxnew;
                                particle[jp[p1]].vy = vynew;
                                particle[jq[q1]].vx = vxnew;
                                particle[jq[q1]].vy = vynew;
                                
                                if (equalize_charge == 1) 
                                {
                                    particle[jp[p1]].charge = particle[i].charge;
                                    particle[jq[q1]].charge = particle[i].charge;
                                }
                                else if (equalize_charge == 2)
                                {
                                    particle[jp[p1]].coulomb = 0;
                                    particle[jq[q1]].coulomb = 0;
                                }
                            }
                            
                            /* update molecule data */
                            mp = molecule[moli].npartners;
                            if (mp < NMAXPARTNERMOLECULES)
                            {
                                newpartner = 1;
                                for (pp=0; pp<mp; pp++)
                                    if (molecule[moli].partner[pp] == molk) newpartner = 0;
                                if (newpartner)
                                {
//                                     printf("adding molecule %i as %ith partner of molecule %i\n", molk, mp, moli);
                                    molecule[moli].partner[mp] = molk;
                                    molecule[moli].connection_type[mp] = particle[i].type;
                                    molecule[moli].npartners++;
                                }
                            }
                            mp = molecule[molk].npartners;
                            if (mp < NMAXPARTNERMOLECULES)
                            {
                                newpartner = 1;
                                for (pp=0; pp<mp; pp++)
                                    if (molecule[molk].partner[pp] == moli) newpartner = 0;
                                if (newpartner)
                                {
//                                     printf("adding molecule %i as %ith partner of molecule %i\n", moli, mp, molk);
                                    molecule[molk].partner[mp] = moli;
                                    molecule[molk].connection_type[mp] = particle[k].type;
                                    molecule[molk].npartners++;
                                }
                            }
                            
                            /* merge secondary partners, experimental */
                            if (SECONDARY_PAIRING) 
                            {
                                np2 = particle[k].npartners;
                                
                                for (p1 = 0; p1 < np2; p1++)
//                                 if (np2 > 1)
                                {
                                    n = particle[k].partner[p1];
                                    x = particle[i].xc - particle[n].xc;
                                    y = particle[i].yc - particle[n].yc;
                                    /* deal with periodic boundary conditions */
                                    if (x > 0.5*(XMAX - XMIN)) x -= (XMAX - XMIN);
                                    else if (x < 0.5*(XMIN - XMAX)) x += (XMAX - XMIN);
                                    if (y > 0.5*(YMAX - YMIN)) y -= (YMAX - YMIN);
                                    else if (y < 0.5*(YMIN - YMAX)) y += (YMAX - YMIN);
                                    distance = module2(x, y);
                                    
                                    if (distance > particle[i].radius)
                                    {
                                        np = particle[i].npartners;
                                        particle[i].npartners++;
                                        particle[i].partner[np] = n;
                                        particle[i].partner_eqd[np] = distance;    
                    
                                        nq = particle[n].npartners;
                                        particle[n].npartners++;
                                        particle[n].partner[nq] = i;
                                        particle[n].partner_eqd[nq] = distance;
                                    }
                                }
                            }
                        }
                    }
                }
                    
                if (reaction)
                {
                    /* TEST */
                    if ((RD_REACTION == CHEM_DNA_ENZYME)||(RD_REACTION == CHEM_DNA_ENZYME_REPAIR))
                    {
                        /* consider merged marticles as reactive */
                        if ((particle[i].added == 1)&&(particle[i].type >= 3)&&(particle[i].type <= 6))
                        {
                            particle[i].paired = 1;
                            particle[i].partner_molecule = molk;
                            for (p1 = 0; p1<particle[i].npartners; p1++)
                            {
                                n = particle[i].partner[p1];
                                particle[n].paired = 1;
                                particle[n].partner_molecule = molk;
                            }
                        }
                        
                        if ((particle[k].added == 1)&&(type2 >= 3)&&(type2 <= 6))
                        {
                            particle[k].paired = 1;
                            particle[k].partner_molecule = moli;
                            for (p1 = 0; p1<particle[k].npartners; p1++)
                            {
                                n = particle[k].partner[p1];
                                particle[n].paired = 1;
                                particle[n].partner_molecule = moli;
                            }
                        }
                    }
                    
                    collisions[ncollisions].x = particle[i].xc;
                    collisions[ncollisions].y = particle[i].yc;
                    collisions[ncollisions].time = COLLISION_TIME;
                    collisions[ncollisions].color = 0.0;
                    
                    if (ncollisions < 2*NMAXCOLLISIONS - 1) ncollisions++;
                    else printf("Too many collisions\n");
                }
            }
        }
        p++;
    }
    return(ncollisions);
}


void translate_cluster(int j, t_cluster cluster[NMAXCIRCLES], t_particle particle[NMAXCIRCLES], double dx, double dy)
/* translate a cluster and all partcles it contains */
{
    int k, p, np;
        
    np = cluster[j].nparticles;
    
    if (np == 1)
    {
        p = cluster[j].particle[0];
        particle[p].xc += dx;
        particle[p].yc += dy;
        cluster[j].xg = particle[p].xc;
        cluster[j].yg = particle[p].yc;
    }
    else 
    {
        for (k=0; k<np; k++)
        {
            p = cluster[j].particle[k];
            particle[p].xc += dx;
            particle[p].yc += dy;
        }
        cluster[j].xg += dx;
        cluster[j].yg += dy;
    }
}

void rotate_cluster(int j, t_cluster cluster[NMAXCIRCLES], t_particle particle[NMAXCIRCLES], double angle)
/* translate a cluster and all partcles it contains */
{
    int k, p, np;
    double ca, sa, x, y;
        
    np = cluster[j].nparticles;
    
    if (np == 1)
    {
        p = cluster[j].particle[0];
        particle[p].angle += angle;
        cluster[j].angle = particle[p].angle;
    }
    else 
    {
        ca = cos(angle);
        sa = sin(angle);
        for (k=0; k<np; k++)
        {
            p = cluster[j].particle[k];
            particle[p].angle += angle;
            x = particle[p].xc - cluster[j].xg;
            y = particle[p].yc - cluster[j].yg;
            particle[p].xc = cluster[j].xg + ca*x - sa*y;
            particle[p].yc = cluster[j].yg + sa*x + ca*y;
        }
        cluster[j].angle += angle;
    }
}

void rotate_cluster_around_particle(int j, int i, t_cluster cluster[NMAXCIRCLES], t_particle particle[NMAXCIRCLES], double angle)
/* translate a cluster and all partcles it contains around center of particle i */
{
    int k, p, np;
    double ca, sa, x, y;
        
    np = cluster[j].nparticles;
    
    if (np == 1)
    {
        p = cluster[j].particle[0];
        particle[p].angle += angle;
        cluster[j].angle = particle[p].angle;
    }
    else 
    {
        ca = cos(angle);
        sa = sin(angle);
        for (k=0; k<np; k++)
        {
            p = cluster[j].particle[k];
            particle[p].angle += angle;
            x = particle[p].xc - particle[i].xc;
            y = particle[p].yc - particle[i].yc;
            particle[p].xc = particle[i].xc + ca*x - sa*y;
            particle[p].yc = particle[i].yc + sa*x + ca*y;
        }
        cluster[j].angle += angle;
    }
}

void translate_and_rotate_cluster(int j, t_cluster cluster[NMAXCIRCLES], t_particle particle[NMAXCIRCLES], double dx, double dy, double angle)
/* translate and rotate a cluster and all partcles it contains */
{
    int k, p, np;
    double ca, sa, x, y;
    
    np = cluster[j].nparticles;
    
    if (np == 1)
    {
        p = cluster[j].particle[0];
        particle[p].xc += dx;
        particle[p].yc += dy;
        particle[p].angle += angle;
        cluster[j].xg = particle[p].xc;
        cluster[j].yg = particle[p].yc;
        cluster[j].angle = particle[p].angle;
    }
    else
    {
        ca = cos(angle);
        sa = sin(angle);
        for (k=0; k<np; k++)
        {
            p = cluster[j].particle[k];
            x = particle[p].xc - cluster[j].xg;
            y = particle[p].yc - cluster[j].yg;
            particle[p].xc = cluster[j].xg + ca*x - sa*y + dx;
            particle[p].yc = cluster[j].yg + sa*x + ca*y + dy;
        }
        cluster[j].xg += dx;
        cluster[j].yg += dy;
        cluster[j].angle += angle;
    }
}


int merge_clusters(int i, int j, int pi, int pj, t_particle particle[NMAXCIRCLES], t_cluster cluster[NMAXCIRCLES], int adjust_angles, short int flip)
/* merge clusters i and j to form new cluster il */
/* il is number of largest cluster aming i and j */
/* pi and pj are the numbers of the contact particles */
/* returns number of new (largest) cluster il */
{
    int p, q, q0, c0, is, il, ns, nl, ps, pl;
    double mtot, dangle, d2, newangle;
    static double alpha, alpha2;
    static int first = 1, nmergers = 0;
    
    if (first)
    {
        alpha = PI/(double)NPOLY;
        alpha2 = 2.0*alpha;
        first = 0;
    }
    
    if (i == j) return(i);
    
    nmergers++;
    
    /* find largest cluster */
    if (cluster[i].nparticles >= cluster[j].nparticles)
    {
        il = i;
        is = j;
        pl = pi;
        ps = pj;
    }
    else 
    {
        il = j;
        is = i;
        pl = pj;
        ps = pi;
    }
    ns = cluster[is].nparticles;
    nl = cluster[il].nparticles;
    
    if (ns + nl >= NMAXPARTINCLUSTER) 
    {
        printf("Warning: NMAXPARTINCLUSTER is too small\n");
        printf("Try increasing to %i\n", cluster[i].nparticles + cluster[j].nparticles);
        return(0);
    }
    
    printf("Merging clusters %i and %i, having %i and %i particles\n", il, is, nl, ns);
    fprintf(lj_log, "Merging clusters %i and %i, having %i and %i particles\n", il, is, nl, ns);
    
    fprintf(lj_log, "Large cluster %i:\n", il);
    for (p=0; p<nl; p++) 
    {
        q = cluster[il].particle[p];
        fprintf(lj_log, "Particle %i (%i), angle %.3lg, flip %i\n", p, q, particle[q].angle*180.0/PI, particle[q].flip); 
    }
    fprintf(lj_log, "Reacting particle %i\n\n", pl);
    
//     fprintf(lj_log, "Cluster %i:\n", is);
//     for (p=0; p<ns; p++) 
//     {
//         q = cluster[is].particle[p];
//         fprintf(lj_log, "Particle %i (%i), angle %.3lg, flip %i\n", p, q, particle[q].angle*180.0/PI, particle[q].flip); 
//     }
//     fprintf(lj_log, "Reacting particle %i\n", ps);
    
    cluster[is].active = 0;
    
    /* update list of particles */
    cluster[il].nparticles = ns + nl;
    c0 = particle[cluster[il].particle[0]].cluster_color;
    for (p=0; p<ns; p++) 
    {
        q = cluster[is].particle[p];
        cluster[il].particle[nl+p] = q;
        particle[q].cluster = il;
        particle[q].cluster_color = c0;
        particle[q].collision = nmergers + 1;
    }
    
    /* update cluster size for P_CLUSTER_SIZE color scheme */
    for (p=0; p<ns+nl; p++) 
    {
        q = cluster[il].particle[p];
        particle[q].cluster_size = ns + nl;    
    }
    
    if ((NPOLY%2 == 1)&&(flip)) for (p=0; p<ns; p++)
    {
        q = cluster[il].particle[nl + p];
        particle[q].flip = 1 - particle[q].flip;
    }
        
    fprintf(lj_log, "Small cluster %i:\n", is);
    for (p=0; p<ns; p++) 
    {
        q = cluster[is].particle[p];
        fprintf(lj_log, "Particle %i (%i), angle %.3lg, flip %i\n", p, q, particle[q].angle*180.0/PI, particle[q].flip); 
    }
    fprintf(lj_log, "Reacting particle %i\n\n", ps);
    
    
    /* adjust angles */
    /* use rotate_cluster ? */
    if (adjust_angles)
    {
        for (p=0; p<nl; p++)
        {
            q = cluster[il].particle[p];
//             printf("p = %i, angle = %.5lg\n", p, particle[q].angle*180.0/PI); 
            while (particle[q].angle >= alpha2) particle[q].angle -= alpha2; 
            while (particle[q].angle < 0.0) particle[q].angle += alpha2; 
        }
        q = cluster[il].particle[0];
        newangle =  particle[q].angle;
        cluster[il].angle = newangle;
//         if (NPOLY%2==1) newangle += alpha;
        for (p=nl; p<ns+nl; p++) 
        {
            q = cluster[il].particle[p];
            particle[q].angle = newangle;
//             if ((NPOLY%2==1)&&(particle[q].flip != particle[ps].flip)) particle[q].angle += PI;
        }
//         if (NPOLY%2==1) for (p=nl; p<ns+nl; p++)
//         {
//             q = cluster[il].particle[p];
//             if (particle[q].flip != particle[ps].flip) particle[q].angle += alpha;
//         }
        
        if (NPOLY%2==1) for (p=0; p<ns+nl; p++)
        {
            q = cluster[il].particle[p];
            if (particle[q].flip == particle[0].flip) particle[q].angle = newangle;
            else particle[q].angle = newangle + alpha;
        }
    }
    
    /* total mass */
    mtot = cluster[is].mass + cluster[il].mass;
    
    /* distance between centers of gravity squared */
    d2 = (cluster[i].xg - cluster[j].xg)*(cluster[i].xg - cluster[j].xg);
    d2 += (cluster[i].yg - cluster[j].yg)*(cluster[i].yg - cluster[j].yg);
    
    /* update center of gravity */
    cluster[il].xg = (cluster[is].xg + cluster[il].xg)/mtot;
    cluster[il].yg = (cluster[is].yg + cluster[il].yg)/mtot;
    
    /* update moment of intertia (using parallel axis thm) */
    cluster[il].inertia_moment += cluster[is].inertia_moment;
    cluster[il].inertia_moment += d2*cluster[il].mass*cluster[is].mass/mtot;
    cluster[il].inertia_moment_inv = 1.0/cluster[il].inertia_moment;
    
    /* update total mass */
    cluster[il].mass = mtot;
    cluster[il].mass_inv = 1.0/mtot;
    
    
    fprintf(lj_log, "Merged cluster %i:\n", is);
    for (p=0; p<ns+nl; p++) 
    {
        q = cluster[il].particle[p];
        fprintf(lj_log, "Particle %i (%i), angle %.3lg, flip %i\n", p, q, particle[q].angle*180.0/PI, particle[q].flip); 
    }
    fprintf(lj_log, "Reacting particle %i\n\n", ps);
    
//     short int active;           /* has value 1 if cluster is active */
//     short int thermostat;       /* has value 1 if cluster is coupled to thermostat */
//     double xg, yg;              /* center of gravity */
//     double vx, vy;              /* velocity of center of gravity */
//     double angle;               /* orientation of cluster */
//     double omega;               /* angular velocity of cluster */
//     double mass, mass_inv;      /* mass of cluster and its inverse */
//     double inertia_moment, inertia_moment_inv;   /* moment of inertia */
//     double fx, fy, torque;      /* force and torque */
//     double energy, emean;       /* energy and averaged energy */
//     double dirmean;             /* time-averaged direction */
//     int nparticles;             /* number of particles in cluster */
//     int particle[NMAXPARTINCLUSTER];    /* list of particles in cluster */    
    return(il);
}


int chem_multi_glue_polygon(int i, int type2, int maxpartners, t_particle particle[NMAXCIRCLES], t_molecule molecule[NMAXCIRCLES], t_cluster cluster[NMAXCIRCLES], t_collision *collisions, int ncollisions, double reaction_prob, int mergeclusters)
/* simplified version of chem_multi_glue_molecule without triangle and similar */
/* conditions, but taking polygon geometry into account */
{
    int k, p, q, p1, q1, p2, np, nq, n, m, closeby = 0, reaction = 0, jp[NMAXPARTNERS], jq[NMAXPARTNERS], jtot[NMAXPARTNERS], r, kk, different, np2, moli, molk, mp, newpartner, pp, type1, poly_condition, nangle, norient;
    double distance, angle, move_factor, xnew, ynew, vxnew, vynew, m1, m2, mr1, mr2, deltav, x, y, ri, rk, alpha, alpha2, ca, rel_angle, rorient, delta, rorient_new, anglei_new, deltax, deltay, x1, y1; 
    
    type1 = particle[i].type;
    np = particle[i].npartners;
    moli = particle[i].molecule;
    if (np > maxpartners) return(ncollisions);
    
    m1 = 1.0/particle[i].mass_inv;
    for (p=0; p<np; p++) 
    {
        jp[p] = particle[i].partner[p];
        m1 += 1.0/particle[jp[p]].mass_inv;
    }
    alpha = PI/(double)NPOLY;
    alpha2 = 2.0*alpha;
    ca = cos(alpha);
    ri = particle[i].radius*ca;
    rk = particle[k].radius*ca;
    
    p = 0;
    while ((p<particle[i].hash_nneighb)&&(!reaction))
    {
        k = particle[i].hashneighbour[p];
        nq = particle[k].npartners;
        molk = particle[k].molecule;
        
        if (np + nq > maxpartners) return(0); 
        
        if ((particle[k].active)&&(particle[k].type == type2)&&(nq <= maxpartners))
        {
            distance  = module2(particle[i].deltax[p], particle[i].deltay[p]);
            rel_angle  = argument(particle[i].deltax[p], particle[i].deltay[p]);
            poly_condition = (distance < REACTION_DIST*(ri+rk));
            
            rorient = particle[k].angle - particle[i].angle;
            if (NPOLY%2 == 1) rorient -= alpha;
            norient = (int)(rorient/alpha2 + 0.5);
            delta = rorient - norient*alpha2;
            poly_condition *= (vabs(delta) < DELTAMAX);
            
            if ((poly_condition)&&((double)rand()/RAND_MAX < reaction_prob))
            {
                /* find nearest orientation facing particle k */
                rorient_new = norient*alpha2;
                nangle = (int)((particle[i].angle - rel_angle)/alpha2);
                anglei_new = rel_angle + (double)(2*nangle+1)*alpha;
            
                m2 = 1.0/particle[k].mass_inv;
                for (q=0; q<nq; q++) 
                {
                    jq[q] = particle[k].partner[q];
                    m2 += 1.0/particle[jq[q]].mass_inv;
                }
                
                mr1 = m1/(m1 + m2);
                mr2 = 1.0 - mr1;   
                
                if ((np == 0)&&(nq == 0))
                {
                    reaction = 1;
                    printf("Merging molecule %i with particle %i\n", i, k);
                
                    distance = (ri+rk)*PAIR_DRATIO;
//                     particle[i].angle = rel_angle - PI/(double)NPOLY;
                    particle[i].angle = anglei_new;
                    particle[k].angle = particle[i].angle + rorient_new;
                    if (NPOLY%2 == 1) particle[k].angle += alpha;
                    
                    particle[k].xc = particle[i].xc + distance*cos(rel_angle);
                    particle[k].yc = particle[i].yc + distance*sin(rel_angle);

                    particle[i].npartners++;
                    particle[i].partner[0] = k;
                    particle[i].partner_eqd[0] = distance;  
                    particle[i].partner_eqa[0] = rorient_new;
                    if (NPOLY%2 == 1) particle[i].partner_eqa[0] += alpha;
                    
                    particle[k].npartners++;
                    particle[k].partner[0] = i;
                    particle[k].partner_eqd[0] = distance;
                    particle[k].partner_eqa[0] = -rorient_new;
                    if (NPOLY%2 == 1)  particle[k].partner_eqa[0] -= alpha;
                    
                    /* equalize speeds */
                    vxnew = mr1*particle[i].vx + mr2*particle[k].vx;
                    vynew = mr1*particle[i].vy + mr2*particle[k].vy;
                    
                    particle[i].vx = vxnew;
                    particle[i].vy = vynew;
                    particle[k].vx = vxnew;
                    particle[k].vy = vynew;
                    
                    if (mergeclusters)
                        merge_clusters(particle[i].cluster, particle[k].cluster, i, k, particle, cluster, 0, 0);
                }
                else if (np+nq+2 < NMAXPARTNERS)
                {
                    /* check if i and k belong to different molecules */
                    different = 1;
                    for (p1=0; p1<np; p1++) 
                        if (k == jp[p1]) different = 0;
                    for (q1=0; q1<nq; q1++) 
                        if (i == jq[q1]) different = 0;
                    
                    if (different)
                    {
                        deltav = module2(particle[i].vx - particle[k].vx, particle[i].vy - particle[k].vy);
//                         printf("Delta v is %.3lg\n", deltav);
                    
                        if (deltav < DELTAVMAX) 
                        {
                            reaction = 1;
                            printf("Pairing clusters containing particles %i and %i\n", i, k);
                            printf("np = %i, nq = %i\n", np, nq);
                            fprintf(lj_log, "Pairing clusters containing particles %i and %i\n", i, k);
                            fprintf(lj_log, "np = %i, nq = %i\n", np, nq);
                        
                            /* TODO: fix angles */
                            
                            distance = (ri+rk)*PAIR_DRATIO;
                            particle[i].angle = anglei_new;
                            for (p1=0; p1<np; p1++) particle[jp[p1]].angle = anglei_new;
                            particle[k].angle = particle[i].angle + rorient_new;
                            if (NPOLY%2 == 1) particle[k].angle += alpha;
                            
                            /* translate/rotate molecule containing particle k */
                            deltax = -particle[k].xc + particle[i].xc + distance*cos(rel_angle);
                            deltay = -particle[k].yc + particle[i].yc + distance*sin(rel_angle);
                            particle[k].xc += deltax;
                            particle[k].yc += deltay;
                            for (q1=0; q1<nq; q1++)
                            {
                                particle[jq[q1]].xc += deltax;
                                particle[jq[q1]].yc += deltay;
                                particle[jq[q1]].angle = particle[i].angle + rorient_new;
                                if (NPOLY%2 == 1) particle[jq[q1]].angle += alpha;
                            }
                            
                            particle[i].npartners++;
                            particle[i].partner[np] = k;
                            particle[i].partner_eqd[np] = distance;  
                            particle[i].partner_eqa[np] = rorient_new;
                            if (NPOLY%2 == 1) particle[i].partner_eqa[np] += alpha;
                            for (p1=0; p1<np; p1++) 
                            {
                                particle[jp[p1]].partner_eqa[p1] = rorient_new;
                                if (NPOLY%2 == 1) particle[i].partner_eqa[np] += alpha;
                            }
                    
                            particle[k].npartners++;
                            particle[k].partner[nq] = i;
                            particle[k].partner_eqd[nq] = distance;
                            particle[k].partner_eqa[nq] = -rorient_new;
                            if (NPOLY%2 == 1)  particle[k].partner_eqa[nq] -= alpha;

                            
                            /* equalize speeds */
                            vxnew = mr1*particle[i].vx + mr2*particle[k].vx;
                            vynew = mr1*particle[i].vy + mr2*particle[k].vy;
                    
                            particle[i].vx = vxnew;
                            particle[i].vy = vynew;
                            particle[k].vx = vxnew;
                            particle[k].vy = vynew;
                            
                            /* TODO: fix angles */
                            for (p1=0; p1<np; p1++) for (q1=0; q1<nq; q1++)
                            {
                                particle[jp[p1]].vx = vxnew;
                                particle[jp[p1]].vy = vynew;
                                particle[jq[q1]].vx = vxnew;
                                particle[jq[q1]].vy = vynew;
                            }
                            
                            if (mergeclusters)
                                merge_clusters(particle[i].cluster, particle[k].cluster, i, k, particle, cluster, 0, 0);
                        }
                    }
                }
                    
                if (reaction)
                {
                    collisions[ncollisions].x = particle[i].xc;
                    collisions[ncollisions].y = particle[i].yc;
                    collisions[ncollisions].time = COLLISION_TIME;
                    collisions[ncollisions].color = 0.0;
                    
                    if (ncollisions < 2*NMAXCOLLISIONS - 1) ncollisions++;
                    else printf("Too many collisions\n");
                }
            }
        }
        p++;
    }
    return(ncollisions);
}


void repair_cluster(int cl, t_particle particle[NMAXCIRCLES], t_cluster cluster[NMAXCIRCLES], int nsteps, int kill_overlays)
/* repair cluster by letting positions converge to equilibrium position */
{
    int t, j, k, p, q, np, nq, m;
    double dist, eqdist, fx, fy, dx, dy, r, phi, newphi, newx, newy, eqdist2;
    static double alpha, alpha2, ca, sa, kspring;
    static int first = 1;
    
    if (first)
    {
        alpha = PI/(double)NPOLY;
        alpha2 = 2.0*alpha; 
        ca = cos(alpha);
        sa = sin(alpha);
        kspring = 5000.0;
        first = 0;
    }
    
    for (t=0; t<nsteps; t++)
    {
        for (p=0; p<cluster[cl].nparticles; p++)
        {
            q = cluster[cl].particle[p];
            fx = 0.0;
            fy = 0.0;
            for (j=0; j<particle[q].hash_nneighb; j++)
            {
                k = particle[q].hashneighbour[j];
                if (particle[k].cluster == cl)
                {
                    dx = particle[q].deltax[j];
                    dy = particle[q].deltay[j];
                    dist = module2(dx, dy);
                    eqdist = (particle[k].radius + particle[q].radius)*ca;
                    
                    if (dist < 1.2*eqdist)
                    {
                        if ((NPOLY%2==0)||(particle[q].flip != particle[k].flip))
                        {
                            phi = argument(dx, dy); 
//                         printf("phi = %.3lg, angle = %.3lg\t", phi*180.0/PI, particle[q].angle*180.0/PI);
                        
                            m = (int)(0.5 + (double)NPOLY + (phi - particle[q].angle + alpha)/alpha2) - NPOLY;
                            newphi = particle[q].angle - alpha + (double)m*alpha2;
//                         printf("m = %i, newphi = %.3lg\n", m, newphi*180.0/PI);
                        
                            newx = particle[k].xc - eqdist*cos(newphi);
                            newy = particle[k].yc - eqdist*sin(newphi);
                        
//                         printf("(x,y) = (%.3lg, %.3lg) -> (%.3lg, %.3lg)\n", particle[q].xc, particle[q].yc, newx, newy); 
                            dx = newx - particle[q].xc;
                            dy = newy - particle[q].yc;
                            r = module2(dx, dy);
                        
//                         printf("Particle %i distance to target = %.3lg\n", k, r); 
//                         if ((t == nsteps-1)&&(q == 470)&&(k == 472)) fprintf(lj_log, "Particle %i relative distance to target %i = %.5lg\n", q, k, r/eqdist); 
                        
                            fx += r*dx;
                            fy += r*dy;
                        }
                        else
                        {
                            eqdist2 = 1.0*(particle[k].radius + particle[q].radius)*sa;
                            if (dist < eqdist2)
                            {
                                fx -= 0.1*(eqdist2 - dist)*dx/dist;
                                fy -= 0.1*(eqdist2 - dist)*dy/dist;
                            }
                        }
                    }
                    
                    if ((kill_overlays)&&(dist < REPAIR_MIN_DIST*eqdist))
                    {
                        nq = particle[q].npartners;
                        np = particle[k].npartners;
                        if (nq == 1) particle[q].active = 0;
                        if (np == 1) particle[k].active = 0;
//                         if (np > nq) particle[q].active = 0;
//                         else particle[k].active = 0;
                        
                    }
                }
            }
            particle[q].xc += kspring*fx*DT_PARTICLE;
            particle[q].yc += kspring*fy*DT_PARTICLE;
        }
    }
}


int chem_multi_glue_polygon2(int i, int maxpartners, t_particle particle[NMAXCIRCLES], t_molecule molecule[NMAXCIRCLES], t_cluster cluster[NMAXCIRCLES], t_collision *collisions, int ncollisions, double reaction_prob, int mergeclusters, int singlecluster)
/* simplified version of chem_multi_glue_molecule without triangle and similar */
/* conditions, but taking polygon geometry into account */
/* version 2, for CHEM_POLYGON_CLUSTER interaction */
{
    int k, p, q, p1, q1, p2, np, nq, n, m, closeby = 0, reaction = 0, jp[NMAXPARTNERS], jq[NMAXPARTNERS], jtot[NMAXPARTNERS], r, kk, different, np2, moli, molk, mp, newpartner, pp, type1, poly_condition, nangle, norient, cl1, cl2, cp, cq, newcluster, nmin, clmin, clmax, newreaction;
    double distance, angle, move_factor, xnew, ynew, vxnew, vynew, m1, m2, mr1, mr2, deltav, x, y, ri, rk, alpha, alpha2, ca, rel_angle, rorient, delta, rorient_new, anglei_new, deltax, deltay, x1, y1, aratio, delta2, eqdistance, dx, dy, dalpha1, dalpha2, ssize1, ssize2; 
    static short int first = 1;
    static int scluster;
    
    if ((singlecluster)&&(first))
    {
        scluster = 0;
        while (cluster[scluster].active == 0) scluster = rand()%ncircles;
        cluster[scluster].selected = 1;
        particle[scluster].collision = 1;
        printf("Selected cluster %i as single cluster\n", scluster);
        first = 0;
    }
    
    type1 = particle[i].type;
    np = particle[i].npartners;
    if (np > maxpartners) return(ncollisions);
    
    cl1 = particle[i].cluster;
    cp = cluster[cl1].nparticles;
    m1 = cluster[cl1].mass;

    alpha = PI/(double)NPOLY;
    alpha2 = 2.0*alpha;
    ca = cos(alpha);
    ri = particle[i].radius*ca;
    rk = particle[k].radius*ca;
    
    p = 0;
    while ((p<particle[i].hash_nneighb)&&(!reaction))
    {
        k = particle[i].hashneighbour[p];
        nq = particle[k].npartners;
        cl2 = particle[k].cluster;
        
        if (cluster[cl1].nparticles < cluster[cl2].nparticles) 
        {
            nmin = np;
            clmin = cluster[cl1].nparticles;
            clmax = cluster[cl2].nparticles;
        }
        else 
        {
            nmin = nq;
            clmin = cluster[cl2].nparticles;
            clmax = cluster[cl1].nparticles;
        }
            
//         if (np + nq + 2 > maxpartners) return(ncollisions); 
        
        if ((!reaction)&&(particle[k].active)&&(np+nq+1 < maxpartners)&&(cl2!=cl1)&&(nmin <= SMALL_NP_MAXSIZE)&&(clmin <= SMALL_CLUSTER_MAXSIZE))
        {
            cq = cluster[cl2].nparticles;
//             poly_condition = (cl1!=cl2);
            
            /* test distance */
            distance  = module2(particle[i].deltax[p], particle[i].deltay[p]);
            rel_angle  = argument(particle[i].deltax[p], particle[i].deltay[p]);
            poly_condition = (distance < REACTION_DIST*(ri+rk));
            
            /* test relative angle */
            rorient = particle[k].angle - particle[i].angle;
            if (NPOLY%2 == 1) rorient -= alpha;
            
            while (rorient > alpha) rorient -= alpha2;
            while (rorient < -alpha) rorient += alpha2;
            poly_condition *= (vabs(rorient) < DELTAMAX);
            
//             norient = (int)(rorient/alpha2 + 0.5 + 2.0*(double)NPOLY);
//             delta = rorient - (double)(norient-2*NPOLY)*alpha2;
//             poly_condition *= (vabs(delta) < DELTAMAX);
            
            /* test orientation */
            aratio = 0.5*vabs(rel_angle)/alpha;
            delta2 = aratio - (double)((int)aratio);
            poly_condition *= (vabs(delta2 - 0.5) < 0.5*DELTAMAX*alpha);
            
            /* option singlecluster: only allow merger if one of the clusters is selected cluster */
            if (singlecluster) poly_condition *= ((cluster[cl1].selected)||(cluster[cl2].selected)||(clmax <= NOTSELECTED_CLUSTER_MAXSIZE));
            
            newreaction = 0;
            
            if ((poly_condition)&&((double)rand()/RAND_MAX < reaction_prob))
            {
                /* find nearest orientation facing particle k */
                rorient_new = norient*alpha2;
                nangle = (int)((particle[i].angle - rel_angle)/alpha2);
                anglei_new = rel_angle + (double)(2*nangle+1)*alpha;
            
                m2 = cluster[cl2].mass;
                    
                mr1 = m1/(m1 + m2);
                mr2 = 1.0 - mr1;   
                
                if ((np == 0)&&(nq == 0))
                {
                    reaction = 1;
                    newreaction = 1;
                    printf("Merging molecule %i with particle %i\n", i, k);
                
                    eqdistance = ri + rk;
                    particle[i].angle = anglei_new;
                    particle[k].angle = particle[i].angle + rorient_new;
                    if (NPOLY%2 == 1) particle[k].angle += alpha;
                                        
                    dx = 0.5*(distance - eqdistance)*cos(rel_angle);
                    dy = 0.5*(distance - eqdistance)*sin(rel_angle);
                    
                    translate_cluster(cl1, cluster, particle, dx, dy);
                    translate_cluster(cl2, cluster, particle, -dx, -dy);

                    particle[i].npartners++;
                    particle[i].partner[0] = k;
                    
                    particle[k].npartners++;
                    particle[k].partner[0] = i;
                    
                    /* equalize speeds */
                    vxnew = mr1*particle[i].vx + mr2*particle[k].vx;
                    vynew = mr1*particle[i].vy + mr2*particle[k].vy;
                    
                    particle[i].vx = vxnew;
                    particle[i].vy = vynew;
                    particle[k].vx = vxnew;
                    particle[k].vy = vynew;
                    
                    if (mergeclusters)
                        newcluster = merge_clusters(particle[i].cluster, particle[k].cluster, i, k, particle, cluster, 0, (particle[i].flip == particle[k].flip));
                }
//                 else if ((np+nq+2 < NMAXPARTNERS)&&(vabs(delta) < 0.5*DELTAMAX))
                else if ((np+nq+2 < NMAXPARTNERS)&&(cp+cq <= CLUSTER_MAXSIZE))
                {
                    deltav = module2(particle[i].vx - particle[k].vx, particle[i].vy - particle[k].vy);
//                         printf("Delta v is %.3lg\n", deltav);
                    
                    if (deltav < DELTAVMAX) 
                    {
                        reaction = 1;
                        newreaction = 1;
                        printf("Pairing clusters containing particles %i and %i\n", i, k);
                        printf("np = %i, nq = %i\n", np, nq);
                        fprintf(lj_log, "Pairing clusters containing particles %i and %i\n", i, k);
                        fprintf(lj_log, "np = %i, nq = %i\n", np, nq);
                        
                        eqdistance = ri + rk;
                        
                        dx = 0.5*(distance - eqdistance)*cos(rel_angle);
                        dy = 0.5*(distance - eqdistance)*sin(rel_angle);
                        dalpha1 = anglei_new - particle[i].angle;
                        dalpha2 = anglei_new + rorient_new - particle[k].angle;
                        if (NPOLY%2 == 1) dalpha2 += alpha;
                        /* round to closest multiple of alpha2 */
                        while (dalpha1 > 0.5*DELTAMAX) dalpha1 -= alpha2;
                        while (dalpha1 < -0.5*DELTAMAX) dalpha1 += alpha2;
                        while (dalpha2 > 0.5*DELTAMAX) dalpha2 -= alpha2;
                        while (dalpha2 < -0.5*DELTAMAX) dalpha2 += alpha2;
                        
                        translate_cluster(cl1, cluster, particle, dx, dy);
                        translate_cluster(cl2, cluster, particle, -dx, -dy);

                        ssize1 = 0.02*sqrt((double)cluster[cl1].nparticles);
                        ssize2 = 0.02*sqrt((double)cluster[cl2].nparticles);
                        
                        if (vabs(dalpha1) + ssize1 < 0.5*DELTAMAX)
                            rotate_cluster_around_particle(cl1, i, cluster, particle, dalpha1);
                        else printf("Did not rotate cluster of size %.3lg\n", ssize1);
                        if (vabs(dalpha2) + ssize2 < 0.5*DELTAMAX)
                            rotate_cluster_around_particle(cl2, k, cluster, particle, dalpha2);
                        else printf("Did not rotate cluster of size %.3lg\n", ssize2);
                            
                        particle[i].npartners++;
                        particle[i].partner[np] = k;
                            
                        particle[k].npartners++;
                        particle[k].partner[nq] = i;
                            
                        /* equalize speeds */
                        vxnew = mr1*particle[i].vx + mr2*particle[k].vx;
                        vynew = mr1*particle[i].vy + mr2*particle[k].vy;
                    
                        for (p1=0; p1<cluster[cl1].nparticles; p1++)
                        {
                            q = cluster[cl1].particle[p1];
                            particle[q].vx = vxnew;
                            particle[q].vy = vynew;
                        }
                        for (p1=0; p1<cluster[cl2].nparticles; p1++)
                        {
                            q = cluster[cl2].particle[p1];
                            particle[q].vx = vxnew;
                            particle[q].vy = vynew;
                        }
                            
                        if (mergeclusters)
                        {
                            newcluster = merge_clusters(particle[i].cluster, particle[k].cluster, i, k, particle, cluster, 1, (particle[i].flip == particle[k].flip));
                        }
                    }
                }
                
                if ((newreaction)&&(singlecluster)&&((cluster[cl1].selected)||(cluster[cl2].selected)))
                {
                    scluster = newcluster;
                    printf("[chem_multi_glue_polygon2] cluster 1 = %i, cluster 2 = %i, newcluster = %i\n", cl1, cl2, newcluster);
                    cluster[cl1].selected = 1;
                    cluster[cl2].selected = 1;
//                     cluster[newcluster].selected = 1;
                    printf("[chem_multi_glue_polygon2] Selected cluster: %i\n", scluster);
                }
                
                if (reaction)
                {
                    collisions[ncollisions].x = particle[i].xc;
                    collisions[ncollisions].y = particle[i].yc;
                    collisions[ncollisions].time = COLLISION_TIME;
                    collisions[ncollisions].color = 0.0;
                    
                    if (ncollisions < 2*NMAXCOLLISIONS - 1) ncollisions++;
                    else printf("Too many collisions\n");
                }
            }
        }
//         else if ((nmin > SMALL_NP_MAXSIZE)||(clmin > SMALL_CLUSTER_MAXSIZE))
//         {
//             printf("No merger - Parameters: nmin = %i, clmin = %i\n", nmin, clmin);
//         }
        p++;
    }
    return(ncollisions);
}


int chem_split_molecule(int i, t_particle particle[NMAXCIRCLES], t_molecule molecule[NMAXCIRCLES], t_collision *collisions, int ncollisions)
/* split molecule containing particle i from all other molecules */
{
    int np, nmolp, mol, n, m, p, q, r, mm;
    int m_table[NMAXPARTNERMOLECULES];
    short int reaction = 0;
    
    np = particle[i].npartners;
    
    if (np==0) return(ncollisions);
    
    mol = particle[i].molecule;
    nmolp = molecule[mol].npartners; 
    
    printf("Molecule %i has %i partners: ", mol, nmolp); 
    fprintf(lj_log, "[chem_split_molecule] Molecule %i has %i partners: ", mol, nmolp); 
    for (mm=0; mm<nmolp; mm++) 
    {
        m_table[mm] = molecule[mol].partner[mm];
        printf("%i ", m_table[mm]);
        fprintf(lj_log, "%i ", m_table[mm]);
    }
    printf("\n");
    fprintf(lj_log, "\n"); 
    
    for (p=0; p<np; p++)
    {
        n = particle[i].partner[p];
        if (particle[n].molecule != mol)
        {
            dissociate_particles(i, n, p, particle);
            reaction = 1;
        }
    }
    
    np = particle[i].npartners;
    for (p=0; p<np; p++)
    {
        for (q=0; q<np; q++)
        {
            n = particle[i].partner[p];
            for (r=0; r<particle[n].npartners; r++)
            {
                m = particle[n].partner[r];
                if (particle[m].molecule != mol)
                {
                    dissociate_particles(n, m, r, particle);
                    reaction = 1;
                }
            }
        }
    }
       
    if (reaction)
    {
        update_single_molecule_data(mol, particle, molecule);
        for (mm=0; mm<nmolp; mm++) 
            update_single_molecule_data(m_table[mm], particle, molecule);
        
        printf("[chem_split_molecule] - dissociating molecule %i\n", mol);
        fprintf(lj_log, "[chem_split_molecule] - dissociating molecule %i\n", mol);
        collisions[ncollisions].x = particle[i].xc;
        collisions[ncollisions].y = particle[i].yc;
        collisions[ncollisions].time = COLLISION_TIME;
        collisions[ncollisions].color = type_hue(particle[i].type);
            
        if (ncollisions < 2*NMAXCOLLISIONS - 1) ncollisions++;
        else printf("Too many collisions\n");    
    }
    
    return(ncollisions);
}

int chem_split_molecule_if_nearby(int i, int type, t_particle particle[NMAXCIRCLES], t_molecule molecule[NMAXCIRCLES], t_collision *collisions, int ncollisions, double reaction_dist, double reaction_prob)
/* split molecules containing i if i is close to particle of type type */
{
    int p, k, mol;
    double distance;
    short int reaction = 0;
    
    mol = particle[i].molecule;
    p = 0;
    while ((p<particle[i].hash_nneighb)&&(!reaction))
    {
        k = particle[i].hashneighbour[p];
        if ((particle[k].type == type)&&(particle[k].molecule != mol))
        {
            distance = module2(particle[i].xc - particle[k].xc, particle[i].yc - particle[k].yc);
            if ((distance < reaction_dist)&&((double)rand()/RAND_MAX < reaction_prob))
            {
                chem_split_molecule(i, particle, molecule, collisions, ncollisions);
                reaction = 1;
            }
        }
        p++;
    }
    
    if (reaction)
    {
        printf("dissociating molecule %i\n", mol);
        collisions[ncollisions].x = particle[i].xc;
        collisions[ncollisions].y = particle[i].yc;
        collisions[ncollisions].time = COLLISION_TIME;
        collisions[ncollisions].color = type_hue(particle[i].type);
            
        if (ncollisions < 2*NMAXCOLLISIONS - 1) ncollisions++;
        else printf("Too many collisions\n");    
    }
    
    return(ncollisions);
}

int chem_local_split_molecule_if_nearby(int i, int type, int enzyme_type, t_particle particle[NMAXCIRCLES], t_molecule molecule[NMAXCIRCLES], t_collision *collisions, int ncollisions, double reaction_dist, double reaction_prob)
/* split atoms in molecule containing i from type type, if i is close to particle of type enzyme_type */
{
    int p, k, moli, molk, mol_diss, q, n, q1, m, added, mp, pp;
    double distance;
    short int reaction = 0, reaction1 = 0;
    
    moli = particle[i].molecule;
    p = 0;
    while ((p<particle[i].hash_nneighb)&&(!reaction))
    {
        reaction1 = 0;
        k = particle[i].hashneighbour[p];
        molk = particle[k].molecule;
        
        added = particle[i].added;
//         added = 0;
        if ((added == 0)&&(particle[k].type == enzyme_type)&&(molk != moli))
        {
            distance = module2(particle[i].xc - particle[k].xc, particle[i].yc - particle[k].yc);
            if ((distance < reaction_dist)&&((double)rand()/RAND_MAX < reaction_prob))
            {
                if (particle[i].coulomb < 5) particle[i].coulomb = 5;
                particle[i].reactive = 0;
                for (q=0; q<particle[i].npartners; q++)
                {
                    n = particle[i].partner[q];
                    if (particle[n].added == 0)
                    {
                        if (particle[n].type == type)
                        {
                            mol_diss = particle[n].molecule;
                            dissociate_particles(i, n, q, particle);
                            if (particle[n].coulomb < 5) particle[n].coulomb = 5;
                            reaction = 1;
                            reaction1 = 1;
                        }
                        else if (particle[n].type == particle[i].type)
                        {
                            for (q1=0; q1<particle[n].npartners; q1++)
                            {
                                m = particle[n].partner[q1];
                                if ((particle[m].type == type)&&(particle[m].added == 0))
                                {
                                    mol_diss = particle[m].molecule;
                                    dissociate_particles(n, m, q1, particle);
                                    if (particle[n].coulomb < 5) particle[n].coulomb = 5;
                                    if (particle[m].coulomb < 5) particle[m].coulomb = 5;
                                    particle[n].reactive = 0;
                                    particle[m].reactive = 0;
                                    reaction = 1;
                                    reaction1 = 1;
                                }
                            }
                        }
                    }
                }
                
                /* TEST */
                if ((RD_REACTION == CHEM_DNA_ENZYME)||(RD_REACTION == CHEM_DNA_ENZYME_REPAIR))
                {
                    if ((particle[i].type >= 3)&&(particle[i].type <= 6))
                    {
                        particle[i].paired = 0;
                        for (q1 = 0; q1 <= particle[i].npartners; q1++)
                        {
                            n = particle[i].partner[q1];
                            particle[n].paired = 0;
                        }
                        
                        particle[k].paired = 0;
                        for (q1 = 0; q1 <= particle[k].npartners; q1++)
                        {
                            n = particle[k].partner[q1];
                            particle[n].paired = 0;
                        }
                    }
                }   
            }
            
            /* update molecule data */
            if ((reaction1)&&(moli != mol_diss))
            {
                printf("Splitting molecules %i and %i\n", moli, mol_diss);
                
                mp = molecule[moli].npartners;
                pp = 0;
                while (molecule[moli].partner[pp] != mol_diss) pp++;
                if (pp < mp)
                {
                    while (pp < mp-1) 
                    {
                        molecule[moli].partner[pp] = molecule[moli].partner[pp+1];
                        molecule[moli].connection_type[pp] = molecule[moli].connection_type[pp+1];
                        pp++;
                    }
                    molecule[moli].npartners--;
                }
                
                mp = molecule[mol_diss].npartners;
                pp = 0;
                while (molecule[mol_diss].partner[pp] != moli) pp++;
                if (pp < mp)
                {
                    while (pp < mp-1) 
                    {
                        molecule[mol_diss].partner[pp] = molecule[mol_diss].partner[pp+1];
                        molecule[mol_diss].connection_type[pp] = molecule[mol_diss].connection_type[pp+1];
                        pp++;
                    }
                    molecule[mol_diss].npartners--;
                }
//                 sleep(3);
            }
        }
        p++;
    }
    
    if (reaction)
    {
        printf("dissociating molecule %i\n", moli);
        collisions[ncollisions].x = particle[i].xc;
        collisions[ncollisions].y = particle[i].yc;
        collisions[ncollisions].time = COLLISION_TIME;
        collisions[ncollisions].color = type_hue(particle[i].type);
            
        if (ncollisions < 2*NMAXCOLLISIONS - 1) ncollisions++;
        else printf("Too many collisions\n");    
    }
    
    return(ncollisions);
}

void change_type_proportion(t_particle particle[NMAXCIRCLES], double prop)
/* change proportion of particles of types 1 or 2 */
{
    int i, n0=0, n1=0, n2=0, ntot, nc=0, nmod, counter=0, cmax=1000;
    
    for (i=0; i<ncircles; i++) if (particle[i].active)
    {
        if (particle[i].type == 0) n0++;
        if (particle[i].type == 1) n1++;
        if (particle[i].type == 2) n2++;
    }
    
    ntot = n0 + n1 + n2;
    n0 += n1;
    
    if ((double)(n0+1) < prop*(double)ntot)
    {
        nmod = (int)((double)ntot*prop - (double)n0);
        if (nmod > 0) while ((nc<nmod)&&(counter<cmax))
        {   
            i = rand()%ncircles;
            if ((particle[i].type == 2)&&(particle[i].active))
            {
                particle[i].type = 1;
                particle[i].radius = MU;
                particle[i].interaction = INTERACTION;
                particle[i].eq_dist = EQUILIBRIUM_DIST;
                particle[i].spin_range = SPIN_RANGE;
                particle[i].spin_freq = SPIN_INTER_FREQUENCY;
                particle[i].mass_inv = 1.0/PARTICLE_MASS;
                particle[i].inertia_moment_inv = 1.0/PARTICLE_INERTIA_MOMENT;
                particle[i].charge = CHARGE;
                nc++;
            }
            counter++;
        }
    }
    else if ((double)(n0-1) > prop*(double)ntot)
    {
        nmod = (int)((double)n0 - (double)ntot*prop);
        if (nmod > 0) while ((nc<nmod)&&(counter<cmax))
        {   
            i = rand()%ncircles;
            if ((particle[i].type <= 1)&&(particle[i].active)) 
            {
                particle[i].type = 2;
                particle[i].radius = MU_B;
                particle[i].interaction = INTERACTION_B;
                particle[i].eq_dist = EQUILIBRIUM_DIST_B;
                particle[i].spin_range = SPIN_RANGE_B;
                particle[i].spin_freq = SPIN_INTER_FREQUENCY_B;
                particle[i].mass_inv = 1.0/PARTICLE_MASS_B;
                particle[i].inertia_moment_inv = 1.0/PARTICLE_INERTIA_MOMENT_B;
                particle[i].charge = CHARGE_B;
                nc++;
            }
            counter++;
        }
    }
}

int find_partners(t_particle *particle, int position, t_particle* *pstack, int *stacksize, int cluster)
/* look for active partners of particle[position], returns difference in partners */
{
    int k, nopen = 0, n;
    
    if (!pstack[position]->active) return(-1);
    
    pstack[position]->cluster = cluster;
    
    for (k=0; k<pstack[position]->npartners; k++)
    {
        n = pstack[position]->partner[k];
        if ((particle[n].active)&&(particle[n].cluster != cluster))
        {
            particle[n].tested = 1;
            particle[n].cluster = cluster;
            particle[n].cactive = 1;
            nopen++;
            if (*stacksize < ncircles)
            {
                (*stacksize)++;
                pstack[*stacksize-1] = &particle[n];
            }
        }
    }

    if (nopen == 0) 
    {
        pstack[position]->cactive = 0;
        return(-1);
    }
    else return(nopen);
}

int update_cluster_color(t_particle particle[NMAXCIRCLES])
/* update the colors of clusters */
{
    int i, k, p, position, nactive, stacksize, nclusters = 0;
    t_particle **cstack;
    
    cstack = (t_particle* *)malloc(ncircles*sizeof(struct t_particle *));
    
    #pragma omp parallel for private(i)
    for (i=0; i<ncircles; i++) particle[i].tested = 0;
    
    for (i=0; i<ncircles; i++) if ((particle[i].active)&&(!particle[i].tested))
    {
        position = 0;
        stacksize = 1;
        nactive = 1;
        cstack[0] = &particle[i];
        
        particle[i].cactive = 1;
        particle[i].tested = 1;
        
        /* do a depth first search starting in i */
        while (nactive > 0) 
        {
            while (!cstack[position]->cactive) position++;
            if (position == stacksize) position = 0;
            
            /* find an active partner in stack */
            nactive += find_partners(particle, position, cstack, &stacksize, particle[i].cluster);
        }
        nclusters++;
    }
    
    free(cstack);
    return(nclusters);
}


int repair_dna(t_particle particle[NMAXCIRCLES], t_molecule molecule[NMAXCIRCLES], t_collision *collisions, int ncollisions)
/* repair unwanted connections in DNA, for CHEM_DNA_ENZYME_REPAIR reaction type */
{
    int i, mol, moli, j, molj, p, q, p1, q1, molp, molpp, molq, molqq, r, type, type1, type2, type3, type4, np, nq, pp, qq, mp, split, common_partner, smol1, smol2, npartners;
    int p_table[NMAXPARTNERMOLECULES], q_table[NMAXPARTNERMOLECULES];
    double dist;
    
    printf("Repairing DNA\n"); 
    /* Prevent single base from connecting to two different molecules */
    for (moli=0; moli<nmolecules; moli++) if (molecule[moli].added)
    {
        split = 0;
        /* search for double connections */
        for (j=0; j<molecule[moli].npartners; j++)
        {
            molj = molecule[moli].partner[j];
            type1 = molecule[moli].connection_type[j];
            if (type1 >= 3)
            {
                for (p=0; p<molecule[moli].npartners; p++)
                {
                    molp = molecule[moli].partner[p];
                    type2 = molecule[moli].connection_type[p];
                    if ((type2 >= 3)&&(molp != molj)) 
                    {
                        smol1 = molj;
                        smol2 = molp;
                        split = 1;
                    }
                }
            }
        }
        if (split)
        {
            printf("Splitting molecule %i from molecules %i and %i\n\n\n", moli, smol1, smol2);
                                        
                    for (p=0; p<molecule[moli].nparticles; p++)
                    {
                        p1 = molecule[moli].particle[p];
                        type = particle[p1].type;
                        if (type >= 3)
                        {
                            i = p1;
                            for (q=0; q<particle[p1].npartners; q++)
                            {
                                q1 = particle[p1].partner[q];
                                if (particle[q1].molecule == molj)
                                    dissociate_particles(p1, q1, q, particle);
                            }
                        }
                    }
                    for (q=0; q<molecule[smol1].nparticles; q++)
                    {
                        q1 = molecule[smol1].particle[q];
                        type = particle[q1].type;
                        if (type >= 3)
                            for (p=0; p<particle[q1].npartners; p++)
                            {
                                p1 = particle[q1].partner[p];
                                if (particle[p1].molecule == moli)
                                    dissociate_particles(q1, p1, p, particle);
                            }
                    }
                    for (q=0; q<molecule[smol2].nparticles; q++)
                    {
                        q1 = molecule[smol2].particle[q];
                        type = particle[q1].type;
                        if (type >= 3)
                            for (p=0; p<particle[q1].npartners; p++)
                            {
                                p1 = particle[q1].partner[p];
                                if (particle[p1].molecule == moli)
                                    dissociate_particles(q1, p1, p, particle);
                            }
                    }
                    
                    /* update molecule data */
                    mp = molecule[moli].npartners;
                    pp = 0;
                    while (molecule[moli].partner[pp] != smol1) pp++;
                    if (pp < mp)
                    {
                        while (pp < mp-1) 
                        {
                            molecule[moli].partner[pp] = molecule[moli].partner[pp+1];
                            molecule[moli].connection_type[pp] = molecule[moli].connection_type[pp+1];
                            pp++;
                        }
                        molecule[moli].npartners--;
                    }
                    pp = 0;
                    while (molecule[moli].partner[pp] != smol2) pp++;
                    if (pp < mp)
                    {
                        while (pp < mp-1) 
                        {
                            molecule[moli].partner[pp] = molecule[moli].partner[pp+1];
                            molecule[moli].connection_type[pp] = molecule[moli].connection_type[pp+1];
                            pp++;
                        }
                        molecule[moli].npartners--;
                    }
                    
                    mp = molecule[smol1].npartners;
                    pp = 0;
                    while (molecule[smol1].partner[pp] != moli) pp++;
                    if (pp < mp)
                    {
                        while (pp < mp-1) 
                        {
                            molecule[smol1].partner[pp] = molecule[smol1].partner[pp+1];
                            molecule[smol1].connection_type[pp] = molecule[smol1].connection_type[pp+1];
                            pp++;
                        }
                        molecule[smol1].npartners--;
                    }
                    mp = molecule[smol2].npartners;
                    pp = 0;
                    while (molecule[smol2].partner[pp] != moli) pp++;
                    if (pp < mp)
                    {
                        while (pp < mp-1) 
                        {
                            molecule[smol2].partner[pp] = molecule[smol2].partner[pp+1];
                            molecule[smol2].connection_type[pp] = molecule[smol2].connection_type[pp+1];
                            pp++;
                        }
                        molecule[smol2].npartners--;
                    }
            
                    printf("dissociating molecule %i\n", moli);
                    collisions[ncollisions].x = particle[i].xc;
                    collisions[ncollisions].y = particle[i].yc;
                    collisions[ncollisions].time = COLLISION_TIME;
                    collisions[ncollisions].color = type_hue(particle[i].type);
            
                    if (ncollisions < 2*NMAXCOLLISIONS - 1) ncollisions++;
                    else printf("Too many collisions\n"); 
        }
    }
    
    /* Test for skipped bases */
    /* if moli and molj are paired molecules, check whether they form a "square" */
    /* i.e. a partner of a partner of moli is a partner of molj */
    /* taking the connection types into account */
    for (moli=0; moli<nmolecules; moli++)  if (molecule[moli].added)
        for (j=0; j<molecule[moli].npartners; j++) 
        {
            molj = molecule[moli].partner[j];
            type1 = molecule[moli].connection_type[j];
            type2 = 3 - type1;
            if ((molecule[molj].added)&&((type1 == 1)||(type1 == 2)))   /* i-j connection is a backbone */
            {
                /* find base-neighbours for moli (normally there should be at most one) */
                np = 0;
                for (p=0; p<molecule[moli].npartners; p++)
                {
                    molp = molecule[moli].partner[p];
                    type3 = molecule[moli].connection_type[p];
                    if (type3 >= 3)
                    {
                        p_table[np] = molp;
                        np++;
                    }
                }
                
                /* find base-neighbours for molj (normally there should be at most one) */
                nq = 0;
                for (q=0; q<molecule[molj].npartners; q++)
                {
                    molq = molecule[molj].partner[q];
                    type3 = molecule[molj].connection_type[q];
//                     if ((molq != molj)&&(type3 >= 3))
                    if (type3 >= 3)
                    {
                        q_table[nq] = molq;
                        nq++;
                    }
                }
                
                /* test for split, split occurs if each molecule has a base-neighbour, 
                 and they are not neighbours */
                if ((np == 0)||(nq == 0)) split = 0;
                else
                {
                    split = 1;
                    
                    printf("Base partners of molecule %i(%i): ", moli, type1);
                    for (p=0; p<np; p++) printf("%i ", p_table[p]);
                    printf("\n");
                    printf("Base partners of molecule %i(%i): ", molj, type2);
                    for (q=0; q<nq; q++) printf("%i ", q_table[q]);
                    printf("\n");
                
                    for (p=0; p<np; p++)
                    {
                        molp = p_table[p];
                        for (q=0; q<nq; q++)
                        {
                            molq = q_table[q];
                            
                            printf("Checking base partners %i and %i\n", molp, molq);
                            printf("type1 = %i, type2 = %i\n", type1, type2);
                            printf("Molecule %i has %i partners: ", molp, molecule[molp].npartners);
                            for (pp = 0; pp < molecule[molp].npartners; pp++) 
                                printf("%i(%i) ", molecule[molp].partner[pp], molecule[molp].connection_type[pp]);
                            printf("\n");
                            printf("Molecule %i has %i partners: ", molq, molecule[molq].npartners);
                            for (pp = 0; pp < molecule[molq].npartners; pp++) 
                                printf("%i(%i) ", molecule[molq].partner[pp], molecule[molq].connection_type[pp]);
                            printf("\n");
                            
                            for (pp = 0; pp < molecule[molp].npartners; pp++)
                            {
                                molpp = molecule[molp].partner[pp];
                                type4 = molecule[molp].connection_type[pp];
                                printf("Testing molecule %i(%i)\n", molpp, type4);
                                if ((type4 == type2)&&(molpp == molq))
                                    split = 0;
                            }
                        }
                    }
                }
                
                if (split)      /* break backbone link between i and j */ 
                {
                    printf("Splitting molecules %i and %i\n\n\n", moli, molj);
                                        
                    for (p=0; p<molecule[moli].nparticles; p++)
                    {
                        p1 = molecule[moli].particle[p];
                        type = particle[p1].type;
                        if (type == type1)
                            i = p1;
                            for (q=0; q<particle[p1].npartners; q++)
                            {
                                q1 = particle[p1].partner[q];
                                if (particle[q1].molecule == molj)
                                    dissociate_particles(p1, q1, q, particle);
                            }
                    }
                    for (q=0; q<molecule[molj].nparticles; q++)
                    {
                        q1 = molecule[molj].particle[q];
                        type = particle[q1].type;
                        if (type == type2)
                            for (p=0; p<particle[q1].npartners; p++)
                            {
                                p1 = particle[q1].partner[p];
                                if (particle[p1].molecule == moli)
                                    dissociate_particles(q1, p1, p, particle);
                            }
                    }
                    
                    
//                         printf("Splitting molecules %i and %i\n", moli, mol_diss);
                    /* update molecule data */
                    mp = molecule[moli].npartners;
                    pp = 0;
                    while (molecule[moli].partner[pp] != molj) pp++;
                    if (pp < mp)
                    {
                        while (pp < mp-1) 
                        {
                            molecule[moli].partner[pp] = molecule[moli].partner[pp+1];
                            molecule[moli].connection_type[pp] = molecule[moli].connection_type[pp+1];
                            pp++;
                        }
                        molecule[moli].npartners--;
                    }
                
                    mp = molecule[molj].npartners;
                    pp = 0;
                    while (molecule[molj].partner[pp] != moli) pp++;
                    if (pp < mp)
                    {
                        while (pp < mp-1) 
                        {
                            molecule[molj].partner[pp] = molecule[molj].partner[pp+1];
                            molecule[molj].connection_type[pp] = molecule[molj].connection_type[pp+1];
                            pp++;
                        }
                        molecule[molj].npartners--;
                    }
            
                    printf("dissociating molecule %i\n", moli);
                    collisions[ncollisions].x = particle[i].xc;
                    collisions[ncollisions].y = particle[i].yc;
                    collisions[ncollisions].time = COLLISION_TIME;
                    collisions[ncollisions].color = type_hue(particle[i].type);
            
                    if (ncollisions < 2*NMAXCOLLISIONS - 1) ncollisions++;
                    else printf("Too many collisions\n");    
                }
            }
        }
        
        return(ncollisions);
}
     
int check_dna_pairing(t_particle particle[NMAXCIRCLES], t_molecule molecule[NMAXCIRCLES], t_collision *collisions, int ncollisions)
/* repair unwanted connections in DNA, for CHEM_DNA_ENZYME_REPAIR reaction type */
{
    int i, moli, j, type1, molj, p, pp, q, qq, deltatype;
    short int split;
    double dist;
    
    printf("Checking DNA pairing\n"); 
//     fprintf(lj_log, "Checking DNA pairing\n"); 
    /* Test for wrong connections */
    /* On rare occasions, it happens that a base pair tries to connect, but the connection fails */
    /* However, the nucleotides are still recorded as connected */
    /* This part tests for such falsely recorded connections based on the distance, and removes them */
    for (moli=0; moli<nmolecules; moli++) if (molecule[moli].added)
    {
        for (j=0; j<molecule[moli].npartners; j++)
        {
            split = 0;
            molj = molecule[moli].partner[j];
            type1 = molecule[moli].connection_type[j];
            
            /* search for wrongly recorded connections */
            for (p=0; p<molecule[moli].nparticles; p++)
            {
                pp = molecule[moli].particle[p];
                if (particle[pp].type == type1)
                {
                    for (q=0; q<particle[p].npartners; q++)
                    {
                        qq = particle[p].partner[q];
                        deltatype = particle[qq].type - type1;
                        if (deltatype < 0) deltatype *= -1;
//                         if ((particle[qq].molecule == molj)&&(particle[qq].type >= 1))
                        if ((particle[qq].molecule == molj)&&(deltatype == 1))
                        {
                            dist = module2(particle[pp].xc - particle[qq].xc, particle[pp].yc - particle[qq].yc);
//                             if (dist > 10.0*MU) 
//                             if (dist > 15.0*MU) 
                            if (dist > 20.0*MU) 
                            {
                                split = 1;
                                i = pp;
                                fprintf(lj_log, "\n\n\n Time = %i\n", frame_time);  
                                fprintf(lj_log, "Particles %i(%i) in molecule %i and %i(%i) in molecule %i are at distance %.3lg\n", pp, particle[pp].type, moli, qq, particle[qq].type, molj, dist);
                            }
                        }
                    }
                }
            }
        
            if (split)
            {
                printf("Splitting molecule %i from molecule %i because of wrong connection\n", moli, molj);
                fprintf(lj_log, "Splitting molecule %i from molecule %i because of wrong connection\n", moli, molj);
                
                ncollisions = chem_split_molecule(i, particle, molecule, collisions, ncollisions);
            }
        }
    }
    
    return(ncollisions);
}


int update_types(t_particle particle[NMAXCIRCLES], t_molecule molecule[NMAXCIRCLES], t_cluster cluster[NMAXCIRCLES], t_collision *collisions, int ncollisions, int *particle_numbers, int time, double *delta_e)
/* update the types in case of reaction-diffusion equation */
{
    int i, j, k, n, n3, n4, p, type, atype, btype, oldncollisions, delta_n;
    short int reacted;
    double distance, rnd, p1, p2, reac_dist;
    static double inv_masses[RD_TYPES+1], radii[RD_TYPES+1];
    static int first = 1, repair_counter = 0;
    
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
        case (CHEM_H2O_H_OH):
        {
            for (i=0; i<ncircles; i++) 
            {
                if ((particle[i].active)&&(particle[i].type <= 1))
                {
                    if (particle[i].npartners == 2)
                        ncollisions = chem_dissociate_molecule(i, particle, collisions, ncollisions, DISSOCIATION_PROB);
                    else if (particle[i].npartners == 1)
                        ncollisions = chem_merge_molecule(i, 2, 0, particle, collisions, ncollisions, REACTION_PROB);
                }
            }
            printf("%i collisions\n", ncollisions);
            delta_n = ncollisions - oldncollisions; 
            printf("delta_n = %i\n", delta_n);
            return(ncollisions);
        }
        case (CHEM_2H2O_H3O_OH):
        {
            for (i=0; i<ncircles; i++) 
            {
                if ((particle[i].active)&&(particle[i].type <= 1))
                {
                    if (particle[i].npartners == 2)
                    {
                        oldncollisions = ncollisions;
                        ncollisions = chem_dissociate_molecule(i, particle, collisions, ncollisions, DISSOCIATION_PROB);
                        if (ncollisions == oldncollisions)
                            ncollisions = chem_merge_molecule(i, 2, 0, particle, collisions, ncollisions, REACTION_PROB);                       
                    }
                    else if (particle[i].npartners == 1)
                        ncollisions = chem_merge_molecule(i, 2, 0, particle, collisions, ncollisions, REACTION_PROB);
                    else if (particle[i].npartners == 3)
                        ncollisions = chem_dissociate_molecule(i, particle, collisions, ncollisions, 2.0*DISSOCIATION_PROB);
                }
            }
            printf("%i collisions\n", ncollisions);
            delta_n = ncollisions - oldncollisions; 
            printf("delta_n = %i\n", delta_n);
            return(ncollisions);
        }
        case (CHEM_AGGREGATION):
        {
            for (i=0; i<ncircles; i++) if (particle[i].active)
            {
                if (particle[i].npartners < AGREGMAX) for (k=0; k<3; k++)
                {
                    ncollisions = chem_multi_glue_molecule(i, k, AGREGMAX, 0, 0, 0, 0, particle, molecule, collisions, ncollisions, REACTION_PROB);
                }
                
                /* decouple particles with several partners from thermostat */
                if (particle[i].npartners >= AGREG_DECOUPLE) particle[i].thermostat = 0;
            }
            
            /* update cluster color scheme */
            if ((PLOT == P_CLUSTER)||(PLOT_B == P_CLUSTER))
                update_cluster_color(particle);
            
            printf("%i collisions\n", ncollisions);
            delta_n = ncollisions - oldncollisions; 
            printf("delta_n = %i\n", delta_n);
            return(ncollisions);
        }
        case (CHEM_AGGREGATION_CHARGE):
        {
            for (i=0; i<ncircles; i++) if (particle[i].active)
            {
                if (particle[i].npartners < AGREGMAX) for (k=0; k<3; k++)
                {
                    ncollisions = chem_multi_glue_molecule(i, k, AGREGMAX, 0, 0, 1, 0, particle, molecule, collisions, ncollisions, REACTION_PROB);
                }
                
                /* decouple particles with several partners from thermostat */
                if (particle[i].npartners >= AGREG_DECOUPLE) particle[i].thermostat = 0;
            }
            
            /* update cluster color scheme */
            if ((PLOT == P_CLUSTER)||(PLOT_B == P_CLUSTER))
                update_cluster_color(particle);
            
            printf("%i collisions\n", ncollisions);
            delta_n = ncollisions - oldncollisions; 
            printf("delta_n = %i\n", delta_n);
            return(ncollisions);
        }
        case (CHEM_AGGREGATION_NNEIGH):
        {
            for (i=0; i<ncircles; i++) if (particle[i].active)
            {
                if (particle[i].npartners < AGREGMAX) for (k=0; k<3; k++)
                {
                    ncollisions = chem_multi_glue_molecule(i, k, AGREGMAX, 1, 0, 1, 0, particle, molecule, collisions, ncollisions, REACTION_PROB);
                }
                
                /* decouple particles with several partners from thermostat */
                if (particle[i].npartners >= AGREG_DECOUPLE) particle[i].thermostat = 0;
            }
            
            /* update cluster color scheme */
            if ((PLOT == P_CLUSTER)||(PLOT_B == P_CLUSTER))
                update_cluster_color(particle);
            
            printf("%i collisions\n", ncollisions);
            delta_n = ncollisions - oldncollisions; 
            printf("delta_n = %i\n", delta_n);
            return(ncollisions);
        }
        case (CHEM_DNA):
        {
            for (i=0; i<ncircles; i++) if (particle[i].active)
            {
                if (particle[i].npartners < AGREGMAX) switch (particle[i].type) 
                {
                    case (2): 
                    {
                        ncollisions = chem_multi_glue_molecule(i, 2, AGREGMAX, 0, 0, 0, 0, particle, molecule, collisions, ncollisions, REACTION_PROB);
                        break;
                    }
                    case (3): 
                    {
                        ncollisions = chem_multi_glue_molecule(i, 4, AGREGMAX, 0, 0, 0, 0, particle, molecule, collisions, ncollisions, REACTION_PROB);
                        break;
                    }
                    case (4): 
                    {
                        ncollisions = chem_multi_glue_molecule(i, 3, AGREGMAX, 0, 0, 0, 0, particle, molecule, collisions, ncollisions, REACTION_PROB);
                        break;
                    }
                    case (5): 
                    {
                        ncollisions = chem_multi_glue_molecule(i, 6, AGREGMAX, 0, 0, 0, 0, particle, molecule, collisions, ncollisions, REACTION_PROB);
                        break;
                    }
                    case (6): 
                    {
                        ncollisions = chem_multi_glue_molecule(i, 5, AGREGMAX, 0, 0, 0, 0, particle, molecule, collisions, ncollisions, REACTION_PROB);
                        break;
                    }
                }
                
                /* decouple particles with several partners from thermostat */
                if (particle[i].npartners >= AGREG_DECOUPLE) particle[i].thermostat = 0;
            }
            
            /* update cluster color scheme */
            if ((PLOT == P_CLUSTER)||(PLOT_B == P_CLUSTER))
                update_cluster_color(particle);
            
            printf("%i collisions\n", ncollisions);
            delta_n = ncollisions - oldncollisions; 
            printf("delta_n = %i\n", delta_n);
            return(ncollisions);
        }
        case (CHEM_DNA_ALT):
        {
            for (i=0; i<ncircles; i++) if (particle[i].active)
            {
                if (particle[i].npartners < AGREGMAX) 
                {
                    atype = particle[i].type;
                    if (atype%2 == 1) btype = atype+1;
                    else btype = atype-1;
                    if (atype > 0)
                        ncollisions = chem_multi_glue_molecule(i, btype, AGREGMAX, 0, 0, 0, 0, particle, molecule, collisions, ncollisions, REACTION_PROB);
                }
                
                /* decouple particles with several partners from thermostat */
                if (particle[i].npartners >= AGREG_DECOUPLE) particle[i].thermostat = 0;
            }
            
            /* update cluster color scheme */
//             if ((PLOT == P_CLUSTER)||(PLOT_B == P_CLUSTER))
            update_cluster_color(particle);
            
            printf("%i collisions\n", ncollisions);
            delta_n = ncollisions - oldncollisions; 
            printf("delta_n = %i\n", delta_n);
            return(ncollisions);
        }
        case (CHEM_DNA_DOUBLE):
        {
            for (i=0; i<ncircles; i++) if (particle[i].active)
            {
                if (particle[i].npartners < AGREGMAX) 
                {
                    atype = particle[i].type;
                    if (atype%2 == 1) btype = atype+1;
                    else btype = atype-1;
                    if (atype > 0)
                        ncollisions = chem_multi_glue_molecule(i, btype, AGREGMAX, 0, 0, 0, 0, particle, molecule, collisions, ncollisions, REACTION_PROB);
                }
                
                /* decouple particles with several partners from thermostat */
                if (particle[i].npartners >= AGREG_DECOUPLE) particle[i].thermostat = 0;
            }
            
            /* update cluster color scheme */
            if ((PLOT == P_CLUSTER)||(PLOT_B == P_CLUSTER))
                update_cluster_color(particle);
            
            printf("%i collisions\n", ncollisions);
            delta_n = ncollisions - oldncollisions; 
            printf("delta_n = %i\n", delta_n);
            return(ncollisions);
        }
        case (CHEM_DNA_DSPLIT):
        {
            for (i=0; i<ncircles; i++) if (particle[i].active)
            {
                if (particle[i].npartners < AGREGMAX) 
                {
                    atype = particle[i].type;
                    if (atype%2 == 1) btype = atype+1;
                    else btype = atype-1;
                    if (atype > 0)
                        ncollisions = chem_multi_glue_molecule(i, btype, AGREGMAX, 3, 0, 0, 0, particle, molecule, collisions, ncollisions, REACTION_PROB);
                }
                
                /* split molecule with a certain probability */
                if ((double)rand()/RAND_MAX < DISSOCIATION_PROB)
                    ncollisions = chem_split_molecule(i, particle, molecule, collisions, ncollisions);
                
                /* decouple particles with several partners from thermostat */
                if (particle[i].npartners >= AGREG_DECOUPLE) particle[i].thermostat = 0;
            }
            
            /* update cluster color scheme */
            if ((PLOT == P_CLUSTER)||(PLOT_B == P_CLUSTER))
                update_cluster_color(particle);
            
            printf("%i collisions\n", ncollisions);
            delta_n = ncollisions - oldncollisions; 
            printf("delta_n = %i\n", delta_n);
            return(ncollisions);
        }
        case (CHEM_DNA_BASE_SPLIT):
        {
            for (i=0; i<ncircles; i++) if (particle[i].active)
            {
                if (particle[i].npartners < AGREGMAX) 
                {
                    atype = particle[i].type;
                    if (atype%2 == 1) btype = atype+1;
                    else btype = atype-1;
                    if (atype > 0)
                        ncollisions = chem_multi_glue_molecule(i, btype, AGREGMAX, 3, 0, 0, 0, particle, molecule, collisions, ncollisions, REACTION_PROB);
                }
                
                /* split molecule with a certain probability */
                if ((double)rand()/RAND_MAX < DISSOCIATION_PROB)
                    ncollisions = chem_split_molecule(i, particle, molecule, collisions, ncollisions);
                
                reac_dist = 2.5*MU;
                
                switch (particle[i].type) {
                    case (3):
                    {
                        ncollisions = chem_split_molecule_if_nearby(i, 3, particle, molecule, collisions, ncollisions, reac_dist, REACTION_PROB);
                        ncollisions = chem_split_molecule_if_nearby(i, 5, particle, molecule, collisions, ncollisions, reac_dist, REACTION_PROB);
                        ncollisions = chem_split_molecule_if_nearby(i, 6, particle, molecule, collisions, ncollisions, reac_dist, REACTION_PROB);
                        break;
                    }
                    case (4):
                    {
                        ncollisions = chem_split_molecule_if_nearby(i, 4, particle, molecule, collisions, ncollisions, reac_dist, REACTION_PROB);
                        ncollisions = chem_split_molecule_if_nearby(i, 5, particle, molecule, collisions, ncollisions, reac_dist, REACTION_PROB);
                        ncollisions = chem_split_molecule_if_nearby(i, 6, particle, molecule, collisions, ncollisions, reac_dist, REACTION_PROB);
                        break;
                    }
                    case (5):
                    {
                        ncollisions = chem_split_molecule_if_nearby(i, 3, particle, molecule, collisions, ncollisions, reac_dist, REACTION_PROB);
                        ncollisions = chem_split_molecule_if_nearby(i, 3, particle, molecule, collisions, ncollisions, reac_dist, REACTION_PROB);
                        ncollisions = chem_split_molecule_if_nearby(i, 5, particle, molecule, collisions, ncollisions, reac_dist, REACTION_PROB);
                        break;
                    }
                    case (6):
                    {
                        ncollisions = chem_split_molecule_if_nearby(i, 3, particle, molecule, collisions, ncollisions, reac_dist, REACTION_PROB);
                        ncollisions = chem_split_molecule_if_nearby(i, 4, particle, molecule, collisions, ncollisions, reac_dist, REACTION_PROB);
                        ncollisions = chem_split_molecule_if_nearby(i, 6, particle, molecule, collisions, ncollisions, reac_dist, REACTION_PROB);
                        break;
                    }
                }
                
                /* decouple particles with several partners from thermostat */
                if (particle[i].npartners >= AGREG_DECOUPLE) particle[i].thermostat = 0;
            }
            
            /* update cluster color scheme */
            if ((PLOT == P_CLUSTER)||(PLOT_B == P_CLUSTER))
                update_cluster_color(particle);
            
            printf("%i collisions\n", ncollisions);
            delta_n = ncollisions - oldncollisions; 
            printf("delta_n = %i\n", delta_n);
            return(ncollisions);
        }
        case (CHEM_DNA_ENZYME):
        {
            for (i=0; i<ncircles; i++) if (particle[i].active)
            {
                if (particle[i].npartners < AGREGMAX) 
                {
                    atype = particle[i].type;
                    if (atype%2 == 1) btype = atype+1;
                    else btype = atype-1;
                    if (atype > 0)
                        ncollisions = chem_multi_glue_molecule(i, btype, AGREGMAX, 3, 0, 0, 0, particle, molecule, collisions, ncollisions, REACTION_PROB);
                }
                
                /* split molecule with a certain probability */
                /* TODO update molecules after split */
                if ((double)rand()/RAND_MAX < DISSOCIATION_PROB)
                    ncollisions = chem_split_molecule(i, particle, molecule, collisions, ncollisions);
                
                reac_dist = 2.5*MU;
                
                switch (particle[i].type) {
                    case (3):
                    {
                        ncollisions = chem_local_split_molecule_if_nearby(i, 4, 7, particle, molecule, collisions, ncollisions, reac_dist, REACTION_PROB);
                        break;
                    }
                    case (4):
                    {
                        ncollisions = chem_local_split_molecule_if_nearby(i, 3, 7, particle, molecule, collisions, ncollisions, reac_dist, REACTION_PROB);
                        break;
                    }
                    case (5):
                    {
                        ncollisions = chem_local_split_molecule_if_nearby(i, 6, 7, particle, molecule, collisions, ncollisions, reac_dist, REACTION_PROB);
                        break;
                    }
                    case (6):
                    {
                        ncollisions = chem_local_split_molecule_if_nearby(i, 5, 7, particle, molecule, collisions, ncollisions, reac_dist, REACTION_PROB);
                        break;
                    }
                    case (7):
                    {
                        if ((double)rand()/RAND_MAX < KILLING_PROB) particle[i].active = 0;
                        break;
                    }
                }
                
                /* decouple particles with several partners from thermostat */
                if (particle[i].npartners >= AGREG_DECOUPLE) particle[i].thermostat = 0;
            }
            
            /* update cluster color scheme */
            if ((PLOT == P_CLUSTER)||(PLOT_B == P_CLUSTER))
                update_cluster_color(particle);
            
            /* for debugging */
//             for (i=0; i<nmolecules; i++)
//                 if (molecule[i].npartners > 0)
//                 {
//                     printf("Molecule %i has %i partners: ", i, molecule[i].npartners);
//                     for (j=0; j<molecule[i].npartners; j++)
//                         printf("%i ", molecule[i].partner[j]);
//                     printf("\n");
//                 }
//             sleep(3);
            
            printf("%i collisions\n", ncollisions);
            delta_n = ncollisions - oldncollisions; 
            printf("delta_n = %i\n", delta_n);
            return(ncollisions);
        }
        case (CHEM_DNA_ENZYME_REPAIR):
        {
            for (i=0; i<ncircles; i++) if (particle[i].active)
            {
                if (particle[i].npartners < AGREGMAX) 
                {
                    atype = particle[i].type;
                    if (atype%2 == 1) btype = atype+1;
                    else btype = atype-1;
                    if (atype > 0)  /* TEST */
                        ncollisions = chem_multi_glue_molecule(i, btype, AGREGMAX, 3, 0, 0, 0, particle, molecule, collisions, ncollisions, REACTION_PROB);
//                         ncollisions = chem_multi_glue_molecule(i, btype, AGREGMAX, 3, 0, 0, 2, particle, molecule, collisions, ncollisions, REACTION_PROB);
                }
                
                /* split molecule with a certain probability */
                /* TODO update molecules after split */
                if ((double)rand()/RAND_MAX < DISSOCIATION_PROB)
                    ncollisions = chem_split_molecule(i, particle, molecule, collisions, ncollisions);
                
                reac_dist = 2.5*MU;
                
                switch (particle[i].type) {
                    case (3):
                    {
                        ncollisions = chem_local_split_molecule_if_nearby(i, 4, 7, particle, molecule, collisions, ncollisions, reac_dist, REACTION_PROB);
                        break;
                    }
                    case (4):
                    {
                        ncollisions = chem_local_split_molecule_if_nearby(i, 3, 7, particle, molecule, collisions, ncollisions, reac_dist, REACTION_PROB);
                        break;
                    }
                    case (5):
                    {
                        ncollisions = chem_local_split_molecule_if_nearby(i, 6, 7, particle, molecule, collisions, ncollisions, reac_dist, REACTION_PROB);
                        break;
                    }
                    case (6):
                    {
                        ncollisions = chem_local_split_molecule_if_nearby(i, 5, 7, particle, molecule, collisions, ncollisions, reac_dist, REACTION_PROB);
                        break;
                    }
                    case (7):
                    {
                        if ((double)rand()/RAND_MAX < KILLING_PROB) particle[i].active = 0;
                        break;
                    }
                }
                
                /* decouple particles with several partners from thermostat */
                if (particle[i].npartners >= AGREG_DECOUPLE) particle[i].thermostat = 0;
            }
            
//             ncollisions = repair_dna(particle, molecule, collisions, ncollisions);
            
            repair_counter++;
            if (repair_counter >= 2)
            {
                ncollisions = repair_dna(particle, molecule, collisions, ncollisions);
//                 ncollisions = check_dna_pairing(particle, molecule, collisions, ncollisions);
                repair_counter = 0;
            }
            
            /* update cluster color scheme */
            if ((PLOT == P_CLUSTER)||(PLOT_B == P_CLUSTER))
                update_cluster_color(particle);
            
                        /* for debugging */
//             for (i=0; i<nmolecules; i++)
//                 if (molecule[i].npartners > 0)
//                 {
//                     printf("Molecule %i has %i partners: ", i, molecule[i].npartners);
//                     for (j=0; j<molecule[i].npartners; j++)
//                         printf("%i[%i] ", molecule[i].partner[j], molecule[i].connection_type[j]);
//                     printf("\n");
//                 }
//             sleep(3);


            printf("%i collisions\n", ncollisions);
            delta_n = ncollisions - oldncollisions; 
            printf("delta_n = %i\n", delta_n);
            return(ncollisions);
        }
        case (CHEM_POLYGON_AGGREGATION):
        {
            for (i=0; i<ncircles; i++) if (particle[i].active)
            {
                if (particle[i].npartners < AGREGMAX) for (k=0; k<3; k++)
                {
                    ncollisions = chem_multi_glue_polygon(i, k, AGREGMAX, particle, molecule, cluster, collisions, ncollisions, REACTION_PROB, 0);
                }
                
                /* decouple particles with several partners from thermostat */
                if (particle[i].npartners >= AGREG_DECOUPLE) particle[i].thermostat = 0;
            }
            
            /* update cluster color scheme */
            if ((PLOT == P_CLUSTER)||(PLOT_B == P_CLUSTER))
                update_cluster_color(particle);
            
            printf("%i collisions\n", ncollisions);
            delta_n = ncollisions - oldncollisions; 
            printf("delta_n = %i\n", delta_n);
            return(ncollisions);
        }
        case (CHEM_POLYGON_CLUSTER):
        {
            for (i=0; i<ncircles; i++) if (particle[i].active)
            {
                if (particle[i].npartners < NPOLY) 
                {
                    ncollisions = chem_multi_glue_polygon2(i, NPOLY, particle, molecule, cluster, collisions, ncollisions, REACTION_PROB, 1, 0);
                }
                
                /* decouple particles with several partners from thermostat */
                if (particle[i].npartners >= AGREG_DECOUPLE) particle[i].thermostat = 0;
            }
            
            /* repair clusters */
            if (!REPAIR_CLUSTERS) for (i=0; i<ncircles; i++) 
                if ((cluster[i].active)&&(cluster[i].nparticles >= 2))
                    repair_cluster(i, particle, cluster, 1000, 0);
                        
            printf("%i collisions\n", ncollisions);
            delta_n = ncollisions - oldncollisions; 
            printf("delta_n = %i\n", delta_n);
            return(ncollisions);
        }case (CHEM_POLYGON_ONECLUSTER):
        {
            i = 0;
            reacted = 0;
            while ((i<ncircles)&&(!reacted)) 
            {
                if (particle[i].active)
                {
                    if (particle[i].npartners < NPOLY) 
                    {
                        oldncollisions = ncollisions;
                        ncollisions = chem_multi_glue_polygon2(i, NPOLY, particle, molecule, cluster, collisions, ncollisions, REACTION_PROB, 1, 1);
                        if (ncollisions > oldncollisions) reacted = 1;
                    }
                
                    /* decouple particles with several partners from thermostat */
                    if (particle[i].npartners >= AGREG_DECOUPLE) particle[i].thermostat = 0;
                }
                i++;
            }
            
            /* repair clusters */
            if (!REPAIR_CLUSTERS) for (i=0; i<ncircles; i++) 
                if ((cluster[i].active)&&(cluster[i].nparticles >= 2))
                    repair_cluster(i, particle, cluster, 1000, 1);
                        
            printf("%i collisions\n", ncollisions);
            delta_n = ncollisions - oldncollisions; 
            printf("delta_n = %i\n", delta_n);
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
    static int second = 2, gshift = INITIAL_TIME + NSTEPS, nmax, ygrad, largewin;
    
    if (second > 0)
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
        
        second--;
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

int init_cluster_config(t_particle particle[NMAXCIRCLES], t_cluster cluster[NMAXCIRCLES])
/* initialize the clusters for option CLUSTER_PARTICLES */
/* returns number of active clusters */
{
    int i, j, tmp, nclusters = 0;
    
    for (i=0; i<ncircles; i++) /*if (particle[i].active)*/
    {
        nclusters++;
        cluster[i].active = particle[i].active;
        cluster[i].thermostat = 1;
        cluster[i].selected = 0;
        cluster[i].xg = particle[i].xc;
        cluster[i].yg = particle[i].yc;
        cluster[i].angle = particle[i].angle;
        cluster[i].vx = particle[i].vx;
        cluster[i].vy = particle[i].vy;
        cluster[i].omega = particle[i].omega;
        cluster[i].mass_inv = particle[i].mass_inv;
        cluster[i].mass = 1.0/particle[i].mass_inv;
        cluster[i].inertia_moment_inv = particle[i].inertia_moment_inv;
        cluster[i].inertia_moment = 1.0/particle[i].inertia_moment_inv;
        cluster[i].fx = particle[i].fx;
        cluster[i].fy = particle[i].fy;
        cluster[i].torque = particle[i].torque;
        cluster[i].energy = particle[i].energy;
        cluster[i].emean = particle[i].emean;
        
        cluster[i].nparticles = 1;
        cluster[i].particle[0] = i;
//         cluster[i].angle_ref = 0;
        
        particle[i].cluster = i;
        particle[i].cluster_color = i;
        particle[i].cluster_size = 1;
        particle[i].flip = 0;
    }
    
    /* randomize cluster number using Fisher-Yates algorithm */
    for (i=0; i<ncircles-1; i++)
    {
        j = i + rand()%(ncircles-i);
        tmp = particle[i].cluster_color;
        particle[i].cluster_color = particle[j].cluster_color;
        particle[j].cluster_color = tmp;
    }
//     for (i=0; i<ncircles; i++)
//     {
//         j = cluster[i].particle[0];
//         particle[j].cluster = i;
//     }
//     
//     for (i=0; i<ncircles; i++)
//     {
//         fprintf(lj_log, "Particle %i -> cluster %i -> particle %i\n", i, particle[i].cluster, cluster[particle[i].cluster].particle[0]);
//     }
    
    return(nclusters);
}



void compute_cluster_force(t_cluster cluster[NMAXCIRCLES], t_particle particle[NMAXCIRCLES])
/* compute force and torque on clusters from force and torque on their particles */
{
    int i, p, j;
    double fx, fy, torque, energy;
    
    for (i=0; i<ncircles; i++) if (cluster[i].active)
    {
        fx = 0.0; 
        fy = 0.0; 
        torque = 0.0;
        for (p=0; p<cluster[i].nparticles; p++)
        {
            j = cluster[i].particle[p];
            fx += particle[j].fx;
            fy += particle[j].fy;
            torque += particle[j].torque;
            torque += (particle[j].xc - cluster[i].xg)*particle[j].fy;
            torque -= (particle[j].yc - cluster[i].yg)*particle[j].fx;
        }
        cluster[i].fx = fx;
        cluster[i].fy = fy;
        cluster[i].torque = torque;
    }
}


void update_conveyor_belts(t_segment segment[NMAXSEGMENTS], t_belt belt[NMAXBELTS])
/* move shovels of conveyor belt */
{
    int b, shovel, firstseg, seg;
    double position, length, width, angle, newpos, deltapos, beltlength, shift, beta;
    
    for (b = 0; b < nbelts; b++)
    {
        length = belt[b].length;
        width = belt[b].width;
        angle = belt[b].angle;
        beltlength = 2.0*length + DPI*width;
        for (shovel = 0; shovel < belt[b].nshovels; shovel++)
        {
            position = belt[b].shovel_pos[shovel];
            deltapos = belt[b].speed*DT_PARTICLE;
            newpos = position + deltapos;
            firstseg = belt[b].shovel_segment[shovel];
            if (newpos < length)
            {
                shift = deltapos;
                for (seg = firstseg; seg < firstseg+6; seg++)
                    translate_one_segment(segment, seg, shift*belt[b].tx, shift*belt[b].ty);
            }
            else if (newpos < length + deltapos)
            {
                shift = length - position;
                beta = (newpos - length)/width;
                for (seg = firstseg; seg < firstseg+6; seg++)
                {
                    translate_one_segment(segment, seg, shift*belt[b].tx, shift*belt[b].ty);
                    rotate_one_segment(segment, seg, -beta, belt[b].x2, belt[b].y2);
                }
            }
            else if (newpos < length + PI*width)
            {
                beta = deltapos/width;
                for (seg = firstseg; seg < firstseg+6; seg++)
                {
                    rotate_one_segment(segment, seg, -beta, belt[b].x2, belt[b].y2);
                }
            }
            else if (newpos < length + PI*width + deltapos)
            {
                beta = (length + PI*width - position)/width;
                shift = newpos - length - PI*width;
                for (seg = firstseg; seg < firstseg+6; seg++)
                {
                    rotate_one_segment(segment, seg, -beta, belt[b].x2, belt[b].y2);
                    translate_one_segment(segment, seg, -shift*belt[b].tx, -shift*belt[b].ty);
                }
            }
            else if (newpos < 2.0*length + PI*width)
            {
                shift = deltapos;
                for (seg = firstseg; seg < firstseg+6; seg++)
                    translate_one_segment(segment, seg, -shift*belt[b].tx, -shift*belt[b].ty);
            }
            else if (newpos < 2.0*length + PI*width + deltapos)
            {
                shift = 2.0*length + PI*width - position;
                beta = (newpos - 2.0*length - PI*width)/width;
                for (seg = firstseg; seg < firstseg+6; seg++)
                {
                    translate_one_segment(segment, seg, -shift*belt[b].tx, -shift*belt[b].ty);
                    rotate_one_segment(segment, seg, -beta, belt[b].x2, belt[b].y2);
                }
            }
            else if (newpos < beltlength)
            {
                beta = deltapos/width;
                for (seg = firstseg; seg < firstseg+6; seg++)
                {
                    rotate_one_segment(segment, seg, -beta, belt[b].x1, belt[b].y1);
                }
            }
            else
            {
                beta = (beltlength - position)/width;
                shift = newpos - beltlength;
                for (seg = firstseg; seg < firstseg+6; seg++)
                {
                    rotate_one_segment(segment, seg, -beta, belt[b].x1, belt[b].y1);
                    translate_one_segment(segment, seg, shift*belt[b].tx, shift*belt[b].ty);
                }
            }
            
            if (newpos > beltlength) newpos -= beltlength;
            belt[b].shovel_pos[shovel] = newpos;
        }
    }
}

void draw_frame(int i, int plot, int bg_color, int ncollisions, int traj_position, 
                int traj_length, 
                int wall, double pressure[N_PRESSURES], double pleft, double pright, 
                int *particle_numbers, short int refresh, 
                t_lj_parameters params, t_particle particle[NMAXCIRCLES], 
                t_cluster cluster[NMAXCIRCLES], 
                t_collision *collisions, t_hashgrid hashgrid[HASHX*HASHY], 
                t_tracer trajectory[TRAJECTORY_LENGTH*N_TRACER_PARTICLES], 
                t_obstacle obstacle[NMAXOBSTACLES], t_segment segment[NMAXSEGMENTS], 
                t_group_data *group_speeds, t_group_segments *segment_group, 
                t_belt *conveyor_belt, int *tracer_n)
/* draw a movie frame */
{
    if (TRACER_PARTICLE) draw_trajectory(trajectory, traj_position, traj_length, particle, cluster, tracer_n, plot);
    draw_particles(particle, cluster, plot, params.beta, collisions, ncollisions, bg_color, hashgrid, params);
    draw_container(params.xmincontainer, params.xmaxcontainer, obstacle, segment, conveyor_belt, wall);
    if (PRINT_PARAMETERS) print_parameters(params, PRINT_LEFT, pressure, refresh);
    if (PLOT_SPEEDS) draw_speed_plot(group_speeds, i);
    if (PLOT_TRAJECTORIES) draw_trajectory_plot(group_speeds, i);
    if ((i > INITIAL_TIME)&&(PLOT_PARTICLE_NUMBER)) draw_particle_nb_plot(particle_numbers, i - INITIAL_TIME);
//     if (PLOT_PARTICLE_NUMBER) draw_particle_nb_plot(particle_numbers, i - INITIAL_TIME);
    if (BOUNDARY_COND == BC_EHRENFEST) print_ehrenfest_parameters(particle, pleft, pright);
    else if (PRINT_PARTICLE_NUMBER) 
    {
        if (REACTION_DIFFUSION) print_particle_types_number(particle, RD_TYPES);
        else print_particle_number(ncircles);
    }
    else if (PRINT_PARTICLE_SPEEDS) print_particles_speeds(particle);
    else if (PRINT_SEGMENTS_SPEEDS) print_segment_group_speeds(segment_group);
}

