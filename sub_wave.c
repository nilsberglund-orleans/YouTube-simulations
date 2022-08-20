/*********************/
/* Graphics routines */
/*********************/

#include "colors_waves.c"

int writetiff_new(char *filename, char *description, int x, int y, int width, int height, int compression)
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
  free(image); /* prenvents RAM consumption*/
  TIFFClose(file);
  return 0;
}


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

void blank()
{
    if (BLACK) glClearColor(0.0, 0.0, 0.0, 1.0);
    else glClearColor(1.0, 1.0, 1.0, 1.0);
    glClear(GL_COLOR_BUFFER_BIT);
}


void test_save_frame()  /* some tests with various resolutions */
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
    writetiff(n2, "Wave equation in a planar domain", 0, 0, WINWIDTH, WINHEIGHT, COMPRESSION_LZW);  // works for 1080p -> "-50px" 
    // choose one of the following according to the comment beside.
//     writetiff(n2, "Wave equation in a planar domain", 0, 0, WINWIDTH, WINHEIGHT-40, COMPRESSION_LZW);  
    /* to use with 1080p in drop_billiard.c- probably the best because it's
                                                                                                // generating 1080p image, lighter, and then cropping those 40 pixels to
                                                                                                // avoid the strange band*/
//     writetiff(n2, "Wave equation in a planar domain", 0, 0, WINWIDTH, WINHEIGHT-50, COMPRESSION_LZW);  // works for 1080p -> "-50px" band!!!
    // writetiff(n2, "Wave equation in a planar domain", 0, 0, 1920, 1080-40, COMPRESSION_LZW);           //another perfect 1080p from 1440p in setup
    // writetiff(n2, "Wave equation in a planar domain", -WINWIDTH/8+320, -WINHEIGHT/8+180, WINWIDTH-640, WINHEIGHT-400, COMPRESSION_LZW); // perfect 1040p from 1440p in setup
}

void test_save_frame_counter(int counter)  /* some tests with various resolutions */
/* same as save_frame, but with imposed image number (for option DOUBLE_MOVIE) */
{
  char *name="wave.", n2[100];
  char format[6]=".%05i";
  
    strcpy(n2, name);
    sprintf(strstr(n2,"."), format, counter);
    strcat(n2, ".tif");
    printf(" saving frame %s \n",n2);
    writetiff(n2, "Wave equation in a planar domain", 0, 0, WINWIDTH, WINHEIGHT, COMPRESSION_LZW);  // works for 1080p -> "-50px" 
    
    // choose one of the following according to the comment beside.
//     writetiff(n2, "Wave equation in a planar domain", 0, 0, WINWIDTH, WINHEIGHT-40, COMPRESSION_LZW);  
    /* to use with 1080p in drop_billiard.c- probably the best because it's
                                                                                                // generating 1080p image, lighter, and then cropping those 40 pixels to
                                                                                                // avoid the strange band*/
//     writetiff(n2, "Wave equation in a planar domain", 0, 0, WINWIDTH, WINHEIGHT-50, COMPRESSION_LZW);  // works for 1080p -> "-50px" band!!!
    // writetiff(n2, "Wave equation in a planar domain", 0, 0, 1920, 1080-40, COMPRESSION_LZW);           //another perfect 1080p from 1440p in setup
    // writetiff(n2, "BWave equation in a planar domain", -WINWIDTH/8+320, -WINHEIGHT/8+180, WINWIDTH-640, WINHEIGHT-400, COMPRESSION_LZW); // perfect 1040p from 1440p in setup
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

void save_frame_counter(int counter)
/* same as save_frame, but with imposed image number (for option DOUBLE_MOVIE) */
{
  char *name="wave.", n2[100];
  char format[6]=".%05i";

    strcpy(n2, name);
    sprintf(strstr(n2,"."), format, counter);
    strcat(n2, ".tif");
    printf(" saving frame %s \n",n2);
    writetiff(n2, "Wave equation in a planar domain", 0, 0,
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
 
 
// int in_polygon(double x, double y, double r, int npoly, double apoly)
// /* test whether (x,y) is in regular polygon of npoly sides inscribed in circle of radious r, turned by apoly Pi/2 */
// {
//     int condition = 1, k;
//     double omega, cw, angle; 
//     
//     omega = DPI/((double)npoly);
//     cw = cos(omega*0.5);
//     for (k=0; k<npoly; k++)  
//     {
//         angle = apoly*PID + ((double)k+0.5)*omega;
//         condition = condition*(x*cos(angle) + y*sin(angle) < r*cw);
//     }
//     return(condition);
// }



int in_tpolygon(double x, double y, t_polygon polygon)
/* test whether (x,y) is in polygon */
{
    int condition = 1, k;
    double omega, cw, angle, x1, y1; 
    
    x1 = (x-polygon.xc)/polygon.radius;
    y1 = (y-polygon.yc)/polygon.radius;
    
    /* first test whether point is in circumcircle */
    if (x1*x1 + y1*y1 >= 1.0) return(0);
    
    omega = DPI/((double)polygon.nsides);
    cw = cos(omega*0.5);
    for (k=0; k<polygon.nsides; k++)  
    {
        angle = polygon.angle*PID + ((double)k+0.5)*omega;
        condition = condition*(x1*cos(angle) + y1*sin(angle) < cw);
    }
    return(condition);
}

 /*********************/
/* drawing routines  */
/*********************/

/* The billiard boundary is drawn in (x,y) coordinates               */
/* However for the grid points, we use integer coordinates (i,j)     */
/* GL would allow to always work in (x,y) coordinates but using both */
/* sets of coordinates decreases number of double computations when  */
/* drawing the field                                                 */

void xy_to_ij(double x, double y, int ij[2])
/* convert (x,y) position to (i,j) in table representing wave */
{
    double x1, y1;

    x1 = (x - XMIN)/(XMAX - XMIN);
    y1 = (y - YMIN)/(YMAX - YMIN);

    ij[0] = (int)(x1 * (double)NX);
    ij[1] = (int)(y1 * (double)NY);
}


void xy_to_pos(double x, double y, double pos[2])
/* convert (x,y) position to double-valued position in table representing wave */
{
    double x1, y1;

    x1 = (x - XMIN)/(XMAX - XMIN);
    y1 = (y - YMIN)/(YMAX - YMIN);

    pos[0] = x1 * (double)NX;
    pos[1] = y1 * (double)NY;
}


void ij_to_xy(int i, int j, double xy[2])
/* convert (i,j) position in table representing wave to (x,y) */
{
    double x1, y1;

    xy[0] = XMIN + ((double)i)*(XMAX-XMIN)/((double)NX);
    xy[1] = YMIN + ((double)j)*(YMAX-YMIN)/((double)NY);
}

void erase_area(double x, double y, double dx, double dy)
{
    double pos[2], rgb[3];
    
    hsl_to_rgb(220.0, 0.8, 0.7, rgb);
    
    glColor3f(rgb[0], rgb[1], rgb[2]);
    glBegin(GL_QUADS);
    xy_to_pos(x - dx, y - dy, pos);
    glVertex2d(pos[0], pos[1]);
    xy_to_pos(x + dx, y - dy, pos);
    glVertex2d(pos[0], pos[1]);
    xy_to_pos(x + dx, y + dy, pos);
    glVertex2d(pos[0], pos[1]);
    xy_to_pos(x - dx, y + dy, pos);
    glVertex2d(pos[0], pos[1]);
    glEnd();
}


void erase_area_rgb(double x, double y, double dx, double dy, double rgb[3])
{
    double pos[2];
    
    glColor3f(rgb[0], rgb[1], rgb[2]);
    glBegin(GL_QUADS);
    xy_to_pos(x - dx, y - dy, pos);
    glVertex2d(pos[0], pos[1]);
    xy_to_pos(x + dx, y - dy, pos);
    glVertex2d(pos[0], pos[1]);
    xy_to_pos(x + dx, y + dy, pos);
    glVertex2d(pos[0], pos[1]);
    xy_to_pos(x - dx, y + dy, pos);
    glVertex2d(pos[0], pos[1]);
    glEnd();
}


void erase_area_hsl(double x, double y, double dx, double dy, double h, double s, double l)
{
    double pos[2], rgb[3];
    
    hsl_to_rgb(h, s, l, rgb);
    erase_area_rgb(x, y, dx, dy, rgb);
}

void draw_line(double x1, double y1, double x2, double y2)
{
    double pos[2];
    
    glBegin(GL_LINE_STRIP);
    xy_to_pos(x1, y1, pos);
    glVertex2d(pos[0], pos[1]);
    xy_to_pos(x2, y2, pos);
    glVertex2d(pos[0], pos[1]);
    glEnd();    
}

void draw_rectangle(double x1, double y1, double x2, double y2)
{
    double pos[2];
    
    glBegin(GL_LINE_LOOP);
    xy_to_pos(x1, y1, pos);
    glVertex2d(pos[0], pos[1]);
    xy_to_pos(x2, y1, pos);
    glVertex2d(pos[0], pos[1]);
    xy_to_pos(x2, y2, pos);
    glVertex2d(pos[0], pos[1]);
    xy_to_pos(x1, y2, pos);
    glVertex2d(pos[0], pos[1]);
    glEnd();    
}

void draw_rotated_rectangle(double x1, double y1, double x2, double y2)
{
    double pos[2];
    double xa, ya, xb, yb, xc, yc;
    
    xa = 0.5*(x1 - y2);
    xb = 0.5*(x2 - y1);
    xc = 0.5*(x1 - y1);
    ya = 0.5*(x1 + y1);
    yb = 0.5*(x2 + y2);
    yc = 0.5*(x2 + y1);
    
    glBegin(GL_LINE_LOOP);
    xy_to_pos(xc, ya, pos);
    glVertex2d(pos[0], pos[1]);
    xy_to_pos(xb, yc, pos);
    glVertex2d(pos[0], pos[1]);
    xy_to_pos(xc, yb, pos);
    glVertex2d(pos[0], pos[1]);
    xy_to_pos(xa, yc, pos);
    glVertex2d(pos[0], pos[1]);
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
        xy_to_pos(x + r*cos(alpha), y + r*sin(alpha), pos);
        glVertex2d(pos[0], pos[1]);
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
    xy_to_pos(x, y, pos);
    glVertex2d(pos[0], pos[1]);
    for (i=0; i<=nseg; i++)
    {
        alpha = (double)i*dalpha;
        xy_to_pos(x + r*cos(alpha), y + r*sin(alpha), pos);
        glVertex2d(pos[0], pos[1]);
    }
    
    glEnd();
}

void draw_tpolygon(t_polygon polygon)
{
    int i;
    double pos[2], alpha, dalpha;
    
    dalpha = DPI/(double)polygon.nsides;
    
    glBegin(GL_LINE_LOOP);
    for (i=0; i<=polygon.nsides; i++)
    {
        alpha = PID*polygon.angle + (double)i*dalpha;
        xy_to_pos(polygon.xc + polygon.radius*cos(alpha), polygon.yc + polygon.radius*sin(alpha), pos);
        glVertex2d(pos[0], pos[1]);
    }
    glEnd();
}


int init_circle_config_pattern(t_circle circles[NMAXCIRCLES], int circle_pattern)
/* initialise the arrays circlex, circley, circlerad and circleactive */
/* for billiard shape D_CIRCLES */
{
    int i, j, k, n, ncirc0, n_p_active, ncandidates=5000, naccepted; 
    double dx, dy, p, phi, r, r0, ra[5], sa[5], height, x, y = 0.0, gamma, dpoisson = 3.25*MU, xx[4], yy[4], dr, dphi;
    short int active_poisson[NMAXCIRCLES], far;
    
    switch (circle_pattern) {
        case (C_SQUARE):
        {
            ncircles = NGRIDX*NGRIDY;
            dy = (YMAX - YMIN)/((double)NGRIDY);
            for (i = 0; i < NGRIDX; i++)
                for (j = 0; j < NGRIDY; j++)
                {
                    n = NGRIDY*i + j;
                    circles[n].xc = ((double)(i-NGRIDX/2) + 0.5)*dy;
                    circles[n].yc = YMIN + ((double)j + 0.5)*dy;
                    circles[n].radius = MU;
                    circles[n].active = 1;
                }
            break;
        }
        case (C_HEX):
        {
            ncircles = NGRIDX*(NGRIDY+1);
            dy = (YMAX - YMIN)/((double)NGRIDY);
            dx = dy*0.5*sqrt(3.0);
            for (i = 0; i < NGRIDX; i++)
                for (j = 0; j < NGRIDY+1; j++)
                {
                    n = (NGRIDY+1)*i + j;
                    circles[n].xc = ((double)(i-NGRIDX/2) + 0.5)*dy;   /* is +0.5 needed? */
                    circles[n].yc = YMIN + ((double)j - 0.5)*dy;
                    if ((i+NGRIDX)%2 == 1) circles[n].yc += 0.5*dy;
                    circles[n].radius = MU;
                    /* activate only circles that intersect the domain */
                    if ((circles[n].yc < YMAX + MU)&&(circles[n].yc > YMIN - MU)) circles[n].active = 1;
                    else circles[n].active = 0;
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
                    circles[n].xc = ((double)(i-NGRIDX/2) + 0.5 + 0.5*((double)rand()/RAND_MAX - 0.5))*dy;
                    circles[n].yc = YMIN + 0.5 + ((double)j + 0.5 + 0.5*((double)rand()/RAND_MAX - 0.5))*dy;
                    circles[n].radius = MU*sqrt(1.0 + 0.8*((double)rand()/RAND_MAX - 0.5));
                    circles[n].active = 1;
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
                    circles[n].xc = ((double)(i-NGRIDX/2) + 0.5)*dy;
                    circles[n].yc = YMIN + ((double)j + 0.5)*dy;
                    circles[n].radius = MU;
                    p = (double)rand()/RAND_MAX;
                    if (p < P_PERCOL) circles[n].active = 1;
                    else circles[n].active = 0;
                }
            break;
        }
        case (C_RAND_POISSON):
        {
            ncircles = NPOISSON;
            for (n = 0; n < NPOISSON; n++)
            {
                circles[n].xc = LAMBDA*(2.0*(double)rand()/RAND_MAX - 1.0);
                circles[n].yc = (YMAX - YMIN)*(double)rand()/RAND_MAX + YMIN;
                circles[n].radius = MU;
                circles[n].active = 1;
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
                    circles[n].xc = r*cos(phi);
                    circles[n].yc = r*sin(phi);
                    circles[n].radius = MU;
                    circles[n].active = 1;
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
                    circles[n].xc = r*cos(phi);
                    circles[n].yc = r*sin(phi);
                    circles[n].radius = LAMBDA*ra[j];
                    circles[n].active = 1;
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
                    circles[4*i + j].xc = xx[i];
                    circles[4*i + j].yc = yy[j];
                    
                }
                
            circles[ncircles - 1].xc = X_TARGET;
            circles[ncircles - 1].yc = Y_TARGET;
            
            for (i=0; i<ncircles - 1; i++)
            {
                circles[i].radius = MU;
                circles[i].active = 1;
            }
            
            circles[ncircles - 1].radius = 0.5*MU;
            circles[ncircles - 1].active = 2;
            
            break;
        }        
        case (C_POISSON_DISC):
        {
            printf("Generating Poisson disc sample\n");
            /* generate first circle */
            circles[0].xc = LAMBDA*(2.0*(double)rand()/RAND_MAX - 1.0);
            circles[0].yc = (YMAX - YMIN)*(double)rand()/RAND_MAX + YMIN;
            active_poisson[0] = 1;
//             circles[0].active = 1;
            n_p_active = 1;
            ncircles = 1;
            
            while ((n_p_active > 0)&&(ncircles < NMAXCIRCLES))
            {
                /* randomly select an active circle */
                i = rand()%(ncircles);
                while (!active_poisson[i]) i = rand()%(ncircles);                 
//                 printf("Starting from circle %i at (%.3f,%.3f)\n", i, circles[i].xc, circles[i].yc);
                /* generate new candidates */
                naccepted = 0;
                for (j=0; j<ncandidates; j++)
                {
                    r = dpoisson*(2.0*(double)rand()/RAND_MAX + 1.0);
                    phi = DPI*(double)rand()/RAND_MAX;
                    x = circles[i].xc + r*cos(phi);
                    y = circles[i].yc + r*sin(phi);
//                        printf("Testing new circle at (%.3f,%.3f)\t", x, y);
                    far = 1;
                    for (k=0; k<ncircles; k++) if ((k!=i))
                    {
                        /* new circle is far away from circle k */
                        far = far*((x - circles[k].xc)*(x - circles[k].xc) + (y - circles[k].yc)*(y - circles[k].yc) >=     dpoisson*dpoisson);
                        /* new circle is in domain */
                        far = far*(vabs(x) < LAMBDA)*(y < YMAX)*(y > YMIN);
                    }
                    if (far)    /* accept new circle */
                    {
                        printf("New circle at (%.3f,%.3f) accepted\n", x, y);
                        circles[ncircles].xc = x;
                        circles[ncircles].yc = y;
                        circles[ncircles].radius = MU;
                        circles[ncircles].active = 1;
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
                circles[n].xc = -LAMBDA + n*dx;
                circles[n].yc = y;
                y += height*gamma; 
                if (y > YMAX) y -= height;
                circles[n].radius = MU;
                circles[n].active = 1;
            }
            
            /* test for circles that overlap top or bottom boundary */
            ncirc0 = ncircles;
            for (n=0; n < ncirc0; n++)
            {
                if (circles[n].yc + circles[n].radius > YMAX)
                {
                    circles[ncircles].xc = circles[n].xc;
                    circles[ncircles].yc = circles[n].yc - height;
                    circles[ncircles].radius = MU;
                    circles[ncircles].active = 1;
                    ncircles ++;
                }
                else if (circles[n].yc - circles[n].radius < YMIN)
                {
                    circles[ncircles].xc = circles[n].xc;
                    circles[ncircles].yc = circles[n].yc + height;
                    circles[ncircles].radius = MU;
                    circles[ncircles].active = 1;
                    ncircles ++;
                }
            }
            break;
        }
        case (C_GOLDEN_SPIRAL):
        {
            ncircles = 1;
            circles[0].xc = 0.0;
            circles[0].yc = 0.0;
            
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
                    circles[ncircles].xc = x;
                    circles[ncircles].yc = y;
                    ncircles++;
                }
            }
            
            for (i=0; i<ncircles; i++)
            {
                circles[i].radius = MU;
                /* inactivate circles outside the domain */
                if ((circles[i].yc < YMAX + MU)&&(circles[i].yc > YMIN - MU)) circles[i].active = 1;
//                 printf("i = %i, circlex = %.3lg, circley = %.3lg\n", i, circles[i].xc, circles[i].yc);
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
                    circles[n].xc = ((double)(i-NGRIDX/2) + 0.5)*dy;   /* is +0.5 needed? */
                    circles[n].yc = YMIN + ((double)j - 0.5)*dy;
                    if (((i+NGRIDX)%4 == 2)||((i+NGRIDX)%4 == 3)) circles[n].yc += 0.5*dy;
                    circles[n].radius = MU;
                    /* activate only circles that intersect the domain */
                    if ((circles[n].yc < YMAX + MU)&&(circles[n].yc > YMIN - MU)) circles[n].active = 1;
                    else circles[n].active = 0;
                }
            break;
        }
        case (C_RINGS):
        {
            ncircles = NGRIDX*NGRIDY;
            dphi = DPI/((double)NGRIDX);
            dr = 0.5*LAMBDA/(double)NGRIDY;
            for (i = 0; i < NGRIDX; i++)
                for (j = 0; j < NGRIDY; j++)
                {
                    n = NGRIDY*i + j;
                    phi = (double)i*dphi;
                    r = 0.5*LAMBDA + (double)j*dr;
                    circles[n].xc = r*cos(phi);
                    circles[n].yc = r*sin(phi);
                    circles[n].radius = MU;
                    /* activate only circles that intersect the domain */
                    if ((circles[n].yc < YMAX + MU)&&(circles[n].yc > YMIN - MU)) circles[n].active = 1;
                    else circles[n].active = 0;
                }
            break;
        }
        case (C_RINGS_T):
        {
            ncircles = NGRIDX*NGRIDY;
            dphi = DPI/((double)NGRIDX);
            dr = 0.5*LAMBDA/(double)NGRIDY;
            for (i = 0; i < NGRIDX; i++)
                for (j = 0; j < NGRIDY; j++)
                {
                    n = NGRIDY*i + j;
                    phi = (double)i*dphi;
                    phi += 0.5*(double)j*dphi;
                    r = 0.5*LAMBDA + (double)j*dr;
                    circles[n].xc = r*cos(phi);
                    circles[n].yc = r*sin(phi);
                    circles[n].radius = MU;
                    /* activate only circles that intersect the domain */
                    if ((circles[n].yc < YMAX + MU)&&(circles[n].yc > YMIN - MU)) circles[n].active = 1;
                    else circles[n].active = 0;
                }
            break;
        }
        case (C_RINGS_SPIRAL):
        {
            ncircles = 0;
//             circles[0].xc = 0.5*LAMBDA;
//             circles[0].yc = 0.0;
            
            gamma = (sqrt(5.0) - 1.0)*PI;    /* golden mean times 2Pi */
            phi = 0.0;
            r0 = 0.5*LAMBDA;
            r = r0 + MU;
            
            for (i=0; i<1000; i++) 
            {
                x = r*cos(phi);
                y = r*sin(phi);
                
                phi += gamma;
                r += 0.1*MU*r0/r;
                
                if (x*x + y*y < LAMBDA)
                {
                    circles[ncircles].xc = x;
                    circles[ncircles].yc = y;
                    ncircles++;
                }
            }
            
            for (i=0; i<ncircles; i++)
            {
                circles[i].radius = MU;
                /* inactivate circles outside the domain */
                if ((circles[i].yc < YMAX + MU)&&(circles[i].yc > YMIN - MU)) circles[i].active = 1;
            }
            break;
        }
        case (C_ONE):
        {
            circles[ncircles].xc = 0.0;
            circles[ncircles].yc = 0.0;
            circles[ncircles].radius = MU;
            circles[ncircles].active = 1;
            ncircles += 1;
            break;
        }
        case (C_TWO):   /* used for comparison with cloak */
        {
            circles[ncircles].xc = 0.0;
            circles[ncircles].yc = 0.0;
            circles[ncircles].radius = MU;
            circles[ncircles].active = 2;
            ncircles += 1;

            circles[ncircles].xc = 0.0;
            circles[ncircles].yc = 0.0;
            circles[ncircles].radius = 2.0*MU;
            circles[ncircles].active = 1;
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
    return(ncircles);
}

int init_circle_config(t_circle circles[NMAXCIRCLES])
/* for backward compatibility */
{
    return (init_circle_config_pattern(circles, CIRCLE_PATTERN));
}

int init_polygon_config_pattern(t_polygon polygons[NMAXCIRCLES], int circle_pattern)
/* initialise the polygon configuration, for billiard shape D_CIRCLES */
/* uses init_circle_config, this is where C++ would be more elegant */
{
    int i, ncircles;
    t_circle circle[NMAXCIRCLES];
    
    ncircles = init_circle_config_pattern(circle, circle_pattern);
    for (i=0; i<NMAXCIRCLES; i++)
    {
        polygons[i].xc = circle[i].xc;
        polygons[i].yc = circle[i].yc;
        polygons[i].radius = circle[i].radius;
        polygons[i].active = circle[i].active;
        polygons[i].nsides = NPOLY;
        
        if (RANDOM_POLY_ANGLE) polygons[i].angle = DPI*(double)rand()/RAND_MAX;
        else polygons[i].angle = APOLY;
/*        
        if (i < ncircles) printf("(x,y) = (%.2f, %.2f), r = %.2f, angle = %.2f, sides = %i\n", polygons[i].xc, polygons[i].yc, polygons[i].radius, polygons[i].angle, polygons[i].nsides);*/
    }
    
    /* adjust angles for C_RINGS configuration */
    if ((circle_pattern == C_RINGS)||(circle_pattern == C_RINGS_T)||(circle_pattern == C_RINGS_SPIRAL))
        for (i=0; i<ncircles; i++) if (polygons[i].active)
            polygons[i].angle += argument(polygons[i].xc, polygons[i].yc)/PID;
        
    return(ncircles);
}
    
int init_polygon_config(t_polygon polygons[NMAXCIRCLES])
/* for backward compatibility */
{
    return (init_polygon_config_pattern(polygons, CIRCLE_PATTERN));
}

int axial_symmetry(double z1[2], double z2[2], double z[2], double zprime[2])
/* compute reflection of point z wrt axis through z1 and z2 */
{
    double u[2], r, zdotu, zparallel[2], zperp[2];
    
    /* compute unit vector parallel to z1-z2 */
    u[0] = z2[0] - z1[0];
    u[1] = z2[1] - z1[1];
    r = module2(u[0], u[1]);
    if (r == 0) return(0);      /* z1 and z2 are the same */
    
    u[0] = u[0]/r;
    u[1] = u[1]/r;
//     printf("u = (%.2f, %.2f)\n", u[0], u[1]);
    
    /* projection of z1z on z1z2 */
    zdotu = (z[0] - z1[0])*u[0] + (z[1] - z1[1])*u[1];
    zparallel[0] = zdotu*u[0];
    zparallel[1] = zdotu*u[1];
//     printf("zparallel = (%.2f, %.2f)\n", zparallel[0], zparallel[1]);
    
    /* normal vector to z1z2 */
    zperp[0] = z[0] - z1[0] - zparallel[0];
    zperp[1] = z[1] - z1[1] - zparallel[1];
//     printf("zperp = (%.2f, %.2f)\n", zperp[0], zperp[1]);
    
    /* reflected point */
    zprime[0] = z[0] - 2.0*zperp[0];
    zprime[1] = z[1] - 2.0*zperp[1];
    
    return(1);
}

int axial_symmetry_tvertex(t_vertex z1, t_vertex z2, t_vertex z, t_vertex *zprime)
/* compute reflection of point z wrt axis through z1 and z2 */
{
    double r, zdotu;
    t_vertex u, zparallel, zperp;
    
    /* compute unit vector parallel to z1-z2 */
    u.x = z2.x - z1.x;
    u.y = z2.y - z1.y;
    r = module2(u.x, u.y);
    if (r == 0) return(0);      /* z1 and z2 are the same */
    
    u.x = u.x/r;
    u.y = u.y/r;
    
    /* projection of z1z on z1z2 */
    zdotu = (z.x - z1.x)*u.x + (z.y - z1.y)*u.y;
    zparallel.x = zdotu*u.x;
    zparallel.y = zdotu*u.y;
    
    /* normal vector to z1z2 */
    zperp.x = z.x - z1.x - zparallel.x;
    zperp.y = z.y - z1.y - zparallel.y;
    
    /* reflected point */
    zprime->x = z.x - 2.0*zperp.x;
    zprime->y = z.y - 2.0*zperp.y;
    
    return(1);
}

int compute_tokarsky_coordinates(double xshift, double yshift, double scaling, 
                                     t_vertex polyline[NMAXPOLY])
/* compute positions of vertices of tokarsky room */
{
    int i;
    double pos[2];
    
    polyline[0].x = 0.0;    polyline[0].y = 2.0;
    polyline[1].x = 1.0;    polyline[1].y = 3.0;
    polyline[2].x = 1.0;    polyline[2].y = 4.0;
    polyline[3].x = 2.0;    polyline[3].y = 4.0;
    polyline[4].x = 2.0;    polyline[4].y = 3.0;
    polyline[5].x = 3.0;    polyline[5].y = 3.0;
    polyline[6].x = 3.0;    polyline[6].y = 2.0;
    polyline[7].x = 5.0;    polyline[7].y = 2.0;
    polyline[8].x = 5.0;    polyline[8].y = 3.0;
    polyline[9].x = 6.0;    polyline[9].y = 3.0;

    polyline[10].x = 6.0;    polyline[10].y = 4.0;
    polyline[11].x = 7.0;    polyline[11].y = 3.0;
    polyline[12].x = 8.0;    polyline[12].y = 3.0;
    polyline[13].x = 8.0;    polyline[13].y = 2.0;
    polyline[14].x = 7.0;    polyline[14].y = 2.0;
    polyline[15].x = 7.0;    polyline[15].y = 1.0;
    polyline[16].x = 6.0;    polyline[16].y = 0.0;
    polyline[17].x = 6.0;    polyline[17].y = 1.0;
    polyline[18].x = 5.0;    polyline[18].y = 1.0;
    polyline[19].x = 4.0;    polyline[19].y = 0.0;

    polyline[20].x = 4.0;    polyline[20].y = 1.0;
    polyline[21].x = 3.0;    polyline[21].y = 1.0;
    polyline[22].x = 2.0;    polyline[22].y = 0.0;
    polyline[23].x = 2.0;    polyline[23].y = 1.0;
    polyline[24].x = 1.0;    polyline[24].y = 1.0;
    polyline[25].x = 1.0;    polyline[25].y = 2.0;
            
    for (i=0; i<26; i++)
    {
        polyline[i].x = (polyline[i].x + xshift)*scaling;
        polyline[i].y = (polyline[i].y + yshift)*scaling;
        xy_to_pos(polyline[i].x, polyline[i].y, pos);
        polyline[i].posi = pos[0];
        polyline[i].posj = pos[1];
    }
    return(26);        
}

void compute_isospectral_coordinates(int type, int ishift, double xshift, double yshift, double scaling, 
                                     t_vertex polyline[NMAXPOLY])
/* compute positions of vertices of isospectral billiards */
/* central triangle has coordinates (0,0), (1,0) and (LAMBDA,MU) fed into affine transformation */
/* defined by (xshift - 0.5), (yshift - 0.25) and scaling*/
{
    int i;
    double pos[2];
    
    polyline[ishift].x = (xshift - 0.5)*scaling;
    polyline[ishift].y = (yshift - 0.25)*scaling;
    
    polyline[ishift+1].x = (0.5+xshift)*scaling;
    polyline[ishift+1].y = (yshift - 0.25)*scaling;
    
    polyline[ishift+2].x = (LAMBDA+xshift - 0.5)*scaling;
    polyline[ishift+2].y = (MU+yshift - 0.25)*scaling; 
    
    axial_symmetry_tvertex(polyline[ishift], polyline[ishift+2], polyline[ishift+1], &polyline[ishift+3]);
    axial_symmetry_tvertex(polyline[ishift], polyline[ishift+1], polyline[ishift+2], &polyline[ishift+4]);
    axial_symmetry_tvertex(polyline[ishift+1], polyline[ishift+2], polyline[ishift], &polyline[ishift+5]);

    if (type == 0)
    {
        axial_symmetry_tvertex(polyline[ishift], polyline[ishift+3], polyline[ishift+2], &polyline[ishift+6]);
        axial_symmetry_tvertex(polyline[ishift+1], polyline[ishift+4], polyline[ishift], &polyline[ishift+7]);
        axial_symmetry_tvertex(polyline[ishift+2], polyline[ishift+5], polyline[ishift+1], &polyline[ishift+8]);
    }
    else
    {
        axial_symmetry_tvertex(polyline[ishift+2], polyline[ishift+3], polyline[ishift], &polyline[ishift+6]);
        axial_symmetry_tvertex(polyline[ishift], polyline[ishift+4], polyline[ishift+1], &polyline[ishift+7]);
        axial_symmetry_tvertex(polyline[ishift+1], polyline[ishift+5], polyline[ishift+2], &polyline[ishift+8]);
    }
    
    for (i=ishift; i<ishift+9; i++)
    {
        xy_to_pos(polyline[i].x, polyline[i].y, pos);
        polyline[i].posi = pos[0];
        polyline[i].posj = pos[1];
    }
}


int compute_tokaprime_coordinates(double xshift, t_vertex polyline[NMAXPOLY])
/* compute positions of vertices of Tokarsky room made of 86 triangles */
{
    double ta, tb, a, b, pos[2];
    int i;
    
    polyline[0].x = 0.0;
    polyline[0].y = 1.0;
    
    polyline[1].x = 0.0;
    polyline[1].y = 1.0 - LAMBDA;
    
    ta = tan(0.05*PI);
    tb = tan(0.4*PI);
    
    a = LAMBDA*tb/(ta + tb);
    b = a*ta;
    
    polyline[2].x = b;
    polyline[2].y = 1.0 - a;
    
    axial_symmetry_tvertex(polyline[0], polyline[2], polyline[1], &polyline[3]);
    axial_symmetry_tvertex(polyline[0], polyline[3], polyline[2], &polyline[4]);
    axial_symmetry_tvertex(polyline[3], polyline[4], polyline[0], &polyline[43]);
    
    for (i=4; i<42; i++)
        axial_symmetry_tvertex(polyline[i], polyline[43], polyline[i-1], &polyline[i+1]);
    
    for (i=2; i<44; i++)
    {
        polyline[i+42].x = -polyline[i].x;
        polyline[i+42].y = polyline[i].y;
    }
    
    for (i=0; i<86; i++)
    {
        polyline[i].x += xshift;
        xy_to_pos(polyline[i].x, polyline[i].y, pos);
        polyline[i].posi = pos[0];
        polyline[i].posj = pos[1];
    }
    
    return(86);
}

void compute_homophonic_coordinates(int type, int ishift, double xshift, double yshift, double scaling, 
                                    t_vertex polyline[NMAXPOLY])
/* compute positions of vertices of homophonic billiards */
{
    int i;
    double pos[2];
    
    polyline[ishift].x = (0.5 + xshift)*scaling;
    polyline[ishift].y = (yshift - 0.25)*scaling;
    
    polyline[ishift+1].x = (0.25 + xshift)*scaling;
    polyline[ishift+1].y = (0.25*sqrt(3.0) + yshift - 0.25)*scaling; 
    
    polyline[ishift+2].x = (xshift - 0.5)*scaling;
    polyline[ishift+2].y = (yshift - 0.25)*scaling;
    
    axial_symmetry_tvertex(polyline[ishift+1],  polyline[ishift+2],  polyline[ishift],    &polyline[ishift+3]);
    axial_symmetry_tvertex(polyline[ishift],    polyline[ishift+1],  polyline[ishift+2],  &polyline[ishift+21]);
    axial_symmetry_tvertex(polyline[ishift],    polyline[ishift+21], polyline[ishift+1],  &polyline[ishift+10]);
    axial_symmetry_tvertex(polyline[ishift+10], polyline[ishift+21], polyline[ishift+0],  &polyline[ishift+11]);
    axial_symmetry_tvertex(polyline[ishift+11], polyline[ishift+21], polyline[ishift+10], &polyline[ishift+13]);
    axial_symmetry_tvertex(polyline[ishift+11], polyline[ishift+13], polyline[ishift+21], &polyline[ishift+12]);
    axial_symmetry_tvertex(polyline[ishift+13], polyline[ishift+21], polyline[ishift+11], &polyline[ishift+14]);
    axial_symmetry_tvertex(polyline[ishift+14], polyline[ishift+21], polyline[ishift+13], &polyline[ishift+20]);
    axial_symmetry_tvertex(polyline[ishift+14], polyline[ishift+20], polyline[ishift+21], &polyline[ishift+15]);
    axial_symmetry_tvertex(polyline[ishift+20], polyline[ishift+15], polyline[ishift+14], &polyline[ishift+19]);
    
    if (type == 0)
    {
        axial_symmetry_tvertex(polyline[ishift],   polyline[ishift+2], polyline[ishift+1], &polyline[ishift+8]);
        axial_symmetry_tvertex(polyline[ishift+2], polyline[ishift+8], polyline[ishift+0], &polyline[ishift+7]);
        axial_symmetry_tvertex(polyline[ishift+2], polyline[ishift+7], polyline[ishift+8], &polyline[ishift+5]);
        axial_symmetry_tvertex(polyline[ishift+2], polyline[ishift+5], polyline[ishift+7], &polyline[ishift+4]);
        axial_symmetry_tvertex(polyline[ishift+5], polyline[ishift+7], polyline[ishift+2], &polyline[ishift+6]);
        axial_symmetry_tvertex(polyline[ishift],   polyline[ishift+8], polyline[ishift+2], &polyline[ishift+9]);
        
        axial_symmetry_tvertex(polyline[ishift+15], polyline[ishift+19], polyline[ishift+20], &polyline[ishift+16]);
        axial_symmetry_tvertex(polyline[ishift+16], polyline[ishift+19], polyline[ishift+15], &polyline[ishift+18]);
        axial_symmetry_tvertex(polyline[ishift+16], polyline[ishift+18], polyline[ishift+19], &polyline[ishift+17]);        
    }
    else
    {
        axial_symmetry_tvertex(polyline[ishift+2], polyline[ishift+3], polyline[ishift+1], &polyline[ishift+5]);
        axial_symmetry_tvertex(polyline[ishift+3], polyline[ishift+5], polyline[ishift+2], &polyline[ishift+4]);
        axial_symmetry_tvertex(polyline[ishift+2], polyline[ishift+5], polyline[ishift+3], &polyline[ishift+6]);
        axial_symmetry_tvertex(polyline[ishift+2], polyline[ishift+6], polyline[ishift+5], &polyline[ishift+8]);
        axial_symmetry_tvertex(polyline[ishift+6], polyline[ishift+8], polyline[ishift+2], &polyline[ishift+7]);
        axial_symmetry_tvertex(polyline[ishift+2], polyline[ishift+8], polyline[ishift+6], &polyline[ishift+9]);
        
        axial_symmetry_tvertex(polyline[ishift+10], polyline[ishift+11], polyline[ishift+21], &polyline[ishift+16]);
        axial_symmetry_tvertex(polyline[ishift+11], polyline[ishift+12], polyline[ishift+13], &polyline[ishift+18]);        
        axial_symmetry_tvertex(polyline[ishift+16], polyline[ishift+18], polyline[ishift+11], &polyline[ishift+17]);
    }  
    
    for (i=ishift; i<44; i++)
    {
        xy_to_pos(polyline[i].x, polyline[i].y, pos);
        polyline[i].posi = pos[0];
        polyline[i].posj = pos[1];
    }
}


int compute_vonkoch_coordinates(int depth, t_vertex polyline[NMAXPOLY])
/* compute positions of vertices of von Koch snowflake fractal */
{
    int nsides = 3, i, j, k, l, n, z, ii, jj, quater[MDEPTH], cond;
    short int vkoch[NMAXPOLY], turnright; 
    double angle, length, x, y, pos[2];
    
    for (k=0; k<depth; k++) nsides *= 4;
    ncircles = nsides;
    
    if (nsides > NMAXPOLY)
    {
        printf("NMAXPOLY needs to be increased to %i\n", nsides);
        nsides = NMAXPOLY;
    }

    for (i=0; i<nsides/3; i++)
    {
        /* compute quaternary expansion of i */
        ii = i;
        for (l=0; l<depth; l++)
        {
            quater[l] = ii%4;
            ii = ii - (ii%4);
            ii = ii/4;
        }
                
        /* find first nonzero digit */
        z = 0;
        while ((quater[z] == 0)&&(z<depth)) z++;
            
        /* compute left/right turns */
        if (i==0) vkoch[0] = 0;
        else if (z != depth)
        {   
            if (quater[z] == 2) vkoch[i] = 0;
            else vkoch[i] = 1;
        }
    }
            
    /* compute vertices */
    angle = APOLY*PID + 2.0*PI/3.0;
    x = LAMBDA*cos(APOLY*PID - PI/6.0);
    y = LAMBDA*sin(APOLY*PID - PI/6.0);
    length = 2.0*LAMBDA*sin(PI/3.0);
            
    for (k=0; k<depth; k++) length = length/3.0;
    printf("Length = %.2f\n", length);
            
    for (i=0; i<nsides; i++)
    {
        polyline[i].x = x*MU;
        polyline[i].y = y*MU;
            
        x += length*cos(angle);
        y += length*sin(angle);
                
        turnright = vkoch[i%(nsides/3)+1];
        if (turnright) angle -= PI/3.0;
        else angle += 2.0*PI/3.0;
            
        while (angle > DPI) angle -= DPI;
        while (angle < 0.0) angle += DPI;    
        
        xy_to_pos(x*MU, y*MU, pos);
        polyline[i].posi = pos[0];
        polyline[i].posj = pos[1];
    }
    
    return(nsides);
}


int compute_star_coordinates(t_vertex polyline[NMAXPOLY])
/* compute positions of vertices of star-shaped domain */
{
    int i;
    double alpha, r, x, y, pos[2];
    
    alpha = DPI/(double)NPOLY;
    
    for (i=0; i<NPOLY; i++)
    {
        if (i%2 == 0) r = LAMBDA - MU;
        else r = LAMBDA;
        
        x = r*cos(APOLY*PID + alpha*(double)i);
        y = r*sin(APOLY*PID + alpha*(double)i);
        polyline[i].x = x;
        polyline[i].y = y;
        
        xy_to_pos(x, y, pos);        
        polyline[i].posi = pos[0];
        polyline[i].posj = pos[1];
    }
    
    /* add origin to compute xy_in_billiard */
    polyline[NPOLY].x = 0.0;
    polyline[NPOLY].y = 0.0;
    
    return(NPOLY);
}

int compute_fresnel_coordinates(t_vertex polyline[NMAXPOLY])
/* compute positions of vertices approximating Fresnel lens */
{
    int i;
    double ymax, dy, x, y, x1, pos[2];
    
    ymax = 0.9*LAMBDA;
    dy = 2.0*ymax/(double)NSEG;
    
    if (LAMBDA > 0.0) x = -MU;
    else x = MU;
    polyline[0].x = x;
    polyline[0].y = -ymax;
    xy_to_pos(x, -ymax, pos);        
    polyline[0].posi = pos[0];
    polyline[0].posj = pos[1];
    
    for (i=1; i<NSEG; i++)
    {
        y = -ymax + dy*(double)i;
        x = sqrt(LAMBDA*LAMBDA - y*y) - vabs(LAMBDA);
//         x = sqrt(LAMBDA*LAMBDA - y*y) - LAMBDA*LAMBDA;
        
        while (x <= 0.0) x+= MU;
        if (LAMBDA < 0.0) x = -x;
        
        polyline[i].x = x;
        polyline[i].y = y;
        
        xy_to_pos(x, y, pos);        
        polyline[i].posi = pos[0];
        polyline[i].posj = pos[1];
    }
        
    if (LAMBDA > 0.0) x = -MU;
    else x = MU;
    polyline[NSEG].x = x;
    polyline[NSEG].y = ymax;
    xy_to_pos(x, ymax, pos);
    polyline[NSEG].posi = pos[0];
    polyline[NSEG].posj = pos[1];

    return(NSEG+1);
}

int compute_double_fresnel_coordinates(t_vertex polyline[NMAXPOLY], double xshift)
/* compute positions of vertices approximating two facing Fresnel lenses */
{
    int i;
    double pos[2];
    
    compute_fresnel_coordinates(polyline);
        
    for (i=0; i<=NSEG; i++)
    {
        polyline[i].x -= xshift;
        xy_to_pos(polyline[i].x, polyline[i].y, pos);        
        polyline[i].posi = pos[0];
        polyline[i].posj = pos[1];
        
        polyline[NSEG + 1 + i].x = -polyline[i].x;
        polyline[NSEG + 1 + i].y = polyline[i].y;
        xy_to_pos(polyline[NSEG + 1 + i].x, polyline[NSEG + 1 + i].y, pos);        
        polyline[NSEG + 1 + i].posi = pos[0];
        polyline[NSEG + 1 + i].posj = pos[1];
    }

    return(2*NSEG+2);
}


int compute_noisepanel_coordinates(t_vertex polyline[NMAXPOLY])
/* compute positions of vertices of noise panel */
{
    int i, n, even;
    double ymax, dy, x, y, x1, pos[2];
    
    /* find the leftmost point */
    x = 0.0;
    n = 0;
    while (x > XMIN) 
    {
        x -= LAMBDA;
        n++;
    }
    if (n%2 == 0) even = 1;
    else even = 0;
    
    i = 0;
    while (x <= XMAX + LAMBDA)
    {
        if (even) y = YMIN + 0.1;
        else y = YMIN + 0.1 + MU;
        
        x1 = x;
        if (x1 > XMAX) x1 = XMAX;
        else if (x1 < XMIN) x1 = XMIN;
        
        polyline[i].x = x1;
        polyline[i].y = y;
        
        xy_to_pos(x1, y, pos);        
        polyline[i].posi = pos[0];
        polyline[i].posj = pos[1];
        
        x += LAMBDA;
        even = 1 - even;
        i++;
    }
    n = i;
    for (i=0; i<n; i++)
    {
        polyline[n+i].x = polyline[n-i-1].x;
        polyline[n+i].y = -polyline[n-i-1].y;
        polyline[n+i].posi = polyline[n-i-1].posi;
        polyline[n+i].posj = NY - polyline[n-i-1].posj;
    }
        
    return(2*n);
}

int compute_noisepanel_rect_coordinates(t_vertex polyline[NMAXPOLY])
/* compute positions of vertices of noise panel */
{
    int i, n, even;
    double ymax, dy, x, y, x1, pos[2];
    
    /* find the leftmost point */
    x = -NPWIDTH;
    n = 0;
    while (x > XMIN) 
    {
        x -= LAMBDA;
        n++;
    }
    if (n%2 == 0) even = 1;
    else even = 0;
    
    i = 0;
    while (x <= 0.0)
    {
        if (even) y = YMIN + 0.1;
        else y = YMIN + 0.1 + MU;
        
        x1 = x;
        if (x1 > XMAX) x1 = XMAX;
        else if (x1 < XMIN) x1 = XMIN;
        
        polyline[i].x = x1;
        polyline[i].y = y;
        
        xy_to_pos(x1, y, pos);        
        polyline[i].posi = pos[0];
        polyline[i].posj = pos[1];
        
        x += LAMBDA;
        even = 1 - even;
        i++;
    }
    n = i;
    for (i=0; i<n; i++)
    {
        polyline[n+i].x = polyline[n-i-1].x;
        polyline[n+i].y = -polyline[n-i-1].y;
        polyline[n+i].posi = polyline[n-i-1].posi;
        polyline[n+i].posj = NY - polyline[n-i-1].posj;
    }
        
    return(2*n);
}

int compute_qrd_coordinates(t_vertex polyline[NMAXPOLY])
/* compute positions of quadratic noise diffuser */
{
    int n = 0, b, k, k1, kmin, kmax;
    double x, y, x1, y1 = YMIN, pos[2];
    
    kmin = (int)(XMIN/LAMBDA) - 2;
    kmax = (int)(XMAX/LAMBDA) + 2;
    
    for (b = -1; b <= 1; b+= 2)
    {
        if (b == 1) y1 = YMAX;
        for (k = kmin; k < kmax; k++)
        {
            x = LAMBDA*((double)(k) - 0.5);
            k1 = (k*k) % 13;
            if (b == -1) y = YMIN + (MU/13.0)*(14.0 - (double)k1);
            else y = YMAX - (MU/13.0)*(14.0 - (double)k1);
        
            polyline[n].x = x;
            polyline[n].y = y1;
            xy_to_pos(x, y1, pos);        
            polyline[n].posi = pos[0];
            polyline[n].posj = pos[1];
            n++;

            polyline[n].x = x;
            polyline[n].y = y;
            xy_to_pos(x, y, pos);        
            polyline[n].posi = pos[0];
            polyline[n].posj = pos[1];
            n++;
        
            y1 = y;
        }
    }
    
    return(n);
}

int init_polyline(int depth, t_vertex polyline[NMAXPOLY])
/* initialise variable polyline, for certain polygonal domain shapes */
{
    switch (B_DOMAIN) {
        case (D_TOKARSKY):
        {
            return(compute_tokarsky_coordinates(-4.0, -2.0, (XMAX - XMIN)/8.4, polyline));
        }
        case (D_TOKA_PRIME):
        {
            return(compute_tokaprime_coordinates(-MU, polyline));
        }
        case (D_ISOSPECTRAL):
        {
            compute_isospectral_coordinates(0, 0, ISO_XSHIFT_LEFT, ISO_YSHIFT_LEFT, ISO_SCALE, polyline);
            compute_isospectral_coordinates(1, 9, ISO_XSHIFT_RIGHT, ISO_YSHIFT_RIGHT, ISO_SCALE, polyline);
            return(18);
        }
        case (D_HOMOPHONIC):
        {
            compute_homophonic_coordinates(0, 0, ISO_XSHIFT_LEFT, ISO_YSHIFT_LEFT, ISO_SCALE, polyline);
            compute_homophonic_coordinates(1, 22, ISO_XSHIFT_RIGHT, ISO_YSHIFT_RIGHT, ISO_SCALE, polyline);
            return(44);
        }
        case (D_VONKOCH):
        {
            return(compute_vonkoch_coordinates(depth, polyline));
        }
        case (D_VONKOCH_HEATED):
        {
            return(compute_vonkoch_coordinates(depth, polyline));
        }
        case (D_STAR):
        {
            return(compute_star_coordinates(polyline));
        }
        case (D_FRESNEL):
        {
            return(compute_fresnel_coordinates(polyline));
        }
        case (D_DOUBLE_FRESNEL):
        {
            return(compute_double_fresnel_coordinates(polyline, LAMBDA));
        }
        case (D_NOISEPANEL):
        {
            return(compute_noisepanel_coordinates(polyline));
        }
        case (D_NOISEPANEL_RECT):
        {
            return(compute_noisepanel_rect_coordinates(polyline));
        }
        case (D_QRD):
        {
            return(compute_qrd_coordinates(polyline));
        }
        default:
        {
            return(0);
        }
    }
}

void isospectral_initial_point(double x, double y, double left[2], double right[2])
/* compute initial coordinates in isospectral billiards */
{
    left[0] = (x + ISO_XSHIFT_LEFT)*ISO_SCALE;
    left[1] = (y + ISO_YSHIFT_LEFT)*ISO_SCALE;
    right[0] = (x + ISO_XSHIFT_RIGHT)*ISO_SCALE;
    right[1] = (y + ISO_YSHIFT_RIGHT)*ISO_SCALE;    
}

void homophonic_initial_point(double xleft, double yleft, double xright, double yright, double left[2], double right[2])
/* compute initial coordinates in isospectral billiards */
{
    left[0] = (xleft + ISO_XSHIFT_LEFT)*ISO_SCALE;
    left[1] = (yleft + ISO_YSHIFT_LEFT)*ISO_SCALE;
    right[0] = (xright + ISO_XSHIFT_RIGHT)*ISO_SCALE;
    right[1] = (yright + ISO_YSHIFT_RIGHT)*ISO_SCALE;    
}

int xy_in_triangle(double x, double y, double z1[2], double z2[2], double z3[2])
/* returns 1 iff (x,y) is inside the triangle with vertices z1, z2, z3 */
{
    double v1, v2, v3;

    /* compute wedge products */
    v1 = (z2[0] - z1[0])*(y - z1[1]) - (z2[1] - z1[1])*(x - z1[0]);
    v2 = (z3[0] - z2[0])*(y - z2[1]) - (z3[1] - z2[1])*(x - z2[0]);
    v3 = (z1[0] - z3[0])*(y - z3[1]) - (z1[1] - z3[1])*(x - z3[0]);
    
    if ((v1 >= 0.0)&&(v2 >= 0.0)&&(v3 >= 0.0)) return(1);
    else return(0);
}

int xy_in_triangle_tvertex(double x, double y, t_vertex z1, t_vertex z2, t_vertex z3)
/* returns 1 iff (x,y) is inside the triangle with vertices z1, z2, z3 */
{
    double v1, v2, v3;

    /* compute wedge products */
    v1 = (z2.x - z1.x)*(y - z1.y) - (z2.y - z1.y)*(x - z1.x);
    v2 = (z3.x - z2.x)*(y - z2.y) - (z3.y - z2.y)*(x - z2.x);
    v3 = (z1.x - z3.x)*(y - z3.y) - (z1.y - z3.y)*(x - z3.x);
    
    if ((v1 >= 0.0)&&(v2 >= 0.0)&&(v3 >= 0.0)) return(1);
    else return(0);
}


int xy_in_billiard_single_domain(double x, double y, int b_domain, int ncirc, t_circle *circles)
/* returns 1 if (x,y) represents a point in the billiard */
{
    double l2, r2, r2mu, omega, b, c, angle, z, x1, y1, x2, y2, u, v, u1, v1, dx, dy, width, alpha, s, a;
    int i, j, k, k1, k2, condition = 0, m;
    static int first = 1, nsides;

    switch (b_domain) {
        case (D_NOTHING):
        {
            return(1);
            break;
        }
        case (D_RECTANGLE):
        {
            if ((vabs(x) <LAMBDA)&&(vabs(y) < 1.0)) return(1);
            else return(0);
            break;
        }
        case (D_ELLIPSE):
        {
            if (x*x/(LAMBDA*LAMBDA) + y*y < 1.0) return(1);
            else return(0);
            break;
        }
        case (D_STADIUM):
        {
            if ((x > -0.5*LAMBDA)&&(x < 0.5*LAMBDA)&&(y > -1.0)&&(y < 1.0)) return(1);
            else if (module2(x+0.5*LAMBDA, y) < 1.0) return(1);
            else if (module2(x-0.5*LAMBDA, y) < 1.0) return(1);
            else return(0);
            break;
        }
        case (D_SINAI):
        {
            if (x*x + y*y > LAMBDA*LAMBDA) return(1);
            else return(0);
            break;
        }
        case (D_DIAMOND):
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
        case (D_TRIANGLE):
        {
            if ((x>-LAMBDA)&&(y>-1.0)&&(LAMBDA*y+x<0.0)) return(1);
            else return(0);
            break;
        }
        case (D_FLAT):
        {
            if (y > -LAMBDA) return(1);
            else return(0);
            break;
        }
        case (D_ANNULUS):
        {
            l2 = LAMBDA*LAMBDA;
            r2 = x*x + y*y;
            if ((r2 > l2)&&(r2 < 1.0)) return(1);
            else return(0);
        }
        case (D_POLYGON):
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
        case (D_YOUNG):
        {
            if ((x < -MU)||(x > MU)) return(1);
            else if ((vabs(y-LAMBDA) < MU)||(vabs(y+LAMBDA) < MU)) return (1);
            else return(0);
        }
        case (D_GRATING):
        {
            k1 = -(int)((-YMIN)/LAMBDA);
            k2 = (int)(YMAX/LAMBDA);
            condition = 1;
            for (i=k1; i<= k2; i++)
            {
                z = (double)i*LAMBDA;
                condition = condition*(x*x + (y-z)*(y-z) > MU*MU);
            }
//             printf("x = %.3lg, y = %.3lg, k1 = %i, k2 = %i, condition = %i\n", x, y, k1, k2, condition);
            return(condition);
        }
        case (D_EHRENFEST):
        {
            return(((x-1.0)*(x-1.0) + y*y < LAMBDA*LAMBDA)||((x+1.0)*(x+1.0) + y*y < LAMBDA*LAMBDA)||((vabs(x) < 1.0)&&(vabs(y) < MU)));
        }
        case (D_DISK_GRID):
        {
            dy = (YMAX - YMIN)/((double)NGRIDY);
            for (i = -NGRIDX/2; i < NGRIDX/2; i++)
                for (j = 0; j < NGRIDY; j++)
                {
                    x1 = ((double)i + 0.5)*dy;
                    y1 = YMIN + ((double)j + 0.5)*dy;
                    if ((x-x1)*(x-x1) + (y-y1)*(y-y1) < MU*MU) return(0); 
                }
            return(1);
        }
        case (D_DISK_HEX):
        {
            dy = (YMAX - YMIN)/((double)NGRIDY);
            dx = dy*0.5*sqrt(3.0);
            for (i = -NGRIDX/2; i < NGRIDX/2; i++)
                for (j = -1; j < NGRIDY; j++)
                {
                    x1 = ((double)i + 0.5)*dy;
                    y1 = YMIN + ((double)j + 0.5)*dy;
                    if ((i+NGRIDX)%2 == 1) y1 += 0.5*dy;
                    if ((x-x1)*(x-x1) + (y-y1)*(y-y1) < MU*MU) return(0); 
                }
            return(1);
        }
        case (D_PARABOLA):
        {
            return(x > 0.25*y*y/LAMBDA - LAMBDA);
        }
        case (D_TWO_PARABOLAS):
        {
            x1 = 0.25*y*y/MU - MU - LAMBDA;
            x2 = -x1;
            width = 0.25*MU;
            if (width > 0.2) width = 0.2;
            if (vabs(y) > 1.5*MU) return(1);
            else if ((x < x1 - width)||(x > x2 + width)) return(1);
            else if ((x > x1)&&(x < x2)) return(1);
            else return(0);
        }
        case (D_FOUR_PARABOLAS):
        {
            x1 = MU + LAMBDA - 0.25*y*y/MU;
            y1 = MU + LAMBDA - 0.25*x*x/MU;
            return((vabs(x) < x1)&&(vabs(y) < y1));
        }
        case (D_POLY_PARABOLAS):
        {
            condition = 1;
            omega = DPI/((double)NPOLY);
            for (k=0; k<NPOLY; k++)  
            {
                angle = APOLY*PID + (k+0.5)*omega;
                x1 = x*cos(angle) + y*sin(angle);
                y1 = -x*sin(angle) + y*cos(angle);
                condition = condition*(x1 < LAMBDA + MU - 0.25*y1*y1/MU);
            }
            return(condition);
        }
        case (D_PENROSE):
        {
            c = sqrt(LAMBDA*LAMBDA - (1.0 - MU)*(1.0 - MU));
            width = 0.1*MU;
            x1 = vabs(x);
            y1 = vabs(y);
            /* sides */
            if (vabs(x) >= LAMBDA) return(0);
            /* upper and lower ellipse */
            else if ((vabs(y) >= MU)&&(x*x/(LAMBDA*LAMBDA) + (y1-MU)*(y1-MU)/((1.0-MU)*(1.0-MU)) >= 1.0)) return(0);
            /* small ellipses */
            else if ((vabs(x) <= c)&&(4.0*(x1-c)*(x1-c)/(MU*MU) + y*y/(MU*MU) <= 1.0)) return(0);
            /* straight parts */
            else if ((vabs(x) >= c)&&(vabs(y) <= width)) return(0);
            else return(1);
        }
        case (D_HYPERBOLA):
        {
            b = MU*sqrt(1.0 + x*x/(LAMBDA*LAMBDA - MU*MU)); 
            if (y > 1.02*b) return(1);
            else if (y < 0.98*b) return (1);
            else return(0);
        }
        case (D_TOKARSKY):
        {
            x1 = 4.0 + x/(XMAX - XMIN)*8.4;
            y1 = 2.0 + y/(XMAX - XMIN)*8.4;
            if ((x1 <= 0.0)||(x1 >= 8.0)) return(0);
            else if (x1 < 1.0)
            {
                if (y1 <= 2.0) return(0);
                else if (y1 >= x1 + 2.0) return(0);
                else return(1);
            }
            else if (x1 < 2.0)
            {
                if (y1 <= 1.0) return(0);
                else if (y1 >= 4.0) return(0);
                else return(1);
            }
            else if (x1 < 3.0)
            {
                if (y1 <= x1 - 2.0) return(0);
                else if (y1 >= 3.0) return(0);
                else return(1);
            }
            else if (x1 < 4.0)
            {
                if (y1 <= 1.0) return(0);
                else if (y1 >= 2.0) return(0);
                else return(1);
            }
            else if (x1 < 5.0)
            {
                if (y1 <= x1 - 4.0) return(0);
                else if (y1 >= 2.0) return(0);
                else return(1);
            }
            else if (x1 < 6.0)
            {
                if (y1 <= 1.0) return(0);
                else if (y1 >= 3.0) return(0);
                else return(1);
            }
            else if (x1 < 7.0)
            {
                if (y1 <= x1 - 6.0) return(0);
                else if (y1 >= 10.0 - x1) return(0);
                else return(1);
            }
            else
            {
                if (y1 <= 2.0) return(0);
                else if (y1 >= 3.0) return(0);
                else return(1);
            }
        }
        case (D_TOKA_PRIME):
        {
//             x1 = vabs(x);
            if (x + MU > 0.0) x1 = x;
            else x1 = -2.0*MU - x;
            
            condition = xy_in_triangle_tvertex(x1, y, polyline[0], polyline[1], polyline[2]);
            condition += xy_in_triangle_tvertex(x1, y, polyline[0], polyline[2], polyline[3]);
            condition += xy_in_triangle_tvertex(x1, y, polyline[i], polyline[3], polyline[4]);
            
            for (i=3; i<42; i++)
                condition += xy_in_triangle_tvertex(x1, y, polyline[i], polyline[43], polyline[i+1]);
            
            condition += xy_in_triangle_tvertex(x1, y, polyline[42], polyline[43], polyline[3]);
            return(condition >= 1);
        }
        case (D_ISOSPECTRAL):
        {
            /* 1st triangle */
            condition  = xy_in_triangle_tvertex(x, y, polyline[0], polyline[1], polyline[2]);
            condition += xy_in_triangle_tvertex(x, y, polyline[0], polyline[4], polyline[1]);
            condition += xy_in_triangle_tvertex(x, y, polyline[1], polyline[5], polyline[2]);
            condition += xy_in_triangle_tvertex(x, y, polyline[0], polyline[2], polyline[3]);
            condition += xy_in_triangle_tvertex(x, y, polyline[1], polyline[4], polyline[7]);
            condition += xy_in_triangle_tvertex(x, y, polyline[2], polyline[5], polyline[8]);
            condition += xy_in_triangle_tvertex(x, y, polyline[0], polyline[3], polyline[6]);

            /* 2nd triangle */
            condition += xy_in_triangle_tvertex(x, y,  polyline[9], polyline[10], polyline[11]);
            condition += xy_in_triangle_tvertex(x, y,  polyline[9], polyline[13], polyline[10]);
            condition += xy_in_triangle_tvertex(x, y, polyline[10], polyline[14], polyline[11]);
            condition += xy_in_triangle_tvertex(x, y,  polyline[9], polyline[11], polyline[12]);
            condition += xy_in_triangle_tvertex(x, y,  polyline[9], polyline[16], polyline[13]);
            condition += xy_in_triangle_tvertex(x, y, polyline[10], polyline[17], polyline[14]);
            condition += xy_in_triangle_tvertex(x, y, polyline[11], polyline[15], polyline[12]);
            return(condition >= 1);
        }
        case (D_HOMOPHONIC):
        {
            /* conditions could be summarised in larger triangles, but this is to keep */
            /* the option of using triangles with other angles than 30-60-90 */
            
            /* 1st triangle */
            condition = xy_in_triangle_tvertex(x, y, polyline[2], polyline[0], polyline[1]);
            condition += xy_in_triangle_tvertex(x, y, polyline[2], polyline[1], polyline[3]);
            condition += xy_in_triangle_tvertex(x, y, polyline[0], polyline[21], polyline[1]);
            condition += xy_in_triangle_tvertex(x, y, polyline[0], polyline[10], polyline[21]);
            condition += xy_in_triangle_tvertex(x, y, polyline[10], polyline[11], polyline[21]);
            condition += xy_in_triangle_tvertex(x, y, polyline[11], polyline[13], polyline[21]);
            condition += xy_in_triangle_tvertex(x, y, polyline[11], polyline[12], polyline[13]);
            condition += xy_in_triangle_tvertex(x, y, polyline[13], polyline[14], polyline[21]);
            condition += xy_in_triangle_tvertex(x, y, polyline[14], polyline[20], polyline[21]);
            condition += xy_in_triangle_tvertex(x, y, polyline[14], polyline[15], polyline[20]);
            condition += xy_in_triangle_tvertex(x, y, polyline[15], polyline[19], polyline[20]);

            condition += xy_in_triangle_tvertex(x, y, polyline[2], polyline[4], polyline[5]);
            condition += xy_in_triangle_tvertex(x, y, polyline[2], polyline[5], polyline[7]);
            condition += xy_in_triangle_tvertex(x, y, polyline[5], polyline[6], polyline[7]);
            condition += xy_in_triangle_tvertex(x, y, polyline[2], polyline[7], polyline[8]);
            condition += xy_in_triangle_tvertex(x, y, polyline[2], polyline[8], polyline[0]);
            condition += xy_in_triangle_tvertex(x, y, polyline[0], polyline[8], polyline[9]);
            condition += xy_in_triangle_tvertex(x, y, polyline[0], polyline[9], polyline[10]);

            condition += xy_in_triangle_tvertex(x, y, polyline[15], polyline[16], polyline[19]);
            condition += xy_in_triangle_tvertex(x, y, polyline[16], polyline[17], polyline[18]);
            condition += xy_in_triangle_tvertex(x, y, polyline[16], polyline[18], polyline[19]);

            /* 2nd triangle */
            condition += xy_in_triangle_tvertex(x, y, polyline[22+2], polyline[22+0], polyline[22+1]);
            condition += xy_in_triangle_tvertex(x, y, polyline[22+2], polyline[22+1], polyline[22+3]);
            condition += xy_in_triangle_tvertex(x, y, polyline[22+0], polyline[22+21], polyline[22+1]);
            condition += xy_in_triangle_tvertex(x, y, polyline[22+0], polyline[22+10], polyline[22+21]);
            condition += xy_in_triangle_tvertex(x, y, polyline[22+10], polyline[22+11], polyline[22+21]);
            condition += xy_in_triangle_tvertex(x, y, polyline[22+11], polyline[22+13], polyline[22+21]);
            condition += xy_in_triangle_tvertex(x, y, polyline[22+11], polyline[22+12], polyline[22+13]);
            condition += xy_in_triangle_tvertex(x, y, polyline[22+13], polyline[22+14], polyline[22+21]);
            condition += xy_in_triangle_tvertex(x, y, polyline[22+14], polyline[22+20], polyline[22+21]);
            condition += xy_in_triangle_tvertex(x, y, polyline[22+14], polyline[22+15], polyline[22+20]);
            condition += xy_in_triangle_tvertex(x, y, polyline[22+15], polyline[22+19], polyline[22+20]);
            
            condition += xy_in_triangle_tvertex(x, y, polyline[22+2], polyline[22+3], polyline[22+5]);
            condition += xy_in_triangle_tvertex(x, y, polyline[22+3], polyline[22+4], polyline[22+5]);
            condition += xy_in_triangle_tvertex(x, y, polyline[22+2], polyline[22+5], polyline[22+6]);
            condition += xy_in_triangle_tvertex(x, y, polyline[22+2], polyline[22+6], polyline[22+8]);
            condition += xy_in_triangle_tvertex(x, y, polyline[22+2], polyline[22+8], polyline[22+9]);
            condition += xy_in_triangle_tvertex(x, y, polyline[22+6], polyline[22+7], polyline[22+8]);

            condition += xy_in_triangle_tvertex(x, y, polyline[22+11], polyline[22+10], polyline[22+16]);
            condition += xy_in_triangle_tvertex(x, y, polyline[22+11], polyline[22+16], polyline[22+18]);
            condition += xy_in_triangle_tvertex(x, y, polyline[22+11], polyline[22+18], polyline[22+12]);
            condition += xy_in_triangle_tvertex(x, y, polyline[22+16], polyline[22+17], polyline[22+18]);
            
            return(condition >= 1);
        }
        case (D_CIRCLES):
        {
            for (i = 0; i < ncirc; i++)
                if (circles[i].active) 
                {
                    x1 = circles[i].xc;
                    y1 = circles[i].yc;
                    r2 = circles[i].radius*circles[i].radius;
                    if ((x-x1)*(x-x1) + (y-y1)*(y-y1) < r2) return(0); 
                }
            return(1);
        }
        case (D_CIRCLES_IN_RECT):   /* returns 2 inside circles, 0 outside rectangle */
        {
            for (i = 0; i < ncirc; i++)
                if (circles[i].active) 
                {
                    x1 = circles[i].xc;
                    y1 = circles[i].yc;
                    r2 = circles[i].radius*circles[i].radius;
                    if ((x-x1)*(x-x1) + (y-y1)*(y-y1) < r2) return(2); 
                }
            if ((vabs(x) >= LAMBDA)||(vabs(y) >= 1.0)) return(0);
            else return(1);
        }
        case (D_POLYGONS):
        {
            for (i = 0; i < ncirc; i++) 
                if ((polygons[i].active)&&(in_tpolygon(x, y, polygons[i]))) return(0);
            return(1);
        }
        case (D_VONKOCH):
        {
            condition = xy_in_triangle_tvertex(x, y, polyline[0], polyline[npolyline/3], polyline[2*npolyline/3]);
            m = 1;
            k = 1;
            for (i = 0; i < MDEPTH; i++)
            {
                m = m*4;
                for (j = 0; j < npolyline/m; j++)
                    condition += xy_in_triangle_tvertex(x, y, polyline[j*m + k], polyline[j*m + 2*k], polyline[j*m + 3*k]);
                k = k*4;
            }
            return(condition >= 1);
        }
        case (D_STAR):
        {
            condition = xy_in_triangle_tvertex(x, y, polyline[NPOLY], polyline[NPOLY-1], polyline[0]);
            for (i = 0; i < NPOLY-1; i++)
                condition += xy_in_triangle_tvertex(x, y, polyline[NPOLY], polyline[i], polyline[i+1]);
            return(condition >= 1);
        }
        case (D_FRESNEL):
        {
            if (vabs(y) > 0.9*vabs(LAMBDA)) return(1);
            if (vabs(x) > MU) return(1);
            
            x1 = sqrt(LAMBDA*LAMBDA - y*y) - vabs(LAMBDA);
            while (x1 <= 0.0) x1 += MU;
            if (LAMBDA > 0.0)
            {
                if (x < x1) return(0);
                else return(1);
            }
            else 
            {
                x1 = -x1;
                if (x > x1) return(0);
                else return(1);
            }
        }
        case (D_DOUBLE_FRESNEL):
        {
            if (vabs(y) > 0.9*vabs(LAMBDA)) return(1);
            if (LAMBDA > 0.0)
            {
                if (vabs(x) > LAMBDA + MU) return(1);
            
                x1 = sqrt(LAMBDA*LAMBDA - y*y) - LAMBDA;
                while (x1 <= 0.0) x1 += MU;
                x1 -= LAMBDA;
                if (vabs(x) > -x1) return(0);
                else return(1);
            }
            else
            {
                if (vabs(x) < -LAMBDA - MU) return(1);
            
                x1 = sqrt(LAMBDA*LAMBDA - y*y) + LAMBDA;
                while (x1 <= 0.0) x1 += MU;
                x1 -= LAMBDA;
                if (vabs(x) > x1) return(1);
                else return(0);
            }
        }
        case (D_NOISEPANEL):
        {
            x1 = vabs(x);
            while (x1 > 2.0*LAMBDA) x1 -= 2.0*LAMBDA;
            if (x1 <= LAMBDA) y1 = 0.1 + MU*x1/LAMBDA;
            else y1 = 0.1 + 2.0*MU - MU*x1/LAMBDA;
            return((y > YMIN + y1)&&(y < YMAX - y1));
        }
        case (D_NOISEPANEL_RECT):
        {
            x1 = -x;
            if (x1 > NPWIDTH)
            {
                while (x1 > 2.0*LAMBDA) x1 -= 2.0*LAMBDA;
                if (x1 <= LAMBDA) y1 = 0.1 + MU*x1/LAMBDA;
                else y1 = 0.1 + 2.0*MU - MU*x1/LAMBDA;
                return((y > YMIN + y1)&&(y < YMAX - y1)&&(x > XMIN + 0.1));
            }
            else if (x > NPWIDTH)
            {
                return((vabs(y) < YMAX - 0.1)&&(x < XMAX - 0.1));
            }
            else return(0);
        }
        case (D_QRD):
        {
            x1 = vabs(x)/LAMBDA;
            k = (int)(x1 + 0.5);
            k1 = (k*k) % 13;
            y1 = (MU/13.0)*(14.0 - (double)k1);
            return ((y > YMIN + y1)&&(y < YMAX - y1));
        }
        case (D_QRD_ASYM):
        {
            if (y > 0.0)
            {   
                x1 = vabs(x)/LAMBDA;
                k = (int)(x1 + 0.5);
                k1 = (k*k) % 13;
                y1 = (MU/13.0)*(14.0 - (double)k1);
                return (y < YMAX - y1);
            }
            else
            {   
                x1 = vabs(x + 1.0)/LAMBDA;
                k = (int)(x1 + 0.5);
                k1 = (k*k) % 17;
                y1 = (MU/17.0)*(18.0 - (double)k1);
                return (y > YMIN + y1);
            }
        }
        case (D_CIRCLE_SEGMENT):
        {
            if (vabs(y) > 0.9*vabs(LAMBDA)) return(1);
            
            y1 = 0.9*LAMBDA;
            x1 = sqrt(LAMBDA*LAMBDA - y1*y1) - vabs(LAMBDA) + MU;
            if ((LAMBDA > 0.0)&&(x < x1)) return(1);
            else if ((LAMBDA < 0.0)&&(x > -x1)) return(1);
            
            x1 = sqrt(LAMBDA*LAMBDA - y*y) - vabs(LAMBDA) + MU;
            if (LAMBDA > 0.0)
            {
                if (x < x1) return(0);
                else return(1);
            }
            else
            {
                if (x > -x1) return(0);
                else return(1);                
            }
        }
        case (D_GROOVE):
        {
            s = 0.85*LAMBDA;
            a = 0.5*LAMBDA;
            x1 = x - XMIN - (double)((int)((x - XMIN)/LAMBDA))*LAMBDA;
            if (x1 < a) return (y > YMIN + LAMBDA);
            else return (y > YMIN + LAMBDA + s);
        }
        case (D_FABRY_PEROT):
        {
            return(vabs(x - y*LAMBDA/YMAX) > 0.5*MU);
        }
        case (D_LSHAPE):
        {
            if (vabs(x) > LAMBDA) return(0);
            else if (vabs(y) > 1.0) return(0);
            else if ((x > 0.0)&&(y > 0.0)) return(0);
            else return(1);
        }
        case (D_MENGER):       
        {
            x1 = 0.5*(x+1.0);
            y1 = 0.5*(y+1.0);
            for (k=0; k<MDEPTH; k++)
            {
                x1 = x1*(double)MRATIO;
                y1 = y1*(double)MRATIO;
                if ((vabs(x)<1.0)&&(vabs(y)<1.0)&&(((int)x1 % MRATIO)==MRATIO/2)&&(((int)y1 % MRATIO)==MRATIO/2)) return(0);
            }
            return(1);
        }
        case (D_JULIA_INT):  
        {
            u = x/JULIA_SCALE;
            v = y/JULIA_SCALE;
            i = 0;
            while ((i<MANDELLEVEL)&&(u*u+v*v < 1000.0*MANDELLIMIT))
            {
                u1 = u*u - v*v + julia_x;
                v = 2.0*u*v + julia_y;
                u = u1;
                i++;
            }
            if (u*u + v*v < MANDELLIMIT) return(1);
            else return(0);
        }
        case (D_MENGER_ROTATED):       
        {
            x2 = 1.0*(x + y);
            y2 = 1.0*(x - y);
            if ((vabs(x2) < 1.0)&&(vabs(y2) < 1.0))
            {
                x1 = 0.5*(x2 + 1.0);
                y1 = 0.5*(y2 + 1.0);
                for (k=0; k<MDEPTH; k++)
                {
                    x1 = x1*(double)MRATIO;
                    y1 = y1*(double)MRATIO;
                    if ((vabs(x)<1.0)&&(vabs(y)<1.0)&&(((int)x1 % MRATIO)==MRATIO/2)&&(((int)y1 % MRATIO)==MRATIO/2)) return(0);
                }                
            }
            return(1);
        }
        case (D_ANNULUS_HEATED):      /* returns 2 if in inner circle */ 
        {
            l2 = LAMBDA*LAMBDA;
            r2 = x*x + y*y;
            r2mu = (x-MU)*(x-MU) + y*y;
            if ((r2mu > l2)&&(r2 < 1.0)) return(1);
            else if (r2mu <= l2) return(2);
            else return (0);
        }
        case (D_MENGER_HEATED): 
        {
            if ((vabs(x) >= 1.0)||(vabs(y) >= 1.0)) return(0);
            else
            {
                x1 = 0.5*(x+1.0);
                y1 = 0.5*(y+1.0);
                for (k=0; k<MDEPTH; k++)
                {
                    x1 = x1*(double)MRATIO;
                    y1 = y1*(double)MRATIO;
                    if ((((int)x1 % MRATIO)==MRATIO/2)&&(((int)y1 % MRATIO)==MRATIO/2)) return(k+2);
                }
                return(1);
            }
        }
        case (D_MENGER_H_OPEN):       /* returns 2 if in inner circle */ 
        {
            x1 = 0.5*(x+1.0);
            y1 = 0.5*(y+1.0);
            for (k=0; k<MDEPTH; k++)
            {
                x1 = x1*(double)MRATIO;
                y1 = y1*(double)MRATIO;
                if ((vabs(x)<1.0)&&(vabs(y)<1.0)&&(((int)x1 % MRATIO)==MRATIO/2)&&(((int)y1 % MRATIO)==MRATIO/2)) return(k+2);
            }
            return(1);
        }
        case (D_MANDELBROT):  
        {
            u = 0.0;
            v = 0.0;
            i = 0;
            while ((i<MANDELLEVEL)&&(u*u+v*v < 1000.0*MANDELLIMIT))
            {
                u1 = u*u - v*v + x;
                v = 2.0*u*v + y;
                u = u1;
                i++;
                /* old version used */
                /* u1 = u*u - v*v - x; */
                /* v = 2.0*u*v - y; */
            }
            if (u*u + v*v < MANDELLIMIT) return(0);
            else if ((x-0.5)*(x-0.5)/3.0 + y*y/1.0 > 1.2) return(2);
            else return(1);
        }
        case (D_MANDELBROT_CIRCLE):  
        {
            u = 0.0;
            v = 0.0;
            i = 0;
            while ((i<MANDELLEVEL)&&(u*u+v*v < 1000.0*MANDELLIMIT))
            {
                u1 = u*u - v*v + x;
                v = 2.0*u*v + y;
                u = u1;
                i++;
            }
            if (u*u + v*v < MANDELLIMIT) return(0);
            else if ((x-LAMBDA)*(x-LAMBDA) + (y-0.5)*(y-0.5) < MU*MU) return(2);
            else return(1);
        }
        case (D_JULIA):  
        {
            u = x/JULIA_SCALE;
            v = y/JULIA_SCALE;
            i = 0;
            while ((i<MANDELLEVEL)&&(u*u+v*v < 1000.0*MANDELLIMIT))
            {
                u1 = u*u - v*v + julia_x;
                v = 2.0*u*v + julia_y;
                u = u1;
                i++;
//                 printf("x = %.5lg y = %.5lg i = %i r2 = %.5lg\n", x, y, i, u*u+v*v);
            }
//             printf("i = %i x = %.5lg y = %.5lg r2 = %.5lg\n", i, x, y, u*u+v*v);
            if (u*u + v*v < MANDELLIMIT) return(0);
            else if (x*x/3.0 + y*y/1.0 > 1.2) return(2);
//             else if ((vabs(x) > XMAX - 0.01)||(vabs(y) > YMAX - 0.01)) return(2);
            else return(1);
        }
        case (D_VONKOCH_HEATED):
        {
            if (x*x + y*y > LAMBDA*LAMBDA) return(2);

            x1 = x;
            y1 = y;
            condition = xy_in_triangle_tvertex(x1, y1, polyline[0], polyline[npolyline/3], polyline[2*npolyline/3]);
            m = 1;
            k = 1;
            for (i = 0; i < MDEPTH; i++)
            {
                m = m*4;
                for (j = 0; j < npolyline/m; j++)
                    condition += xy_in_triangle_tvertex(x1, y1, polyline[j*m + k], polyline[j*m + 2*k], polyline[j*m + 3*k]);
                k = k*4;
            }
            if (condition > 0) return(0);
            else return(1);
        }
        default:
        {
            printf("Function ij_in_billiard not defined for this billiard \n");
            return(0);
        }
    }
}

int xy_in_billiard(double x, double y)
/* returns 1 if (x,y) represents a point in the billiard */
{
    if (COMPARISON)
    {
        if (y > 0.0) return (xy_in_billiard_single_domain(x, y, B_DOMAIN, ncircles, circles));
        else return (xy_in_billiard_single_domain(x, y, B_DOMAIN_B, ncircles_b, circles_b));
    }
    else return (xy_in_billiard_single_domain(x, y, B_DOMAIN, ncircles, circles));
}

int ij_in_billiard(int i, int j)
/* returns 1 if (i,j) represents a point in the billiard */
{
    double xy[2];

    ij_to_xy(i, j, xy);

    return(xy_in_billiard(xy[0], xy[1]));
}

void tvertex_lineto(t_vertex z)
/* draws boundary segments of isospectral billiard */
{
    glVertex2d(z.posi, z.posj);
}


void draw_billiard(int fade, double fade_value)      /* draws the billiard boundary */
{
    double x0, x, y, x1, y1, dx, dy, phi, r = 0.01, pos[2], pos1[2], alpha, dphi, omega, z, l, width, a, b, c, ymax;
    int i, j, k, k1, k2, mr2;
    static int first = 1, nsides;

    if (fade)
    {
        if (BLACK) glColor3f(fade_value, fade_value, fade_value);
        else glColor3f(1.0 - fade_value, 1.0 - fade_value, 1.0 - fade_value);        
    }
    else
    {
        if (BLACK) glColor3f(1.0, 1.0, 1.0);
        else glColor3f(0.0, 0.0, 0.0);
    }
    glLineWidth(BOUNDARY_WIDTH);

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
                if (fade) glColor3f(0.3*fade_value, 0.3*fade_value, 0.3*fade_value);
                else glColor3f(0.3, 0.3, 0.3);
                x0 = sqrt(LAMBDA*LAMBDA-1.0);

                glLineWidth(2);
                glEnable(GL_LINE_SMOOTH);
                
                draw_circle(x0, 0.0, r, NSEG);
                draw_circle(-x0, 0.0, r, NSEG);
            }
            break;
        }
        case (D_STADIUM):
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
        case (D_SINAI):
        {
            draw_circle(0.0, 0.0, LAMBDA, NSEG);
            break;
        }
        case (D_DIAMOND):
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
            draw_circle(0.0, 0.0, LAMBDA, NSEG);
            draw_circle(0.0, 0.0, 1.0, NSEG);
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
        case (D_YOUNG):
        {
            glBegin(GL_LINE_STRIP);
            xy_to_pos(-MU, YMIN, pos);
            glVertex2d(pos[0], pos[1]);
            xy_to_pos(-MU, -LAMBDA-MU, pos);
            glVertex2d(pos[0], pos[1]);
            xy_to_pos(MU, -LAMBDA-MU, pos);
            glVertex2d(pos[0], pos[1]);
            xy_to_pos(MU, YMIN, pos);
            glVertex2d(pos[0], pos[1]);
            glEnd();
            
            glBegin(GL_LINE_STRIP);
            xy_to_pos(-MU, YMAX, pos);
            glVertex2d(pos[0], pos[1]);
            xy_to_pos(-MU, LAMBDA+MU, pos);
            glVertex2d(pos[0], pos[1]);
            xy_to_pos(MU, LAMBDA+MU, pos);
            glVertex2d(pos[0], pos[1]);
            xy_to_pos(MU, YMAX, pos);
            glVertex2d(pos[0], pos[1]);
            glEnd();

            glBegin(GL_LINE_LOOP);
            xy_to_pos(-MU, -LAMBDA+MU, pos);
            glVertex2d(pos[0], pos[1]);
            xy_to_pos(-MU, LAMBDA-MU, pos);
            glVertex2d(pos[0], pos[1]);
            xy_to_pos(MU, LAMBDA-MU, pos);
            glVertex2d(pos[0], pos[1]);
            xy_to_pos(MU, -LAMBDA+MU, pos);
            glVertex2d(pos[0], pos[1]);
            glEnd();
            break;
        }
        case (D_GRATING):
        {
            k1 = -(int)(-YMIN/LAMBDA);
            k2 = (int)(YMAX/LAMBDA);
            for (i=k1; i<= k2; i++)
            {
                z = (double)i*LAMBDA;
                draw_circle(0.0, z, MU, NSEG);
            }
            break;
        }
        case (D_EHRENFEST):
        {
            alpha = asin(MU/LAMBDA);
            x0 = 1.0 - sqrt(LAMBDA*LAMBDA - MU*MU);
            dphi = 2.0*(PI-alpha)/((double)NSEG); 
            glBegin(GL_LINE_LOOP);
            for (i=0; i<=NSEG; i++)
            {
                phi = -PI + alpha + (double)i*dphi;
                x = 1.0 + LAMBDA*cos(phi);
                y = LAMBDA*sin(phi);
                xy_to_pos(x, y, pos);
                glVertex2d(pos[0], pos[1]);
            }
            for (i=0; i<=NSEG; i++)
            {
                phi = alpha + (double)i*dphi;
                x = -1.0 + LAMBDA*cos(phi);
                y = LAMBDA*sin(phi);
                xy_to_pos(x, y, pos);
                glVertex2d(pos[0], pos[1]);
            }
            glEnd ();
            break;
        }
        case (D_DISK_GRID):
        {
            glLineWidth(2);
            for (i = -NGRIDX/2; i < NGRIDX/2; i++)
                for (j = 0; j < NGRIDY; j++)
                {
                    dy = (YMAX - YMIN)/((double)NGRIDY);
                    dx = dy*0.5*sqrt(3.0);
                    x1 = ((double)i + 0.5)*dy;
                    y1 = YMIN + ((double)j + 0.5)*dy;
                    draw_circle(x1, y1, MU, NSEG);
                }
            break;
        }
        case (D_DISK_HEX):
        {
            glLineWidth(2);
            for (i = -NGRIDX/2; i < NGRIDX/2; i++)
                for (j = -1; j < NGRIDY; j++)
                {
                    dy = (YMAX - YMIN)/((double)NGRIDY);
                    x1 = ((double)i + 0.5)*dy;
                    y1 = YMIN + ((double)j + 0.5)*dy;
                    if ((i+NGRIDX)%2 == 1) y1 += 0.5*dy;
                    draw_circle(x1, y1, MU, NSEG);
                }
            break;
        }
        case (D_PARABOLA):
        {
            dy = (YMAX - YMIN)/(double)NSEG;
            glBegin(GL_LINE_STRIP);
            
            for (i = 0; i < NSEG+1; i++) 
            {
                y = YMIN + dy*(double)i;
                x = 0.25*y*y/LAMBDA - LAMBDA;
                xy_to_pos(x, y, pos);
                glVertex2d(pos[0], pos[1]);
            }
            glEnd ();
            
            if (FOCI)
            {
                glColor3f(0.3, 0.3, 0.3);
                draw_circle(0.0, 0.0, r, NSEG);
            }
            break;
        }
        case (D_TWO_PARABOLAS):
        {
            dy = 3.0*MU/(double)NSEG;
            width = 0.25*MU;
            if (width > 0.2) width = 0.2;
            glBegin(GL_LINE_LOOP);
            for (i = 0; i < NSEG+1; i++) 
            {
                y = -1.5*MU + dy*(double)i;
                x = 0.25*y*y/MU - MU - LAMBDA;
                xy_to_pos(x, y, pos);
                glVertex2d(pos[0], pos[1]);
            }
            for (i = 0; i < NSEG+1; i++) 
            {
                y = 1.5*MU - dy*(double)i;
                x = 0.25*y*y/MU - (MU + width) - LAMBDA;
                xy_to_pos(x, y, pos);
                glVertex2d(pos[0], pos[1]);
            }
            glEnd ();
            
            glBegin(GL_LINE_LOOP);
            for (i = 0; i < NSEG+1; i++) 
            {
                y = -1.5*MU + dy*(double)i;
                x = LAMBDA + MU - 0.25*y*y/MU;
                xy_to_pos(x, y, pos);
                glVertex2d(pos[0], pos[1]);
            }
            for (i = 0; i < NSEG+1; i++) 
            {
                y = 1.5*MU - dy*(double)i;
                x = LAMBDA + (MU + width) - 0.25*y*y/MU;
                xy_to_pos(x, y, pos);
                glVertex2d(pos[0], pos[1]);
            }
            glEnd ();
            
            if (FOCI)
            {
                glColor3f(0.3, 0.3, 0.3);
                draw_circle(-LAMBDA, 0.0, r, NSEG);
                draw_circle(LAMBDA, 0.0, r, NSEG);
            }

            break;
        }
        case (D_FOUR_PARABOLAS):
        {
            x1 = 2.0*(sqrt(MU*(2.0*MU + LAMBDA)) - MU);
            
            dy = 2.0*x1/(double)NSEG; 
            glBegin(GL_LINE_LOOP);
            for (i = 0; i < NSEG+1; i++) 
            {
                y = -x1 + dy*(double)i;
                x = MU + LAMBDA - 0.25*y*y/MU;
                xy_to_pos(x, y, pos);
                glVertex2d(pos[0], pos[1]);
            }
            for (i = 0; i < NSEG+1; i++) 
            {
                x = x1 - dy*(double)i;
                y = MU + LAMBDA - 0.25*x*x/MU;
                xy_to_pos(x, y, pos);
                glVertex2d(pos[0], pos[1]);
            }
            for (i = 0; i < NSEG+1; i++) 
            {
                y = x1 - dy*(double)i;
                x = -MU - LAMBDA + 0.25*y*y/MU;
                xy_to_pos(x, y, pos);
                glVertex2d(pos[0], pos[1]);
            }
            for (i = 0; i < NSEG+1; i++) 
            {
                x = -x1 + dy*(double)i;
                y = -MU - LAMBDA + 0.25*x*x/MU;
                xy_to_pos(x, y, pos);
                glVertex2d(pos[0], pos[1]);
            }
            glEnd ();
            
            if (FOCI)
            {
                glColor3f(0.3, 0.3, 0.3);
                draw_circle(-LAMBDA, 0.0, r, NSEG);
                draw_circle(LAMBDA, 0.0, r, NSEG);
                draw_circle(0.0, -LAMBDA, r, NSEG);
                draw_circle(0.0, LAMBDA, r, NSEG);
            }
            
            break;
        }
        case (D_POLY_PARABOLAS):
        {
            omega = PI/((double)NPOLY);
            a = 0.25/MU;
            b = 1.0/tan(omega);
            c = LAMBDA + MU;
            ymax = (-b + sqrt(b*b + 4.0*a*c))/(2.0*a);
            dy = 2.0*ymax/(double)NSEG; 
            
//             printf("a = %.3lg, b = %.3lg, ymax = %.3lg\n", a, b,ymax);
            glBegin(GL_LINE_LOOP);
            for (k=0; k<NPOLY; k++)  
            {
                alpha = APOLY*PID + (2.0*(double)k+1.0)*omega;
                for (i = 0; i < NSEG+1; i++) 
                {
                    y1 = -ymax + dy*(double)i;
                    x1 = MU + LAMBDA - 0.25*y1*y1/MU;
                    x = x1*cos(alpha) - y1*sin(alpha);
                    y = x1*sin(alpha) + y1*cos(alpha);
                    xy_to_pos(x, y, pos);
                    glVertex2d(pos[0], pos[1]);
                }
            }
            glEnd ();
            
            if (FOCI)
            {
                glColor3f(0.3, 0.3, 0.3);
                for (k=0; k<NPOLY; k++) 
                {
                    alpha = APOLY*PID + (2.0*(double)k+1.0)*omega;
                    draw_circle(LAMBDA*cos(alpha), LAMBDA*sin(alpha), r, NSEG);
                }
            }
            
            break;
        }
        case (D_PENROSE):
        {
            c = sqrt(LAMBDA*LAMBDA - (1.0 - MU)*(1.0 - MU));
            width = 0.1*MU;
            x1 = vabs(x);
            y1 = vabs(y);
            dphi = PI/(double)NSEG;
            
            glBegin(GL_LINE_LOOP);
            /* upper half ellipse */
            for (i=0; i<=NSEG; i++)
            {
                phi = (double)i*dphi;
                x = LAMBDA*cos(phi);
                y = MU + (1.0-MU)*sin(phi);
                xy_to_pos(x, y, pos);
                glVertex2d(pos[0], pos[1]);
            }
            
            /* straight parts */
            xy_to_pos(-LAMBDA, width, pos);
            glVertex2d(pos[0], pos[1]);
            xy_to_pos(-c, width, pos);
            glVertex2d(pos[0], pos[1]);
            xy_to_pos(-c, MU, pos);
            glVertex2d(pos[0], pos[1]);
            
            /* left half ellipse */
            for (i=0; i<=NSEG; i++)
            {
                phi = (double)i*dphi;
                x = -c + 0.5*MU*sin(phi);
                y = MU*cos(phi);
                xy_to_pos(x, y, pos);
                glVertex2d(pos[0], pos[1]);
            }
            
            /* straight parts */
            xy_to_pos(-c, -width, pos);
            glVertex2d(pos[0], pos[1]);
            xy_to_pos(-LAMBDA, -width, pos);
            glVertex2d(pos[0], pos[1]);
            xy_to_pos(-LAMBDA, -MU, pos);
            glVertex2d(pos[0], pos[1]);
            
            /* lower half ellipse */
            for (i=0; i<=NSEG; i++)
            {
                phi = (double)i*dphi;
                x = -LAMBDA*cos(phi);
                y = -MU - (1.0-MU)*sin(phi);
                xy_to_pos(x, y, pos);
                glVertex2d(pos[0], pos[1]);
            }
            
            /* straight parts */
            xy_to_pos(LAMBDA, -width, pos);
            glVertex2d(pos[0], pos[1]);
            xy_to_pos(c, -width, pos);
            glVertex2d(pos[0], pos[1]);
            xy_to_pos(c, -MU, pos);
            glVertex2d(pos[0], pos[1]);
            
            /* right half ellipse */
            for (i=0; i<=NSEG; i++)
            {
                phi = (double)i*dphi;
                x = c - 0.5*MU*sin(phi);
                y = -MU*cos(phi);
                xy_to_pos(x, y, pos);
                glVertex2d(pos[0], pos[1]);
            }
            
            /* straight parts */
            xy_to_pos(c, width, pos);
            glVertex2d(pos[0], pos[1]);
            xy_to_pos(LAMBDA, width, pos);
            glVertex2d(pos[0], pos[1]);
            xy_to_pos(LAMBDA, MU, pos);
            glVertex2d(pos[0], pos[1]);
            
            glEnd ();
            break; 
        }
        case (D_HYPERBOLA):
        {
            dx = (XMAX - XMIN)/(double)NSEG;
            glBegin(GL_LINE_STRIP);
            for (i = 0; i < NSEG+1; i++) 
            {
                x = XMIN + dx*(double)i;
                y = MU*1.02*sqrt(1.0 + x*x/(LAMBDA*LAMBDA - MU*MU));
                xy_to_pos(x, y, pos);
                glVertex2d(pos[0], pos[1]);
            }
            glEnd ();
            glBegin(GL_LINE_STRIP);
            for (i = 0; i < NSEG+1; i++) 
            {
                x = XMIN + dx*(double)i;
                y = MU*0.98*sqrt(1.0 + x*x/(LAMBDA*LAMBDA - MU*MU));
                xy_to_pos(x, y, pos);
                glVertex2d(pos[0], pos[1]);
            }
            glEnd ();
            
            if (FOCI)
            {
                glColor3f(0.3, 0.3, 0.3);
                draw_circle(0.0, LAMBDA, r, NSEG);
                draw_circle(0.0, -LAMBDA, r, NSEG);
            }
            break;
        }
        case (D_TOKARSKY):
        {
            glBegin(GL_LINE_LOOP);
            for (i=0; i<npolyline; i++) tvertex_lineto(polyline[i]);
            glEnd();
            if (FOCI)
            {
                x = (XMAX - XMIN)/4.2;
                glColor3f(0.3, 0.3, 0.3);
                draw_circle(x, 0.0, r, NSEG);
                draw_circle(-x, 0.0, r, NSEG);
            }
            break;
        }
        case (D_TOKA_PRIME):
        {
            glBegin(GL_LINE_LOOP);
            tvertex_lineto(polyline[0]);
            for (i=4; i<43; i++) tvertex_lineto(polyline[i]);
            tvertex_lineto(polyline[3]);
            tvertex_lineto(polyline[2]);
            tvertex_lineto(polyline[1]);

            tvertex_lineto(polyline[44]);
            tvertex_lineto(polyline[45]);
            for (i=84; i>45; i--) tvertex_lineto(polyline[i]);
            glEnd();
            
            /* inner lines */ 
//             glLineWidth(BOUNDARY_WIDTH/2);
            glLineWidth(1);
            glColor3f(0.75, 0.75, 0.75);
            glBegin(GL_LINE_STRIP);
            tvertex_lineto(polyline[0]);
            tvertex_lineto(polyline[1]);
            tvertex_lineto(polyline[2]);
            tvertex_lineto(polyline[0]);
            tvertex_lineto(polyline[3]);
            tvertex_lineto(polyline[4]);
            glEnd();
            
            glBegin(GL_LINE_STRIP);
            tvertex_lineto(polyline[0]);
            tvertex_lineto(polyline[44]);
            tvertex_lineto(polyline[45]);
            tvertex_lineto(polyline[0]);
            tvertex_lineto(polyline[46]);
            tvertex_lineto(polyline[45]);
            glEnd();
            
            for (i=3; i<43; i++)
            {
                glBegin(GL_LINE_STRIP);
                tvertex_lineto(polyline[i]);
                tvertex_lineto(polyline[43]);
                glEnd();
                glBegin(GL_LINE_STRIP);
                tvertex_lineto(polyline[i+42]);
                tvertex_lineto(polyline[85]);
                glEnd();
            }
            
            break;
        }
        case (D_ISOSPECTRAL):
        {
            /* 1st triangle */
            glBegin(GL_LINE_LOOP);
            tvertex_lineto(polyline[0]);
            tvertex_lineto(polyline[4]);
            tvertex_lineto(polyline[7]);
            tvertex_lineto(polyline[1]);
            tvertex_lineto(polyline[5]);
            tvertex_lineto(polyline[8]);
            tvertex_lineto(polyline[2]);
            tvertex_lineto(polyline[3]);
            tvertex_lineto(polyline[6]);
            glEnd();
            
            /* inner lines */
            glBegin(GL_LINE_LOOP);
            tvertex_lineto(polyline[0]);
            tvertex_lineto(polyline[1]);
            tvertex_lineto(polyline[2]);
            tvertex_lineto(polyline[0]);
            tvertex_lineto(polyline[3]);
            tvertex_lineto(polyline[2]);
            tvertex_lineto(polyline[5]);
            tvertex_lineto(polyline[1]);
            tvertex_lineto(polyline[4]);
            glEnd();
            
            /* 2nd triangle */
            glBegin(GL_LINE_LOOP);
            tvertex_lineto( polyline[9]);
            tvertex_lineto(polyline[16]);
            tvertex_lineto(polyline[13]);
            tvertex_lineto(polyline[10]);
            tvertex_lineto(polyline[17]);
            tvertex_lineto(polyline[14]);
            tvertex_lineto(polyline[11]);
            tvertex_lineto(polyline[15]);
            tvertex_lineto(polyline[12]);
            glEnd();
            
            /* inner lines */
            glBegin(GL_LINE_LOOP);
            tvertex_lineto( polyline[9]);
            tvertex_lineto(polyline[10]);
            tvertex_lineto(polyline[11]);
            tvertex_lineto( polyline[9]);
            tvertex_lineto(polyline[13]);
            tvertex_lineto(polyline[10]);
            tvertex_lineto(polyline[14]);
            tvertex_lineto(polyline[11]);
            tvertex_lineto(polyline[12]);
            glEnd();
            break;
        }
        case (D_HOMOPHONIC):
        {
            /* 1st triangle */
            glBegin(GL_LINE_LOOP);
            tvertex_lineto(polyline[1]);
            tvertex_lineto(polyline[3]);
            tvertex_lineto(polyline[4]);
            tvertex_lineto(polyline[5]);
            tvertex_lineto(polyline[6]);
            tvertex_lineto(polyline[8]);
            tvertex_lineto(polyline[9]);
            tvertex_lineto(polyline[10]);
            tvertex_lineto(polyline[12]);
            tvertex_lineto(polyline[13]);
            tvertex_lineto(polyline[15]);
            tvertex_lineto(polyline[16]);
            tvertex_lineto(polyline[17]);
            tvertex_lineto(polyline[18]);
            tvertex_lineto(polyline[20]);
            glEnd();
            
            /* inner lines */
            glLineWidth(BOUNDARY_WIDTH/2);
            glBegin(GL_LINE_STRIP);
            tvertex_lineto(polyline[9]);
            tvertex_lineto(polyline[1]);
            tvertex_lineto(polyline[2]);
            tvertex_lineto(polyline[5]);
            tvertex_lineto(polyline[7]);
            tvertex_lineto(polyline[2]);
            tvertex_lineto(polyline[8]);
            tvertex_lineto(polyline[21]);
            tvertex_lineto(polyline[10]);
            tvertex_lineto(polyline[2]);
            tvertex_lineto(polyline[21]);
            tvertex_lineto(polyline[11]);
            tvertex_lineto(polyline[13]);
            tvertex_lineto(polyline[21]);
            tvertex_lineto(polyline[14]);
            tvertex_lineto(polyline[20]);
            tvertex_lineto(polyline[15]);
            tvertex_lineto(polyline[19]);
            tvertex_lineto(polyline[16]);
            tvertex_lineto(polyline[18]);
            glEnd();
            
            /* 2nd triangle */
            glLineWidth(BOUNDARY_WIDTH);
            glBegin(GL_LINE_LOOP);
            tvertex_lineto(polyline[22+10]);
            tvertex_lineto(polyline[22+16]);
            tvertex_lineto(polyline[22+17]);
            tvertex_lineto(polyline[22+18]);
            tvertex_lineto(polyline[22+12]);
            tvertex_lineto(polyline[22+13]);
            tvertex_lineto(polyline[22+15]);
            tvertex_lineto(polyline[22+19]);
            tvertex_lineto(polyline[22+20]);
            tvertex_lineto(polyline[22+1]);
            tvertex_lineto(polyline[22+4]);
            tvertex_lineto(polyline[22+5]);
            tvertex_lineto(polyline[22+7]);
            tvertex_lineto(polyline[22+8]);
            tvertex_lineto(polyline[22+9]);
            glEnd();
            
            /* inner lines */
            glLineWidth(BOUNDARY_WIDTH/2);
            glBegin(GL_LINE_STRIP);
            tvertex_lineto(polyline[22+2]);
            tvertex_lineto(polyline[22+6]);
            tvertex_lineto(polyline[22+8]);
            tvertex_lineto(polyline[22+2]);
            tvertex_lineto(polyline[22+5]);
            tvertex_lineto(polyline[22+3]);
            tvertex_lineto(polyline[22+2]);
            tvertex_lineto(polyline[22+1]);
            tvertex_lineto(polyline[22+0]);
            tvertex_lineto(polyline[22+21]);
            tvertex_lineto(polyline[22+18]);
            tvertex_lineto(polyline[22+16]);
            tvertex_lineto(polyline[22+13]);
            tvertex_lineto(polyline[22+21]);
            tvertex_lineto(polyline[22+10]);
            tvertex_lineto(polyline[22+12]);
            tvertex_lineto(polyline[22+21]);
            tvertex_lineto(polyline[22+14]);
            tvertex_lineto(polyline[22+20]);
            tvertex_lineto(polyline[22+15]);
            glEnd();
            break;
        }
        case (D_VONKOCH):
        {
            glLineWidth(BOUNDARY_WIDTH/2);
            glBegin(GL_LINE_LOOP);
            for (i=0; i<npolyline; i++) tvertex_lineto(polyline[i]);
            glEnd();
            break;
        }
        case (D_STAR):
        {
            glLineWidth(BOUNDARY_WIDTH);
            glBegin(GL_LINE_LOOP);
            for (i=0; i<npolyline; i++) tvertex_lineto(polyline[i]);
            glEnd();
            break;
        }
        case (D_FRESNEL):
        {
            glLineWidth(BOUNDARY_WIDTH);
            glBegin(GL_LINE_LOOP);
            for (i=0; i<npolyline; i++) tvertex_lineto(polyline[i]);
            glEnd();
            break;
        }
        case (D_DOUBLE_FRESNEL):
        {
            glLineWidth(BOUNDARY_WIDTH);
            glBegin(GL_LINE_LOOP);
            for (i=0; i<npolyline/2; i++) tvertex_lineto(polyline[i]);
            glEnd();
            glBegin(GL_LINE_LOOP);
            for (i=npolyline/2; i<npolyline; i++) tvertex_lineto(polyline[i]);
            glEnd();
            break;
        }
        case (D_NOISEPANEL):
        {
            glLineWidth(BOUNDARY_WIDTH);
            glBegin(GL_LINE_STRIP);
            for (i=0; i<npolyline; i++) tvertex_lineto(polyline[i]);
            glEnd();
            break;
        }
        case (D_NOISEPANEL_RECT):
        {
            glLineWidth(BOUNDARY_WIDTH);
            glBegin(GL_LINE_STRIP);
            for (i=0; i<npolyline; i++) tvertex_lineto(polyline[i]);
            glEnd();
            break;
        }
        case (D_QRD):
        {
            glLineWidth(BOUNDARY_WIDTH);
            glBegin(GL_LINE_STRIP);
            for (i=0; i<npolyline/2; i++) tvertex_lineto(polyline[i]);
            glEnd();
            glBegin(GL_LINE_STRIP);
            for (i=npolyline/2; i<npolyline; i++) tvertex_lineto(polyline[i]);
            glEnd();
            break;
        }
        case (D_FABRY_PEROT):
        {
            glLineWidth(BOUNDARY_WIDTH);
            draw_line(-LAMBDA - 0.5*MU, YMIN, LAMBDA - 0.5*MU, YMAX);
            draw_line(-LAMBDA + 0.5*MU, YMIN, LAMBDA + 0.5*MU, YMAX);
            break;
        }
        case (D_CIRCLE_SEGMENT):
        {
            glLineWidth(BOUNDARY_WIDTH);
            glBegin(GL_LINE_LOOP);
            for (i=0; i<NSEG; i++) 
            {
                y = -0.9*LAMBDA + (double)i*1.8*LAMBDA/(double)NSEG;
                if (LAMBDA > 0.0) x = sqrt(LAMBDA*LAMBDA - y*y) - LAMBDA + MU;
                else x = -sqrt(LAMBDA*LAMBDA - y*y) - LAMBDA - MU;
                xy_to_pos(x, y, pos);
                glVertex2d(pos[0], pos[1]);
            }
            y = 0.9*LAMBDA;
            if (LAMBDA > 0.0) x = sqrt(LAMBDA*LAMBDA - y*y) - LAMBDA + MU;
            else x = -sqrt(LAMBDA*LAMBDA - y*y) - LAMBDA - MU;
            xy_to_pos(x, y, pos);
            glVertex2d(pos[0], pos[1]);
            glEnd();
            break;
        }
        case (D_CIRCLES):
        {
            glLineWidth(BOUNDARY_WIDTH);
            for (i = 0; i < ncircles; i++) 
                if (circles[i].active) draw_circle(circles[i].xc, circles[i].yc, circles[i].radius, NSEG);
            break;
        }
        case (D_CIRCLES_IN_RECT):
        {
            glLineWidth(BOUNDARY_WIDTH);
            for (i = 0; i < ncircles; i++) 
                if (circles[i].active) draw_circle(circles[i].xc, circles[i].yc, circles[i].radius, NSEG);
            draw_rectangle(-LAMBDA, -1.0, LAMBDA, 1.0);
            if ((FOCI)&&(CIRCLE_PATTERN == C_LASER))
            {
                glColor3f(0.3, 0.3, 0.3);
                draw_circle(X_SHOOTER, Y_SHOOTER, r, NSEG);
            }
            break;
        }
        case (D_POLYGONS):
        {
            glLineWidth(BOUNDARY_WIDTH);
            for (i = 0; i < ncircles; i++) 
                if (polygons[i].active) draw_tpolygon(polygons[i]);
            break;
        }
        case (D_MENGER):
        {
            glLineWidth(3);
//             draw_rectangle(XMIN, -1.0, XMAX, 1.0);
            
            /* level 1 */
            if (MDEPTH > 0)
            {
                glLineWidth(2);
                x = 1.0/((double)MRATIO);
                draw_rectangle(x, x, -x, -x);
            }
            
            /* level 2 */
            if (MDEPTH > 1)
            {
                glLineWidth(1);
                mr2 = MRATIO*MRATIO;
                l = 2.0/((double)mr2);
                
                for (i=0; i<MRATIO; i++)
                    for (j=0; j<MRATIO; j++)
                        if ((i!=MRATIO/2)||(j!=MRATIO/2))
                        {
                            x = -1.0 - 0.5*l + 2.0*((double)i + 0.5)/((double)MRATIO);
                            y = -1.0 - 0.5*l + 2.0*((double)j + 0.5)/((double)MRATIO);
                            draw_rectangle(x, y, x+l, y+l);
                        }
            }
            
            /* level 3 */
            if (MDEPTH > 2)
            {
                glLineWidth(1);
                l = 2.0/((double)(mr2*MRATIO));
                
                for (i=0; i<mr2; i++)
                    for (j=0; j<mr2; j++)
                        if ( (((i%MRATIO!=MRATIO/2))||(j%MRATIO!=MRATIO/2)) && (((i/MRATIO!=MRATIO/2))||(j/MRATIO!=MRATIO/2)) )
                        {
                            x = -1.0 - 0.5*l + 2.0*((double)i + 0.5)/((double)mr2);
                            y = -1.0 - 0.5*l + 2.0*((double)j + 0.5)/((double)mr2);
                            draw_rectangle(x, y, x+l, y+l);
                        }
            }
            
            break;
        }        
        case (D_JULIA_INT):
        {
            /* Do nothing */
            break;
        }
        case (D_MENGER_ROTATED):
        {
            glLineWidth(3);
//             draw_rectangle(XMIN, -1.0, XMAX, 1.0);

            /* level 1 */
            if (MDEPTH > 0)
            {
                glLineWidth(2);
                x = 1.0/((double)MRATIO);
                draw_rotated_rectangle(x, x, -x, -x);
            }
            
            /* level 2 */
            if (MDEPTH > 1)
            {
                glLineWidth(1);
                mr2 = MRATIO*MRATIO;
                l = 2.0/((double)mr2);
                
                for (i=0; i<MRATIO; i++)
                    for (j=0; j<MRATIO; j++)
                        if ((i!=MRATIO/2)||(j!=MRATIO/2))
                        {
                            x = -1.0 - 0.5*l + 2.0*((double)i + 0.5)/((double)MRATIO);
                            y = -1.0 - 0.5*l + 2.0*((double)j + 0.5)/((double)MRATIO);
                            draw_rotated_rectangle(x, y, x+l, y+l);
                        }
            }
            
            /* level 3 */
            if (MDEPTH > 2)
            {
                glLineWidth(1);
                l = 2.0/((double)(mr2*MRATIO));
                
                for (i=0; i<mr2; i++)
                    for (j=0; j<mr2; j++)
                        if ( (((i%MRATIO!=MRATIO/2))||(j%MRATIO!=MRATIO/2)) && (((i/MRATIO!=MRATIO/2))||(j/MRATIO!=MRATIO/2)) )
                        {
                            x = -1.0 - 0.5*l + 2.0*((double)i + 0.5)/((double)mr2);
                            y = -1.0 - 0.5*l + 2.0*((double)j + 0.5)/((double)mr2);
                            draw_rotated_rectangle(x, y, x+l, y+l);
                        }
            }
            
            break;
        }        
        case (D_ANNULUS_HEATED):
        {
            draw_circle(MU, 0.0, LAMBDA, NSEG);
            draw_circle(0.0, 0.0, 1.0, NSEG);
            break;
        }
        case (D_MENGER_HEATED):
        {
            glLineWidth(3);
            draw_rectangle(-1.0, -1.0, 1.0, 1.0);
            
            /* level 1 */
            if (MDEPTH > 0)
            {
                glLineWidth(2);
                x = 1.0/((double)MRATIO);
                draw_rectangle(x, x, -x, -x);
            }
            
            /* level 2 */
            if (MDEPTH > 1)
            {
                glLineWidth(1);
                mr2 = MRATIO*MRATIO;
                l = 2.0/((double)mr2);
                
                for (i=0; i<MRATIO; i++)
                    for (j=0; j<MRATIO; j++)
                        if ((i!=MRATIO/2)||(j!=MRATIO/2))
                        {
                            x = -1.0 - 0.5*l + 2.0*((double)i + 0.5)/((double)MRATIO);
                            y = -1.0 - 0.5*l + 2.0*((double)j + 0.5)/((double)MRATIO);
                            draw_rectangle(x, y, x+l, y+l);
                        }
            }
            
            /* level 3 */
            if (MDEPTH > 2)
            {
                glLineWidth(1);
                l = 2.0/((double)(mr2*MRATIO));
                
                for (i=0; i<mr2; i++)
                    for (j=0; j<mr2; j++)
                        if ( (((i%MRATIO!=MRATIO/2))||(j%MRATIO!=MRATIO/2)) && (((i/MRATIO!=MRATIO/2))||(j/MRATIO!=MRATIO/2)) )
                        {
                            x = -1.0 - 0.5*l + 2.0*((double)i + 0.5)/((double)mr2);
                            y = -1.0 - 0.5*l + 2.0*((double)j + 0.5)/((double)mr2);
                            draw_rectangle(x, y, x+l, y+l);
                        }
            }
            
            break;
        }
        case (D_MENGER_H_OPEN):
        {
            glLineWidth(3);
//             draw_rectangle(XMIN, -1.0, XMAX, 1.0);
            
            /* level 1 */
            if (MDEPTH > 0)
            {
                glLineWidth(2);
                x = 1.0/((double)MRATIO);
                draw_rectangle(x, x, -x, -x);
            }
            
            /* level 2 */
            if (MDEPTH > 1)
            {
                glLineWidth(1);
                mr2 = MRATIO*MRATIO;
                l = 2.0/((double)mr2);
                
                for (i=0; i<MRATIO; i++)
                    for (j=0; j<MRATIO; j++)
                        if ((i!=MRATIO/2)||(j!=MRATIO/2))
                        {
                            x = -1.0 - 0.5*l + 2.0*((double)i + 0.5)/((double)MRATIO);
                            y = -1.0 - 0.5*l + 2.0*((double)j + 0.5)/((double)MRATIO);
                            draw_rectangle(x, y, x+l, y+l);
                        }
            }
            
            /* level 3 */
            if (MDEPTH > 2)
            {
                glLineWidth(1);
                l = 2.0/((double)(mr2*MRATIO));
                
                for (i=0; i<mr2; i++)
                    for (j=0; j<mr2; j++)
                        if ( (((i%MRATIO!=MRATIO/2))||(j%MRATIO!=MRATIO/2)) && (((i/MRATIO!=MRATIO/2))||(j/MRATIO!=MRATIO/2)) )
                        {
                            x = -1.0 - 0.5*l + 2.0*((double)i + 0.5)/((double)mr2);
                            y = -1.0 - 0.5*l + 2.0*((double)j + 0.5)/((double)mr2);
                            draw_rectangle(x, y, x+l, y+l);
                        }
            }
            
            break;
        }        
        case (D_MANDELBROT):
        {
            /* Do nothing */
            break;
        }
        case (D_MANDELBROT_CIRCLE):
        {
            /* Do nothing */
            break;
        }
        case (D_JULIA):
        {
            /* Do nothing */
            break;
        }
        case (D_VONKOCH_HEATED):
        {
            glLineWidth(BOUNDARY_WIDTH/2);
            glBegin(GL_LINE_LOOP);
            for (i=0; i<npolyline; i++) glVertex2d(polyline[i].posi, polyline[i].posj);
            glEnd();
            break;
        }
        case (D_NOTHING):
        {
            break;
        }   
        default:
        {
            printf("Function draw_billiard not defined for this billiard \n");
        }
    }
}

void draw_color_scheme(double x1, double y1, double x2, double y2, int plot, double min, double max)
{
    int j, k, ij_botleft[2], ij_topright[2], imin, imax, jmin, jmax;
    double y, dy, dy_e, rgb[3], value, lum, amp;
    
    xy_to_ij(x1, y1, ij_botleft);
    xy_to_ij(x2, y2, ij_topright);
    
    rgb[0] = 0.0;   rgb[1] = 0.0;   rgb[2] = 0.0;
    erase_area_rgb(0.5*(x1 + x2), x2 - x1, 0.5*(y1 + y2), y2 - y1, rgb);

    if (ROTATE_COLOR_SCHEME)
    {
        jmin = ij_botleft[0];
        jmax = ij_topright[0];
        imin = ij_botleft[1];
        imax = ij_topright[1];    
    }
    else
    {
        imin = ij_botleft[0];
        imax = ij_topright[0];
        jmin = ij_botleft[1];
        jmax = ij_topright[1];    
    }
        
        
    glBegin(GL_QUADS);
    dy = (max - min)/((double)(jmax - jmin));
    dy_e = max/((double)(jmax - jmin));
    
    for (j = jmin; j < jmax; j++)
    {
        switch (plot) {
            case (P_AMPLITUDE):
            {
                value = min + 1.0*dy*(double)(j - jmin);
                color_scheme(COLOR_SCHEME, value, 1.0, 1, rgb);
                break;
            }
            case (P_ENERGY):
            {
                value = dy_e*(double)(j - jmin)*100.0/E_SCALE;
                if (COLOR_PALETTE >= COL_TURBO) color_scheme_asym(COLOR_SCHEME, value, 1.0, 1, rgb);
                else color_scheme(COLOR_SCHEME, value, 1.0, 1, rgb);
                break;
            }
            case (P_MEAN_ENERGY):
            {
                value = dy_e*(double)(j - jmin)*100.0/E_SCALE;
                if (COLOR_PALETTE >= COL_TURBO) color_scheme_asym(COLOR_SCHEME, value, 1.0, 1, rgb);
                else color_scheme(COLOR_SCHEME, value, 1.0, 1, rgb);
                break;
            }
            case (P_LOG_ENERGY):
            {
                value = LOG_SHIFT + LOG_SCALE*log(dy_e*(double)(j - jmin)*100.0/E_SCALE);
//                 if (value <= 0.0) value = 0.0;
                color_scheme(COLOR_SCHEME, value, 1.0, 1, rgb);
                break;
            }
            case (P_LOG_MEAN_ENERGY):
            {
                value = LOG_SHIFT + LOG_SCALE*log(dy_e*(double)(j - jmin)*100.0/E_SCALE);
//                 if (value <= 0.0) value = 0.0;
                color_scheme(COLOR_SCHEME, value, 1.0, 1, rgb);
                break;
            }
            case (P_PHASE):
            {
                value = min + 1.0*dy*(double)(j - jmin);
//                 lum = (color_amplitude(value, 1.0, 1))*0.5;
//                 if (lum < 0.0) lum = 0.0;
//                 hsl_to_rgb(value*360.0, 0.9, 0.5, rgb);
//                 color_scheme(COLOR_SCHEME, value, 1.0, 1, rgb);
//                 amp = color_amplitude_linear(value, 1.0, 1);
                amp = 0.5*color_amplitude_linear(value, 1.0, 1);
                while (amp > 1.0) amp -= 2.0;
                while (amp < -1.0) amp += 2.0;
                amp_to_rgb(0.5*(1.0 + amp), rgb);
                break;
            }
        }
        glColor3f(rgb[0], rgb[1], rgb[2]);
        if (ROTATE_COLOR_SCHEME)
        {
            glVertex2i(j, imin);
            glVertex2i(j, imax);
            glVertex2i(j+1, imax);
            glVertex2i(j+1, imin);            
        }
        else
        {
            glVertex2i(imin, j);
            glVertex2i(imax, j);
            glVertex2i(imax, j+1);
            glVertex2i(imin, j+1);
        }
    }
    glEnd ();
    
    glColor3f(1.0, 1.0, 1.0);
    glLineWidth(BOUNDARY_WIDTH);
    draw_rectangle(x1, y1, x2, y2);
}

void draw_color_scheme_palette(double x1, double y1, double x2, double y2, int plot, double min, double max, int palette)
{
    int j, k, ij_botleft[2], ij_topright[2], imin, imax, jmin, jmax;
    double y, dy, dy_e, rgb[3], value, lum, amp;
    
    xy_to_ij(x1, y1, ij_botleft);
    xy_to_ij(x2, y2, ij_topright);
    
    rgb[0] = 0.0;   rgb[1] = 0.0;   rgb[2] = 0.0;
//     erase_area_rgb(0.5*(x1 + x2), x2 - x1, 0.5*(y1 + y2), y2 - y1, rgb);

    if (ROTATE_COLOR_SCHEME)
    {
        jmin = ij_botleft[0];
        jmax = ij_topright[0];
        imin = ij_botleft[1];
        imax = ij_topright[1];    
    }
    else
    {
        imin = ij_botleft[0];
        imax = ij_topright[0];
        jmin = ij_botleft[1];
        jmax = ij_topright[1];    
    }
        
        
    glBegin(GL_QUADS);
    dy = (max - min)/((double)(jmax - jmin));
    dy_e = max/((double)(jmax - jmin));
    
    for (j = jmin; j < jmax; j++)
    {
        switch (plot) {
            case (P_AMPLITUDE):
            {
                value = min + 1.0*dy*(double)(j - jmin);
                color_scheme_palette(COLOR_SCHEME, palette, value, 1.0, 1, rgb);
                break;
            }
            case (P_ENERGY):
            {
                value = dy_e*(double)(j - jmin)*100.0/E_SCALE;
                if (COLOR_PALETTE >= COL_TURBO) color_scheme_asym_palette(COLOR_SCHEME, palette, value, 1.0, 1, rgb);
                else color_scheme_palette(COLOR_SCHEME, palette, value, 1.0, 1, rgb);
                break;
            }
            case (P_MEAN_ENERGY):
            {
                value = dy_e*(double)(j - jmin)*100.0/E_SCALE;
                if (COLOR_PALETTE >= COL_TURBO) color_scheme_asym_palette(COLOR_SCHEME, palette, value, 1.0, 1, rgb);
                else color_scheme_palette(COLOR_SCHEME, palette, value, 1.0, 1, rgb);
                break;
            }
            case (P_LOG_ENERGY):
            {
                value = LOG_SHIFT + LOG_SCALE*log(dy_e*(double)(j - jmin)*100.0/E_SCALE);
//                 if (value <= 0.0) value = 0.0;
                color_scheme_palette(COLOR_SCHEME, palette, value, 1.0, 1, rgb);
                break;
            }
            case (P_LOG_MEAN_ENERGY):
            {
                value = LOG_SHIFT + LOG_SCALE*log(dy_e*(double)(j - jmin)*100.0/E_SCALE);
//                 if (value <= 0.0) value = 0.0;
                color_scheme_palette(COLOR_SCHEME, palette, value, 1.0, 1, rgb);
                break;
            }
            case (P_PHASE):
            {
                value = min + 1.0*dy*(double)(j - jmin);
//                 lum = (color_amplitude(value, 1.0, 1))*0.5;
//                 if (lum < 0.0) lum = 0.0;
//                 hsl_to_rgb(value*360.0, 0.9, 0.5, rgb);
//                 color_scheme(COLOR_SCHEME, value, 1.0, 1, rgb);
//                 amp = color_amplitude_linear(value, 1.0, 1);
                amp = 0.5*color_amplitude_linear(value, 1.0, 1);
                while (amp > 1.0) amp -= 2.0;
                while (amp < -1.0) amp += 2.0;
                amp_to_rgb(0.5*(1.0 + amp), rgb);
                break;
            }
        }
        glColor3f(rgb[0], rgb[1], rgb[2]);
        if (ROTATE_COLOR_SCHEME)
        {
            glVertex2i(j, imin);
            glVertex2i(j, imax);
            glVertex2i(j+1, imax);
            glVertex2i(j+1, imin);            
        }
        else
        {
            glVertex2i(imin, j);
            glVertex2i(imax, j);
            glVertex2i(imax, j+1);
            glVertex2i(imin, j+1);
        }
    }
    glEnd ();
    
    glColor3f(1.0, 1.0, 1.0);
    glLineWidth(BOUNDARY_WIDTH);
    draw_rectangle(x1, y1, x2, y2);
}

void draw_color_scheme_palette_fade(double x1, double y1, double x2, double y2, int plot, double min, double max, int palette, int fade, double fade_value)
{
    int j, k, ij_botleft[2], ij_topright[2], imin, imax, jmin, jmax;
    double y, dy, dy_e, rgb[3], value, lum, amp;
    
    xy_to_ij(x1, y1, ij_botleft);
    xy_to_ij(x2, y2, ij_topright);
    
    rgb[0] = 0.0;   rgb[1] = 0.0;   rgb[2] = 0.0;
//     erase_area_rgb(0.5*(x1 + x2), x2 - x1, 0.5*(y1 + y2), y2 - y1, rgb);

    if (ROTATE_COLOR_SCHEME)
    {
        jmin = ij_botleft[0];
        jmax = ij_topright[0];
        imin = ij_botleft[1];
        imax = ij_topright[1];    
    }
    else
    {
        imin = ij_botleft[0];
        imax = ij_topright[0];
        jmin = ij_botleft[1];
        jmax = ij_topright[1];    
    }
        
        
    glBegin(GL_QUADS);
    dy = (max - min)/((double)(jmax - jmin));
    dy_e = max/((double)(jmax - jmin));
    
    for (j = jmin; j < jmax; j++)
    {
        switch (plot) {
            case (P_AMPLITUDE):
            {
                value = min + 1.0*dy*(double)(j - jmin);
                color_scheme_palette(COLOR_SCHEME, palette, value, 1.0, 1, rgb);
                break;
            }
            case (P_ENERGY):
            {
                value = dy_e*(double)(j - jmin)*100.0/E_SCALE;
                if (COLOR_PALETTE >= COL_TURBO) color_scheme_asym_palette(COLOR_SCHEME, palette, value, 1.0, 1, rgb);
                else color_scheme_palette(COLOR_SCHEME, palette, value, 1.0, 1, rgb);
                break;
            }
            case (P_MEAN_ENERGY):
            {
                value = dy_e*(double)(j - jmin)*100.0/E_SCALE;
                if (COLOR_PALETTE >= COL_TURBO) color_scheme_asym_palette(COLOR_SCHEME, palette, value, 1.0, 1, rgb);
                else color_scheme_palette(COLOR_SCHEME, palette, value, 1.0, 1, rgb);
                break;
            }
            case (P_LOG_ENERGY):
            {
                value = LOG_SHIFT + LOG_SCALE*log(dy_e*(double)(j - jmin)*100.0/E_SCALE);
//                 if (value <= 0.0) value = 0.0;
                color_scheme_palette(COLOR_SCHEME, palette, value, 1.0, 1, rgb);
                break;
            }
            case (P_LOG_MEAN_ENERGY):
            {
                value = LOG_SHIFT + LOG_SCALE*log(dy_e*(double)(j - jmin)*100.0/E_SCALE);
//                 if (value <= 0.0) value = 0.0;
                color_scheme_palette(COLOR_SCHEME, palette, value, 1.0, 1, rgb);
                break;
            }
            case (P_PHASE):
            {
                value = min + 1.0*dy*(double)(j - jmin);
//                 lum = (color_amplitude(value, 1.0, 1))*0.5;
//                 if (lum < 0.0) lum = 0.0;
//                 hsl_to_rgb(value*360.0, 0.9, 0.5, rgb);
//                 color_scheme(COLOR_SCHEME, value, 1.0, 1, rgb);
//                 amp = color_amplitude_linear(value, 1.0, 1);
                amp = 0.5*color_amplitude_linear(value, 1.0, 1);
                while (amp > 1.0) amp -= 2.0;
                while (amp < -1.0) amp += 2.0;
                amp_to_rgb(0.5*(1.0 + amp), rgb);
                break;
            }
        }
        if (fade) for (k=0; k<3; k++) rgb[k] *= fade_value;
        glColor3f(rgb[0], rgb[1], rgb[2]);
        if (ROTATE_COLOR_SCHEME)
        {
            glVertex2i(j, imin);
            glVertex2i(j, imax);
            glVertex2i(j+1, imax);
            glVertex2i(j+1, imin);            
        }
        else
        {
            glVertex2i(imin, j);
            glVertex2i(imax, j);
            glVertex2i(imax, j+1);
            glVertex2i(imin, j+1);
        }
    }
    glEnd ();
    
    if (fade) glColor3f(fade_value, fade_value, fade_value);
    else glColor3f(1.0, 1.0, 1.0);
    glLineWidth(BOUNDARY_WIDTH);
    draw_rectangle(x1, y1, x2, y2);
}


void print_speed(double speed, int fade, double fade_value)
{
    char message[100];
    double y = YMAX - 0.1, pos[2];
    static double xleftbox, xlefttext;
    static int first = 1;
    
    if (first)
    {
        xleftbox = XMIN + 0.3;
        xlefttext = xleftbox - 0.45;
        first = 0;
    }
    
    erase_area_hsl(xleftbox, y + 0.025, 0.22, 0.05, 0.0, 0.9, 0.0);
    if (fade) glColor3f(fade_value, fade_value, fade_value);
    else glColor3f(1.0, 1.0, 1.0);
    xy_to_pos(xlefttext + 0.28, y, pos);
    sprintf(message, "Mach %.3lg", speed);
    write_text(pos[0], pos[1], message);
}


void init_laplacian_coords(t_laplacian laplace[NX*NY], double phi[NX*NY])
/* compute coordinates of neighbours to compute Laplacian */
{
    int i, j, iplus, iminus, i1, i2, i3, j1, j2, j3, ij[2];;
    
    printf("Initialising Laplacian table\n");
    
    /* Laplacian in the bulk */
    #pragma omp parallel for private(i,j)
    for (i=1; i<NX-1; i++){
        for (j=1; j<NY-1; j++){
            laplace[i*NY+j].nneighb = 4;
            laplace[i*NY+j].nghb[0] = &phi[(i+1)*NY+j];
            laplace[i*NY+j].nghb[1] = &phi[(i-1)*NY+j];
            laplace[i*NY+j].nghb[2] = &phi[i*NY+j+1];
            laplace[i*NY+j].nghb[3] = &phi[i*NY+j-1];
        }
    }
    
    switch (B_COND) {
        case (BC_DIRICHLET):
        {
            /* left boundary */
            #pragma omp parallel for private(j)
            for (j=1; j<NY-1; j++)
            {
                laplace[j].nneighb = 3;
                laplace[j].nghb[0] = &phi[NY+j];
                laplace[j].nghb[1] = &phi[j+1];
                laplace[j].nghb[2] = &phi[j-1];
            }
            /* right boundary */
            #pragma omp parallel for private(j)
            for (j=1; j<NY-1; j++)
            {
                laplace[(NX-1)*NY+j].nneighb = 3;
                laplace[(NX-1)*NY+j].nghb[0] = &phi[(NX-2)*NY+j];
                laplace[(NX-1)*NY+j].nghb[1] = &phi[(NX-1)*NY+j+1];
                laplace[(NX-1)*NY+j].nghb[2] = &phi[(NX-1)*NY+j-1];
            }
            
            /* top boundary */
            #pragma omp parallel for private(i,iplus,iminus)
            for (i=0; i<NX; i++)
            {
                iplus = i+1;   if (iplus == NX) iplus = NX-1;
                iminus = i-1;  if (iminus == -1) iminus = 0;
                
                laplace[i*NY+NY-1].nneighb = 3;
                laplace[i*NY+NY-1].nghb[0] = &phi[iplus*NY+NY-1];
                laplace[i*NY+NY-1].nghb[1] = &phi[iminus*NY+NY-1];
                laplace[i*NY+NY-1].nghb[2] = &phi[i*NY+NY-2];
            }
            
            /* bottom boundary */
            #pragma omp parallel for private(i,iplus,iminus)
            for (i=0; i<NX; i++)
            {
                iplus = i+1;   if (iplus == NX) iplus = NX-1;
                iminus = i-1;  if (iminus == -1) iminus = 0;
                
                laplace[i*NY].nneighb = 3;
                laplace[i*NY].nghb[0] = &phi[iplus*NY];
                laplace[i*NY].nghb[1] = &phi[iminus*NY];
                laplace[i*NY].nghb[2] = &phi[i*NY+1];
            }
            break;
        }
        case (BC_PERIODIC):
        {
            /* left boundary */
            #pragma omp parallel for private(j)
            for (j=1; j<NY-1; j++)
            {
                laplace[j].nneighb = 4;
                laplace[j].nghb[0] = &phi[NY+j];
                laplace[j].nghb[1] = &phi[(NX-1)*NY+j];
                laplace[j].nghb[2] = &phi[j+1];
                laplace[j].nghb[3] = &phi[j-1];
            }
            
            /* right boundary */
            #pragma omp parallel for private(j)
            for (j=1; j<NY-1; j++)
            {
                laplace[(NX-1)*NY+j].nneighb = 4;
                laplace[(NX-1)*NY+j].nghb[0] = &phi[(NX-2)*NY+j];
                laplace[(NX-1)*NY+j].nghb[1] = &phi[j];
                laplace[(NX-1)*NY+j].nghb[2] = &phi[(NX-1)*NY+j+1];
                laplace[(NX-1)*NY+j].nghb[3] = &phi[(NX-1)*NY+j-1];
            }
            
            /* top boundary */
            #pragma omp parallel for private(i,iplus,iminus)
            for (i=0; i<NX; i++)
            {
                iplus = (i+1) % NX;
                iminus = (i-1) % NX;
                if (iminus < 0) iminus += NX;
                    
                laplace[i*NY+NY-1].nneighb = 4;
                laplace[i*NY+NY-1].nghb[0] = &phi[iplus*NY+NY-1];
                laplace[i*NY+NY-1].nghb[1] = &phi[iminus*NY+NY-1];
                laplace[i*NY+NY-1].nghb[2] = &phi[i*NY+NY-2];
                laplace[i*NY+NY-1].nghb[3] = &phi[i*NY];
            }
            
            /* bottom boundary */
            #pragma omp parallel for private(i,iplus,iminus)
            for (i=0; i<NX; i++)
            {
                iplus = (i+1) % NX;
                iminus = (i-1) % NX;
                if (iminus < 0) iminus += NX;
                    
                laplace[i*NY].nneighb = 4;
                laplace[i*NY].nghb[0] = &phi[iplus*NY];
                laplace[i*NY].nghb[1] = &phi[iminus*NY];
                laplace[i*NY].nghb[2] = &phi[i*NY+1];
                laplace[i*NY].nghb[3] = &phi[i*NY+NY-1];
            }
            break;
        }
        case (BC_ABSORBING):
        {
            /* left boundary */
            #pragma omp parallel for private(j)
            for (j=1; j<NY-1; j++)
            {
                laplace[j].nneighb = 3;
                laplace[j].nghb[0] = &phi[NY+j];
                laplace[j].nghb[1] = &phi[j+1];
                laplace[j].nghb[2] = &phi[j-1];
            }
            
            /* right boundary */
            #pragma omp parallel for private(j)
            for (j=1; j<NY-1; j++)
            {
                laplace[(NX-1)*NY+j].nneighb = 3;
                laplace[(NX-1)*NY+j].nghb[0] = &phi[(NX-2)*NY+j];
                laplace[(NX-1)*NY+j].nghb[1] = &phi[(NX-1)*NY+j+1];
                laplace[(NX-1)*NY+j].nghb[2] = &phi[(NX-1)*NY+j-1];
            }
            
            /* top boundary */
            #pragma omp parallel for private(i,iplus,iminus)
            for (i=0; i<NX; i++)
            {
                iplus = (i+1);   if (iplus == NX) iplus = NX-1;
                iminus = (i-1);  if (iminus == -1) iminus = 0;
                    
                laplace[i*NY+NY-1].nneighb = 3;
                laplace[i*NY+NY-1].nghb[0] = &phi[iplus*NY+NY-1];
                laplace[i*NY+NY-1].nghb[1] = &phi[iminus*NY+NY-1];
                laplace[i*NY+NY-1].nghb[2] = &phi[i*NY+NY-2];
            }
            
            /* bottom boundary */
            #pragma omp parallel for private(i,iplus,iminus)
            for (i=0; i<NX; i++)
            {
                iplus = (i+1);   if (iplus == NX) iplus = NX-1;
                iminus = (i-1);  if (iminus == -1) iminus = 0;
                    
                laplace[i*NY].nneighb = 3;
                laplace[i*NY].nghb[0] = &phi[iplus*NY];
                laplace[i*NY].nghb[1] = &phi[iminus*NY];
                laplace[i*NY].nghb[2] = &phi[i*NY+1];
            }
            break;            
        }
        case (BC_VPER_HABS):
        {
            /* left boundary */
            #pragma omp parallel for private(j)
            for (j=1; j<NY-1; j++)
            {
                laplace[j].nneighb = 3;
                laplace[j].nghb[0] = &phi[NY+j];
                laplace[j].nghb[1] = &phi[j+1];
                laplace[j].nghb[2] = &phi[j-1];
            }
            
            /* right boundary */
            #pragma omp parallel for private(j)
            for (j=1; j<NY-1; j++)
            {
                laplace[(NX-1)*NY+j].nneighb = 3;
                laplace[(NX-1)*NY+j].nghb[0] = &phi[(NX-2)*NY+j];
                laplace[(NX-1)*NY+j].nghb[1] = &phi[(NX-1)*NY+j+1];
                laplace[(NX-1)*NY+j].nghb[2] = &phi[(NX-1)*NY+j-1];
            }
            
            /* top boundary */
            #pragma omp parallel for private(i,iplus,iminus)
            for (i=0; i<NX; i++)
            {
                iplus = (i+1);   if (iplus == NX) iplus = NX-1;
                iminus = (i-1);  if (iminus == -1) iminus = 0;
                    
                laplace[i*NY+NY-1].nneighb = 4;
                laplace[i*NY+NY-1].nghb[0] = &phi[iplus*NY+NY-1];
                laplace[i*NY+NY-1].nghb[1] = &phi[iminus*NY+NY-1];
                laplace[i*NY+NY-1].nghb[2] = &phi[i*NY+NY-2];
                laplace[i*NY+NY-1].nghb[3] = &phi[i*NY];
            }
            
            /* bottom boundary */
            #pragma omp parallel for private(i,iplus,iminus)
            for (i=0; i<NX; i++)
            {
                iplus = (i+1);   if (iplus == NX) iplus = NX-1;
                iminus = (i-1);  if (iminus == -1) iminus = 0;
                    
                laplace[i*NY].nneighb = 4;
                laplace[i*NY].nghb[0] = &phi[iplus*NY];
                laplace[i*NY].nghb[1] = &phi[iminus*NY];
                laplace[i*NY].nghb[2] = &phi[i*NY+1];
                laplace[i*NY].nghb[3] = &phi[i*NY+NY-1];
            }
            break;            
        }
        case (BC_LSHAPE):
        {
            /* boundaries */
            xy_to_ij(-LAMBDA, -1.0, ij);
            i1 = ij[0] + 1;     j1 = ij[1] + 1;
            xy_to_ij(0.0, 0.0, ij);
            i2 = ij[0] - 1;     j2 = ij[1] - 1;
            xy_to_ij(LAMBDA, 1.0, ij);
            i3 = ij[0] - 1;     j3 = ij[1] - 1;
            
            printf("L shape corners (%i,%i), (%i,%i), (%i,%i)\n", i1, j1, i2, j2, i3, j3);
    
            /* left boundary */
            #pragma omp parallel for private(j)
            for (j=j1+1; j<j2; j++)
            {
                laplace[i1*NY+j].nneighb = 4;
                laplace[i1*NY+j].nghb[0] = &phi[(i1+1)*NY+j];
                laplace[i1*NY+j].nghb[1] = &phi[(i3)*NY+j];
                laplace[i1*NY+j].nghb[2] = &phi[i1*NY+j+1];
                laplace[i1*NY+j].nghb[3] = &phi[i1*NY+j-1];
            }
            #pragma omp parallel for private(j)
            for (j=j2; j<j3-1; j++)
            {
                laplace[i1*NY+j].nneighb = 4;
                laplace[i1*NY+j].nghb[0] = &phi[(i1+1)*NY+j];
                laplace[i1*NY+j].nghb[1] = &phi[(i2-1)*NY+j];
                laplace[i1*NY+j].nghb[2] = &phi[i1*NY+j+1];
                laplace[i1*NY+j].nghb[3] = &phi[i1*NY+j-1];
            }
            
            /* right boundary */
            #pragma omp parallel for private(j)
            for (j=j1+1; j<j2; j++)
            {
                laplace[(i3)*NY+j].nneighb = 4;
                laplace[(i3)*NY+j].nghb[0] = &phi[i1*NY+j];
                laplace[(i3)*NY+j].nghb[1] = &phi[(i3-1)*NY+j];
                laplace[(i3)*NY+j].nghb[2] = &phi[(i3)*NY+j+1];
                laplace[(i3)*NY+j].nghb[3] = &phi[(i3)*NY+j-1];
            }
            #pragma omp parallel for private(j)
            for (j=j2; j<j3-1; j++)
            {
                laplace[(i2)*NY+j].nneighb = 4;
                laplace[(i2)*NY+j].nghb[0] = &phi[i1*NY+j];
                laplace[(i2)*NY+j].nghb[1] = &phi[(i2-1)*NY+j];
                laplace[(i2)*NY+j].nghb[2] = &phi[(i2)*NY+j+1];
                laplace[(i2)*NY+j].nghb[3] = &phi[(i2)*NY+j-1];
            }
            
            /* top boundary */
            #pragma omp parallel for private(i)
            for (i=i1; i<i2; i++)
            {
                laplace[i*NY+j3].nneighb = 4;
                laplace[i*NY+j3].nghb[0] = &phi[(i+1)*NY+j3];
                laplace[i*NY+j3].nghb[1] = &phi[(i-1)*NY+j3];
                laplace[i*NY+j3].nghb[2] = &phi[i*NY+j1];
                laplace[i*NY+j3].nghb[3] = &phi[i*NY+j3-1];
            }
            #pragma omp parallel for private(i)
            for (i=i2; i<i3; i++)
            {
                laplace[i*NY+j2].nneighb = 4;
                laplace[i*NY+j2].nghb[0] = &phi[(i+1)*NY+j2];
                laplace[i*NY+j2].nghb[1] = &phi[(i-1)*NY+j2];
                laplace[i*NY+j2].nghb[2] = &phi[i*NY+j1];
                laplace[i*NY+j2].nghb[3] = &phi[i*NY+j2-1];
            }
            
            /* bottom boundary */
            #pragma omp parallel for private(i)
            for (i=i1; i<i2; i++)
            {
                laplace[i*NY+j1].nneighb = 4;
                laplace[i*NY+j1].nghb[0] = &phi[(i+1)*NY+j1];
                laplace[i*NY+j1].nghb[1] = &phi[(i-1)*NY+j1];
                laplace[i*NY+j1].nghb[2] = &phi[i*NY+j1+1];
                laplace[i*NY+j1].nghb[3] = &phi[i*NY+j3];
            }
            #pragma omp parallel for private(i)
            for (i=i2; i<i3; i++)
            {
                laplace[i*NY+j1].nneighb = 4;
                laplace[i*NY+j1].nghb[0] = &phi[(i+1)*NY+j1];
                laplace[i*NY+j1].nghb[1] = &phi[(i-1)*NY+j1];
                laplace[i*NY+j1].nghb[2] = &phi[i*NY+j1+1];
                laplace[i*NY+j1].nghb[3] = &phi[i*NY+j2];
            }
            
            /* corners */
            laplace[i1*NY+j1].nneighb = 4;
            laplace[i1*NY+j1].nghb[0] = &phi[(i1+1)*NY+j1];
            laplace[i1*NY+j1].nghb[1] = &phi[(i3-1)*NY+j1];
            laplace[i1*NY+j1].nghb[2] = &phi[i1*NY+j1+1];
            laplace[i1*NY+j1].nghb[3] = &phi[i1*NY+j3-1];

            laplace[(i3)*NY+j1].nneighb = 4;
            laplace[(i3)*NY+j1].nghb[0] = &phi[i1*NY+j1];
            laplace[(i3)*NY+j1].nghb[1] = &phi[(i3-1)*NY+j1];
            laplace[(i3)*NY+j1].nghb[2] = &phi[(i3)*NY+j1+1];
            laplace[(i3)*NY+j1].nghb[3] = &phi[(i3)*NY+j3];

            laplace[i1*NY+j3].nneighb = 4;
            laplace[i1*NY+j3].nghb[0] = &phi[(i1+1)*NY+j3];
            laplace[i1*NY+j3].nghb[1] = &phi[(i3)*NY+j3];
            laplace[i1*NY+j3].nghb[2] = &phi[i1*NY+j1];
            laplace[i1*NY+j3].nghb[3] = &phi[i1*NY+j3-1];

            break;
        }
    }
}


void compute_laplacian(double phi[NX*NY], t_laplacian laplace[NX*NY], double delta[NX*NY], short int xy_in[NX*NY])
/* compute the discretized Laplacian of phi */
{
    int i, j, k, n;
    
    /* in the bulk */
    if (B_COND == BC_LSHAPE)
    {
        #pragma omp parallel for private(i,j,k)
        for (i=1; i<NX-1; i++)
            for (j=1; j<NY-1; j++)
                if (xy_in[i*NY+j])
                {
                    delta[i*NY+j] = -4.0*phi[i*NY+j];
                    for (k=0; k<4; k++)
                        delta[i*NY+j] += *(laplace[i*NY+j].nghb[k]);                
                }
    }
    else
    {
        #pragma omp parallel for private(i,j,k)
        for (i=1; i<NX-1; i++)
            for (j=1; j<NY-1; j++)
                if (xy_in[i*NY+j])
                {
                    delta[i*NY+j] = phi[(i+1)*NY+j] + phi[(i-1)*NY+j] + phi[i*NY+j+1] + phi[i*NY+j-1] - 4.0*phi[i*NY+j];
                }
    }
            
    /* top and bottom boundaries */
    #pragma omp parallel for private(i,k,n)
    for (i=0; i<NX; i++)
    {
        if (xy_in[i*NY])
        {
            n = laplace[i*NY].nneighb;
            delta[i*NY] = -(double)n*phi[i*NY];
            for (k=0; k<n; k++)
                delta[i*NY] += *(laplace[i*NY].nghb[k]);
        }
        if (xy_in[i*NY+NY-1])
        {
            n = laplace[i*NY+NY-1].nneighb;
            delta[i*NY+NY-1] = -(double)n*phi[i*NY+NY-1];
            for (k=0; k<n; k++)
                delta[i*NY+NY-1] += *(laplace[i*NY+NY-1].nghb[k]);
        }
    }
    
    /* left and right boundaries */
    #pragma omp parallel for private(j,k,n)
    for (j=1; j<NY-1; j++)
    {
        if (xy_in[j])
        {
            n = laplace[j].nneighb;
            delta[j] = -(double)n*phi[j];
            for (k=0; k<n; k++)
                delta[j] += *(laplace[j].nghb[k]);
        }
         if (xy_in[(NX-1)*NY+j])
        {
            n = laplace[(NX-1)*NY+j].nneighb;
            delta[(NX-1)*NY+j] = -(double)n*phi[(NX-1)*NY+j];
            for (k=0; k<n; k++)
                delta[(NX-1)*NY+j] += *(laplace[(NX-1)*NY+j].nghb[k]);
        }
    }
}

double oscillating_bc(int time)
{
    double t, phase, a, envelope, omega;
    
    switch (OSCILLATION_SCHEDULE)
    {
        case (OSC_PERIODIC): 
        {
            return(AMPLITUDE*cos((double)time*OMEGA)*exp(-(double)time*DAMPING));
        }
        case (OSC_SLOWING):
        {
            a = 0.0025;
            t = (double)time*OMEGA;
            phase = t - a*t*t;
//             if (time%1000 == 0) printf("time = %i, phase = %.4lg\n", time, phase);
            return(AMPLITUDE*cos(phase)*exp(-phase*DAMPING));
        }
        case (OSC_WAVE_PACKET):
        {
            t = (double)time*OMEGA;
//             a = 10.0;
            a = 0.02/OMEGA;
            phase = AMPLITUDE*cos(t);
            envelope = exp(-(t-0.2)*(t-0.2)/(a*a))*sqrt(DPI/a);
            return(phase*envelope);
        }
        case (OSC_CHIRP):
        {
//             a = 0.25;
            t = (double)time*OMEGA;
            phase = t + ACHIRP*t*t;
            return(AMPLITUDE*sin(phase)*exp(-phase*DAMPING));
        }
    }
}

