/*********************/
/* Graphics routines */
/*********************/

#include "colors_waves.c"
#define TIFF_FREE_PERIOD 1

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
  free(image); /* prenvents RAM consumption*/
  TIFFClose(file);
  return 0;
}


int writetiff(char *filename, char *description, int x, int y, int width, int height, int compression)
{
  TIFF *file;
  GLubyte *image, *p;
  int i;
  static int counter = 0;

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
  /* added 9/9/22 and removed again, since it produces an unwanted "band" on the right */
  /* readded 5/11/22 */
  if (SAVE_MEMORY) free(image); /* prevents RAM consumption*/
//   {
//       counter++; 
//       if (counter%TIFF_FREE_PERIOD == 0)
//       {
//         free(image); /* prevents RAM consumption*/
//         counter = 0;
//       }
//   }
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

void init_hres(int res)		/* initialisation of window in higher resolution */
{
    glLineWidth(3);

    glClearColor(0.0, 0.0, 0.0, 1.0);
    glClear(GL_COLOR_BUFFER_BIT);

//     glOrtho(XMIN, XMAX, YMIN, YMAX , -1.0, 1.0);
    glOrtho(0.0, res*NX, 0.0, res*NY, -1.0, 1.0);
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
 
 double iabs(int i)     /* absolute value */
 {
	int res;

	if (i<0) res = -i;
	else res = i;
	return(i);
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


void xy_to_ij_safe(double x, double y, int ij[2])
/* convert (x,y) position to (i,j) in table representing wave, making sure (i,j) are between 0 and NX or NY */
{
    double x1, y1;

    x1 = (x - XMIN)/(XMAX - XMIN);
    y1 = (y - YMIN)/(YMAX - YMIN);

    ij[0] = (int)(x1 * (double)NX);
    ij[1] = (int)(y1 * (double)NY);
    
    if (ij[0] < 0) ij[0] = 0;
    if (ij[0] > NX-1) ij[0] = NX-1;
    if (ij[1] < 0) ij[1] = 0;
    if (ij[1] > NY-1) ij[1] = NY-1;
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


void xy_to_pos_hres(double x, double y, double pos[2])
/* convert (x,y) position to double-valued position in table representing wave */
{
    double x1, y1;

    x1 = (x - XMIN)/(XMAX - XMIN);
    y1 = (y - YMIN)/(YMAX - YMIN);

    pos[0] = x1*(double)(HRES*NX);
    pos[1] = y1*(double)(HRES*NY);
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

void erase_area_ij(int imin, int jmin, int imax, int jmax)
{
    glColor3f(0.0, 0.0, 0.0);
    glBegin(GL_QUADS);
    glVertex2i(imin, jmin);
    glVertex2i(imax, jmin);
    glVertex2i(imax, jmax);
    glVertex2i(imin, jmax);
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

void draw_vertex(double x, double y)
{
    double pos[2];
    
    xy_to_pos(x, y, pos);
    glVertex2d(pos[0], pos[1]);
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

void draw_line_hres(double x1, double y1, double x2, double y2)
{
    double pos[2];
    
    glBegin(GL_LINE_STRIP);
    xy_to_pos_hres(x1, y1, pos);
    glVertex2d(pos[0], pos[1]);
    xy_to_pos_hres(x2, y2, pos);
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

void draw_filled_rectangle(double x1, double y1, double x2, double y2)
{
    double pos[2];
    
    glBegin(GL_TRIANGLE_FAN);
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

void draw_circle_arc(double x, double y, double r, double angle1, double dangle, int nseg)
{
    int i;
    double pos[2], alpha, dalpha;
    
    dalpha = dangle/(double)nseg;
    
    glBegin(GL_LINE_STRIP);
    for (i=0; i<=nseg; i++)
    {
        alpha = angle1 + (double)i*dalpha;
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

double michelson_schedule(int i)
{
    double t;

    t = (double)(i - INITIAL_TIME)/(double)NSTEPS;
    return(4.0*t*WALL_WIDTH);
//     return(2.0*t*WALL_WIDTH);
}   


int generate_poisson_discs(t_circle circles[NMAXCIRCLES], double xmin, double xmax, double ymin, double ymax, double dpoisson)
/* generate a Poisson disc sample in a given rectangle */
{
    int i, j, k, n_p_active, ncandidates=5000, naccepted; 
    double r, phi, x, y;
    short int active_poisson[NMAXCIRCLES], far;
    
    printf("Generating Poisson disc sample\n");
    /* generate first circle */
    circles[0].xc = (xmax - xmin)*(double)rand()/RAND_MAX + xmin;
    circles[0].yc = (ymax - ymin)*(double)rand()/RAND_MAX + ymin;
    active_poisson[0] = 1;
    n_p_active = 1;
    ncircles = 1;
            
    while ((n_p_active > 0)&&(ncircles < NMAXCIRCLES))
    {
        /* randomly select an active circle */
        i = rand()%(ncircles);
        while (!active_poisson[i]) i = rand()%(ncircles);                 
        /* generate new candidates */
        naccepted = 0;
        for (j=0; j<ncandidates; j++)
        {
            r = dpoisson*(2.0*(double)rand()/RAND_MAX + 1.0);
            phi = DPI*(double)rand()/RAND_MAX;
            x = circles[i].xc + r*cos(phi);
            y = circles[i].yc + r*sin(phi);
            far = 1;
            for (k=0; k<ncircles; k++) if ((k!=i))
            {
                /* new circle is far away from circle k */
                far = far*((x - circles[k].xc)*(x - circles[k].xc) + (y - circles[k].yc)*(y - circles[k].yc) >= dpoisson*dpoisson);
                /* new circle is in domain */
                far = far*(x < xmax)*(x > xmin)*(y < ymax)*(y > ymin);
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
        }
        if (naccepted == 0)    /* inactivate circle i */ 
        {
            active_poisson[i] = 0;
            n_p_active--;
        }
        printf("%i active circles\n", n_p_active);
    }
            
    printf("Generated %i circles\n", ncircles);
    return(ncircles);
}


int init_circle_config_pattern(t_circle circles[NMAXCIRCLES], int circle_pattern)
/* initialise the arrays circlex, circley, circlerad and circleactive */
/* for billiard shape D_CIRCLES */
{
    int i, j, k, n, ncirc0, n_p_active, ncandidates=5000, naccepted; 
    double dx, dy, p, phi, r, r0, ra[5], sa[5], height, x, y = 0.0, gamma, dpoisson = PDISC_FACTOR*MU, xx[4], yy[4], dr, dphi;
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
                    circles[n].yc = YMIN + ((double)j + 0.5 + 0.5*((double)rand()/RAND_MAX - 0.5))*dy;
//                     circles[n].yc = YMIN + 0.5 + ((double)j + 0.5 + 0.5*((double)rand()/RAND_MAX - 0.5))*dy;
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
            
            for (i=0; i<2000; i++) 
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
            dr = (1.0 - RADIUS_FACTOR)*LAMBDA/(double)NGRIDY;
            for (i = 0; i < NGRIDX; i++)
                for (j = 0; j < NGRIDY; j++)
                {
                    n = NGRIDY*i + j;
                    phi = (double)i*dphi;
                    r = RADIUS_FACTOR*LAMBDA + (double)j*dr;
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
            dr = (1.0 - RADIUS_FACTOR)*LAMBDA/(double)NGRIDY;
            for (i = 0; i < NGRIDX; i++)
                for (j = 0; j < NGRIDY; j++)
                {
                    n = NGRIDY*i + j;
                    phi = (double)i*dphi;
                    phi += 0.5*(double)j*dphi;
                    r = RADIUS_FACTOR*LAMBDA + (double)j*dr;
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
//             r0 = 0.5*LAMBDA;
            r0 = RADIUS_FACTOR*LAMBDA;
            r = r0 + MU;
            
            for (i=0; i<(int)(1000.0/RADIUS_FACTOR); i++) 
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
        case (C_RINGS_POISSONDISC):
        {
            ncircles = generate_poisson_discs(circles, YMIN, YMAX, YMIN, YMAX, PDISC_FACTOR*MU);
            for (i=0; i<ncircles; i++)
            {
                circles[i].radius = MU;
                /* inactivate circles outside the domain */
                r = module2(circles[i].xc, circles[i].yc);
                if ((r < LAMBDA)&&(r > RADIUS_FACTOR*LAMBDA)) circles[i].active = 1;
                else circles[i].active = 0;
            }
            break;
        }
        case (C_RINGS_LOGSPIRAL):
        {
            ncircles = NGRIDX*NGRIDY;
            dr = (1.0 - RADIUS_FACTOR)*LAMBDA/(double)NGRIDY;
            dphi = DPI/((double)NGRIDX);
            for (i = 0; i < NGRIDX; i++)
                for (j = 0; j < NGRIDY; j++)
                {
                    n = NGRIDY*i + j;
                    r = RADIUS_FACTOR*LAMBDA + (double)j*dr;
                    phi = (double)i*dphi;
                    phi += log(r/(RADIUS_FACTOR*LAMBDA));
                    circles[n].xc = r*cos(phi);
                    circles[n].yc = r*sin(phi);
                    circles[n].radius = MU;
                    /* activate only circles that intersect the domain */
                    if ((circles[n].yc < YMAX + MU)&&(circles[n].yc > YMIN - MU)) circles[n].active = 1;
                    else circles[n].active = 0;
                }
            break;
        }
        case (C_HEX_BOTTOM):
        {
            ncircles = NGRIDX*(NGRIDY+1);
            dy = (- YMIN)/((double)NGRIDY);
//             dx = dy*0.5*sqrt(3.0);
            dx = (XMAX - XMIN)/((double)NGRIDX);
            for (i = 0; i < NGRIDX; i++)
                for (j = 0; j < NGRIDY+1; j++)
                {
                    n = (NGRIDY+1)*i + j;
                    circles[n].xc = ((double)(i-NGRIDX/2) + 0.5)*dx;   /* is +0.5 needed? */
                    circles[n].yc = YMIN + ((double)j - 0.5)*dy;
                    if ((i+NGRIDX)%2 == 1) circles[n].yc += 0.5*dy;
                    circles[n].radius = MU;
                    /* activate only circles that intersect the domain */
                    if ((circles[n].yc < YMAX + MU)&&(circles[n].yc > YMIN - MU)) circles[n].active = 1;
                    else circles[n].active = 0;
                }
            break;
        }
        case (C_HEX_BOTTOM2):
        {
            ncircles = NGRIDX*(NGRIDY+1);
            dy = (- YMIN)/((double)NGRIDY);
            dx = 2.0*LAMBDA/((double)NGRIDX);
            for (i = 0; i < NGRIDX; i++)
                for (j = 0; j < NGRIDY+1; j++)
                {
                    n = (NGRIDY+1)*i + j;
                    circles[n].xc = ((double)(i-NGRIDX/2) + 0.5)*dx;   /* is +0.5 needed? */
                    circles[n].yc = YMIN + ((double)j - 0.5)*dy;
                    if ((i+NGRIDX)%2 == 1) circles[n].yc += 0.5*dy;
                    circles[n].radius = MU;
                    /* activate only circles that intersect the domain */
                    if ((circles[n].yc < YMAX + MU)&&(circles[n].yc > YMIN - MU)) circles[n].active = 1;
                    else circles[n].active = 0;
                }
            break;
        }
        case (C_SQUARE_BOTTOM):
        {
            ncircles = NGRIDX*(NGRIDY+1);
            dy = (- YMIN)/((double)NGRIDY);
            dx = 2.0*LAMBDA/((double)NGRIDX);
            for (i = 0; i < NGRIDX; i++)
                for (j = 0; j < NGRIDY+1; j++)
                {
                    n = (NGRIDY+1)*i + j;
                    circles[n].xc = ((double)(i-NGRIDX/2) + 0.5)*dx;  
                    circles[n].yc = YMIN + ((double)j - 0.5)*dy;
                    circles[n].radius = MU;
                    /* activate only circles that intersect the domain */
                    if ((circles[n].yc < YMAX + MU)&&(circles[n].yc > YMIN - MU)) circles[n].active = 1;
                    else circles[n].active = 0;
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

int compute_maze_coordinates(t_rectangle polyrect[NMAXPOLY], int type)
/* compute positions of maze */
{
    t_maze* maze;
    int i, j, n, nsides = 0, ropening;
    double dx, dy, x1, y1, x0, padding = 0.02, pos[2], width = MAZE_WIDTH;
    
    maze = (t_maze *)malloc(NXMAZE*NYMAZE*sizeof(t_maze));
    
    ropening = (NYMAZE+1)/2;
    
    init_maze(maze);
    
    /* move the entrance for maze type with two channels */
    if (type == 2)
    {
        n = nmaze(0, ropening-1);
        maze[n].west = 0;
        n = nmaze(0, ropening);
        maze[n].west = 1;
    }
            
    /* build walls of maze */
//     x0 = LAMBDA - 1.0;
    dx = (YMAX - YMIN - 2.0*padding)/(double)(NXMAZE);
    dy = (YMAX - YMIN - 2.0*padding)/(double)(NYMAZE);
            
    for (i=0; i<NXMAZE; i++)
        for (j=0; j<NYMAZE; j++)
        {
            n = nmaze(i, j);
            x1 = YMIN + padding + (double)i*dx + MAZE_XSHIFT;
            y1 = YMIN + padding + (double)j*dy;
            
            if (((i>0)||(j!=ropening))&&(maze[n].west)) 
            {
                polyrect[nsides].x1 = x1 - width;
                polyrect[nsides].y1 = y1 - width;
                polyrect[nsides].x2 = x1 + width;
                polyrect[nsides].y2 = y1 + width + dy;
                nsides++;
            }
            
            if (maze[n].south) 
            {
                polyrect[nsides].x1 = x1 - width;
                polyrect[nsides].y1 = y1 - width;
                polyrect[nsides].x2 = x1 + width + dx;
                polyrect[nsides].y2 = y1 + width;
                nsides++;
            }
        }
    
    /* top side of maze */
    polyrect[nsides].x1 = YMIN + padding + MAZE_XSHIFT;
    polyrect[nsides].y1 = YMAX - padding - width;
    polyrect[nsides].x2 = YMAX - padding + MAZE_XSHIFT;
    polyrect[nsides].y2 = YMAX - padding + width;
    nsides++;
    
    /* right side of maze */
    y1 = YMIN + padding + dy*((double)ropening);
    x1 = YMAX - padding + MAZE_XSHIFT;
    polyrect[nsides].x1 = x1 - width;
    polyrect[nsides].y1 = YMIN - 1.0;
    polyrect[nsides].x2 = x1 + width;
    polyrect[nsides].y2 = y1 - dy;
    nsides++;
            
    polyrect[nsides].x1 = x1 - width;
    polyrect[nsides].y1 = y1;
    polyrect[nsides].x2 = x1 + width;
    polyrect[nsides].y2 = YMAX + 1.0;
    nsides++;
    
    /* left side of maze */
    x1 = YMIN + padding + MAZE_XSHIFT;
    polyrect[nsides].x1 = x1 - width;
    polyrect[nsides].y1 = YMIN - 1.0;
    polyrect[nsides].x2 = x1 + width;
    polyrect[nsides].y2 = YMIN + padding;
    nsides++;
            
    polyrect[nsides].x1 = x1 - width;
    polyrect[nsides].y1 = YMAX - padding;
    polyrect[nsides].x2 = x1 + width;
    polyrect[nsides].y2 = YMAX + 1.0;
    nsides++;
    
    if (type == 1)     /* maze with closed sides */
    {
        polyrect[nsides].x1 = XMIN - 0.5*width;
        polyrect[nsides].y1 = YMIN - 0.5*width;
        polyrect[nsides].x2 = XMIN + 0.5*width;
        polyrect[nsides].y2 = YMAX + 0.5*width;
        nsides++;
        
        polyrect[nsides].x1 = XMIN - 0.5*width;
        polyrect[nsides].y1 = YMIN - 0.5*width;
        polyrect[nsides].x2 = x1 + 0.5*width;
        polyrect[nsides].y2 = YMIN + 0.5*width;
        nsides++;

        polyrect[nsides].x1 = XMIN - 0.5*width;
        polyrect[nsides].y1 = YMAX - 0.5*width;
        polyrect[nsides].x2 = x1 + 0.5*width;
        polyrect[nsides].y2 = YMAX + 0.5*width;
        nsides++;
    }
    
    else if (type == 2)   /* maze with channels */
    {
        /* right channel */
        y1 = YMIN + padding + dy*((double)ropening);
        x1 = YMAX - padding + MAZE_XSHIFT;
        polyrect[nsides].x1 = x1 - 0.5*width;
        polyrect[nsides].y1 = YMIN - padding;
        polyrect[nsides].x2 = XMAX + padding;
        polyrect[nsides].y2 = y1 - dy + width;
        nsides++;
    
        polyrect[nsides].x1 = x1 - 0.5*width;
        polyrect[nsides].y1 = y1 - width;
        polyrect[nsides].x2 = XMAX + padding;
        polyrect[nsides].y2 = YMAX + padding;
        nsides++;
        
        /* left channel */
        x1 = YMIN + padding + MAZE_XSHIFT;
        polyrect[nsides].x1 = XMIN - padding;
        polyrect[nsides].y1 = YMIN - padding;
        polyrect[nsides].x2 = x1 + 0.5*width;
        polyrect[nsides].y2 = y1 - dy + width;
        nsides++;
    
        polyrect[nsides].x1 = XMIN - padding;
        polyrect[nsides].y1 = y1 - width;
        polyrect[nsides].x2 = x1 + 0.5*width;
        polyrect[nsides].y2 = YMAX + padding;
        nsides++;
        
    }
    
    for (i=0; i<nsides; i++)
    {
        xy_to_pos(polyrect[i].x1, polyrect[i].y1, pos);        
        polyrect[i].posi1 = pos[0];
        polyrect[i].posj1 = pos[1];
        xy_to_pos(polyrect[i].x2, polyrect[i].y2, pos);        
        polyrect[i].posi2 = pos[0];
        polyrect[i].posj2 = pos[1];        
    }
    
    free(maze);
    return(nsides);
}

int compute_circular_maze_coordinates(t_rect_rotated polyrectrot[NMAXPOLY], t_arc polyarc[NMAXPOLY], int *npolyrect_rot, int *npolyarc)
/* compute positions of circular maze */
{
    int nblocks, block, i, j, n, p, q, np, na;
    double rmin, rmax, angle, r, dr, phi, dphi, ww, width = 0.02; 
    t_maze* maze;
    
    maze = (t_maze *)malloc(NXMAZE*NYMAZE*sizeof(t_maze));
    
    init_circular_maze(maze);
    
    np = 0;
    na = 0;
    
    /* build walls of maze */
    nblocks = NYMAZE/NXMAZE;
    rmin = 0.15;
    rmax = 1.0;
    angle = DPI/(double)nblocks;
        
    dr = (rmax - rmin)/(double)(NXMAZE);
    
    /* add straight walls */
    for (block = 0; block < nblocks; block++)
    {
        dphi = angle;
        
        /* first circle */
        n = nmaze(0, block*NXMAZE);
        r = rmin - 0.5*width;
        phi = (double)block*angle;
            
        if (maze[n].south)
        {
            polyrectrot[np].x1 = r*cos(phi) + MAZE_XSHIFT;
            polyrectrot[np].y1 = r*sin(phi);
            polyrectrot[np].x2 = (r+dr+width)*cos(phi) + MAZE_XSHIFT;
            polyrectrot[np].y2 = (r+dr+width)*sin(phi);
            polyrectrot[np].width = width;
            np++;
        }
                
        /* second circle */
        r = rmin + dr - 0.5*width;
        dphi *= 0.5;
        for (q=0; q<2; q++)
        {
            n = nmaze(1, block*NXMAZE + q);
            phi = (double)(block)*angle + (double)q*dphi;
            
            if (maze[n].south)
            {
                polyrectrot[np].x1 = r*cos(phi) + MAZE_XSHIFT;
                polyrectrot[np].y1 = r*sin(phi);
                polyrectrot[np].x2 = (r+dr+width)*cos(phi) + MAZE_XSHIFT;
                polyrectrot[np].y2 = (r+dr+width)*sin(phi);
                polyrectrot[np].width = width;
                np++;
            }
        }
                
        /* other circles */
        ww = 2;
        i = 2;
        while (ww < NXMAZE)
        {
            dphi *= 0.5;
            for (p = 0; p < ww; p++)
            {
                r = rmin + (double)i*dr - 0.5*width;
//                 printf("Segment, i = %i, dphi = %.2lg, r = %.2lg\n", i, dphi, r);
                for (q = 0; q < 2*ww; q++)
                {
                    j = block*NXMAZE + q;
                    n = nmaze(i,j);
                    phi = (double)(block)*angle + (double)q*dphi;
                    
                    if (maze[n].south)
                    {
                        polyrectrot[np].x1 = r*cos(phi) + MAZE_XSHIFT;
                        polyrectrot[np].y1 = r*sin(phi);
                        polyrectrot[np].x2 = (r+dr+width)*cos(phi) + MAZE_XSHIFT;
                        polyrectrot[np].y2 = (r+dr+width)*sin(phi);
                        polyrectrot[np].width = width;
                        np++;
                    }
                }
                i++;
            }
            ww *= 2;
        }
                
    }
    
    /* add circular arcs */
    for (block = 0; block < nblocks; block++)
    {
        dphi = angle;
        
        /* first circle */
        n = nmaze(0, block*NXMAZE);
        r = rmin;
        phi = (double)block*angle;
        
        if ((block > 0)&&(maze[n].west))
        {
            polyarc[na].xc = MAZE_XSHIFT;
            polyarc[na].yc = 0.0;
            polyarc[na].r = r;
            polyarc[na].angle1 = phi;
            polyarc[na].dangle = dphi;
            polyarc[na].width = width;
            na++;
        }
                
        /* second circle */
        r = rmin + dr;
        dphi *= 0.5;
        for (q=0; q<2; q++)
        {
            n = nmaze(1, block*NXMAZE + q);
            phi = (double)(block)*angle + (double)q*dphi;
            
            if (maze[n].west)
            {
                polyarc[na].xc = MAZE_XSHIFT;
                polyarc[na].yc = 0.0;
                polyarc[na].r = r;
                polyarc[na].angle1 = phi;
                polyarc[na].dangle = dphi;
                polyarc[na].width = width;
                na++;
            }
        }
                
        /* other circles */
        ww = 2;
        i = 2;
        while (ww < NXMAZE)
        {
            dphi *= 0.5;
            for (p = 0; p < ww; p++)
            {
                r = rmin + (double)i*dr;
                printf("Circle, i = %i, dphi = %.2lg, r = %.2lg\n", i, dphi, r);
                for (q = 0; q < 2*ww; q++)
                {
                    j = block*NXMAZE + q;
                    n = nmaze(i,j);
                    phi = (double)(block)*angle + (double)q*dphi;
                    
                    if (maze[n].west)
                    {
                        polyarc[na].xc = MAZE_XSHIFT;
                        polyarc[na].yc = 0.0;
                        polyarc[na].r = r;
                        polyarc[na].angle1 = phi;
                        polyarc[na].dangle = dphi;
                        polyarc[na].width = width;
                        na++;
                    }
                }
                i++;
            }
            ww *= 2;
        }
    }
    
    /* outer boundary of maze */
    polyarc[na].xc = MAZE_XSHIFT;
    polyarc[na].yc = 0.0;
    polyarc[na].r = rmax;
    polyarc[na].angle1 = dphi;
    polyarc[na].dangle = DPI - dphi;
    polyarc[na].width = width;
    na++;
    
    *npolyrect_rot = np;
    *npolyarc = na;
    
    free(maze);
}

int compute_interior_maze_coordinates(t_rectangle polyrect[NMAXPOLY], int type)
/* compute positions of complement of maze */
{
    t_maze* maze;
    int i, j, n, nsides = 0, ropening;
    double dx, dy, x1, y1, x0, padding = 0.02, pos[2], width = MAZE_WIDTH;
    
    /* TODO */
    
    maze = (t_maze *)malloc(NXMAZE*NYMAZE*sizeof(t_maze));
    
    ropening = (NYMAZE+1)/2;
    
    init_maze(maze);
    
    /* move the entrance for maze type with two channels */
    if (type == 2)
    {
        n = nmaze(0, ropening-1);
        maze[n].west = 0;
        n = nmaze(0, ropening);
        maze[n].west = 1;
    }
            
    /* build walls of maze */
//     x0 = LAMBDA - 1.0;
    dx = (YMAX - YMIN - 2.0*padding)/(double)(NXMAZE);
    dy = (YMAX - YMIN - 2.0*padding)/(double)(NYMAZE);
            
    for (i=0; i<NXMAZE; i++)
        for (j=0; j<NYMAZE; j++)
        {
            n = nmaze(i, j);
            x1 = YMIN + padding + (double)i*dx + MAZE_XSHIFT;
            y1 = YMIN + padding + (double)j*dy;
            
            if (!maze[n].west) 
            {
                polyrect[nsides].x1 = x1 - 0.5*dx - width;
                polyrect[nsides].y1 = y1 + 0.5*dx - width;
                polyrect[nsides].x2 = x1 + 0.5*dx + width;
                polyrect[nsides].y2 = y1 + 0.5*dx + width;
                nsides++;
            }
            
            if (!maze[n].south) 
            {
                polyrect[nsides].x1 = x1 + 0.5*dx - width;
                polyrect[nsides].y1 = y1 - 0.5*dx - width;
                polyrect[nsides].x2 = x1 + 0.5*dx + width;
                polyrect[nsides].y2 = y1 + 0.5*dx + width;
                nsides++;
            }
        }
    
    
    /* right channel of maze */
    y1 = YMIN + padding + dy*((double)ropening);
    x1 = YMAX - padding + MAZE_XSHIFT;
    polyrect[nsides].x1 = x1 - 0.5*dx - width;
    polyrect[nsides].y1 = y1 - 0.5*dy - width;
    polyrect[nsides].x2 = XMAX;
    polyrect[nsides].y2 = y1 - 0.5*dy + width;
    nsides++;
                
    /* left channel of maze */
    x1 = YMIN + padding + MAZE_XSHIFT;
    polyrect[nsides].x1 = XMIN;
    polyrect[nsides].y1 = y1 - 0.5*dy - width;
    polyrect[nsides].x2 = x1 + 0.5*dx + width;
    polyrect[nsides].y2 = y1 - 0.5*dy + width;
    nsides++;
    
    for (i=0; i<nsides; i++)
    {
        xy_to_pos(polyrect[i].x1, polyrect[i].y1, pos);        
        polyrect[i].posi1 = pos[0];
        polyrect[i].posj1 = pos[1];
        xy_to_pos(polyrect[i].x2, polyrect[i].y2, pos);        
        polyrect[i].posi2 = pos[0];
        polyrect[i].posj2 = pos[1];        
    }
    
    free(maze);
    return(nsides);
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

int init_polyrect(t_rectangle polyrect[NMAXPOLY])
/* initialise variable polyrect, for certain polygonal domain shapes */
{
    switch (B_DOMAIN) {
        case (D_MAZE):
        {
            return(compute_maze_coordinates(polyrect, 0));
        }
        case (D_MAZE_CLOSED):
        {
            return(compute_maze_coordinates(polyrect, 1));
        }
        case (D_MAZE_CHANNELS):
        {
            return(compute_maze_coordinates(polyrect, 2));
        }
        default:
        {
            if ((ADD_POTENTIAL)&&(POTENTIAL == POT_MAZE)) return(compute_maze_coordinates(polyrect, 1));
            return(0);
        }
    }
}

int init_polyrect_euler(t_rectangle polyrect[NMAXPOLY], int domain)
/* initialise variable polyrect, for certain polygonal domain shapes */
{
    switch (domain) {
        case (D_MAZE):
        {
            return(compute_maze_coordinates(polyrect, 0));
        }
        case (D_MAZE_CLOSED):
        {
            return(compute_maze_coordinates(polyrect, 1));
        }
        case (D_MAZE_CHANNELS):
        {
            return(compute_maze_coordinates(polyrect, 2));
        }
        case (D_MAZE_CHANNELS_INT):
        {
            return(compute_interior_maze_coordinates(polyrect, 2));
        }
     }
}


void init_polyrect_arc(t_rect_rotated polyrectrot[NMAXPOLY], t_arc polyarc[NMAXPOLY], int *npolyrect, int *npolyarc)
/* initialise variables polyrectrot and polyarc, for certain domain shapes */
{
    switch (B_DOMAIN) {
        case (D_MAZE_CIRCULAR):
        {
            compute_circular_maze_coordinates(polyrectrot, polyarc, npolyrect, npolyarc);
            break;
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


int xy_in_polyrect(double x, double y, t_rectangle rectangle)
/* returns 1 if (x,y) is in rectangle */
{
    double x1, y1, x2, y2;
    
    if (rectangle.x1 < rectangle.x2) 
    {
        x1 = rectangle.x1;
        x2 = rectangle.x2;
    }
    else 
    {
        x1 = rectangle.x2;
        x2 = rectangle.x1;        
    }
    if (rectangle.y1 < rectangle.y2) 
    {
        y1 = rectangle.y1;
        y2 = rectangle.y2;
    }
    else 
    {
        y1 = rectangle.y2;
        y2 = rectangle.y1;        
    }
    if (x < x1) return(0);
    if (x > x2) return(0);
    if (y < y1) return(0);
    if (y > y2) return(0);
    return(1);
}


int ij_in_polyrect(double i, double j, t_rectangle rectangle)
/* returns 1 if (x,y) is in rectangle */
{
    int i1, i2, j1, j2;
    
    if (rectangle.posi1 < rectangle.posi2) 
    {
        i1 = rectangle.posi1;
        i2 = rectangle.posi2;
    }
    else 
    {
        i1 = rectangle.posi2;
        i2 = rectangle.posi1;        
    }
    if (rectangle.posj1 < rectangle.posj2) 
    {
        j1 = rectangle.posj1;
        j2 = rectangle.posj2;
    }
    else 
    {
        j1 = rectangle.posj2;
        j2 = rectangle.posj1;        
    }
    if (i < i1) return(0);
    if (i > i2) return(0);
    if (j < j1) return(0);
    if (j > j2) return(0);
    return(1);
}

int xy_in_rectrotated(double x, double y, t_rect_rotated rectrot)
/* returns 1 if (x,y) is in rectangle */
{
    double l, u1, u2, v1, v2, pscal, h2;
    
    l = module2(rectrot.x2 - rectrot.x1, rectrot.y2 - rectrot.y1);
    if (l == 0.0) return(0);
    
    /* unit vector along axis */
    u1 = (rectrot.x2 - rectrot.x1)/l;
    u2 = (rectrot.y2 - rectrot.y1)/l;
    
    /* vector from one extremity to (x,y) */
    v1 = x - rectrot.x1;
    v2 = y - rectrot.y1;
    
    /* inner product */
    pscal = u1*v1 + u2*v2;
    if (pscal < 0.0) return(0);
    if (pscal > l) return(0);
    
    h2 = v1*v1 + v2*v2 - pscal*pscal;
    return(4.0*h2 <= rectrot.width*rectrot.width);
}

int xy_in_arc(double x, double y, t_arc arc)
/* returns 1 if (x,y) is in arc */
{
    double rho, phi, alpha;
    
    rho = module2(x - arc.xc, y - arc.yc);
    
    if (vabs(rho - arc.r) > 0.5*arc.width) return(0);
    
    phi = argument(x - arc.xc, y - arc.yc);
    
    alpha = phi - arc.angle1;
    while (alpha < 0.0) alpha += DPI;
    while (alpha > DPI) alpha -= DPI;
    
    return(alpha <= arc.dangle);
}


int rc_hyp(double x, double y)
/* xy_in for D_RITCHEY_CHRETIEN_HYPERBOLIC domain */
{
    static int first = 1;
    static double m, d, dprime, b, r1, r2, e1, e2, a1, b1, a2, b2, x01, x02;
    double y1, f, f1, f2;
    
    if (first)
    {
        m = LAMBDA;     /* secondary magnification */
        d = 2.5;        /* distance between mirrors */
        b = 3.0;        /* distance between secondary mirror and effective focal point */
        
        f = m*d + b;
        f1 = f/m;       /* focal distance of primary mirror */
        dprime = f1 - d;    /* distance between secondary mirror and primary focal point */
        f2 = m*dprime/(m+1); /* focal distance of secondary mirror */
        
        r1 = 2.0*f/m;      /* radius of curvature of primary mirror */
        r2 = 2.0*b/(m - 1.0);    /* radius of curvature of secondary mirror */
        e1 = sqrt(1.0 + 2.0*b/(m*m*m*d));    /* eccentricity of primary mirror */ 
        e2 = sqrt(1.0 + 2.0*(m*(2.0*m - 1.0) + b/d)/((m-1)*(m-1)*(m-1)));    
                                            /* eccentricity of secondary mirror */
        a1 = f1/e1;
        b1 = sqrt(a1*r1);       /* semi-axes of primary mirror */
        x01 = 0.5*d + a1;       /* center of primary hyperbola */
        
        a2 = f2/e2;
        b2 = sqrt(a2*r2);       /* semi-axes of secondary mirror */
        x02 = -0.5*d + a2;       /* center of secondary hyperbola */
                
        first = 0;
    }
    
    
    y1 = vabs(y);  
            
    if (x > 0.0)    /* primary mirror */   
    {
        if (y1 < MU) return(1);
        if (x > 0.5*d + WALL_WIDTH) return(1);
        return((x-x01)*(x-x01) > a1*a1*(1.0 + y*y/(b1*b1)));
    }
    else    /* secondary mirror */
    {
        if (y1 > 0.3) return(1);
//         if (x < -0.5*d - WALL_WIDTH) return(1);
        return((x > x02)||((x-x02)*(x-x02) < a2*a2*(1.0 + y*y/(b2*b2))));
    }
}


int init_xyin_from_image(short int * xy_in[NX])
/* initialize table xy_in from an image file */
{
    FILE *image_file;
    int nx, ny, maxrgb, i, j, k, ii, jj, nmaxpixels = 2000000, rgbtot, scan, rgbval;
    int *rgb_values;
    double scalex, scaley;

    image_file = fopen("PHOTONSband.ppm", "r");
    scan = fscanf(image_file,"%i %i\n", &nx, &ny);
    scan = fscanf(image_file,"%i\n", &maxrgb);
    
    rgb_values = (int *)malloc(3*nmaxpixels*sizeof(int));
    
    scalex = (double)nx/(double)NX;
    scaley = (double)ny/(double)NY;
    
    if (nx*ny > nmaxpixels)
    {
        printf("Image too large, increase nmaxpixels in init_xyin_from_image()\n");
        exit(0);
    }
    
    /* read rgb values */
    for (j=0; j<ny; j++)
        for (i=0; i<nx; i++)
            for (k=0; k<3; k++)
                {
                    scan = fscanf(image_file,"%i\n", &rgbval);
                    rgb_values[3*(j*nx+i)+k] = rgbval;
                }
                
    for (i=0; i<NX; i++)
        for (j=0; j<NY; j++)
        {
            ii = nx/2 + (int)((double)(i-NX/2)*scalex);
            jj = ny/2 - (int)((double)(j-NY/2)*scalex);
            if (ii > nx-1) ii = nx-1;
            if (ii < 0) ii = 0;
            if (jj > ny-1) jj = ny-1;
            if (jj < 0) jj = 0;
            k = 3*(jj*nx+ii);
            rgbtot = rgb_values[k] + rgb_values[k+1] + rgb_values[k+2];
            if (rgbtot > 3*maxrgb/2) xy_in[i][j] = 1;
            else xy_in[i][j] = 0;
        }   
    
    fclose(image_file);
    free(rgb_values);
    return(1);
}


int xy_in_billiard_single_domain(double x, double y, int b_domain, int ncirc, t_circle *circles)
/* returns 1 if (x,y) represents a point in the billiard */
{
    double l2, r1, r2, r2mu, omega, b, c, d, angle, z, x1, y1, x2, y2, y3, u, v, u1, v1, dx, dy, width, alpha, s, a, r, height, ca, sa, l, ht, xshift, zz[9][2], x5, x6, f, fp, h1, cb, sb, c1, c2;
    int i, j, k, k1, k2, condition = 0, m;
    static int first = 1, nsides;
    static double h, hh, ra, rb, ll, salpha;

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
        case (D_EXT_ELLIPSE):
        {
            if (x*x/(LAMBDA*LAMBDA) + y*y/(MU*MU) > 1.0) return(1);
            else return(0);
            break;
        }
        case (D_EXT_ELLIPSE_CURVED):
        {
            y1 = y + 0.4*x*x;
            if (x*x/(LAMBDA*LAMBDA) + y1*y1/(MU*MU) > 1.0) return(1);
            else return(0);
            break;
        }
        case (D_EXT_ELLIPSE_CURVED_BDRY):
        {
            if (y > YMAX - 0.05) return(0);
            if (y < YMIN + 0.05) return(0);
            y1 = y + 0.4*x*x;
            if (x*x/(LAMBDA*LAMBDA) + y1*y1/(MU*MU) > 1.0) return(1);
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
                angle = APOLY*PID + ((double)k+0.5)*omega;
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
        case (D_WAVEGUIDE):
        {
            x1 = XMIN + MU;
            x2 = XMAX - 2.0*MU - 1.5*LAMBDA;
            y1 = 0.5*LAMBDA;
            y2 = 1.5*LAMBDA;
            if (x < x1) return(0);
            if (x > x2 + 1.5*LAMBDA) return(0);
            if (vabs(y) > y2) return(0);
            if (x < x2) return(vabs(y) >= y1);
            r = module2(x-x2, y);
            if (r < 0.5*LAMBDA) return(0);
            if (r > 1.5*LAMBDA) return(0);
            return(1);
        }
        case (D_WAVEGUIDE_W):
        {
            x1 = vabs(x);
            width = LAMBDA - 2.0*MU;
            height = 0.5*MU;
            if (x1 > 2.0*LAMBDA - MU) return(0);
            if (vabs(y) > MU + width + height) return(0);
            if (y >= height)
            {
                r = module2(x1, y-height);
                if ((r > MU)&&(r < MU + width)) return(1);
                if (x1 > LAMBDA + MU) return(1);
                return(0);
            }
            if (y <= -height)
            {
                r = module2(x1-LAMBDA, y+height);
                if ((r > MU)&&(r < MU + width)) return(1);
                return(0);
            }
            if (x1 > LAMBDA + MU) return(1);
            if ((x1 > MU)&&(x1 < MU + width)) return(1);
            return(0);
        }   
        case (D_WAVEGUIDES_W):
        {
            x1 = vabs(x);
            
            /* upper fiber */
            width = LAMBDA - 2.0*MU;
            height = MU + width;
            if (y > MU + width + height) return(0);
            if (y >= height)
            {
                r = module2(x1, y-height);
                if ((r > MU)&&(r < MU + width)) return(1);
                if ((y < height + MU)&&(x1 > LAMBDA + MU)&&(x1 < LAMBDA + MU + width)) return(1);
                return(0);
            }
            else
            {
                r = module2(x1-LAMBDA, y-height);
                if ((r > MU)&&(r < MU + width)) return(1);
            }
            
            /* lower fiber */
            height = -MU - width;
            width = LAMBDA - 2.0*MU_B;
            if (y < -MU_B - width + height) return(0);
            if (y >= +height)
            {
                r = module2(x1, y-height);
                if ((r > MU_B)&&(r < MU_B + width)) return(1);
                if ((y < height + MU_B)&&(x1 > LAMBDA + MU_B)&&(x1 < LAMBDA + MU_B + width)) return(1);
                return(0);
            }
            else
            {
                r = module2(x1-LAMBDA, y-height);
                if ((r > MU_B)&&(r < MU_B + width)) return(1);
            }
            
            return(0);
        }   
        case (D_WAVEGUIDES_COUPLED):
        {
            x1 = vabs(x)*2.0/XMAX;
            y1 = vabs(y);
            
            width = 0.03;
            b = 1.0*(1.0 - 2.0*MU - width);
            a = MU + 0.5*width - b;
            x6 = x1*x1*x1;
            x5 = x6*x1*x1;
            x6 = x5*x1;
            
            f = a + b*(1.0 + 2.0*x6)/(1.0 + x6);
            fp = 6.0*b*x5/((1.0 + x6)*(1.0 + x6));
            
            return(vabs(f - y1) < MU*sqrt(1.0 + fp*fp));
        }
        case (D_WAVEGUIDES_COUPLED_N):
        {
            x1 = vabs(x);
            y1 = vabs(y);
            
            width = 0.01;
            b = 1.0*(1.0 - 2.0*MU - width);
            a = MU + 0.5*width - b;
            c = 15.0;
            x6 = x1*x1*x1;
            x5 = x6*x1*x1;
            x6 = x5*x1;
            
            f = a + b*(1.0 + 2.0*c*x6)/(1.0 + c*x6);
            fp = 6.0*b*c*x5/((1.0 + c*x6)*(1.0 + c*x6));
            
            if (vabs(f - y1) < 0.9*MU*sqrt(1.0 + fp*fp)) return(1);
            if (vabs(f - y1) < MU*sqrt(1.0 + fp*fp)) return(2);
            return(0);
        }
        case (D_WAVEGUIDE_S):
        {
            a = 0.5*(1.0 - MU);
            if (x < 0.0)
            {
                x = -x;
                y = -y;
            }
            if (x < LAMBDA)
            {
                if (vabs(y - 2.0*a) < MU) return(1);
                if (vabs(y) < MU) return(1);
                if (vabs(y + 2.0*a) < MU) return(1);
                return(0);
            }
            r = module2(x-LAMBDA, y-a);
            return(vabs(r - a) < MU);
        }
        case (D_WAVEGUIDE_S_SHORT):
        {
            a = -1.0 + MU;
            if (y < 0.0)
            {
                x = -x;
                y = -y;
            }
            if (x < -LAMBDA)
            {
                if (vabs(y - 1.0 + MU) > MU) return(0);
                if (vabs(y - 1.0 + MU) > 0.8*MU) return(2);
                return(1);
            }
            r = module2(x-a, y);
            if (vabs(r - 1.0 + MU) < 0.8*MU) return(1);
            if (vabs(r - 1.0 + MU) < MU) return(2);
            return(0);
        }
        case (D_WAVEGUIDE_BADSPLICE):
        {
            if (x < 0.0) 
            {
                y1 = vabs(vabs(y) - 0.5);
            }
            else
            {
                if (y > 0.0) y1 = vabs(y - 0.5 - MU);
                else y1 = vabs(y + 0.5 - MU_B);
            }
            if (y1 < 0.9*LAMBDA) return(1);
            if (y1 < LAMBDA) return(2);
            return(0);
        }
        case (D_MAZE):
        {
            for (i=0; i<npolyrect; i++)
                if ((x > polyrect[i].x1)&&(x < polyrect[i].x2)&&(y > polyrect[i].y1)&&(y < polyrect[i].y2)) return(0);
            return(1);
        }
        case (D_MAZE_CLOSED):
        {
            for (i=0; i<npolyrect; i++)
                if ((x > polyrect[i].x1)&&(x < polyrect[i].x2)&&(y > polyrect[i].y1)&&(y < polyrect[i].y2)) return(0);
            return(1);
        }
        case (D_MAZE_CHANNELS):
        {
            for (i=0; i<npolyrect; i++)
                if ((x > polyrect[i].x1)&&(x < polyrect[i].x2)&&(y > polyrect[i].y1)&&(y < polyrect[i].y2)) return(0);
            return(1);
        }
        case (D_MAZE_CIRCULAR):
        {
            for (i=0; i<npolyrect_rot; i++)
                if (xy_in_rectrotated(x, y, polyrectrot[i])) return(0); 
            for (i=0; i<npolyarc; i++)
                if (xy_in_arc(x, y, polyarc[i])) return(0); 
            return(1);
        }
        case (D_CHESSBOARD):
        {
            i = (int)(vabs(x)/LAMBDA + 0.5);
            j = (int)(vabs(y)/LAMBDA + 0.5);
            if ((i+j)%2 == 0) return(1);
            else return(0);
        }
        case (D_TRIANGLE_TILES):
        {
            if (first)
            {
                h = LAMBDA/(2.0*sqrt(3.0));
                hh = h*3.0;
                first = 0;
            }
            i = (int)((y + h)/hh + 10.0);
            y1 = sin(DPI/3.0)*x - 0.5*y;
            j = (int)((y1 + h)/hh + 10.0);
            y1 = sin(-DPI/3.0)*x -0.5*y;
            k = (int)((y1 + h)/hh + 10.0);
            if ((i+j+k)%2 == 0) return(1);
            else return(0);
        }
        case (D_HEX_TILES):
        {
            if (first)
            {
                ra = -1.0/sqrt(3.0);
                rb = -2.0*ra;
                first = 0;
            }
            x1 = (x + ra*y)/LAMBDA + 30.0;
            y1 = rb*y/LAMBDA + 30.0;
            
            x1 = x1 - (double)(3*(int)(x1/3.0));
            y1 = y1 - (double)(3*(int)(y1/3.0));
            
            if ((x1 > 2.0)&&(y1 < 1.0)) return(1);
            if ((x1 < 1.0)&&(y1 > 2.0)) return(1);
            if (x1 + y1 < 1.0) return(1);
            if (x1 + y1 > 5.0) return(1);
            return(0);
        }
        case (D_FUNNELS):
        {
            y1 = y;
            if (y > 0.5*YMAX) y1 -= YMAX;
            if (y < -0.5*YMAX) y1 += YMAX;
            y1 = vabs(y1 - MU*x);
            y1 = vabs(y1 - 0.5*YMAX)*2.0/YMAX;
            if (y1 > 0.25*(1.0 + LAMBDA + x*x)) return(0);
            return(1);
        }
        case (D_ONE_FUNNEL):
        {
            y1 = vabs(y);
            if (y1 > MU + LAMBDA*x*x) return(0);
            return(1);
        }
        case (D_LENSES_RING):
        {
            if (first)
            {
                salpha = DPI/(double)NPOLY;
                h = LAMBDA*tan(PI/(double)NPOLY);
                if (h < MU) ll = sqrt(MU*MU - h*h);
                else ll = 0.0;
                first = 0;
            }
            for (i=0; i<NPOLY; i++)
            {
                ca = cos((double)i*salpha + APOLY*PID);
                sa = sin((double)i*salpha + APOLY*PID);
                x1 = x*ca + y*sa;
                y1 = -x*sa + y*ca;
                if ((module2(x1 - LAMBDA - ll, y1) < MU)&&(module2(x1 - LAMBDA + ll, y1) < MU)) return(0);
            }
            return(1);
        }
        case (D_LENS):
        {
            if (vabs(y) > MU) return(1);
            a = LAMBDA - 0.25;
            if (x > 0) return ((x+a)*(x+a) + y*y > LAMBDA*LAMBDA);
            return ((x-a)*(x-a) + y*y > LAMBDA*LAMBDA);
        }
        case (D_LENS_WALL):
        {
            if (vabs(y) > MU) return(1);
            a = LAMBDA - 0.25;
            if (x > 0) return ((x+a)*(x+a) + y*y > LAMBDA*LAMBDA);
            return ((x-a)*(x-a) + y*y > LAMBDA*LAMBDA);
        }
        case (D_TWO_LENSES_WALL):
        {
            width = 0.2;
            a = 1.9;
            b = 1.9 - 2.0*LAMBDA + 2.0*width;
            x1 = vabs(x);
            if ((x1 < WALL_WIDTH)&&(vabs(y) > MU)) return(2);
            if (module2(x1 - a, y) > LAMBDA) return(1);
            if (module2(x1 - b, y) > LAMBDA) return(1);
            return(0);
        }
        case (D_TWO_LENSES_OBSTACLE):
        {
            width = 0.2;
            a = 1.9;
            b = 1.9 - 2.0*LAMBDA + 2.0*width;
            x1 = vabs(x);
            if (module2(x1 - a, y) > LAMBDA) return(1);
            if (module2(x1 - b, y) > LAMBDA) return(1);
            return(0);
        }
        case (D_FRESNEL_ZONE_PLATE):
        {
            a = 6.25;
            b = 0.125;
            if (vabs(x + LAMBDA) >= MU) return(1);
            return(cos(DPI*(a*y*y + b)) < 0.0);
        }
        case (D_FRESNEL_ZONE_PLATE_INV):
        {
            a = 6.25;
            if (vabs(x + LAMBDA) >= MU) return(1);
            return(cos(DPI*a*y*y) > 0.0);
        }
        case (D_LENS_ROTATED):
        {
            ca = cos(APOLY*PID);
            sa = sin(APOLY*PID);
            x1 = x*ca + y*sa;
            y1 = -x*sa + y*ca;
            if (vabs(y1) > MU) return(1);
            a = LAMBDA - 0.25;
            if (x1 > 0) return ((x1+a)*(x1+a) + y1*y1 > LAMBDA*LAMBDA);
            return ((x1-a)*(x1-a) + y1*y1 > LAMBDA*LAMBDA);
        }
        case (D_LENS_CONCAVE):
        {
            width = 0.1;
            x1 = vabs(x);
            if (vabs(y) > MU) return(1);
            return((x1 - LAMBDA - width)*(x1 - LAMBDA - width) + y*y < LAMBDA*LAMBDA);
        }
        case (D_LENS_CONVEX_CONCAVE):
        {
            if (vabs(y) > MU) return(1);
            xshift = 0.6;
            if (x > 0.0)
            {
                x1 = x - xshift;
                a = LAMBDA - 0.25;
                if (x1 > 0) return ((x1+a)*(x1+a) + y*y > LAMBDA*LAMBDA);
                return ((x1-a)*(x1-a) + y*y > LAMBDA*LAMBDA);
            }
            else
            {
                width = 0.1;
                x1 = vabs(x + xshift);
                return((x1 - LAMBDA - width)*(x1 - LAMBDA - width) + y*y < LAMBDA*LAMBDA);
            }
        }
        case (D_TREE):
        {
            /* upper triangle */
            h = 1.5*LAMBDA;
            r = 1.5*LAMBDA;
            x1 = vabs(x);
            y1 = y - h;
            zz[0][0] = r*cos(PI/6.0);   zz[0][1] = -0.5*r;
            zz[1][0] = 0.0;             zz[1][1] = r;
            zz[2][0] = 0.0;             zz[2][1] = -0.5*r;
            
            if (y1 > 0.95*zz[1][1]) return(1);
            
            /* middle triangle */
            h = 0.25*LAMBDA;
            r = 2.0*LAMBDA;
            y2 = y - h;
            zz[3][0] = r*cos(PI/6.0);   zz[3][1] = -0.5*r;
            zz[4][0] = 0.0;             zz[4][1] = r;
            zz[5][0] = 0.0;             zz[5][1] = -0.5*r;
            
            /* bottom triangle */
            h = -1.5*LAMBDA;
            r = 3.0*LAMBDA;
            y3 = y - h;
            zz[6][0] = r*cos(PI/6.0);   zz[6][1] = -0.5*r;
            zz[7][0] = 0.0;             zz[7][1] = r;
            zz[8][0] = 0.0;             zz[8][1] = -0.5*r;
            
            if (xy_in_triangle(x1, y1, zz[0], zz[1], zz[2])) return(0);
            if (xy_in_triangle(x1, y2, zz[3], zz[4], zz[5])) return(0);
            if (xy_in_triangle(x1, y3, zz[6], zz[7], zz[8])) return(0);
            
            if (xy_in_triangle(0.9*x1, 0.9*y1, zz[0], zz[1], zz[2])) return(1);
            if (xy_in_triangle(0.9*x1, 0.9*y2, zz[3], zz[4], zz[5])) return(1);
            if (xy_in_triangle(0.9*x1, 0.9*y3, zz[6], zz[7], zz[8])) return(1);
            
            if (xy_in_triangle(0.8*x1, 0.8*y1, zz[0], zz[1], zz[2])) return(2);
            if (xy_in_triangle(0.8*x1, 0.8*y2, zz[3], zz[4], zz[5])) return(2);
            if (xy_in_triangle(0.8*x1, 0.8*y3, zz[6], zz[7], zz[8])) return(2);

            return(1);
        }
        case (D_MICHELSON):
        {
            if ((vabs(x) > LAMBDA)&&(vabs(y) > LAMBDA)) return(2);
            r = 0.5*sqrt(2.0);
            h = YMAX - 0.5*WALL_WIDTH;
            x1 = (x + y)*r;
            y1 = (-x + y)*r;
            if ((vabs(y1) < 0.5*WALL_WIDTH)&&(vabs(x1) < LAMBDA*sqrt(2.0))) return(1);
            if ((vabs(x - h) < 0.5*WALL_WIDTH)&&(vabs(y) < LAMBDA)) return(2);
            if ((vabs(x) < LAMBDA)&&(vabs(y - h) < 0.5*WALL_WIDTH)) return(2);
            return(0);
        }
        case (D_MICHELSON_MOVING):
        {
            if ((vabs(x) > LAMBDA)&&(vabs(y) > LAMBDA)) return(2);
            r = 0.5*sqrt(2.0);
            h = YMAX - MU - 0.5*WALL_WIDTH;
            h1 = YMAX - MU - 0.5*WALL_WIDTH + michelson_position;
            x1 = (x + y)*r;
            y1 = (-x + y)*r;
            if ((vabs(y1) < 0.5*WALL_WIDTH)&&(vabs(x1) < LAMBDA*sqrt(2.0))) return(1);
            if ((vabs(x - h1) < 0.5*WALL_WIDTH)&&(vabs(y) < LAMBDA)) return(2);
            if ((vabs(x) < LAMBDA)&&(vabs(y - h) < 0.5*WALL_WIDTH)) return(2);
            return(0);
        }
        case (D_RITCHEY_CHRETIEN_SPHERICAL):
        {
            /* LAMBDA is magnification M */
            d = 2.5;
            b = 3.5;
            r1 = 2.0*(d + b/LAMBDA);
            r2 = 6.0/(LAMBDA - 1.0);
            x1 = 0.5*d - r1;
            x2 = -0.5*d - r2;
            y1 = vabs(y);
            if (x < 0.0)
            {
                if (y1 > 0.3) return(1);
                if (x < -0.5*d - WALL_WIDTH) return(1);
                return(module2(x-x2, y1) > r2);
            }
            else
            {
                if (y1 < MU) return(1);
                if (x > 0.5*d + WALL_WIDTH) return(1);
                return(module2(x-x1, y1) < r1);
            }
        }
        case (D_RITCHEY_CHRETIEN_HYPERBOLIC):
        {
            return(rc_hyp(x, y));
        }
        case (D_GRADIENT_INDEX_LENS):
        {
            /* use IOR_GRADIENT_INDEX_LENS for index of refraction control */
            return(1);
        }
        case (D_IMAGE):
        {
            /* Not needed, uses option XYIN_INITIALISED */
            return(1);
        }
        case (D_MAGNETRON):
        {
            r = module2(x,y);
            if (r > LAMBDA) return(1);
            if (r < 0.33*LAMBDA) return(1);
            if ((vabs(y) < WALL_WIDTH)&&(x > -0.66*LAMBDA)) return(1);
            if ((vabs(x) < WALL_WIDTH)&&(vabs(y) < 0.66*LAMBDA)) return(1);
            if ((vabs(x-y) < WALL_WIDTH*sqrt(2.0))&&(r < 0.66*LAMBDA)) return(1);
            if ((vabs(x+y) < WALL_WIDTH*sqrt(2.0))&&(r < 0.66*LAMBDA)) return(1);
            for (k=0; k<8; k++)
            {
                x1 = 0.66*cos((double)k*0.5*PID);
                y1 = 0.66*sin((double)k*0.5*PID);
                r = module2(x - x1, y - y1);
                if (r < MU) return(1);
            }
            return(0);
        }
        case (D_MAGNETRON_CATHODE):
        {
            r = module2(x,y);
            if (r < WALL_WIDTH) return(0);
            if (r > LAMBDA) return(1);
            if (r < 0.33*LAMBDA) return(1);
            if ((vabs(y) < WALL_WIDTH)&&(x > -0.66*LAMBDA)) return(1);
            if ((vabs(x) < WALL_WIDTH)&&(vabs(y) < 0.66*LAMBDA)) return(1);
            if ((vabs(x-y) < WALL_WIDTH*sqrt(2.0))&&(r < 0.66*LAMBDA)) return(1);
            if ((vabs(x+y) < WALL_WIDTH*sqrt(2.0))&&(r < 0.66*LAMBDA)) return(1);
            for (k=0; k<8; k++)
            {
                x1 = 0.66*cos((double)k*0.5*PID);
                y1 = 0.66*sin((double)k*0.5*PID);
                r = module2(x - x1, y - y1);
                if (r < MU) return(1);
            }
            return(0);
        }
        case (D_TWOCIRCLES):
        {
            x1 = 0.25*(LAMBDA + 5.0*MU);
            x2 = -0.75*(LAMBDA + 2.0*MU);
            if (module2(x - x1, y) < LAMBDA) return(1);
            if (module2(x - x2, y) < MU) return(1);
            if ((x < x1)&&(x > x2)&&(vabs(y) < WALL_WIDTH)) return(1);
            return(0);
        }
        case (D_POLYCIRCLES):
        {
            if (module2(x, y) < LAMBDA) return(1);
            for (k=0; k<NPOLY; k++)
            {
                angle = APOLY*PID + (double)k*DPI/(double)NPOLY;
                ca = cos(angle);
                sa = sin(angle);
                x1 = (LAMBDA + 2.0*MU)*ca;
                y1 = (LAMBDA + 2.0*MU)*sa;
                if (module2(x - x1, y - y1) < MU) return(1);
                x1 = x*ca + y*sa;
                y1 = -x*sa + y*ca;
                if ((x1 > 0.0)&&(x1 < LAMBDA + 2.0*MU)&&(vabs(y1) < WALL_WIDTH)) return(1);
            }
            return(0);
        }
        case (D_POLYCIRCLES_ANGLED):
        {
            if (module2(x, y) < LAMBDA) return(1);
            alpha = -0.5*PID;
            cb = cos(alpha);
            sb = sin(alpha);
            for (k=0; k<NPOLY; k++)
            {
                angle = APOLY*PID + (double)k*DPI/(double)NPOLY;
                ca = cos(angle);
                sa = sin(angle);
                x1 = x*ca + y*sa;
                y1 = -x*sa + y*ca;
                x2 = LAMBDA + 2.0*MU*cb;
                y2 = 2.0*MU*sb;
                if (module2(x1 - x2, y1 - y2) < MU) return(1);
                c1 = -sb*LAMBDA + cb*WALL_WIDTH;
                c2 = -sb*LAMBDA - cb*WALL_WIDTH;
                b = -sb*x1 + cb*y1;
                if ((x1 > 0.5*LAMBDA)&&(x1 < LAMBDA+MU)&&(b < c1)&&(b > c2)) return(1);
            }
            return(0);
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


void hex_transfo(double u, double v, double *x, double *y)
/* linear transformation of plane used for hex tiles */
{
    static double ra, rb;
    static int first = 1;
    
    if (first)
    {
        ra = 0.5;
        rb = 0.5*sqrt(3.0);
        first = 0;
    }   
    
    *x = u + ra*v;
    *y = rb*v;
}

void draw_rc_hyp()
/* draw D_RITCHEY_CHRETIEN_HYPERBOLIC domain */
{
    static int first = 1;
    static double m, d, dprime, b, r1, r2, e1, e2, a1, b1, a2, b2, x01, x02, dy;
    double x, y, f, f1, f2;
    
    if (first)
    {
        m = LAMBDA;     /* secondary magnification */
        d = 2.5;        /* distance between mirrors */
        b = 3.0;        /* distance between secondary mirror and effective focal point */
        
        f = m*d + b;
        f1 = f/m;       /* focal distance of primary mirror */
        dprime = f1 - d;    /* distance between secondary mirror and primary focal point */
        f2 = m*dprime/(m+1); /* focal distance of secondary mirror */
        
        r1 = 2.0*f/m;      /* radius of curvature of primary mirror */
        r2 = 2.0*b/(m - 1.0);    /* radius of curvature of secondary mirror */
        e1 = sqrt(1.0 + 2.0*b/(m*m*m*d));    /* eccentricity of primary mirror */ 
        e2 = sqrt(1.0 + 2.0*(m*(2.0*m - 1.0) + b/d)/((m-1)*(m-1)*(m-1)));    

        a1 = f1/e1;
        b1 = sqrt(a1*r1);       /* semi-axes of primary mirror */
        x01 = 0.5*d + a1;       /* center of primary hyperbola */
        
        a2 = f2/e2;
        b2 = sqrt(a2*r2);       /* semi-axes of secondary mirror */
        x02 = -0.5*d + a2;       /* center of secondary hyperbola */
                
        dy = (YMAX - YMIN)/(double)NSEG;
        first = 0;
    }
    
    /* primary mirror */
    glBegin(GL_LINE_STRIP);
    draw_vertex(0.5*d + WALL_WIDTH, YMAX);
    draw_vertex(0.5*d + WALL_WIDTH, MU);
    for (y = MU; y < YMAX; y += dy)
    {
        x = x01 - a1*sqrt(1.0 + y*y/(b1*b1));
        draw_vertex(x, y);
    }
    glEnd();
           
    glBegin(GL_LINE_STRIP);
    draw_vertex(0.5*d + WALL_WIDTH, YMIN);
    draw_vertex(0.5*d + WALL_WIDTH, -MU);
    for (y = -MU; y > YMIN; y -= dy)
    {
        x = x01 - a1*sqrt(1.0 + y*y/(b1*b1));
        draw_vertex(x, y);
    }
    glEnd();
    
    /* secondary mirror */
    glBegin(GL_LINE_STRIP);
    draw_vertex(XMIN, -0.3);
    for (y = -0.3; y < 0.3; y += dy)
    {
        x = x02 - a2*sqrt(1.0 + y*y/(b2*b2));
        draw_vertex(x, y);
    }
    draw_vertex(XMIN, 0.3);
    glEnd();
    
    /* focal plane */
    glBegin(GL_LINE_STRIP);
    draw_vertex(b - 0.5*d, YMIN);
    draw_vertex(b - 0.5*d, YMAX);
    glEnd();
}


void draw_billiard(int fade, double fade_value)      /* draws the billiard boundary */
{
    double x0, y0, x, y, x1, y1, x2, y2, dx, dy, phi, r = 0.01, pos[2], pos1[2], alpha, alpha2, dphi, omega, z, l, width, a, b, c, d, r1, r2, ymax, height, xmin, xmax, ca, sa, xshift, x5, x6, f, fp, xratio, w;
    int i, j, k, k1, k2, mr2, ntiles;
    static int first = 1, nsides;
    static double h, hh, sqr3, ll, salpha, arcangle;
    char message[100];

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
        case (D_EXT_ELLIPSE):
        {
            glBegin(GL_LINE_LOOP);
            for (i=0; i<=NSEG; i++)
            {
                phi = (double)i*DPI/(double)NSEG;
                x = LAMBDA*cos(phi);
                y = MU*sin(phi);
                xy_to_pos(x, y, pos);
                glVertex2d(pos[0], pos[1]);
            }
            glEnd ();

            /* draw foci */
            if (FOCI)
            {
                if (fade) glColor3f(0.3*fade_value, 0.3*fade_value, 0.3*fade_value);
                else glColor3f(0.3, 0.3, 0.3);
                x0 = sqrt(LAMBDA*LAMBDA-MU*MU);

                glLineWidth(2);
                glEnable(GL_LINE_SMOOTH);
                
                draw_circle(x0, 0.0, r, NSEG);
                draw_circle(-x0, 0.0, r, NSEG);
            }
            break;
        }
        case (D_EXT_ELLIPSE_CURVED):
        {
            glBegin(GL_LINE_LOOP);
            for (i=0; i<=NSEG; i++)
            {
                phi = (double)i*DPI/(double)NSEG;
                x = LAMBDA*cos(phi);
                y = MU*sin(phi) - 0.4*x*x;
                xy_to_pos(x, y, pos);
                glVertex2d(pos[0], pos[1]);
            }
            glEnd ();

            break;
        }
        case (D_EXT_ELLIPSE_CURVED_BDRY):
        {
            glBegin(GL_LINE_LOOP);
            for (i=0; i<=NSEG; i++)
            {
                phi = (double)i*DPI/(double)NSEG;
                x = LAMBDA*cos(phi);
                y = MU*sin(phi) - 0.4*x*x;
                xy_to_pos(x, y, pos);
                glVertex2d(pos[0], pos[1]);
            }
            glEnd ();
            
            draw_line(XMIN, YMAX - 0.05, XMAX, YMAX - 0.05);
            draw_line(XMIN, YMIN + 0.05, XMAX, YMIN + 0.05);

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
        case (D_WAVEGUIDE):
        {   
            x1 = XMIN + MU;
            x2 = XMAX - 2.0*MU - 1.5*LAMBDA;
            y1 = 0.5*LAMBDA;
            y2 = 1.5*LAMBDA;
             
            glLineWidth(BOUNDARY_WIDTH);
            draw_line(x1, y1, x2, y1);
            draw_line(x1, y1, x1, y2);
            draw_line(x1, y2, x2, y2);
            draw_line(x1, -y1, x2, -y1);
            draw_line(x1, -y1, x1, -y2);
            draw_line(x1, -y2, x2, -y2);
            
            dphi = PI/(double)NSEG;
            glBegin(GL_LINE_STRIP);
            for (i=0; i<NSEG; i++) 
            {
                phi = -PID + dphi*(double)i;
                x = x2 + y2*cos(phi);
                y = y2*sin(phi);
                xy_to_pos(x, y, pos);
                glVertex2d(pos[0], pos[1]);
            }
            glEnd();
            
            glBegin(GL_LINE_STRIP);
            for (i=0; i<NSEG; i++) 
            {
                phi = -PID + dphi*(double)i;
                x = x2 + y1*cos(phi);
                y = y1*sin(phi);
                xy_to_pos(x, y, pos);
                glVertex2d(pos[0], pos[1]);
            }
            glEnd();
            
            break;
        }
        case (D_WAVEGUIDE_W):
        {
            width = LAMBDA - 2.0*MU;
            height = 0.5*MU;
            
            draw_circle_arc(0.0, height, MU, 0.0, PI, NSEG);
            draw_circle_arc(0.0, height, MU + width, 0.0, PI, NSEG);
            
            draw_circle_arc(LAMBDA, -height, MU, PI, PI, NSEG);
            draw_circle_arc(LAMBDA, -height, MU + width, PI, PI, NSEG);
            
            draw_circle_arc(-LAMBDA, -height, MU, PI, PI, NSEG);
            draw_circle_arc(-LAMBDA, -height, MU + width, PI, PI, NSEG);
            
            draw_line(-2.0*LAMBDA + MU, - height, -2.0*LAMBDA + MU, height + MU + width);
            draw_line(-2.0*LAMBDA + MU, height + MU + width, -2.0*LAMBDA + MU + width, height + MU + width);
            draw_line(-2.0*LAMBDA + MU + width, height + MU + width, -2.0*LAMBDA + MU + width, -height);
            
            draw_line(-MU-width, -height, -MU-width, height);
            draw_line(-MU, -height, -MU, height);
            
            draw_line(2.0*LAMBDA - MU, - height, 2.0*LAMBDA - MU, height + MU + width);
            draw_line(2.0*LAMBDA - MU, height + MU + width, 2.0*LAMBDA - MU - width, height + MU + width);
            draw_line(2.0*LAMBDA - MU - width, height + MU + width, 2.0*LAMBDA - MU - width, -height);
            
            draw_line(MU+width, -height, MU+width, height);
            draw_line(MU, -height, MU, height);
            
            break;
        }
        case (D_WAVEGUIDES_W):
        {
            /* upper waveguide */
            width = LAMBDA - 2.0*MU;
            height = MU + width;

            draw_circle_arc(0.0, height, MU, 0.0, PI, NSEG);
            draw_circle_arc(0.0, height, MU + width, 0.0, PI, NSEG);
            
            draw_circle_arc(LAMBDA, height, MU, PI, PI, NSEG);
            draw_circle_arc(LAMBDA, height, MU + width, PI, PI, NSEG);
            
            draw_circle_arc(-LAMBDA, height, MU, PI, PI, NSEG);
            draw_circle_arc(-LAMBDA, height, MU + width, PI, PI, NSEG);
            
            draw_line(LAMBDA + MU, height, LAMBDA + MU, height + MU);
            draw_line(LAMBDA + MU, height + MU, LAMBDA + MU + width, height + MU);
            draw_line(LAMBDA + MU + width, height + MU, LAMBDA + MU + width, height);
            
            draw_line(-LAMBDA - MU, height, -LAMBDA - MU, height + MU);
            draw_line(-LAMBDA - MU, height + MU, -LAMBDA - MU - width, height + MU);
            draw_line(-LAMBDA - MU - width, height + MU, -LAMBDA - MU - width, height);
            
            /* lower waveguide */
            height = -MU - width;
            width = LAMBDA - 2.0*MU_B;
            
            draw_circle_arc(0.0, height, MU_B, 0.0, PI, NSEG);
            draw_circle_arc(0.0, height, MU_B + width, 0.0, PI, NSEG);
            
            draw_circle_arc(LAMBDA, height, MU_B, PI, PI, NSEG);
            draw_circle_arc(LAMBDA, height, MU_B + width, PI, PI, NSEG);
            
            draw_circle_arc(-LAMBDA, height, MU_B, PI, PI, NSEG);
            draw_circle_arc(-LAMBDA, height, MU_B + width, PI, PI, NSEG);
            
            draw_line(LAMBDA + MU_B, height, LAMBDA + MU_B, height + MU_B);
            draw_line(LAMBDA + MU_B, height + MU_B, LAMBDA + MU_B + width, height + MU_B);
            draw_line(LAMBDA + MU_B + width, height + MU_B, LAMBDA + MU_B + width, height);
            
            draw_line(-LAMBDA - MU_B, height, -LAMBDA - MU_B, height + MU_B);
            draw_line(-LAMBDA - MU_B, height + MU_B, -LAMBDA - MU_B - width, height + MU_B);
            draw_line(-LAMBDA - MU_B - width, height + MU_B, -LAMBDA - MU_B - width, height);
            
            break;
        }
        case (D_WAVEGUIDES_COUPLED):
        {
            width = 0.03;
            xratio = XMAX/2.0;
            b = 1.0*(1.0 - 2.0*MU - width);
            a = MU + 0.5*width - b;
            
            if ((DRAW_WAVE_PROFILE)&&(VERTICAL_WAVE_PROFILE))
                dx = (XMAX - XMIN - 0.06)/(double)NSEG;
            else dx = (XMAX - XMIN + 0.1)/(double)NSEG;
            dx = dx/xratio;
            xmin = XMIN/xratio;
            
            x1 = xmin - 0.1;
            y1 = 1.0;
            for (i=0; i<NSEG; i++)
            {
                x2 = xmin - 0.1 + (double)i*dx;
                x6 = x2*x2*x2;
                x5 = x6*x2*x2;
                x6 = x5*x2;
            
                f = a + b*(1.0 + 2.0*x6)/(1.0 + x6);
                fp = 6.0*b*x5/((1.0 + x6)*(1.0 + x6));
                
                y2 = f + MU*sqrt(1.0 + fp*fp);
                
                draw_line(x1*xratio, y1, x2*xratio, y2);
                x1 = x2;
                y1 = y2;
            }
            
            x1 = xmin - 0.1;
            y1 = 1.0;
            for (i=0; i<NSEG; i++)
            {
                x2 = xmin - 0.1 + (double)i*dx;
                x6 = x2*x2*x2;
                x5 = x6*x2*x2;
                x6 = x5*x2;
            
                f = a + b*(1.0 + 2.0*x6)/(1.0 + x6);
                fp = 6.0*b*x5/((1.0 + x6)*(1.0 + x6));
                
                y2 = f - MU*sqrt(1.0 + fp*fp);
                
                draw_line(x1*xratio, y1, x2*xratio, y2);
                x1 = x2;
                y1 = y2;
            }
            
            x1 = xmin - 0.1;
            y1 = -1.0;
            for (i=0; i<NSEG; i++)
            {
                x2 = xmin - 0.1 + (double)i*dx;
                x6 = x2*x2*x2;
                x5 = x6*x2*x2;
                x6 = x5*x2;
            
                f = a + b*(1.0 + 2.0*x6)/(1.0 + x6);
                fp = 6.0*b*x5/((1.0 + x6)*(1.0 + x6));
                
                y2 = -(f + MU*sqrt(1.0 + fp*fp));
                
                draw_line(x1*xratio, y1, x2*xratio, y2);
                x1 = x2;
                y1 = y2;
            }
            
            x1 = xmin - 0.1;
            y1 = -1.0;
            for (i=0; i<NSEG; i++)
            {
                x2 = xmin - 0.1 + (double)i*dx;
                x6 = x2*x2*x2;
                x5 = x6*x2*x2;
                x6 = x5*x2;
            
                f = a + b*(1.0 + 2.0*x6)/(1.0 + x6);
                fp = 6.0*b*x5/((1.0 + x6)*(1.0 + x6));
                
                y2 = -f + MU*sqrt(1.0 + fp*fp);
                
                draw_line(x1*xratio, y1, x2*xratio, y2);
                x1 = x2;
                y1 = y2;
            }
            
            x2 = LAMBDA/xratio;
            x6 = x2*x2*x2;
            x5 = x6*x2*x2;
            x6 = x5*x2;
            
            f = a + b*(1.0 + 2.0*x6)/(1.0 + x6);
            fp = 6.0*b*x5/((1.0 + x6)*(1.0 + x6));
                
            y1 = f - MU*sqrt(1.0 + fp*fp);
            y2 = f + MU*sqrt(1.0 + fp*fp);
            
            x2 = LAMBDA;
            
            draw_line(x2, YMIN, x2, -y2);
            draw_line(x2, -y1, x2, y1);
            draw_line(x2, y2, x2, YMAX);
            draw_line(-x2, YMIN, -x2, -y2);
            draw_line(-x2, -y1, -x2, y1);
            draw_line(-x2, y2, -x2, YMAX);
            
            break;
        }
        case (D_WAVEGUIDES_COUPLED_N):
        {
            width = 0.01;
            xratio = XMAX/2.0;
            b = 1.0*(1.0 - 2.0*MU - width);
            a = MU + 0.5*width - b;
            c = 15.0;
            
            if ((DRAW_WAVE_PROFILE)&&(VERTICAL_WAVE_PROFILE))
                dx = (XMAX - XMIN - 0.06)/(double)NSEG;
            else dx = (XMAX - XMIN + 0.1)/(double)NSEG;
            dx = dx/xratio;
            xmin = XMIN/xratio;
            
            x1 = xmin - 0.1;
            y1 = 1.0;
            for (i=0; i<NSEG; i++)
            {
                x2 = xmin - 0.1 + (double)i*dx;
                x6 = x2*x2*x2;
                x5 = x6*x2*x2;
                x6 = x5*x2;
            
                f = a + b*(1.0 + 2.0*c*x6)/(1.0 + c*x6);
                fp = 6.0*b*c*x5/((1.0 + c*x6)*(1.0 + c*x6));
                
                y2 = f + MU*sqrt(1.0 + fp*fp);
                
                draw_line(x1*xratio, y1, x2*xratio, y2);
                x1 = x2;
                y1 = y2;
            }
            
            x1 = xmin - 0.1;
            y1 = 1.0;
            for (i=0; i<NSEG; i++)
            {
                x2 = xmin - 0.1 + (double)i*dx;
                x6 = x2*x2*x2;
                x5 = x6*x2*x2;
                x6 = x5*x2;
            
                f = a + b*(1.0 + 2.0*c*x6)/(1.0 + c*x6);
                fp = 6.0*b*c*x5/((1.0 + c*x6)*(1.0 + c*x6));
                
                y2 = f - MU*sqrt(1.0 + fp*fp);
                
                draw_line(x1*xratio, y1, x2*xratio, y2);
                x1 = x2;
                y1 = y2;
            }
            
            x1 = xmin - 0.1;
            y1 = -1.0;
            for (i=0; i<NSEG; i++)
            {
                x2 = xmin - 0.1 + (double)i*dx;
                x6 = x2*x2*x2;
                x5 = x6*x2*x2;
                x6 = x5*x2;
            
                f = a + b*(1.0 + 2.0*c*x6)/(1.0 + c*x6);
                fp = 6.0*b*c*x5/((1.0 + c*x6)*(1.0 + c*x6));
                
                y2 = -(f + MU*sqrt(1.0 + fp*fp));
                
                draw_line(x1*xratio, y1, x2*xratio, y2);
                x1 = x2;
                y1 = y2;
            }
            
            x1 = xmin - 0.1;
            y1 = -1.0;
            for (i=0; i<NSEG; i++)
            {
                x2 = xmin - 0.1 + (double)i*dx;
                x6 = x2*x2*x2;
                x5 = x6*x2*x2;
                x6 = x5*x2;
            
                f = a + b*(1.0 + 2.0*c*x6)/(1.0 + c*x6);
                fp = 6.0*b*c*x5/((1.0 + c*x6)*(1.0 + c*x6));
                
                y2 = -f + MU*sqrt(1.0 + fp*fp);
                
                draw_line(x1*xratio, y1, x2*xratio, y2);
                x1 = x2;
                y1 = y2;
            }
            
            x2 = LAMBDA/xratio;
            x6 = x2*x2*x2;
            x5 = x6*x2*x2;
            x6 = x5*x2;
            
            f = a + b*(1.0 + 2.0*c*x6)/(1.0 + c*x6);
            fp = 6.0*b*c*x5/((1.0 + c*x6)*(1.0 + c*x6));
                
            y1 = f - MU*sqrt(1.0 + fp*fp);
            y2 = f + MU*sqrt(1.0 + fp*fp);
            
            x2 = LAMBDA;
            b = 0.7;
            
            draw_line(x2, YMIN, x2, -y2);
            draw_line(x2, -y1, x2, -b);
            draw_line(x2, -b, XMAX, -b);
            draw_line(XMAX, b, x2, b);
            draw_line(x2, b, x2, y1);
            draw_line(x2, y2, x2, YMAX);
            
            draw_line(-x2, YMIN, -x2, -y2);
            draw_line(-x2, -y1, -x2, -b);
            draw_line(-x2, -b, XMIN, -b);
            draw_line(XMIN, b, -x2, b);
            draw_line(-x2, b, -x2, y1);
            draw_line(-x2, y2, -x2, YMAX);
                        
            break;
        }
        case (D_WAVEGUIDE_S):
        {
            a = 0.5*(1.0 - MU);
            
            draw_circle_arc(LAMBDA, a, a-MU, -PID, PI, NSEG);
            draw_circle_arc(LAMBDA, a, a+MU, -PID, PI, NSEG);
            draw_circle_arc(-LAMBDA, -a, a-MU, PID, PI, NSEG);
            draw_circle_arc(-LAMBDA, -a, a+MU, PID, PI, NSEG);
            
            draw_line(LAMBDA, 1.0, -LAMBDA, 1.0);
            draw_line(-LAMBDA, 1.0, -LAMBDA, 2.0*a-MU);
            draw_line( -LAMBDA, 2.0*a-MU, LAMBDA, 2.0*a-MU);
            
            draw_line(-LAMBDA, MU, LAMBDA, MU);
            draw_line(-LAMBDA, -MU, LAMBDA, -MU);
            
            draw_line(-LAMBDA, -1.0, LAMBDA, -1.0);
            draw_line(LAMBDA, -1.0, LAMBDA, -2.0*a+MU);
            draw_line(LAMBDA, -2.0*a+MU, -LAMBDA, -2.0*a+MU);
            break;
        }
        case (D_WAVEGUIDE_S_SHORT):
        {
            a = -1.0 + MU;
            
            draw_circle_arc(-1.0+MU, 0.0, 1.0-2.0*MU, 0.0, PID, NSEG);
            draw_circle_arc(-1.0+MU, 0.0, 1.0, 0.0, PID, NSEG);
            draw_circle_arc(1.0-MU, 0.0, 1.0-2.0*MU, PI, PID, NSEG);
            draw_circle_arc(1.0-MU, 0.0, 1.0, PI, PID, NSEG);
            
            draw_line(-1.0, 1.0, -1.0+MU, 1.0);
            draw_line(-1.0, 1.0, -1.0, 1.0-2.0*MU);
            draw_line(-1.0, 1.0-2.0*MU, -1.0+MU, 1.0-2.0*MU);
            draw_line(1.0, -1.0, 1.0-MU, -1.0);
            draw_line(1.0, -1.0, 1.0, -1.0+2.0*MU);
            draw_line(1.0, -1.0+2.0*MU, 1.0-MU, -1.0+2.0*MU);
            
            break;
        }
        case (D_WAVEGUIDE_BADSPLICE):
        {
            draw_rectangle(XMIN-1.0, 0.5-LAMBDA, 0.0, 0.5+LAMBDA);
            draw_rectangle(0.0, 0.5-LAMBDA+MU, XMAX + 1.0, 0.5+LAMBDA+MU);
            draw_rectangle(XMIN-1.0, -0.5-LAMBDA, 0.0, -0.5+LAMBDA);
            draw_rectangle(0.0, -0.5-LAMBDA+MU_B, XMAX + 1.0, -0.5+LAMBDA+MU_B);
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
        case (D_MAZE):
        {
            glLineWidth(BOUNDARY_WIDTH);
            if (fade) glColor3f(0.15*fade_value, 0.15*fade_value, 0.15*fade_value);
            else glColor3f(0.15, 0.15, 0.15);
            for (i=0; i<npolyrect; i++)
                draw_filled_rectangle(polyrect[i].x1, polyrect[i].y1, polyrect[i].x2, polyrect[i].y2);
            break;
        }
        case (D_MAZE_CLOSED):
        {
            glLineWidth(BOUNDARY_WIDTH);
            if (fade) glColor3f(0.15*fade_value, 0.15*fade_value, 0.15*fade_value);
            else glColor3f(0.15, 0.15, 0.15);
            for (i=0; i<npolyrect; i++)
                draw_filled_rectangle(polyrect[i].x1, polyrect[i].y1, polyrect[i].x2, polyrect[i].y2);
            break;
        }
        case (D_MAZE_CHANNELS):
        {
            glLineWidth(BOUNDARY_WIDTH);
            if (fade) glColor3f(0.15*fade_value, 0.15*fade_value, 0.15*fade_value);
            else glColor3f(0.15, 0.15, 0.15);
            for (i=0; i<npolyrect; i++)
                draw_filled_rectangle(polyrect[i].x1, polyrect[i].y1, polyrect[i].x2, polyrect[i].y2);
            break;
        }
        case (D_MAZE_CIRCULAR):
        {
            glLineWidth(BOUNDARY_WIDTH);
            if (fade) glColor3f(0.15*fade_value, 0.15*fade_value, 0.15*fade_value);
            else glColor3f(0.15, 0.15, 0.15);
            /* TODO */
            break;
        }
        case (D_CHESSBOARD):
        {
            glLineWidth(BOUNDARY_WIDTH);
            x = 0.5*LAMBDA;
            while (x < XMAX) 
            {
                draw_line(x, YMIN, x, YMAX);
                x += LAMBDA;
            }
            x = -0.5*LAMBDA;
            while (x > XMIN) 
            {
                draw_line(x, YMIN, x, YMAX);
                x -= LAMBDA;
            }
            y = 0.5*LAMBDA;
            while (y < YMAX) 
            {
                draw_line(XMIN, y, XMAX, y);
                y += LAMBDA;
            }
            y = -0.5*LAMBDA;
            while (y > YMIN) 
            {
                draw_line(XMIN, y, XMAX, y);
                y -= LAMBDA;
            }
            
            break;
        }
        case (D_TRIANGLE_TILES):
        {
            if (first)
            {
                h = LAMBDA/(2.0*sqrt(3.0));
                hh = h*3.0;
                sqr3 = sqrt(3.0);
                first = 0;
            }
            glLineWidth(BOUNDARY_WIDTH);
            y = -h;
            while (y < YMAX) 
            {
                draw_line(XMIN, y, XMAX, y);
                y += hh;
            }
            y = -h;
            while (y > YMIN) 
            {
                draw_line(XMIN, y, XMAX, y);
                y -= hh;
            }
            x = -0.5*LAMBDA;
            y = -h;
            while (x < 1.5*XMAX)
            {
                draw_line(x - 10.0, y - 10.0*sqr3, x + 10.0, y + 10.0*sqr3);
                draw_line(x - 10.0, y + 10.0*sqr3, x + 10.0, y - 10.0*sqr3);
                x += LAMBDA;
            }
            x = -0.5*LAMBDA;
            while (x > 1.5*XMIN)
            {
                draw_line(x - 10.0, y - 10.0*sqr3, x + 10.0, y + 10.0*sqr3);
                draw_line(x - 10.0, y + 10.0*sqr3, x + 10.0, y - 10.0*sqr3);
                x -= LAMBDA;
            }
            break;
        }
        case (D_HEX_TILES):
        {
            ntiles = (int)(XMAX/LAMBDA) + 1;
            for (i=-ntiles; i<ntiles; i++)
                for (j=-ntiles; j<ntiles; j++)
                {
                    x0 = 3.0*LAMBDA*(double)i;
                    y0 = 3.0*LAMBDA*(double)j;
                    
                    hex_transfo(x0, y0 + LAMBDA, &x, &y);
                    hex_transfo(x0, y0 + 2.0*LAMBDA, &x1, &y1);
                    draw_line(x, y, x1, y1);
                    hex_transfo(x0 + LAMBDA, y0 + 2.0*LAMBDA, &x2, &y2);
                    draw_line(x1, y1, x2, y2);
                    hex_transfo(x0 + LAMBDA, y0 + 3.0*LAMBDA, &x1, &y1);
                    draw_line(x1, y1, x2, y2);
                    hex_transfo(x0 + 2.0*LAMBDA, y0 + LAMBDA, &x1, &y1);
                    draw_line(x1, y1, x2, y2);
                    hex_transfo(x0 + 3.0*LAMBDA, y0 + LAMBDA, &x2, &y2);
                    draw_line(x1, y1, x2, y2);
                    hex_transfo(x0 + 2.0*LAMBDA, y0, &x2, &y2);
                    draw_line(x1, y1, x2, y2);
                    hex_transfo(x0 + LAMBDA, y0, &x1, &y1);
                    draw_line(x1, y1, x2, y2);
                    draw_line(x1, y1, x, y);
                
                    hex_transfo(x0 + 2.0*LAMBDA, y0 + 3.0*LAMBDA, &x, &y);
                    hex_transfo(x0 + 3.0*LAMBDA, y0 + 2.0*LAMBDA, &x1, &y1);
                    draw_line(x, y, x1, y1);
                    
                }
            break;
            
        }
        case (D_FUNNELS):
        {
            if (LAMBDA < 3.0)
            {
                xmax = sqrt(3.0 - LAMBDA);
                for (j=-2; j<2; j++)
                    for (k=-1; k<=1; k+=2)
                    {
                        for (i=0; i <= NSEG; i++)
                        {
                            x = -xmax + (2.0*xmax)*(double)i/(double)NSEG;
                            y = (double)j*YMAX + MU*x + 0.5*YMAX*(1.0 + 0.25*(double)k*(1.0 + LAMBDA + x*x));
                            if (i > 0) draw_line(x1, y1, x, y);
                            x1 = x;
                            y1 = y;
                        }
                    }
            }
            break;
        }
        case (D_ONE_FUNNEL):
        {
            for (k=-1; k<2; k+=2)
            {
                x1 = XMIN;
                y1 = (double)k*(MU + LAMBDA*x1*x1);
                for (i=0; i<=NSEG; i++)
                {
                    x = XMIN + (XMAX - XMIN)*(double)i/(double)NSEG;
                    y = (double)k*(MU + LAMBDA*x*x);
                    draw_line(x1, y1, x, y);
                    x1 = x;
                    y1 = y;
                }
            }
            break;
        }
        case (D_LENSES_RING):
        {
            if (first)
            {
                salpha = DPI/(double)NPOLY;
                h = LAMBDA*tan(PI/(double)NPOLY);
                if (h < MU) ll = sqrt(MU*MU - h*h);
                else ll = 0.0;
                arcangle = atan(h/ll);
                first = 0;
            }
            for (i=0; i<NPOLY; i++)
            {
                ca = cos((double)i*salpha + APOLY*PID);
                sa = sin((double)i*salpha + APOLY*PID);
                x = ca*(LAMBDA - ll);
                y = sa*(LAMBDA - ll);
                draw_circle_arc(x, y, MU, -arcangle + APOLY*PID + (double)i*salpha, 2.0*arcangle, NSEG);
                x = ca*(LAMBDA + ll);
                y = sa*(LAMBDA + ll);
                draw_circle_arc(x, y, MU, PI-arcangle + APOLY*PID + (double)i*salpha, 2.0*arcangle, NSEG);
            }
            break;
        }
        case (D_LENS):
        {
            a = LAMBDA - 0.25;
            if (first)
            {
                arcangle = asin(MU/LAMBDA);
                ll = LAMBDA*cos(arcangle) - a;
                first = 0;
            }
            
            draw_circle_arc(-a, 0.0, LAMBDA, -arcangle, 2.0*arcangle, NSEG);
            draw_line(ll, MU, -ll, MU);
            draw_circle_arc(a, 0.0, LAMBDA, PI-arcangle, 2.0*arcangle, NSEG);
            draw_line(-ll, -MU, ll, -MU);
            break;
        }
        case (D_LENS_WALL):
        {
            a = LAMBDA - 0.25;
            width = 0.05;
            if (first)
            {
                arcangle = asin(MU/LAMBDA);
                ll = LAMBDA*cos(arcangle) - a;
                first = 0;
            }
            
            draw_circle_arc(-a, 0.0, LAMBDA, -arcangle, 2.0*arcangle, NSEG);
            draw_line(ll, MU, -ll, MU);
            draw_circle_arc(a, 0.0, LAMBDA, PI-arcangle, 2.0*arcangle, NSEG);
            draw_line(-ll, -MU, ll, -MU);
            
            draw_rectangle(-width, MU, width, YMAX + 1.0);
            draw_rectangle(-width, YMIN-1.0, width, -MU);
            break;
        }
        case (D_TWO_LENSES_WALL):
        {
            width = 0.2;
            a = 1.9;
            b = 1.9 - 2.0*LAMBDA + 2.0*width;
            if (first)
            {
                arcangle = acos(1.0 - width/LAMBDA);
                first = 0;
            }
            draw_circle_arc(a, 0.0, LAMBDA, PI-arcangle, 2.0*arcangle, NSEG);
            draw_circle_arc(b, 0.0, LAMBDA, -arcangle, 2.0*arcangle, NSEG);
            draw_circle_arc(-a, 0.0, LAMBDA, -arcangle, 2.0*arcangle, NSEG);
            draw_circle_arc(-b, 0.0, LAMBDA, PI-arcangle, 2.0*arcangle, NSEG);
            
            width = WALL_WIDTH;
            draw_rectangle(-width, MU, width, YMAX + 1.0);
            draw_rectangle(-width, YMIN-1.0, width, -MU);
            break;
        }
        case (D_FRESNEL_ZONE_PLATE):
        {
            a = 6.25;
            b = 0.125;
//             y = sqrt(1.0/(6.0*a));
            y = sqrt((0.25 - b)/a);
            draw_rectangle(-LAMBDA - MU, -y, -LAMBDA + MU, y);
            for (i=1; i<(int)(a*a+1); i++)
            {
                x = sqrt((-0.25 - b + (double)i)/a);
                y = sqrt((0.25 - b + (double)i)/a);
                draw_rectangle(-LAMBDA - MU, x, -LAMBDA + MU, y);
                draw_rectangle(-LAMBDA - MU, -x, -LAMBDA + MU, -y);
            }
            break;
        }
        case (D_FRESNEL_ZONE_PLATE_INV):
        {
            a = 6.25;
            for (i=0; i<2*(int)a; i++)
            {
                x = sqrt((0.25 + (double)i)/a);
                y = sqrt((-0.25 + (double)(i+1))/a);
                draw_rectangle(-LAMBDA - MU, x, -LAMBDA + MU, y);
                draw_rectangle(-LAMBDA - MU, -x, -LAMBDA + MU, -y);
            }
            break;
        }
        case (D_TWO_LENSES_OBSTACLE):
        {
            width = 0.2;
            a = 1.9;
            b = 1.9 - 2.0*LAMBDA + 2.0*width;
            if (first)
            {
                arcangle = acos(1.0 - width/LAMBDA);
                first = 0;
            }
            draw_circle_arc(a, 0.0, LAMBDA, PI-arcangle, 2.0*arcangle, NSEG);
            draw_circle_arc(b, 0.0, LAMBDA, -arcangle, 2.0*arcangle, NSEG);
            draw_circle_arc(-a, 0.0, LAMBDA, -arcangle, 2.0*arcangle, NSEG);
            draw_circle_arc(-b, 0.0, LAMBDA, PI-arcangle, 2.0*arcangle, NSEG);
            
            width = 0.05;
            draw_rectangle(-width, -MU, width, MU);
            break;
        }
        case (D_LENS_ROTATED):
        {
            ca = cos(APOLY*PID);
            sa = sin(APOLY*PID);
            a = LAMBDA - 0.25;
            if (first)
            {
                arcangle = asin(MU/LAMBDA);
                ll = LAMBDA*cos(arcangle) - a;
                first = 0;
            }
            draw_circle_arc(-a*ca, -a*sa, LAMBDA, -arcangle + APOLY*PID, 2.0*arcangle, NSEG);
            draw_line(ll*ca - MU*sa, ll*sa + MU*ca, -ll*ca - MU*sa, -ll*sa + MU*ca);
            draw_circle_arc(a*ca, a*sa, LAMBDA, PI-arcangle + APOLY*PID, 2.0*arcangle, NSEG);
            draw_line(-ll*ca + MU*sa, -ll*sa - MU*ca, ll*ca + MU*sa, ll*sa - MU*ca);
            break;
        }
        case (D_LENS_CONCAVE):
        {
            arcangle = asin(MU/LAMBDA);
            width = 0.1;
            a = LAMBDA + width - sqrt(LAMBDA*LAMBDA - MU*MU);
            draw_circle_arc(LAMBDA + width, 0.0, LAMBDA, PI-arcangle, 2.0*arcangle, NSEG);
            draw_line(a, MU, -a, MU);
            draw_circle_arc(-LAMBDA - width, 0.0, LAMBDA, -arcangle, 2.0*arcangle, NSEG);
            draw_line(a, -MU, -a, -MU);
            break;
        }
        case (D_LENS_CONVEX_CONCAVE):
        {
            xshift = 0.6;
            
            a = LAMBDA - 0.25;
            if (first)
            {
                arcangle = asin(MU/LAMBDA);
                ll = LAMBDA*cos(arcangle) - a;
                first = 0;
            }
            
            draw_circle_arc(xshift-a, 0.0, LAMBDA, -arcangle, 2.0*arcangle, NSEG);
            draw_line(xshift+ll, MU, xshift-ll, MU);
            draw_circle_arc(xshift+a, 0.0, LAMBDA, PI-arcangle, 2.0*arcangle, NSEG);
            draw_line(xshift-ll, -MU, xshift+ll, -MU);
            
            width = 0.1;
            a = LAMBDA + width - sqrt(LAMBDA*LAMBDA - MU*MU);
            draw_circle_arc(LAMBDA + width - xshift, 0.0, LAMBDA, PI-arcangle, 2.0*arcangle, NSEG);
            draw_line(a - xshift, MU, -a - xshift, MU);
            draw_circle_arc(-LAMBDA - width - xshift, 0.0, LAMBDA, -arcangle, 2.0*arcangle, NSEG);
            draw_line(a - xshift, -MU, -a - xshift, -MU);
            break;
        }
        case (D_TREE):
        {
//             h = LAMBDA;
//             r = 2.0*LAMBDA;
//             
//             x1 = r*cos(PI/6.0);
//             y1 = -h;
//             y2 = r;
//             
//             x1 = x1/0.8;
//             y1 = h + y1/0.8;
//             x2 = 0.0;
//             y2 = h + y2/0.8;
//             
//             draw_line(x1, y1, x2, y2);
//             draw_line(x2, y2, -x1, y1);
//             draw_line(-x1, y1, x1, y1);
            
            /* do nothing */
            break;
        }
        case (D_MICHELSON):
        {
            w = 0.25*sqrt(2.0)*WALL_WIDTH;
            h = YMAX - 0.5*WALL_WIDTH;
            glBegin(GL_LINE_LOOP);
            xy_to_pos(LAMBDA + w, LAMBDA - w, pos);
            glVertex2d(pos[0], pos[1]);
            xy_to_pos(LAMBDA - w, LAMBDA + w, pos);
            glVertex2d(pos[0], pos[1]);
            xy_to_pos(-LAMBDA - w, -LAMBDA + w, pos);
            glVertex2d(pos[0], pos[1]);
            xy_to_pos(-LAMBDA + w, -LAMBDA - w, pos);
            glVertex2d(pos[0], pos[1]);
            glEnd();  
            draw_rectangle(h - 0.5*WALL_WIDTH, -LAMBDA, h + 0.5*WALL_WIDTH, LAMBDA);
            draw_rectangle(-LAMBDA, h - 0.5*WALL_WIDTH, LAMBDA, h + 0.5*WALL_WIDTH);
            draw_rectangle(LAMBDA, LAMBDA, XMAX + 1.0, YMAX + 1.0);
            draw_rectangle(LAMBDA, -LAMBDA, XMAX + 1.0, YMIN - 1.0);
            draw_rectangle(-LAMBDA, -LAMBDA, XMIN - 1.0, YMIN - 1.0);
            draw_rectangle(-LAMBDA, LAMBDA, XMIN - 1.0, YMAX + 1.0);
            break;
        }
        case (D_MICHELSON_MOVING):
        {
            w = 0.25*sqrt(2.0)*WALL_WIDTH;
            h = YMAX - MU - 0.5*WALL_WIDTH;
            glBegin(GL_LINE_LOOP);
            xy_to_pos(LAMBDA + w, LAMBDA - w, pos);
            glVertex2d(pos[0], pos[1]);
            xy_to_pos(LAMBDA - w, LAMBDA + w, pos);
            glVertex2d(pos[0], pos[1]);
            xy_to_pos(-LAMBDA - w, -LAMBDA + w, pos);
            glVertex2d(pos[0], pos[1]);
            xy_to_pos(-LAMBDA + w, -LAMBDA - w, pos);
            glVertex2d(pos[0], pos[1]);
            glEnd();  
            draw_rectangle(h - 0.5*WALL_WIDTH + michelson_position, -LAMBDA, h + 0.5*WALL_WIDTH + michelson_position, LAMBDA);
            draw_rectangle(-LAMBDA, h - 0.5*WALL_WIDTH, LAMBDA, h + 0.5*WALL_WIDTH);
            draw_rectangle(LAMBDA, LAMBDA, XMAX + 1.0, YMAX + 1.0);
            draw_rectangle(LAMBDA, -LAMBDA, XMAX + 1.0, YMIN - 1.0);
            draw_rectangle(-LAMBDA, -LAMBDA, XMIN - 1.0, YMIN - 1.0);
            draw_rectangle(-LAMBDA, LAMBDA, XMIN - 1.0, YMAX + 1.0);
            
            x = h - 0.5*WALL_WIDTH + michelson_schedule(INITIAL_TIME);
            draw_line(x, LAMBDA, x, LAMBDA + WALL_WIDTH);
            draw_line(x, -LAMBDA, x, -LAMBDA - WALL_WIDTH);
            
            if (fade) glColor3f(fade_value, fade_value, fade_value);
            else glColor3f(1.0, 1.0, 1.0);
            xy_to_pos(YMAX - MU - 0.1*XMAX, LAMBDA + 0.1*YMAX + 0.1, pos);
//             xy_to_pos(YMAX - MU - 0.3*XMAX, LAMBDA + 0.1*YMAX + 0.1, pos);
            sprintf(message, "Mirror position");
            write_text(pos[0], pos[1], message);
            xy_to_pos(YMAX - MU - 0.1*XMAX, LAMBDA + 0.1*YMAX, pos);
            /* wavelength eyeballed from simulation result, should be improved */
            sprintf(message, "%.3f wavelengths", michelson_position/0.0667);
            write_text(pos[0], pos[1], message);
            break;
        }
        case (D_RITCHEY_CHRETIEN_SPHERICAL):
        {
            d = 2.5;
            b = 3.5;
            r1 = 2.0*(d + b/LAMBDA);
            r2 = 6.0/(LAMBDA - 1.0);
            x1 = 0.5*d - r1;
            x2 = -0.5*d - r2;
            
            alpha = asin(0.3/r2);
            ca = cos(alpha);
            draw_circle_arc(x2, 0.0, r2, -alpha, 2.0*alpha, NSEG);
            draw_line(x2 + r2*ca, 0.3, -0.5*d - WALL_WIDTH, 0.3);
            draw_line(-0.5*d - WALL_WIDTH, 0.3, -0.5*d - WALL_WIDTH, -0.3);
            draw_line(-0.5*d - WALL_WIDTH, -0.3, x2 + r2*ca, -0.3);
            
            alpha = asin(MU/r1);
            ca = cos(alpha);
            draw_circle_arc(x1, 0.0, r1, alpha, 0.5*PID, NSEG);
            draw_circle_arc(x1, 0.0, r1, -alpha, -0.5*PID, NSEG);
            draw_line(x1 + r1*ca, MU, 0.5*d + WALL_WIDTH, MU);
            draw_line(0.5*d + WALL_WIDTH, MU, 0.5*d + WALL_WIDTH, YMAX);
            draw_line(x1 + r1*ca, -MU, 0.5*d + WALL_WIDTH, -MU);
            draw_line(0.5*d + WALL_WIDTH, -MU, 0.5*d + WALL_WIDTH, YMIN);
            break;
        }
        case (D_RITCHEY_CHRETIEN_HYPERBOLIC):
        {
            draw_rc_hyp();
            break;
        }
        case (D_GRADIENT_INDEX_LENS):
        {
            draw_line(-LAMBDA, YMIN, -LAMBDA, YMAX);
            draw_line(LAMBDA, YMIN, LAMBDA, YMAX);
            draw_line(-LAMBDA, -1.0, LAMBDA, -1.0);
            draw_line(-LAMBDA, 1.0, LAMBDA, 1.0);
            /* focal plane */
            draw_line(WAVE_PROFILE_X, YMIN, WAVE_PROFILE_X, YMAX);
            break;
        }
        case (D_IMAGE):
        {
            /* do nothing */
            break;
        }
        case (D_MAGNETRON):
        {
            /* TODO */
            break;
        }
        case (D_MAGNETRON_CATHODE):
        {
            /* TODO */
            break;
        }
        case (D_TWOCIRCLES):
        {
            x1 = 0.25*(LAMBDA + 5.0*MU);
            x2 = -0.75*(LAMBDA + 2.0*MU);
            alpha = asin(WALL_WIDTH/LAMBDA);
            alpha2 = asin(WALL_WIDTH/MU);
            draw_circle_arc(x1, 0.0, LAMBDA, -PI + alpha, DPI - 2.0*alpha, NSEG);
            draw_line(x1 - LAMBDA*cos(alpha), WALL_WIDTH, x2 + MU*cos(alpha2), WALL_WIDTH);
            draw_circle_arc(x2, 0.0, MU, alpha2, DPI - 2.0*alpha2, NSEG);
            draw_line(x2 + MU*cos(alpha2), -WALL_WIDTH, x1 - LAMBDA*cos(alpha), -WALL_WIDTH);
            break; 
        }
        case (D_POLYCIRCLES):
        {
            alpha = asin(WALL_WIDTH/LAMBDA);
            alpha2 = asin(WALL_WIDTH/MU);
            for (k=0; k<NPOLY; k++)
            {
                phi = APOLY*PID + (double)k*DPI/(double)NPOLY;
                x1 = (LAMBDA + 2.0*MU)*cos(phi);
                y1 = (LAMBDA + 2.0*MU)*sin(phi);
                draw_line(LAMBDA*cos(-alpha + phi), LAMBDA*sin(-alpha + phi), x1 - MU*cos(alpha2 + phi), y1 - MU*sin(alpha2 + phi));
                draw_circle_arc(x1, y1, MU, -PI + alpha2 + phi, DPI - 2.0*alpha2, NSEG); 
                draw_line(x1 - MU*cos(-alpha2 + phi), y1 - MU*sin(-alpha2 + phi), LAMBDA*cos(alpha + phi), LAMBDA*sin(alpha + phi)); 
                draw_circle_arc(0.0, 0.0, LAMBDA, alpha + phi, DPI/(double)NPOLY - 2.0*alpha, NSEG); 
            }
            break;
        }
        case (D_POLYCIRCLES_ANGLED):
        {
            alpha = atan(WALL_WIDTH/LAMBDA);
            alpha2 = asin(WALL_WIDTH/(MU*sqrt(2.0)));
            for (k=0; k<NPOLY; k++)
            {
                phi = APOLY*PID + (double)k*DPI/(double)NPOLY;
                x1 = LAMBDA*cos(phi);
                y1 = LAMBDA*sin(phi);
                x2 = x1 + 2.0*MU*cos(phi - 0.5*PID);
                y2 = y1 + 2.0*MU*sin(phi - 0.5*PID);
                draw_line(LAMBDA*cos(phi-alpha), LAMBDA*sin(phi-alpha), 
                          x2 - MU*cos(phi -0.5*PID + alpha2), 
                          y2 - MU*sin(phi -0.5*PID + alpha2));
                draw_circle_arc(x2, y2, MU, PI + phi - 0.5*PID + alpha2, DPI - 2.0*alpha2, NSEG); 
                draw_line(x2 - MU*cos(phi -0.5*PID - alpha2), 
                          y2 - MU*sin(phi -0.5*PID - alpha2), 
                          LAMBDA*cos(phi+alpha), LAMBDA*sin(phi+alpha)); 
                draw_circle_arc(0.0, 0.0, LAMBDA, phi + alpha, DPI/(double)NPOLY - 2.0*alpha, NSEG);
            }
            
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
    double y, dy, dy_e, rgb[3], value, lum, amp, dy_phase;
    
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
    dy_phase = 1.0/((double)(jmax - jmin));
    
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
            case (P_ENERGY_FLUX):
            {
                value = dy_e*(double)(j - jmin)*100.0/E_SCALE;
                if (COLOR_PALETTE >= COL_TURBO) color_scheme_asym_palette(COLOR_SCHEME, COLOR_PALETTE, value, 1.0, 1, rgb);
                else color_scheme_palette(COLOR_SCHEME, COLOR_PALETTE, value, 1.0, 1, rgb);
//                 value = min + 1.0*dy*(double)(j - jmin);
//                 amp = 0.7*color_amplitude_linear(value, 1.0, 1);
//                 while (amp > 1.0) amp -= 2.0;
//                 while (amp < -1.0) amp += 2.0;
//                 amp_to_rgb(0.5*(1.0 + amp), rgb);
                break;
            }
            case (P_TOTAL_ENERGY_FLUX):
            {
//                 value = min + 1.0*dy*(double)(j - jmin);
//                 amp = 0.7*color_amplitude_linear(value, 1.0, 1);
//                 while (amp > 1.0) amp -= 2.0;
//                 while (amp < -1.0) amp += 2.0;
//                 amp_to_rgb(0.5*(1.0 + amp), rgb);
                value = dy_phase*(double)(j - jmin);
                color_scheme_palette(C_ONEDIM_LINEAR, COLOR_PALETTE, value, 1.0, 1, rgb);
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
    double y, dy, dy_e, rgb[3], value, lum, amp, dy_phase;
    
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
    dy_phase = 1.0/((double)(jmax - jmin));
    
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
                value = LOG_SCALE*log(dy_e*(double)(j - jmin)*100.0/E_SCALE);
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
            case (P_ENERGY_FLUX):
            {
                value = dy_e*(double)(j - jmin)*100.0/E_SCALE;
                if (COLOR_PALETTE >= COL_TURBO) color_scheme_asym_palette(COLOR_SCHEME, palette, value, 1.0, 1, rgb);
                else color_scheme_palette(COLOR_SCHEME, palette, value, 1.0, 1, rgb);
//                 value = min + 1.0*dy*(double)(j - jmin);
//                 amp = 0.7*color_amplitude_linear(value, 1.0, 1);
//                 while (amp > 1.0) amp -= 2.0;
//                 while (amp < -1.0) amp += 2.0;
//                 amp_to_rgb(0.5*(1.0 + amp), rgb);
                break;
            }
            case (P_TOTAL_ENERGY_FLUX):
            {
//                 value = min + 1.0*dy*(double)(j - jmin);
//                 amp = 0.7*color_amplitude_linear(value, 1.0, 1);
//                 while (amp > 1.0) amp -= 2.0;
//                 while (amp < -1.0) amp += 2.0;
//                 amp_to_rgb(0.5*(1.0 + amp), rgb);
                value = dy_phase*(double)(j - jmin);
                color_scheme_palette(C_ONEDIM_LINEAR, palette, value, 1.0, 1, rgb);
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
    double y, dy, dy_e, rgb[3], value, lum, amp, dy_phase;
    
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
    dy_phase = 1.0/((double)(jmax - jmin));
    
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
                value = LOG_SCALE*log(dy_e*(double)(j - jmin)*100.0/E_SCALE);
//                 printf("value = %.2lg\n", value);
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
            case (P_ENERGY_FLUX):
            {
                value = dy_e*(double)(j - jmin);
                if (COLOR_PALETTE >= COL_TURBO) color_scheme_asym_palette(COLOR_SCHEME, palette, value, 1.0, 1, rgb);
                else color_scheme_palette(COLOR_SCHEME, palette, value, 1.0, 1, rgb);
//                 value = min + 1.0*dy*(double)(j - jmin);
//                 amp = 0.7*color_amplitude_linear(value, 1.0, 1);
//                 while (amp > 1.0) amp -= 2.0;
//                 while (amp < -1.0) amp += 2.0;
//                 amp_to_rgb(0.5*(1.0 + amp), rgb);
                break;
            }
            case (P_TOTAL_ENERGY_FLUX):
            {
//                 value = min + 1.0*dy*(double)(j - jmin);
//                 amp = 0.7*color_amplitude_linear(value, 1.0, 1);
//                 while (amp > 1.0) amp -= 2.0;
//                 while (amp < -1.0) amp += 2.0;
//                 amp_to_rgb(0.5*(1.0 + amp), rgb);
                value = dy_phase*(double)(j - jmin);
                color_scheme_palette(C_ONEDIM_LINEAR, palette, value, 1.0, 1, rgb);
                break;
            }
            case (P_AVERAGE_ENERGY):
            {
                value = dy_e*(double)(j - jmin)*100.0/E_SCALE;
                if (COLOR_PALETTE >= COL_TURBO) color_scheme_asym_palette(COLOR_SCHEME, palette, value, 1.0, 1, rgb);
                else color_scheme_palette(COLOR_SCHEME, palette, value, 1.0, 1, rgb);
                break;
            }
            case (P_LOG_AVERAGE_ENERGY):
            {
                value = LOG_SCALE*log(dy_e*(double)(j - jmin)*100.0/E_SCALE);
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
            case (Z_EULER_VORTICITY):      
            {
                value = min + 1.0*dy*(double)(j - jmin);
                color_scheme_palette(COLOR_SCHEME, palette, 0.7*value, 1.0, 0, rgb);
                break;
            }
            case (Z_EULER_LOG_VORTICITY):
            {
                value = min + 1.0*dy*(double)(j - jmin);
                color_scheme_palette(COLOR_SCHEME, palette, 0.7*value, 1.0, 0, rgb);
                break;
            }
            case (Z_EULER_VORTICITY_ASYM):      
            {
                value = min + 1.0*dy*(double)(j - jmin);
                color_scheme_palette(COLOR_SCHEME, palette, 0.7*value, 1.0, 0, rgb);
                break;
            }
            case (Z_EULER_LPRESSURE):     
            {
                value = min + 1.0*dy*(double)(j - jmin);
                color_scheme_palette(COLOR_SCHEME, palette, 0.7*value, 1.0, 0, rgb);
                break;
            }
             case (Z_EULER_PRESSURE):     
            {
                value = min + 1.0*dy*(double)(j - jmin);
                color_scheme_palette(COLOR_SCHEME, palette, 0.7*value, 1.0, 0, rgb);
                break;
            }
            case (Z_EULER_DENSITY):      
            {
                value = min + 1.0*dy*(double)(j - jmin);
                color_scheme_palette(COLOR_SCHEME, palette, 0.7*value, 1.0, 0, rgb);
                break;
            }
            case (Z_EULER_SPEED):      
            {
                value = min + 1.0*dy*(double)(j - jmin);
                color_scheme_palette(COLOR_SCHEME, palette, 0.7*value, 1.0, 0, rgb);
                break;
            }
            case (Z_EULERC_VORTICITY):      
            {
                value = min + 1.0*dy*(double)(j - jmin);
                color_scheme_palette(COLOR_SCHEME, palette, 0.7*value, 1.0, 0, rgb);
                break;
            }
             case (P_3D_LOG_MEAN_ENERGY):
            {
                value = LOG_SCALE*dy_e*(double)(j - jmin)*100.0/E_SCALE;
                color_scheme_asym_palette(COLOR_SCHEME, palette, value, 1.0, 1, rgb);
                break;
            }
            default:
            {
                value = min + 1.0*dy*(double)(j - jmin);
                color_scheme_palette(COLOR_SCHEME, palette, 0.7*value, 1.0, 0, rgb);
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


void draw_circular_color_scheme_palette_fade(double x1, double y1, double radius, int plot, double min, double max, int palette, int fade, double fade_value)
{
    int j, k;
    double x, y, dy, dy_e, dy_phase, rgb[3], value, lum, amp, dphi, pos[2], phi, xy[2], zscale = 0.85;
    
//     glBegin(GL_TRIANGLE_FAN);
    xy_to_pos(x1, y1, xy);
    glVertex2d(xy[0], xy[1]);
    dy = (max - min)/360.0;
    dy_e = max/360.0;
    dy_phase = 1.0/360.0;
    dphi = DPI/360.0;
    
    for (j = 0; j <= 360; j++)
    {
        switch (plot) {
            case (P_AMPLITUDE):
            {
                value = min + 1.0*dy*(double)(j);
                color_scheme_palette(COLOR_SCHEME, palette, value, 1.0, 1, rgb);
                break;
            }
            case (P_ENERGY):
            {
                value = dy_e*(double)(j)*100.0/E_SCALE;
                if (COLOR_PALETTE >= COL_TURBO) color_scheme_asym_palette(COLOR_SCHEME, palette, value, 1.0, 1, rgb);
                else color_scheme_palette(COLOR_SCHEME, palette, value, 1.0, 1, rgb);
                break;
            }
            case (P_MEAN_ENERGY):
            {
                value = dy_e*(double)(j)*100.0/E_SCALE;
                if (COLOR_PALETTE >= COL_TURBO) color_scheme_asym_palette(COLOR_SCHEME, palette, value, 1.0, 1, rgb);
                else color_scheme_palette(COLOR_SCHEME, palette, value, 1.0, 1, rgb);
                break;
            }
            case (P_LOG_ENERGY):
            {
                value = LOG_SCALE*log(dy_e*(double)(j)*100.0/E_SCALE);
                color_scheme_palette(COLOR_SCHEME, palette, value, 1.0, 1, rgb);
                break;
            }
            case (P_LOG_MEAN_ENERGY):
            {
                value = LOG_SHIFT + LOG_SCALE*log(dy_e*(double)(j)*100.0/E_SCALE);
                color_scheme_palette(COLOR_SCHEME, palette, value, 1.0, 1, rgb);
                break;
            }
            case (P_ENERGY_FLUX):
            {
                value = dy_e*(double)(j);
                if (COLOR_PALETTE >= COL_TURBO) color_scheme_asym_palette(COLOR_SCHEME, palette, value, 1.0, 1, rgb);
                else color_scheme_palette(COLOR_SCHEME, palette, value, 1.0, 1, rgb);
                break;
            }
            case (P_TOTAL_ENERGY_FLUX):
            {
                value = dy_phase*(double)(j);
                color_scheme_palette(C_ONEDIM_LINEAR, palette, value, 1.0, 1, rgb);
                break;
            }
            case (P_AVERAGE_ENERGY):
            {
                value = dy_e*(double)(j)*100.0/E_SCALE;
                if (COLOR_PALETTE >= COL_TURBO) color_scheme_asym_palette(COLOR_SCHEME, palette, value, 1.0, 1, rgb);
                else color_scheme_palette(COLOR_SCHEME, palette, value, 1.0, 1, rgb);
                break;
            }
            case (P_LOG_AVERAGE_ENERGY):
            {
                value = LOG_SCALE*log(dy_e*(double)(j)*100.0/E_SCALE);
                color_scheme_palette(COLOR_SCHEME, palette, value, 1.0, 1, rgb);
                break;
            }
            case (P_PHASE):
            {
                value = min + 1.0*dy*(double)(j);
                amp = 0.5*color_amplitude_linear(value, 1.0, 1);
                while (amp > 1.0) amp -= 2.0;
                while (amp < -1.0) amp += 2.0;
                amp_to_rgb(0.5*(1.0 + amp), rgb);
                break;
            }
            case (Z_EULER_VORTICITY):      
            {
                value = min + 1.0*dy*(double)(j);
                color_scheme_palette(COLOR_SCHEME, palette, 0.7*value, 1.0, 0, rgb);
                break;
            }
            case (Z_EULER_LOG_VORTICITY):
            {
                value = min + 1.0*dy*(double)(j);
                color_scheme_palette(COLOR_SCHEME, palette, 0.7*value, 1.0, 0, rgb);
                break;
            }
            case (Z_EULER_VORTICITY_ASYM):      
            {
                value = min + 1.0*dy*(double)(j);
                color_scheme_palette(COLOR_SCHEME, palette, 0.7*value, 1.0, 0, rgb);
                break;
            }
            case (Z_EULER_LPRESSURE):     
            {
                value = min + 1.0*dy*(double)(j);
                color_scheme_palette(COLOR_SCHEME, palette, 0.7*value, 1.0, 0, rgb);
                break;
            }
             case (Z_EULER_PRESSURE):     
            {
                value = min + 1.0*dy*(double)(j);
                color_scheme_palette(COLOR_SCHEME, palette, 0.7*value, 1.0, 0, rgb);
                break;
            }
            case (Z_EULER_DENSITY):      
            {
                value = min + 1.0*dy*(double)(j);
                color_scheme_palette(COLOR_SCHEME, palette, 0.7*value, 1.0, 0, rgb);
                break;
            }
            case (Z_EULER_SPEED):      
            {
                value = min + 1.0*dy*(double)(j);
                color_scheme_palette(COLOR_SCHEME, palette, 0.7*value, 1.0, 0, rgb);
                break;
            }
            case (Z_EULERC_VORTICITY):      
            {
                value = min + 1.0*dy*(double)(j);
                color_scheme_palette(COLOR_SCHEME, palette, 0.7*value, 1.0, 0, rgb);
                break;
            }
            case (Z_SWATER_DIRECTION_SPEED):
            {
                value = dy_phase*(double)(j) - 0.5*PHASE_SHIFT + 0.5;
                if (value > 1.0) value -= 1.0;
                color_scheme_palette(C_ONEDIM_LINEAR, palette, value, 1.0, 1, rgb);
                break;
            }
            default:
            {
                value = min + 1.0*dy*(double)(j);
                color_scheme_palette(COLOR_SCHEME, palette, 0.7*value, 1.0, 0, rgb);
                break;
            }
        }
        if (fade) for (k=0; k<3; k++) rgb[k] *= fade_value;
        glColor3f(rgb[0], rgb[1], rgb[2]);
        
        glBegin(GL_TRIANGLE_FAN);
        xy_to_pos(x1, y1, xy);
        glVertex2d(xy[0], xy[1]);
        xy_to_pos(x1 + radius*cos(dphi*(double)(j-1)), y1 + zscale*radius*sin(dphi*(double)(j-1)), xy);
        glVertex2d(xy[0], xy[1]);
        xy_to_pos(x1 + radius*cos(dphi*(double)j), y1 + zscale*radius*sin(dphi*(double)j), xy);
        glVertex2d(xy[0], xy[1]);
        glEnd ();
    }
    xy_to_pos(x1 + radius*cos(dphi), y1 + radius*sin(dphi), xy);
//     glVertex2d(xy[0], xy[1]);
//     glEnd ();
    
    if (fade) glColor3f(fade_value, fade_value, fade_value);
    else glColor3f(1.0, 1.0, 1.0);
    glLineWidth(BOUNDARY_WIDTH*3/2);
    glEnable(GL_LINE_SMOOTH);
    
    dphi = DPI/(double)NSEG;
    glBegin(GL_LINE_LOOP);
    for (j = 0; j < NSEG; j++)
    {               
        xy_to_pos(x1 + radius*cos(dphi*(double)j), y1 + zscale*radius*sin(dphi*(double)j), xy);
        glVertex2d(xy[0], xy[1]);
        glVertex2d(xy[0], xy[1]);
    }
    glEnd ();
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

void print_frequency(double phase_shift, int fade, double fade_value)
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
    sprintf(message, "Frequency %.2f", 25.0*phase_shift);
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

double oscillating_bc(int time, int j)
{
    int ij[2], jmin, jmax;
    double t, t1, phase, a, envelope, omega, amp, dist2;
    
    switch (OSCILLATION_SCHEDULE)
    {
        case (OSC_PERIODIC): 
        {
            if (OSCIL_LEFT_YSHIFT > 0.0)
                t = (double)time*OMEGA - (double)j*OSCIL_LEFT_YSHIFT/(double)NY;
            else 
                t = (double)time*OMEGA + (double)(NY-j)*OSCIL_LEFT_YSHIFT/(double)NY;
            if (t < 0.0) return(0.0);
            else return(AMPLITUDE*cos(t)*exp(-(double)t*DAMPING));
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
            a = sqrt(INITIAL_VARIANCE)/OMEGA;
            phase = AMPLITUDE*cos(t);
            envelope = exp(-(t-INITIAL_SHIFT)*(t-INITIAL_SHIFT)/(a*a))*sqrt(DPI/a);
            return(phase*envelope);
        }
        case (OSC_WAVE_PACKETS):
        {
            t = (double)time*OMEGA;
            t1 = t - (double)((int)(t/WAVE_PACKET_SHIFT + 0.5))*WAVE_PACKET_SHIFT;
            a = sqrt(INITIAL_VARIANCE)/OMEGA;
            phase = AMPLITUDE*cos(t);
            envelope = exp(-(t1-INITIAL_SHIFT)*(t1-INITIAL_SHIFT)/(a*a))*sqrt(DPI/a);
            return(phase*envelope);
        }
        case (OSC_CHIRP):
        {
//             a = 0.25;
            t = (double)time*OMEGA;
            phase = t + ACHIRP*t*t;
            return(AMPLITUDE*sin(phase)*exp(-phase*DAMPING));
        }
        case (OSC_BEAM):
        {
            xy_to_ij(0.0, OSCIL_YMAX, ij);
            jmax = ij[1];
            jmin = NY - jmax;
            if (j < jmin) dist2 = (double)(jmin-1-j)*(YMAX-YMIN)/(double)NY;                
            else if (j > jmax) dist2 = (double)(j - jmax)*(YMAX-YMIN)/(double)NY;
            else dist2 = 0.0;
            dist2 = 500.0*dist2*dist2;
            t = (double)time*OMEGA;
            amp = AMPLITUDE*cos(t)*exp(-(double)t*DAMPING);
            amp *= exp(-dist2/INITIAL_VARIANCE);
            return(amp);
        }
        case (OSC_BEAM_GAUSSIAN):
        {
            dist2 = (double)(j - NY/2)*(YMAX-YMIN)/(double)NY;  
            dist2 = dist2*dist2;
            t = (double)time*OMEGA;
            amp = AMPLITUDE*cos(t)*exp(-(double)t*DAMPING);
            amp *= exp(-dist2/INITIAL_VARIANCE);
            return(amp);
        }
        case (OSC_BEAM_SINE):
        {
            xy_to_ij(0.0, OSCIL_YMAX, ij);
            jmax = ij[1];
            jmin = NY - jmax;
            dist2 = (double)(j - NY/2)/(double)(jmax - NY/2); 
            if (j > jmax) return(0.0);
            if (j < jmin) return(0.0);
            t = (double)time*OMEGA;
            amp = AMPLITUDE*cos(t)*exp(-(double)t*DAMPING);
            amp *= cos(PID*dist2);
            return(amp);
        }
        case (OSC_BEAM_TWOPERIODS):
        {
            xy_to_ij(0.0, OSCIL_YMAX, ij);
            jmax = ij[1];
            jmin = NY - jmax;
            if (j < jmin) dist2 = (double)(jmin-1-j)*(YMAX-YMIN)/(double)NY;                
            else if (j > jmax) dist2 = (double)(j - jmax)*(YMAX-YMIN)/(double)NY;
            else dist2 = 0.0;
            dist2 = 500.0*dist2*dist2;
            t = (double)time*OMEGA;
            amp = AMPLITUDE*(0.5*cos(t) + 0.5*cos(sqrt(7.0)*t))*exp(-(double)t*DAMPING);
            amp *= exp(-dist2/INITIAL_VARIANCE);
            return(amp);
        }
        case (OSC_TWO_WAVES):
        {
            t = (double)time*OMEGA + (double)(NY-j)*OSCIL_LEFT_YSHIFT/(double)NY;
            t1 = (double)time*OMEGA + (double)(j)*OSCIL_LEFT_YSHIFT/(double)NY;    
            if (j < NY/2) amp = cos(t);
            else amp = cos(t1);
            return(AMPLITUDE*amp);
        }
        case (OSC_TWO_WAVES_ADDED):
        {
            t = (double)time*OMEGA + (double)(NY-j)*OSCIL_LEFT_YSHIFT/(double)NY;
            t1 = (double)time*OMEGA + (double)(j)*OSCIL_LEFT_YSHIFT/(double)NY;    
            amp = 0.0;
            if (t > 0.0) amp += cos(t);
            if (t1 > 0.0) amp += cos(t1);
            return(AMPLITUDE*amp);
        }
    }
}

void init_ior_2d(short int *xy_in[NX], double *tcc_table[NX], double *tgamma_table[NX], double ior_angle)
/* compute variable index of refraction */
/* should be at some point merged with 3D version in suv_wave_3d.c */
{
    int i, j, k, n, inlens, ncircles;
    double courant2 = COURANT*COURANT, courantb2 = COURANTB*COURANTB, lambda1, mu1, courant_med, courant_med2;
    double a, u, v, u1, x, y, xy[2], norm2, speed, r2, c, salpha, h, ll, ca, sa, x1, y1, dx, dy, sum, sigma, x0, y0, rgb[3];
    double xc[NGRIDX*NGRIDY], yc[NGRIDX*NGRIDY], height[NGRIDX*NGRIDY];
    static double xc_stat[NGRIDX*NGRIDY], yc_stat[NGRIDX*NGRIDY], sigma_stat;
    static int first = 1;
    t_circle circles[NMAXCIRCLES];
    
    rgb[0] = 1.0;
    rgb[1] = 1.0;
    rgb[2] = 1.0;
    
    if (VARIABLE_IOR)
    {
        switch (IOR) {
            case (IOR_MANDELBROT):
            {
                #pragma omp parallel for private(i,j)
                for (i=0; i<NX; i++){
                    for (j=0; j<NY; j++){
                        ij_to_xy(i, j, xy);
                        x = xy[0];
                        y = xy[1];
                        u = 0.0;
                        v = 0.0;
                        k = 0;
                        while ((k<MANDELLEVEL)&&(u*u+v*v < 1000.0*MANDELLIMIT))
                        {
                            u1 = u*u - v*v + x;
                            v = 2.0*u*v + y;
                            u = u1;
                            k++;
                        }
                        norm2 = u*u + v*v;
                        if (norm2 < MANDELLIMIT)
                        {
//                             tc[i*NY+j] = COURANT;
                            tcc_table[i][j] = courant2;
                            tgamma_table[i][j] = GAMMA;
                        }
                        else 
                        {
                            speed = 1.0 + MANDEL_IOR_SCALE*log(1.0 + norm2/MANDELLIMIT);
                            if (speed < 0.01) speed = 0.01;
                            tcc_table[i][j] = courant2*speed;
//                             tc[i*NY+j] = COURANT*sqrt(speed);
                            tgamma_table[i][j] = GAMMA;
                        }
                    }
                }
                break;
            }
            case (IOR_MANDELBROT_LIN):
            {
                #pragma omp parallel for private(i,j)
                for (i=0; i<NX; i++){
                    for (j=0; j<NY; j++){
                        ij_to_xy(i, j, xy);
                        x = xy[0];
                        y = xy[1];
                        u = 0.0;
                        v = 0.0;
                        k = 0;
                        while ((k<MANDELLEVEL)&&(u*u+v*v < 1000.0*MANDELLIMIT))
                        {
                            u1 = u*u - v*v + x;
                            v = 2.0*u*v + y;
                            u = u1;
                            k++;
                        }
                        if (k >= MANDELLEVEL)
                        {
//                             tc[i*NY+j] = COURANT;
                            tcc_table[i][j] = courant2;
                            tgamma_table[i][j] = GAMMA;
                        }
                        else 
                        {
                            speed = (double)k/(double)MANDELLEVEL;
                            if (speed < 1.0e-10) speed = 1.0e-10;
                            else if (speed > 10.0) speed = 10.0;
                            tcc_table[i][j] = courant2*speed;
//                             tc[i*NY+j] = COURANT*sqrt(speed);
                            tgamma_table[i][j] = GAMMA;
                        }
                    }
                }
                break;
            }
            case (IOR_EARTH):
            {
                for (i=0; i<NX; i++){
                    for (j=0; j<NY; j++){
                        ij_to_xy(i, j, xy);
                        r2 = xy[0]*xy[0] + xy[1]*xy[1];
                        if (r2 > 1.0) c = 0.0;
                        else if (r2 < 0.25*0.25) c = 0.8*COURANT;
                        else if (r2 < 0.58*0.58) c = COURANT*(0.68 - 0.55*r2);
                        else c = COURANT*(1.3 - 0.9*r2);
//                         tc[i*NY+j] = c;
                        tcc_table[i][j] = c;
                        tgamma_table[i][j] = GAMMA;
                    }
                }
                break;
            }
            case (IOR_EXPLO_LENSING):
            {
                salpha = DPI/(double)NPOLY;
//                 lambda1 = LAMBDA;
//                 mu1 = LAMBDA;
                lambda1 = 0.5*LAMBDA;
                mu1 = 0.5*LAMBDA;
                h = lambda1*tan(PI/(double)NPOLY);
                if (h < mu1) ll = sqrt(mu1*mu1 - h*h);
                else ll = 0.0;
                
//                 #pragma omp parallel for private(i,j)
                for (i=0; i<NX; i++){
                    for (j=0; j<NY; j++) if (xy_in[i*NY+j]) {
                        ij_to_xy(i, j, xy);
                        x = xy[0];
                        y = xy[1];
                        inlens = 0;
                        for (k=0; k<NPOLY; k++)
                        {
                            ca = cos(((double)k+0.5)*salpha + APOLY*PID);
                            sa = sin(((double)k+0.5)*salpha + APOLY*PID);
                            x1 = x*ca + y*sa;
                            y1 = -x*sa + y*ca;
                            if ((module2(x1 - lambda1 - ll, y1) < mu1)&&(module2(x1 - lambda1 + ll, y1) < mu1)) inlens = 1; 
                        }
                        if (inlens) c = COURANTB;
                        else c = COURANT;
//                         tc[i*NY+j] = c;
                        tcc_table[i][j] = c*c;
                        tgamma_table[i][j] = GAMMA;
                    }
                    else
                    {
//                         tc[i*NY+j] = 0.0;
                        tcc_table[i][j] = 0.0;
                        tgamma_table[i][j] = 0.0;
                    }
                }
                break;
            }
            case (IOR_PERIODIC_WELLS):
            {
                dx = (XMAX - XMIN)/(double)NGRIDX;
                dy = (YMAX - YMIN)/(double)NGRIDY;
                sigma = 0.2*dx*dx;
                for (i=0; i<NGRIDX; i++)
                    for (j=0; j<NGRIDY; j++)
                    {
                        
                        n = j*NGRIDX + i;
                        xc[n] = XMIN + dx*((double)i + 0.5);
                        yc[n] = YMIN + dy*((double)j + 0.5);
                        if (j%2 == 1) yc[n] += 0.5*dx;
                    }
                    
//                 #pragma omp parallel for private(i,j)
                for (i=0; i<NX; i++){
                    for (j=0; j<NY; j++){
                        ij_to_xy(i, j, xy);
                        x = xy[0];
                        y = xy[1];
                        sum = 0.0;
                        for (n = 0; n<NGRIDX*NGRIDY; n++)
                        {
                            r2 = (x - xc[n])*(x - xc[n]) + (y - yc[n])*(y - yc[n]);
                            sum += exp(-r2/(sigma));
                        }
//                         sum = tanh(sum);
//                         printf("%.3lg\n", sum);
                        tcc_table[i][j] = COURANT*sum + COURANTB*(1.0-sum);
                        tgamma_table[i][j] = GAMMA;
                    }
                }
                break;
            }
            case (IOR_RANDOM_WELLS):
            {
                dx = (XMAX - XMIN)/(double)NGRIDX;
                dy = (YMAX - YMIN)/(double)NGRIDY;
                sigma = 0.2*dx*dx;
                for (i=0; i<NGRIDX; i++)
                    for (j=0; j<NGRIDY; j++)
                    {
                        
                        n = j*NGRIDX + i;
                        xc[n] = XMIN + dx*((double)i + 0.5 + 0.1*gaussian());
                        yc[n] = YMIN + dy*((double)j + 0.5 + 0.1*gaussian());
//                         if (j%2 == 1) yc[n] += 0.5*dx;
                        height[n] = 0.5 + 0.5*gaussian();
                        if (height[n] > 1.0) height[n] = 1.0;
                        if (height[n] < 0.0) height[n] = 0.0;
                    }
                    
//                 #pragma omp parallel for private(i,j)
                for (i=0; i<NX; i++)
                {
                    if (i%100 == 0) printf("Computing potential for column %i of %i\n", i, NX);
                    for (j=0; j<NY; j++){
                        ij_to_xy(i, j, xy);
                        x = xy[0];
                        y = xy[1];
                        sum = 0.0;
                        for (n = 0; n<NGRIDX*NGRIDY; n++)
                        {
                            r2 = (x - xc[n])*(x - xc[n]) + (y - yc[n])*(y - yc[n]);
                            sum += exp(-r2/(sigma))*height[n];
                        }
                        sum = tanh(sum);
//                         printf("%.3lg\n", sum);
                        tcc_table[i][j] = COURANT*sum + COURANTB*(1.0-sum);
                        tgamma_table[i][j] = GAMMA;
                    }
                }
                break;
            }
            case (IOR_PERIODIC_WELLS_ROTATING):
            {
                if (first)
                {
                    dx = (XMAX - XMIN)/(double)NGRIDX;
                    dy = (YMAX - YMIN)/(double)NGRIDY;
                    sigma_stat = 0.2*dx*dx;
                    for (i=0; i<NGRIDX; i++)
                        for (j=0; j<NGRIDY; j++)
                        {
                            n = j*NGRIDX + i;
                            xc_stat[n] = XMIN + dx*((double)i + 0.5);
                            yc_stat[n] = YMIN + dy*((double)j + 0.5);
                            if (j%2 == 1) yc_stat[n] += 0.5*dx;
                        }
                    first = 0;
                }
                ca = cos(ior_angle);
                sa = sin(ior_angle);
                for (n=0; n<NGRIDX*NGRIDY; n++)
                {
                    xc[n] = xc_stat[n]*ca + yc_stat[n]*sa;
                    yc[n] = -xc_stat[n]*sa + yc_stat[n]*ca;
//                     printf("center[%i] at (%.5lg, %.5lg)\n", n, xc[n], yc[n]);
//                     draw_colored_circle(xc[n], yc[n], 0.05, 10, rgb);
                }
                    
//                 #pragma omp parallel for private(i,j)
                for (i=0; i<NX; i++){
//                     if (i%100 == 0) printf("initializing column %i of %i\n", i, NX);
                    for (j=0; j<NY; j++){
                        ij_to_xy(i, j, xy);
                        x = xy[0];
                        y = xy[1];
                        sum = 0.0;
                        for (n = 0; n<NGRIDX*NGRIDY; n++)
                        {
                            r2 = (x - xc[n])*(x - xc[n]) + (y - yc[n])*(y - yc[n]);
                            sum += exp(-r2/(sigma_stat));
                        }
//                         sum = tanh(sum);
//                         printf("%.3lg\n", sum);
                        tcc_table[i][j] = COURANT*sum + COURANTB*(1.0-sum);
                        tgamma_table[i][j] = GAMMA;
                    }
                }
                break;
            }
            case (IOR_PERIODIC_WELLS_ROTATING_LARGE):
            {
                if (first)
                {
                    dx = (1.5*XMAX - 1.5*XMIN)/(double)NGRIDX;
                    dy = (1.5*YMAX - 1.5*YMIN)/(double)NGRIDY;
                    sigma_stat = 0.2*dx*dx;
                    for (i=0; i<NGRIDX; i++)
                        for (j=0; j<NGRIDY; j++)
                        {
                            n = j*NGRIDX + i;
                            xc_stat[n] = 1.5*XMIN + dx*((double)i + 0.5);
                            yc_stat[n] = 1.5*YMIN + dy*((double)j + 0.5);
                            if (j%2 == 1) yc_stat[n] += 0.5*dx;
                        }
                    first = 0;
                }
                ca = cos(ior_angle);
                sa = sin(ior_angle);
                for (n=0; n<NGRIDX*NGRIDY; n++)
                {
                    xc[n] = xc_stat[n]*ca + yc_stat[n]*sa;
                    yc[n] = -xc_stat[n]*sa + yc_stat[n]*ca;
                }
                    
//                 #pragma omp parallel for private(i,j)
                for (i=0; i<NX; i++){
//                     if (i%100 == 0) printf("initializing column %i of %i\n", i, NX);
                    for (j=0; j<NY; j++){
                        ij_to_xy(i, j, xy);
                        x = xy[0];
                        y = xy[1];
                        sum = 0.0;
                        for (n = 0; n<NGRIDX*NGRIDY; n++)
                        {
                            r2 = (x - xc[n])*(x - xc[n]) + (y - yc[n])*(y - yc[n]);
                            sum += exp(-r2/(sigma_stat));
                        }
                        tcc_table[i][j] = COURANT*sum + COURANTB*(1.0-sum);
                        tgamma_table[i][j] = GAMMA;
                    }
                }
                break;
            }
            case (IOR_POISSON_WELLS):
            {
                ncircles = init_circle_config_pattern(circles, C_POISSON_DISC);
                for (n = 0; n<ncircles; n++) 
                {
                    height[n] = 0.5 + 0.5*gaussian();
                    if (height[n] > 1.0) height[n] = 1.0;
                    if (height[n] < 0.0) height[n] = 0.0;
                }
                
                for (n = 0; n<ncircles; n++) printf("Circle %i at (%.3lg, %.3lg) height %.3lg\n", n, circles[n].xc, circles[n].yc, height[n]);
                
                sigma = 0.2*(XMAX - XMIN)*(YMAX - YMIN)/(double)ncircles;
//                 sigma = MU*MU;
                
                for (i=0; i<NX; i++)
                {
                    if (i%100 == 0) printf("Computing potential for column %i of %i\n", i, NX);
                    for (j=0; j<NY; j++){
                        ij_to_xy(i, j, xy);
                        x = xy[0];
                        y = xy[1];
                        sum = 0.0;
                        for (n = 0; n<ncircles; n++)
                        {
                            r2 = (x - circles[n].xc)*(x - circles[n].xc) + (y - circles[n].yc)*(y - circles[n].yc);
                            sum += exp(-r2/(sigma))*height[n];
                        }
                        sum = tanh(sum);
//                         printf("%.3lg\n", sum);
                        tcc_table[i][j] = COURANT*sum + COURANTB*(1.0-sum);
                        tgamma_table[i][j] = GAMMA;
                    }
                }
                
                break;
            }
            case (IOR_PPP_WELLS):
            {
                ncircles = init_circle_config_pattern(circles, C_RAND_POISSON);
                for (n = 0; n<ncircles; n++) 
                {
                    height[n] = 0.5 + 0.5*gaussian();
                    if (height[n] > 1.0) height[n] = 1.0;
                    if (height[n] < 0.0) height[n] = 0.0;
                }
                
                for (n = 0; n<ncircles; n++) printf("Circle %i at (%.3lg, %.3lg) height %.3lg\n", n, circles[n].xc, circles[n].yc, height[n]);
                
                sigma = 0.2*(XMAX - XMIN)*(YMAX - YMIN)/(double)ncircles;
//                 sigma = MU*MU;
                
                for (i=0; i<NX; i++)
                {
                    if (i%100 == 0) printf("Computing potential for column %i of %i\n", i, NX);
                    for (j=0; j<NY; j++){
                        ij_to_xy(i, j, xy);
                        x = xy[0];
                        y = xy[1];
                        sum = 0.0;
                        for (n = 0; n<ncircles; n++)
                        {
                            r2 = (x - circles[n].xc)*(x - circles[n].xc) + (y - circles[n].yc)*(y - circles[n].yc);
                            sum += exp(-r2/(sigma))*height[n];
                        }
                        sum = tanh(sum);
//                         printf("%.3lg\n", sum);
                        tcc_table[i][j] = COURANT*sum + COURANTB*(1.0-sum);
                        tgamma_table[i][j] = GAMMA;
                    }
                }
                
                break;
            }
            case (IOR_LENS_WALL):
            {
                printf("Initializing IOR_LENS_WALL\n");
                for (i=0; i<NX; i++){
                    printf("i = %i\n", i);
                    for (j=0; j<NY; j++){
                        ij_to_xy(i, j, xy);
                        x = xy[0];
                        y = xy[1];
                        if ((vabs(x) < 0.05)&&(vabs(y) > MU)) 
                        {
                            tcc_table[i][j] = 0.0;
                            tgamma_table[i][j] = GAMMA;
                        }
                        else 
                        {
                            if (xy_in[i][j] != 0) 
                            {
                                tcc_table[i][j] = courant2;
                                tgamma_table[i][j] = GAMMA;
                            }
                            else if (xy_in[i][j] == 0) 
                            {
                                tcc_table[i][j] = courantb2;
                                tgamma_table[i][j] = GAMMAB;
                            }
                            else
                            {
                                tcc_table[i][j] = 0.0;
                                tgamma_table[i][j] = 1.0;
                            }                                
                        }
                    }
                }
                
                break;
            }
            case (IOR_LENS_OBSTACLE):
            {
                printf("Initializing IOR_LENS_OBSTACLE\n");
                for (i=0; i<NX; i++){
                    for (j=0; j<NY; j++){
                        ij_to_xy(i, j, xy);
                        x = xy[0];
                        y = xy[1];
                        if ((vabs(x) < 0.05)&&(vabs(y) < MU))
                        {
                            tcc_table[i][j] = 0.0;
                            tgamma_table[i][j] = GAMMA;
                        }
                        else 
                        {
                            if (xy_in[i][j] != 0) 
                            {
                                tcc_table[i][j] = courant2;
                                tgamma_table[i][j] = GAMMA;
                            }
                            else 
                            {
                                tcc_table[i][j] = courantb2;
                                tgamma_table[i][j] = GAMMAB;
                            }
                        }
                    }
                }
                break;
            }
            case (IOR_LENS_CONCAVE):
            {
                a = LAMBDA + 0.1 - sqrt(LAMBDA*LAMBDA - MU*MU); 
                for (i=0; i<NX; i++){
                    for (j=0; j<NY; j++){
                        ij_to_xy(i, j, xy);
                        x = xy[0];
                        y = xy[1];
                        if ((vabs(x) < a)&&(vabs(y) > MU))
                        {
                            tcc_table[i][j] = 0.0;
                            tgamma_table[i][j] = GAMMA;
                        }
                        else 
                        {
                            if (xy_in[i][j] != 0) 
                            {
                                tcc_table[i][j] = courant2;
                                tgamma_table[i][j] = GAMMA;
                            }
                            else 
                            {
                                tcc_table[i][j] = courantb2;
                                tgamma_table[i][j] = GAMMAB;
                            }
                        }
                    }
                }
                break;
            }
            case (IOR_LENS_CONVEX_CONCAVE):
            {
                a = LAMBDA + 0.1 - sqrt(LAMBDA*LAMBDA - MU*MU) + 0.6; 
                for (i=0; i<NX; i++){
                    for (j=0; j<NY; j++){
                        ij_to_xy(i, j, xy);
                        x = xy[0];
                        y = xy[1];
                        if ((vabs(x) < a)&&(vabs(y) > MU))
                        {
                            tcc_table[i][j] = 0.0;
                            tgamma_table[i][j] = GAMMA;
                        }
                        else 
                        {
                            if (xy_in[i][j] != 0) 
                            {
                                tcc_table[i][j] = courant2;
                                tgamma_table[i][j] = GAMMA;
                            }
                            else 
                            {
                                tcc_table[i][j] = courantb2;
                                tgamma_table[i][j] = GAMMAB;
                            }
                        }
                    }
                }
                break;
            }
            case (IOR_TREE):
            {
                for (i=0; i<NX; i++){
                    for (j=0; j<NY; j++){
                        ij_to_xy(i, j, xy);
                        x = xy[0];
                        y = xy[1];
                        switch (xy_in[i][j]) {
                            case (0):
                            {
                                tcc_table[i][j] = courantb2;
                                tgamma_table[i][j] = GAMMAB;
                                break;
                            }
                            case (1):
                            {
                                tcc_table[i][j] = courant2;
                                tgamma_table[i][j] = GAMMA;
                                break;
                            }
                            case (2):
                            {
                                tcc_table[i][j] = 0.0;
                                tgamma_table[i][j] = GAMMA;
                                break;
                            }
                            case (3):
                            {
                                tcc_table[i][j] = courant2;
                                tgamma_table[i][j] = GAMMA;
                                break;
                            }
                        }
                    }
                }
                break;
            }
            case (IOR_WAVE_GUIDES_COUPLED):
            {
                for (i=0; i<NX; i++){
                    for (j=0; j<NY; j++){
                        ij_to_xy(i, j, xy);
                        x = xy[0];
                        y = xy[1];
                        if (xy_in[i][j] == 1) 
                        {
                            tcc_table[i][j] = courant2;
                            tgamma_table[i][j] = GAMMA;
                        }
                        else if (vabs(x) > LAMBDA)
                        {
                            tcc_table[i][j] = 0.0;
                            tgamma_table[i][j] = 1.0;                            
                        }
                        else
                        {
                            tcc_table[i][j] = courantb2;
                            tgamma_table[i][j] = GAMMAB;
                        }
                    }
                }
                break;
            }
            case (IOR_WAVE_GUIDES_COUPLED_B):
            {
                for (i=0; i<NX; i++){
                    for (j=0; j<NY; j++){
                        ij_to_xy(i, j, xy);
                        x = xy[0];
                        y = xy[1];
                        if (xy_in[i][j] == 1) 
                        {
                            tcc_table[i][j] = courant2;
                            tgamma_table[i][j] = GAMMA;
                        }
                        else if ((vabs(x) > LAMBDA)&&(vabs(y) > 0.7))
                        {
                            tcc_table[i][j] = 0.0;
                            tgamma_table[i][j] = 1.0;                            
                        }
                        else
                        {
                            tcc_table[i][j] = courantb2;
                            tgamma_table[i][j] = GAMMAB;
                        }
                    }
                }
                break;
            }
            case (IOR_WAVE_GUIDE_COATED):
            {
                courant_med = 0.8*COURANT + 0.2*COURANTB;
                courant_med2 = courant_med*courant_med;
                for (i=0; i<NX; i++){
                    for (j=0; j<NY; j++){
//                         ij_to_xy(i, j, xy);
//                         x = xy[0];
//                         y = xy[1];
                        if (xy_in[i][j] == 1) 
                        {
                            tcc_table[i][j] = courant2;
                            tgamma_table[i][j] = GAMMA;
                        }
                        else if (xy_in[i][j] == 2)
                        {
                            tcc_table[i][j] = courant_med2;
                            tgamma_table[i][j] = GAMMA;                            
                        }
                        else
                        {
                            tcc_table[i][j] = courantb2;
                            tgamma_table[i][j] = GAMMAB;
                        }
                    }
                }
                break;
            }
            case (IOR_MICHELSON):
            {
                for (i=0; i<NX; i++){
                    for (j=0; j<NY; j++){
                        if (xy_in[i][j] == 0) 
                        {
                            tcc_table[i][j] = courant2;
                            tgamma_table[i][j] = GAMMA;
                        }
                        else if (xy_in[i][j] == 2)
                        {
                            tcc_table[i][j] = 0.0;
                            tgamma_table[i][j] = 1.0;                            
                        }
                        else
                        {
                            tcc_table[i][j] = courantb2;
                            tgamma_table[i][j] = GAMMAB;
                        }
                    }
                }
                break;
            }
            case (IOR_GRADIENT_INDEX_LENS):
            {
                /* focal distance is f = 1/(4*LAMBDA*n0*MU) */
                /* with n0 = COURANT/COURANTB */
                for (i=0; i<NX; i++){
                    for (j=0; j<NY; j++){
                        ij_to_xy(i, j, xy);
                        x = vabs(xy[0]);
                        y = vabs(xy[1]);
                        if ((x > LAMBDA)) 
                        {
                            tcc_table[i][j] = courant2;
                            tgamma_table[i][j] = GAMMA;
                        }
                        else if (y > 1.0)
                        {
                            tcc_table[i][j] = 0.0;
                            tgamma_table[i][j] = 0.0;                            
                        }
                        else
                        {
                            tcc_table[i][j] = courantb2/(1.0 - MU*y*y);
                            tgamma_table[i][j] = GAMMAB;
                        }
                    }
                }
                break;
            }
            case (IOR_GRADIENT_INDEX_LENS_B):
            {
                /* focal distance is f = 1/(4*LAMBDA*n0*MU) */
                /* with n0 = COURANT/COURANTB */
                for (i=0; i<NX; i++){
                    for (j=0; j<NY; j++){
                        ij_to_xy(i, j, xy);
                        x = vabs(xy[0]);
                        y = vabs(xy[1]);
                        if ((x > LAMBDA)) 
                        {
                            tcc_table[i][j] = courant2;
                            tgamma_table[i][j] = GAMMA;
                        }
                        else if (y > 1.0)
                        {
                            tcc_table[i][j] = 0.0;
                            tgamma_table[i][j] = 0.0;                            
                        }
                        else
                        {
                            speed = 1.0/(1.0 - MU*y*y);
                            speed *= speed;
                            tcc_table[i][j] = courantb2*speed;
                            tgamma_table[i][j] = GAMMAB;
                        }
                    }
                }
                break;
            }
            case (IOR_LINEAR_X_A):
            {
                /* Warning: Depending on COURANT and COURANTB */
                /* this may generate a wave speed of the form |a - bx| */
                /* Use IOR_LINEAR_X_B instead to avoid this */
                for (i=0; i<NX; i++){
                    for (j=0; j<NY; j++){
                        ij_to_xy(i, j, xy);
                        x = xy[0];
                        speed = COURANT*(1.0-x) + COURANTB*x;
                        tcc_table[i][j] = speed*speed;
                        tgamma_table[i][j] = GAMMA;
                    }
                }
                break;
            }
            case (IOR_LINEAR_X_B):
            {
                for (i=0; i<NX; i++){
                    for (j=0; j<NY; j++){
                        ij_to_xy(i, j, xy);
                        x = xy[0];
                        a = (x-XMIN)/(XMAX-XMIN);
                        speed = COURANT*(1.0-a) + COURANTB*a;
                        tcc_table[i][j] = speed*speed;
                        tgamma_table[i][j] = GAMMA;
                    }
                }
                break;
            }
            default:
            {
                for (i=0; i<NX; i++){
                    for (j=0; j<NY; j++){
//                         tc[i*NY+j] = COURANT;
                        tcc_table[i][j] = COURANT;
                        tgamma_table[i][j] = GAMMA;
                    }
                }
            }
        }
    }
    else
    {
        #pragma omp parallel for private(i,j)
        for (i=0; i<NX; i++){
            for (j=0; j<NY; j++){
                if (xy_in[i][j] != 0)
                {
//                     tc[i*NY+j] = COURANT;
                    tcc_table[i][j] = courant2;
                    if (xy_in[i][j] == 1) tgamma_table[i][j] = GAMMA;
                    else tgamma_table[i][j] = GAMMAB;
                }
                else if (TWOSPEEDS)
                {
//                     tc[i*NY+j] = COURANTB;
                    tcc_table[i][j] = courantb2;
                    tgamma_table[i][j] = GAMMAB;
                }
            }
        }
    }
}


double ior_angle_schedule(int i)
/* angle of rotation for variable index of refraction IOR_PERIODIC_WELLS_ROTATING */
{
    return(IOR_TOTAL_TURNS*DPI*(double)i/(double)NSTEPS);
}


int phased_array_schedule(int i)
/* returns time-dependent dephasing in phased array */
{
    int phase;
    
    phase = i*11/NSTEPS;
    
    switch (phase) {
        case (0): return(4);
        case (1): return(3);
        case (2): return(2);
        case (3): return(-2);
        case (4): return(-3);
        case (5): return(-4); 
        case (6): return(-3); 
        case (7): return(-2); 
        case (8): return(2); 
        case (9): return(3); 
        case (10): return(4); 
        default: return(4);
    }
    
}   

void init_wave_packets(t_wave_packet *packet, int radius)
/* initialise table of wave packets */
{
    int i, j, k, ij[2], nx, ny;
    double dx, dy;
    
    printf("Initialising wave packets\n");
    switch (WAVE_PACKET_SOURCE_TYPE) {
        case (WP_RANDOM1): 
        {
            nx = (int)sqrt((double)N_WAVE_PACKETS);
            ny = N_WAVE_PACKETS/nx;
            dx = 0.2*(XMAX - XMIN)/(double)nx;
            dy = 0.4*(YMAX - YMIN)/(double)ny;
            for (i=0; i<N_WAVE_PACKETS; i++)
            {
                j = i/nx;
                k = i%nx;
                packet[i].xc = XMIN + (double)(j+1)*dx + 0.5*dx*(double)rand()/RAND_MAX;
                packet[i].yc = (double)(k-ny/2)*dy + 0.5*dy*(double)rand()/RAND_MAX;
                packet[i].period = 20.0*(1.0 + 0.5*(double)rand()/RAND_MAX);
                packet[i].amp = INITIAL_AMP;
                packet[i].phase = DPI*(double)rand()/RAND_MAX;
                packet[i].var_envelope = 5.0e5;
                packet[i].time_shift = 10 + rand()%200;
                
                xy_to_ij(packet[i].xc, packet[i].yc, ij);
                if(ij[0] <= radius) ij[0] = radius+1;
                packet[i].ix = ij[0];
                packet[i].iy = ij[1]; 
            }
            break;
        }
        case (WP_RANDOM2): 
        {
            for (i=0; i<N_WAVE_PACKETS; i++)
            {
                packet[i].xc = XMIN + 0.15*(XMAX - XMIN)*(double)rand()/RAND_MAX;
                packet[i].yc = 0.4*(YMAX - YMIN)*((double)rand()/RAND_MAX - 0.5);
                packet[i].period = 50.0*(1.0 + 0.5*(double)rand()/RAND_MAX);
                packet[i].amp = INITIAL_AMP;
                packet[i].phase = DPI*(double)rand()/RAND_MAX;
                packet[i].var_envelope = 1500.0 + 500.0*(double)rand()/RAND_MAX;
                packet[i].time_shift = 10 + rand()%200;
                
                xy_to_ij(packet[i].xc, packet[i].yc, ij);
                if(ij[0] <= radius) ij[0] = radius+1;
                packet[i].ix = ij[0];
                packet[i].iy = ij[1]; 
            }
            break;
        }
        case (WP_PAIR):
        {
            for (i=0; i<2; i++)
            {
                packet[i].xc = -1.25;
                packet[i].yc = 0.06;
                if (i==1) packet[i].yc *= -1.0;
                packet[i].period = OSCILLATING_SOURCE_PERIOD;
                packet[i].amp = INITIAL_AMP;
                packet[i].phase = 0.0;
                packet[i].var_envelope = 100000.0;
                packet[i].time_shift = 0;
                
                xy_to_ij(packet[i].xc, packet[i].yc, ij);
                if(ij[0] <= radius) ij[0] = radius+1;
                packet[i].ix = ij[0];
                packet[i].iy = ij[1]; 
            }
            break;
        }
        case (WP_FIVE):
        {
            for (i=0; i<5; i++)
            {
                packet[i].xc = -1.75;
                packet[i].yc = 0.3*(double)(i-2);
                packet[i].period = OSCILLATING_SOURCE_PERIOD;
                packet[i].amp = INITIAL_AMP;
                packet[i].phase = 0.0;
                packet[i].var_envelope = 550.0;
                packet[i].time_shift = ((3*i)%5)*400;
                
                xy_to_ij(packet[i].xc, packet[i].yc, ij);
                if(ij[0] <= radius) ij[0] = radius+1;
                packet[i].ix = ij[0];
                packet[i].iy = ij[1]; 
            }
            break;
        }
    }
}

double wave_packet_height(double t, t_wave_packet packet, int type_envelope)
/* determines height of wave packet at time t */
{
    double cwave, envelope;
    
    cwave = packet.amp*cos(DPI*t/packet.period + packet.phase);
    
    switch (type_envelope) {
        case (WE_SINE): 
        {
            envelope = 0.1 + sin(DPI*t/packet.var_envelope);
            envelope = envelope*envelope;
            break;
        }
        case (WE_CUTOFF): 
        {
            if (t < (double)packet.var_envelope) envelope = 1.0;
            else envelope = 0.0;
            break;
        }
        case (WE_CONSTANT):
        {
            envelope = 1.0;
            break;
        }
    }
    
    return(cwave*envelope);
}

void add_wave_packets_locally(double *phi[NX], double *psi[NX], t_wave_packet *packet, int time, int radius, int add_period, int type_envelope)
/* add some wave packet sources - this local version leads to numerical artifacts */
{
    int i, ij[2], t, j, k, envelope, irad2, rmin2, rmax2;
    double cwave, wave_height, wave_height1, r2, variance;
    
    variance = (double)(radius*radius);
    rmin2 = radius*radius/2;
    rmax2 = radius*radius + 1;
    if (time%add_period == 0) for (i=0; i<N_WAVE_PACKETS; i++)
    {
        t = (double)(time - packet[i].time_shift);
        wave_height = wave_packet_height(t, packet[i], type_envelope);

        for (j=0; j<radius+1; j++) 
            for (k=0; k<radius+1; k++)
            {
                irad2 = j*j + k*k;
                if ((irad2 < rmax2)&&(irad2 > rmin2))
                {
                    r2 = (double)(irad2);
                    wave_height1 = wave_height*exp(-r2/variance);
                    phi[packet[i].ix + j][packet[i].iy + k] = wave_height1;
                    phi[packet[i].ix - j][packet[i].iy + k] = wave_height1;
                    phi[packet[i].ix + j][packet[i].iy - k] = wave_height1;
                    phi[packet[i].ix - j][packet[i].iy - k] = wave_height1;
                }
                else if (irad2 <= rmin2)
                {
                    phi[packet[i].ix + j][packet[i].iy + k] = 0.0;
                    phi[packet[i].ix - j][packet[i].iy + k] = 0.0;
                    phi[packet[i].ix + j][packet[i].iy - k] = 0.0;
                    phi[packet[i].ix - j][packet[i].iy - k] = 0.0;
                }
                
            }
//         printf("Adding wave packet %i with value %.3lg\n", i, wave_height); 
        
//         if (add==0)
        {
//             t -= 1.0/(double)NVID;
            t -= 0.1/(double)NVID;
            wave_height = wave_packet_height(t, packet[i], type_envelope);
            for (j=0; j<radius+1; j++) 
                for (k=0; k<radius+1; k++) 
                {
                    irad2 = j*j + k*k;
                    if ((irad2 < rmax2)&&(irad2 > rmin2))
                    {
                        r2 = (double)(irad2);
                        wave_height1 = wave_height*exp(-r2/variance);
                        psi[packet[i].ix + j][packet[i].iy + k] = wave_height1;
                        psi[packet[i].ix - j][packet[i].iy + k] = wave_height1;
                        psi[packet[i].ix + j][packet[i].iy - k] = wave_height1;
                        psi[packet[i].ix - j][packet[i].iy - k] = wave_height1;
                    }
                    else if (irad2 <= rmin2)
                    {
                        psi[packet[i].ix + j][packet[i].iy + k] = 0.0;
                        psi[packet[i].ix - j][packet[i].iy + k] = 0.0;
                        psi[packet[i].ix + j][packet[i].iy - k] = 0.0;
                        psi[packet[i].ix - j][packet[i].iy - k] = 0.0;
                    }
                }
        }
    }
}

void add_circular_wave_loc(double factor, double x, double y, double *phi[NX], double *psi[NX], short int * xy_in[NX], int xmin, int xmax, int ymin, int ymax)
/* add drop at (x,y) to the field with given prefactor */
{
    int i, j;
    double xy[2], dist2;
    
//     for (i=0; i<NX; i++)
//         for (j=0; j<NY; j++)
    #pragma omp parallel for private(i,j,xy,dist2)
    for (i=xmin; i<xmax; i++)
        for (j=ymin; j<ymax; j++)
        {
            ij_to_xy(i, j, xy);
            dist2 = (xy[0]-x)*(xy[0]-x) + (xy[1]-y)*(xy[1]-y);
            if ((xy_in[i][j])||(TWOSPEEDS)) 
                phi[i][j] += INITIAL_AMP*factor*exp(-dist2/INITIAL_VARIANCE)*cos(-sqrt(dist2)/INITIAL_WAVELENGTH);
        }
}

void add_wave_packets_globally(double *phi[NX], double *psi[NX], short int * xy_in[NX], t_wave_packet *packet, int time, int add_period, int type_envelope)
/* add some wave packet sources */
{
    int i, ij[2];
    double amp, t, omega, wave_height;
    static int xmin, xmax, ymin, ymax, first=1;
    
    if (first)
    {
        xy_to_ij(XMIN, YMIN, ij);
        xmin = ij[0];
        ymin = ij[1];
        xy_to_ij(XMIN + 0.15*(XMAX - XMIN), YMAX, ij);
        xmax = ij[0];
        ymax = ij[1];
        first = 0;
    }
    
    if (time%add_period == 0) for (i=0; i<N_WAVE_PACKETS; i++)
    {
        t = (double)(time - packet[i].time_shift);
        wave_height = wave_packet_height(t, packet[i], type_envelope);
        add_circular_wave_loc(wave_height, packet[i].xc, packet[i].yc, phi, psi, xy_in, xmin, xmax, ymin, ymax);
        
//         omega = DPI/packet[i].period;
//         amp = packet[i].amp*omega*sin(omega*t + packet[i].phase);
//         if (t > (double)packet[i].var_envelope) amp *= exp(-0.1*(t - (double)packet[i].var_envelope));
//         if (t < (double)packet[i].var_envelope + 100.0)
//             add_circular_wave_loc(amp, packet[i].xc, packet[i].yc, phi, psi, xy_in, xmin, xmax, ymin, ymax);
    }
}


void add_wave_packets(double *phi[NX], double *psi[NX], short int * xy_in[NX], t_wave_packet *packet, int time, int radius, int local, int add_period, int type_envelope)
/* add some wave packet sources */
{
    if (local) add_wave_packets_locally(phi, psi, packet, time, radius, add_period, type_envelope);
    else add_wave_packets_globally(phi, psi, xy_in, packet, time, add_period, type_envelope);
}

int old_source_schedule(int i)
{
    int mod;
    
    if (i < 200) return(0);
    mod = i%(10*OSCILLATING_SOURCE_PERIOD);
    if ((mod < 2*OSCILLATING_SOURCE_PERIOD)&&(rand()%3 < 2)) return(1);
    return(0);
}

void old_init_input_signal()
{
    int i;
 
    for (i=0; i<NSTEPS; i++) input_signal[i] = old_source_schedule(i);
}


int init_input_signal()
{
    int i, j, k, mod, tailcounter = 0;
    char text[200] = "..|--/-...|.-|-../..|--/-...|.-|-../", c;

//     char text[200] = "..|.../-|....|.|.-.|./.-|-.|-.--|-...|--- -..|-.--/---|..-|-/-|....|.|.-.|.", c;
//     char text[200] = "..|--/.-|-./..|--|.--.|.-.|---|...-|.|-../..-.|..|-...|.|.-.", c;
//     char text[200] = "...._.-_.--._.--._-.--__-._._.--__-.--_._.-_.-.__..---_-----_..---_....-", c;
            
    for (i=0; i<MESSAGE_INITIAL_TIME; i++) input_signal[i] = 0;
    
    j = 0;
    while (i < NSTEPS)
    {
        while (j < strlen(text))
        {
            c = text[j];
            
            switch(c) {
                case ('-'):
                {
                    for (k=0; k<MESSAGE_LDASH; k++) 
                    {
                        input_signal[i] = 1;
                        i++;
                    }
                    break;
                }
                case ('.'):
                {
                    for (k=0; k<MESSAGE_LDOT; k++) 
                    {
                        input_signal[i] = 1;
                        i++;
                    }
                    break;
                }
                case ('|'):
                {
                    for (k=0; k<MESSAGE_LINTERLETTER; k++) 
                    {
                        input_signal[i] = 0;
                        i++;
                    }
                    break;
                }
                case ('/'):
                {
                    for (k=0; k<MESSAGE_LSPACE; k++) 
                    {
                        input_signal[i] = 0;
                        i++;
                    }
                    break;
                }
            }
            
            for (k=0; k<MESSAGE_LINTERVAL; k++) 
            {
                input_signal[i] = 0;
                i++;
            }
            
            j++;
        }
        
        input_signal[i] = 0;
        i++;
        tailcounter++;
    }
    
    for (i=0; i<NSTEPS; i++) printf("%i ", input_signal[i]);
    
    printf("\n\n j = %i, string length = %i\n", j, (int)strlen(text));
    printf("Tail frames: %i\n", tailcounter); 
}



