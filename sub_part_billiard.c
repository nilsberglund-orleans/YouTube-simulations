#include "colormaps.c"

#define DUMMY_ABSORBING 100000.0  /* dummy value of config[0] for absorbing circles */
#define BOUNDARY_SHIFT 10000.0    /* shift of boundary parametrisation for circles in domain */
#define DUMMY_SIDE_ABS -100000      /* dummy value of returned side for absorbing circles */
#define MAZE_VERTICAL_EXTENSION 1.0e10  /* extension of vertical lines in some mazes */

long int global_time = 0;    /* counter to keep track of global time of simulation */
int nparticles=NPART; 

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
// 	if (alph < 0.0) alph += DPI;
	return(alph);
 }
 
 int polynome(double a, double b, double c, double r[2])
 {
	double delta, rdelta;
	int im = 1;
	
	delta = b*b - 4*a*c;
	if (delta<0.0)
	{
	/*	printf("ca deconne!");*/
		rdelta = 0.0;
		im = 0;
	}
	else rdelta = sqrt(delta);

	r[0] = (-b + rdelta)/(2.0*a);
	r[1] = (-b - rdelta)/(2.0*a);

	return(im);
 }

 double ipow(double x, int n)
 {
    double y;
    int i;
    
    y = x;
    for (i=1; i<n; i++) y *= x;
    
    return(y);
 }
 
/*********************/
/* Graphics routines */
/*********************/

int writetiff(char *filename, char *description, int x, int y, int width, int height, int compression)
{
  TIFF *file;
  GLubyte *image, *p;
  int i;
  static int mem_counter = 0;

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
  
  /* to avoid RAM overflow */
  if (SAVE_MEMORY) free(image);
//   mem_counter++;
//   if (mem_counter >= 12)
//   {
//         free(image);
//         mem_counter = 0;
//   }
  return 0;
}

  
void init()		/* initialisation of window */
{
    glLineWidth(3);
 
    glClearColor(0.0, 0.0, 0.0, 1.0);
    glClear(GL_COLOR_BUFFER_BIT);

    glOrtho(XMIN, XMAX, YMIN, YMAX , -1.0, 1.0);
}



void rgb_color_scheme(int i, double rgb[3]) /* color scheme */
{
    double hue, y, r;
  
    hue = (double)(COLORSHIFT + i*360/NCOLORS);
    r = 0.9;
  
    while (hue < 0.0) hue += 360.0;
    while (hue >= 360.0) hue -= 360.0;
  
    hsl_to_rgb(hue, r, 0.5, rgb);
    /* saturation = r, luminosity = 0.5 */ 
}


void rgb_color_scheme_minmax(int i, double rgb[3]) 
/* color scheme with specified color interval */
{
    double hue, y, r;
    int delta_hue;
    
    delta_hue = COLOR_HUEMAX - COLOR_HUEMIN;
  
    hue = (double)(COLOR_HUEMIN + i*delta_hue/NCOLORS);
    r = 0.9;
  
    while (hue < COLOR_HUEMIN) hue += delta_hue;
    while (hue >= COLOR_HUEMAX) hue -= delta_hue;
  
    hsl_to_rgb(hue, r, 0.5, rgb);
    /* saturation = r, luminosity = 0.5 */ 
}

void rgb_color_scheme_density(int part_number, double rgb[3], double minprop) 
/* color scheme with specified color interval */
{
    int k;
    double hue, y, r, logprop, logmin;
    
    if (part_number < 0) for (k=0; k<3; k++) rgb[k] = 0.25;  /* visited cells */
    else 
    {
        if (part_number<=1) hue = COLOR_HUEMAX;
        else
        {
            logmin = log(minprop*(double)NPART);
            if (logmin > 0.0) logprop = log((double)part_number)/logmin + 0.01;
            else (logprop = 0.0);
        
            if (logprop > 1.0) logprop = 1.0;
            if (logprop < 0.0) logprop = 0.0;
        
            hue = COLOR_HUEMAX + (COLOR_HUEMIN - COLOR_HUEMAX)*logprop;
        }
        r = 0.9;
        hsl_to_rgb(hue, r, 0.5, rgb);
        /* saturation = r, luminosity = 0.5 */
    }
}


void rgb_color_scheme_lum(int i, double lum, double rgb[3]) /* color scheme */
{
    double hue, y, r;
  
    hue = (double)(COLORSHIFT + i*360/NCOLORS);
    r = 0.9;
  
    while (hue < 0.0) hue += 360.0;
    while (hue >= 360.0) hue -= 360.0;
  
    hsl_to_rgb(hue, r, lum, rgb);
    /* saturation = r */ 
}

void rgb_color_scheme_minmax_lum(int i, double lum, double rgb[3]) 
/* color scheme with specified color interval */
{
    double hue, y, r;
    int delta_hue;
    
    delta_hue = COLOR_HUEMAX - COLOR_HUEMIN;
  
    hue = (double)(COLOR_HUEMIN + i*delta_hue/NCOLORS);
    r = 0.9;
  
    while (hue < COLOR_HUEMIN) hue += delta_hue;
    while (hue >= COLOR_HUEMAX) hue -= delta_hue;
  
    hsl_to_rgb(hue, r, lum, rgb);
    /* saturation = r */ 
}

void blank()
{
    double rgb[3];
    
    if (COLOR_OUTSIDE)
    {
        hsl_to_rgb(OUTER_COLOR, 0.9, 0.15, rgb); 
        glClearColor(rgb[0], rgb[1], rgb[2], 1.0);
    }
    else if (BLACK) glClearColor(0.0, 0.0, 0.0, 1.0);
    else glClearColor(1.0, 1.0, 1.0, 1.0);
    glClear(GL_COLOR_BUFFER_BIT);
}


void save_frame()
{
  static int counter = 0;
  char *name="part.", n2[100];
  char format[6]=".%05i";
  
    counter++; 
//     printf (" p2 counter = %d \n",counter);
    strcpy(n2, name);
    sprintf(strstr(n2,"."), format, counter);
    strcat(n2, ".tif");
    printf(" saving frame %s \n",n2);
    writetiff(n2, "Billiard in an ellipse", 0, 0,
         WINWIDTH, WINHEIGHT, COMPRESSION_LZW);

}

void save_frame_counter(int counter)
{
  char *name="part.", n2[100];
  char format[6]=".%05i";
  
    strcpy(n2, name);
    sprintf(strstr(n2,"."), format, counter);
    strcat(n2, ".tif");
    printf(" saving frame %s \n",n2);
    writetiff(n2, "Billiard in an ellipse", 0, 0,
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

void erase_area(double x, double y, double dx, double dy, double rgb[3])
{
    glColor3f(rgb[0], rgb[1], rgb[2]);
    glBegin(GL_QUADS);
    glVertex2d(x - dx, y - dy);
    glVertex2d(x + dx, y - dy);
    glVertex2d(x + dx, y + dy);
    glVertex2d(x - dx, y + dy);
    glEnd();
}

void erase_rectangle(double x1, double y1, double x2, double y2, double rgb[3])
{
    glColor3f(rgb[0], rgb[1], rgb[2]);
    glBegin(GL_QUADS);
    glVertex2d(x1, y1);
    glVertex2d(x2, y1);
    glVertex2d(x2, y2);
    glVertex2d(x1, y2);
    glEnd();
}

void erase_rectangle_outside(double h, double s, double l)
{
    double rgb[3], dx;
    int k;
    
    dx = 0.5*(XMAX - LAMBDA);
    hsl_to_rgb(h, s, l, rgb);
    erase_rectangle(XMIN, YMIN, XMAX, -1.0, rgb);
    erase_rectangle(XMIN, 1.0, XMAX, YMAX, rgb);
    erase_rectangle(XMIN, YMIN, -LAMBDA, YMAX, rgb);
    erase_rectangle(LAMBDA, YMIN, XMAX, YMAX, rgb);
//     erase_area(0.0, 1.1, 2.0, 0.1, rgb);
//     erase_area(0.0, -1.1, 2.0, 0.1, rgb);
//     erase_area(LAMBDA + dx, 0.0, dx, 2.0, rgb);
//     erase_area(-LAMBDA - dx, 0.0, dx, 2.0, rgb);
    
    glColor3f(1.0, 1.0, 1.0);
    glLineWidth(BILLIARD_WIDTH);
    glBegin(GL_LINE_LOOP);    
    glVertex2d(LAMBDA, -1.0);
    glVertex2d(LAMBDA, 1.0);
    glVertex2d(-LAMBDA, 1.0);
    glVertex2d(-LAMBDA, -1.0);
    glEnd();
}

void draw_line(double x1, double y1, double x2, double y2)
{
    glBegin(GL_LINE_STRIP);
    glVertex2d(x1, y1);
    glVertex2d(x2, y2);
    glEnd();    
}

void draw_circle(double x, double y, double r, int nseg)
{
    int i;
    double phi, dphi, x1, y1;
    
    dphi = DPI/(double)nseg;
    
    glEnable(GL_LINE_SMOOTH);
    glBegin(GL_LINE_LOOP);
    for (i=0; i<=nseg; i++)
    {
        phi = (double)i*dphi;
        x1 = x + r*cos(phi);
        y1 = y + r*sin(phi);
        glVertex2d(x1, y1);
    }
    glEnd ();
}

void draw_colored_circle(double x, double y, double r, int nseg, double rgb[3])
{
    int i, ij[2];
    double pos[2], alpha, dalpha;
    
    dalpha = DPI/(double)nseg;
    
//     glLineWidth(2);
    
    glColor3f(rgb[0], rgb[1], rgb[2]);
    glBegin(GL_TRIANGLE_FAN);
    glVertex2d(x,y);
    for (i=0; i<=nseg; i++)
    {
        alpha = (double)i*dalpha;
        glVertex2d(x + r*cos(alpha), y + r*sin(alpha));
    }
    
    glEnd();
}

void draw_colored_octahedron(double x, double y, double r, double rgb[3])
{
    int i, ij[2];
    double pos[2], alpha, dalpha;
    
    dalpha = DPI*0.125;
        
    glColor3f(rgb[0], rgb[1], rgb[2]);
    glBegin(GL_TRIANGLE_FAN);
    glVertex2d(x,y);
    for (i=0; i<=8; i++)
    {
        alpha = ((double)i+0.5)*dalpha;
        glVertex2d(x + r*cos(alpha), y + r*sin(alpha));
    }
    
    glEnd();
}

void draw_colored_rectangle_rgb(double x1, double y1, double x2, double y2, double rgb[3])
{
    glColor3f(rgb[0], rgb[1], rgb[2]);
    
    glBegin(GL_QUADS);
    glVertex2d(x1, y1);
    glVertex2d(x2, y1);
    glVertex2d(x2, y2);
    glVertex2d(x1, y2);
    glEnd();    

    glColor3f(1.0, 1.0, 1.0);
    glBegin(GL_LINE_LOOP);
    glVertex2d(x1, y1);
    glVertex2d(x2, y1);
    glVertex2d(x2, y2);
    glVertex2d(x1, y2);
    glEnd();    
}

void draw_colored_rectangle(double x1, double y1, double x2, double y2, double hue)
{
    double rgb[3];
    
    hsl_to_rgb_palette(hue, 0.9, 0.5, rgb, COLOR_PALETTE);
    draw_colored_rectangle_rgb(x1, y1, x2, y2, rgb);
}



void draw_filled_sector(double x, double y, double rmin, double rmax, double phi1, double dphi, int nseg, double rgb[3])
{
    int i;
    double alpha, dalpha;
    
    dalpha = dphi/(double)nseg;
    
    glColor3f(rgb[0], rgb[1], rgb[2]);
    for (i=0; i<=nseg; i++)
    {
        alpha = phi1 + (double)i*dalpha;
        glBegin(GL_TRIANGLE_FAN);
        glVertex2d(x + rmin*cos(alpha), y + rmin*sin(alpha));
        glVertex2d(x + rmax*cos(alpha), y + rmax*sin(alpha));
        glVertex2d(x + rmax*cos(alpha+dalpha), y + rmax*sin(alpha+dalpha));
        glVertex2d(x + rmin*cos(alpha+dalpha), y + rmin*sin(alpha+dalpha));
        glEnd();
    }
    
}


void draw_initial_condition_circle(double x, double y, double r, int color)
/* draws a colored circle to mark initial condition */
{
    double rgb[3];
    
    rgb_color_scheme(color, rgb);
    draw_colored_circle(x, y, r, NSEG, rgb);
    rgb[0] = 0.0; rgb[1] = 0.0; rgb[2] = 0.0;
    glColor3f(rgb[0], rgb[1], rgb[2]);
    glLineWidth(4);
    draw_circle(x, y, r, NSEG);
}

int in_circle(double x, double y, double r)
/* test whether (x,y) is in circle of radius r */
{
    return(x*x + y*y < r*r);
}

int out_circle(double x, double y, double r)
/* test whether (x,y) is in circle of radius r */
{
    return(x*x + y*y > r*r);
}

int in_polygon(double x, double y, double r, int npoly, double apoly)
/* test whether (x,y) is in regular polygon of npoly sides inscribed in circle of radious r, turned by apoly Pi/2 */
{
    int condition = 1, k;
    double omega, cw, angle; 
    
    omega = DPI/((double)npoly);
    cw = cos(omega*0.5);
    for (k=0; k<npoly; k++)  
    {
        angle = apoly*PID + (k+0.5)*omega;
        condition = condition*(x*cos(angle) + y*sin(angle) < r*cw);
    }
    return(condition);
}

 void compute_flower_parameters(double *omega, double *co, double *so, double *axis1, double *axis2, double *phimax)
 /* compute parameters needed for the flower billiard in terms of LAMBDA and NPOLY */
//  double *omega, *co, *so, *axis1, *axis2, *phimax;
 {
    double omega2, co2, so2, r, a, gamma, axissquare1;
    
    /* various angles */
    *omega = DPI/((double)NPOLY);
    omega2 = PI/((double)NPOLY);
    co2 = cos(omega2);
    so2 = sin(omega2);
    *co = cos(*omega);
    *so = sin(*omega);
//     *co = co2*co2 - so2*so2;
//     *so = 2.0*co2*so2;
    
    /* distance of edge of ellipse to the origin */
    r = LAMBDA*co2/(*co);
    
    a = (r*co2 - *co)*(r*co2 - *co);
    gamma = 0.5*r*r - r*co2*(*co) + 0.5*cos(2.0*(*omega));
    axissquare1 = gamma + sqrt(gamma*gamma + a*(*so)*(*so));
    
    /* semi-minor axis */
    *axis1 = sqrt(axissquare1);
    
    /* semi-major axis */
    *axis2 = sqrt(axissquare1 + (*so)*(*so));
    
    /* max angle in ellipse parametrization */
    *phimax = asin(r*so2/(*axis2));
 }
 
 
void paint_billiard_interior()      /* paints billiard interior, for use before draw_conf */
{
    double x0, x, y, phi, r = 0.01, alpha, dphi, omega, beta2, x2, s, x1, y1, angle, co, so, axis1, axis2, phimax;
    int i, j, k, c;
    
    glLineWidth(4);
    
    glEnable(GL_LINE_SMOOTH);
    
    switch (B_DOMAIN) {
        case (D_POLYGON):
        {
            omega = DPI/((double)NPOLY);
            
            if (PAINT_INT)
            {
                if (BLACK) glColor3f(1.0, 1.0, 1.0);
                else glColor3f(0.0, 0.0, 0.0);
                
                glBegin(GL_TRIANGLE_FAN);
                glVertex2d(0.0, 0.0);
                for (i=0; i<=NPOLY; i++)
                {
                    x = cos(i*omega + APOLY*PID);
                    y = sin(i*omega + APOLY*PID);
                    glVertex2d(x, y);
                    x = cos((i+1)*omega + APOLY*PID);
                    y = sin((i+1)*omega + APOLY*PID);
                    glVertex2d(x, y);
               }
                glEnd();
            }
            break;
        }
        case (D_REULEAUX):
        {
            omega = DPI/((double)NPOLY);
            beta2 = asin(sin(omega*0.5)/LAMBDA);
            if (LAMBDA > 0.0) x2 = cos(omega*0.5) + sqrt(LAMBDA*LAMBDA - sin(omega*0.5)*sin(omega*0.5));
            else x2 = cos(omega*0.5) - sqrt(LAMBDA*LAMBDA - sin(omega*0.5)*sin(omega*0.5));
            
            if (PAINT_INT)
            {
                if (BLACK) glColor3f(1.0, 1.0, 1.0);
                else glColor3f(0.0, 0.0, 0.0);
                
                glBegin(GL_TRIANGLE_FAN);
                glVertex2d(0.0, 0.0);
                for (i=0; i<=NPOLY; i++)
                {
                    for (j=0; j<NSEG; j++)
                {
                    s = 2.0*(((double)j/(double)NSEG)-0.5)*beta2;
                    x1 = x2 - LAMBDA*cos(s);
                    y1 = LAMBDA*sin(s);
                    angle = i*omega + APOLY*PID;
                    x = cos(angle)*x1 - sin(angle)*y1;
                    y = sin(angle)*x1 + cos(angle)*y1;
                    glVertex2d(x, y);
                }
               }
                glEnd();
            }
            break;
        }
        case (D_FLOWER): 
        {
            compute_flower_parameters(&omega, &co, &so, &axis1, &axis2, &phimax);
            if (PAINT_INT)
            {
                if (BLACK) glColor3f(1.0, 1.0, 1.0);
                else glColor3f(0.0, 0.0, 0.0);
                
                glBegin(GL_TRIANGLE_FAN);
                glVertex2d(0.0, 0.0);
                for (i=0; i<=NPOLY; i++)
                {
                    for (j=0; j<NSEG; j++)
                    {
                        s = 2.0*(((double)j/(double)NSEG)-0.5)*phimax;
                        x1 = co + axis1*cos(s);
                        y1 = axis2*sin(s);
                        angle = i*omega + APOLY*PID;
                        x = cos(angle)*x1 - sin(angle)*y1;
                        y = sin(angle)*x1 + cos(angle)*y1;
                        glVertex2d(SCALING_FACTOR*x, SCALING_FACTOR*y);
                    }
               }
                glEnd();
            }
            break;
        }
        default: 
        {
            
        }
    }
}

void init_billiard_color()
/* initialise the color in which the billiard is drawn */
{
    if (PAINT_INT) glColor3f(0.5, 0.5, 0.5);
    else
    {
        if (BLACK) glColor3f(1.0, 1.0, 1.0);
        else glColor3f(0.0, 0.0, 0.0);
    }
}

void draw_billiard()      /* draws the billiard boundary */
{
    double x0, x, y, phi, r = 0.01, alpha, dphi, omega, x1, y1, x2, beta2, angle, s, x2plus, x2minus;
    double omega2, co, so, axis1, axis2, phimax, rgb[3], rgb1[3], a, b, ymax, dy, width, cc;
    int i, j, k, c, color;
    
    init_billiard_color();
    glLineWidth(BILLIARD_WIDTH);
    
    glEnable(GL_LINE_SMOOTH);
    
    switch (B_DOMAIN) {
        case (D_RECTANGLE):
        {
            glBegin(GL_LINE_LOOP);    
            glVertex2d(LAMBDA, -1.0);
            glVertex2d(LAMBDA, 1.0);
            glVertex2d(-LAMBDA, 1.0);
            glVertex2d(-LAMBDA, -1.0);
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
                glVertex2d(x, y);
            }
            glEnd ();
    
            /* draw foci */
            if (FOCI)
            {
                glColor3f(0.3, 0.3, 0.3);
                x0 = sqrt(LAMBDA*LAMBDA-1.0);
    
                glLineWidth(2);
                
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
                glVertex2d(x, y);
            }
            for (i=0; i<=NSEG; i++)
            {
                phi = PID + (double)i*PI/(double)NSEG;
                x = -0.5*LAMBDA + cos(phi);
                y = sin(phi);
                glVertex2d(x, y);
            }
            glEnd();
            
            if (DRAW_CONSTRUCTION_LINES)
            {
                glColor3f(0.5, 0.5, 0.5);
                glBegin(GL_LINE_STRIP);
                glVertex2d(-LAMBDA, -1.0);
                glVertex2d(-LAMBDA, 1.0);
                glEnd();
                glBegin(GL_LINE_STRIP);
                glVertex2d(LAMBDA, -1.0);
                glVertex2d(LAMBDA, 1.0);
                glEnd();
            }
            
            break; 
        }
        case (D_SINAI):
        {
            glColor3f(0.5, 0.5, 0.5);
            glBegin(GL_TRIANGLE_FAN);
            glVertex2d(0.0, 0.0);
            for (i=0; i<=NSEG; i++)
            {
                phi = (double)i*DPI/(double)NSEG;
                x = LAMBDA*cos(phi);
                y = LAMBDA*sin(phi);
                glVertex2d(x, y);
            }
            glEnd();
            
            if (BLACK) glColor3f(1.0, 1.0, 1.0);
            else glColor3f(0.0, 0.0, 0.0);
            
            glBegin(GL_LINE_LOOP);
            for (i=0; i<=NSEG; i++)
            {
                phi = (double)i*DPI/(double)NSEG;
                x = LAMBDA*cos(phi);
                y = LAMBDA*sin(phi);
                glVertex2d(x, y);
            }
            glEnd();
            

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
                glVertex2d(x, y);
            }
            for (i=0; i<=NSEG; i++)
            {
                phi = alpha - PID + (double)i*dphi;
                x = -LAMBDA + r*cos(phi);
                y = LAMBDA + r*sin(phi);
                glVertex2d(x, y);
            }
            for (i=0; i<=NSEG; i++)
            {
                phi = alpha + PI + (double)i*dphi;
                x = LAMBDA + r*cos(phi);
                y = LAMBDA + r*sin(phi);
                glVertex2d(x, y);
            }
            for (i=0; i<=NSEG; i++)
            {
                phi = alpha + PID + (double)i*dphi;
                x = LAMBDA + r*cos(phi);
                y = -LAMBDA + r*sin(phi);
                glVertex2d(x, y);
            }
            glEnd();
            break; 
        }
        case (D_TRIANGLE):
        {
            glBegin(GL_LINE_LOOP);    
            glVertex2d(-LAMBDA, -1.0);
            glVertex2d(LAMBDA, -1.0);
            glVertex2d(-LAMBDA, 1.0);
            glEnd();
            break; 
        }
        case (D_ANNULUS):
        {
            /* color inner circle */
            glColor3f(0.5, 0.5, 0.5);
            glBegin(GL_TRIANGLE_FAN);
            glVertex2d(MU, 0.0);
            for (i=0; i<=NSEG; i++)
            {
                phi = (double)i*DPI/(double)NSEG;
                x = LAMBDA*cos(phi) + MU;
                y = LAMBDA*sin(phi);
                glVertex2d(x, y);
            }
            glEnd();
            
            /* color outer domain */
            glColor3f(0.2, 0.2, 0.2);
            glBegin(GL_TRIANGLE_FAN);
            glVertex2d(XMAX, YMAX);
            for (i=0; i<=NSEG; i++)
            {
                phi = (double)i*PID/(double)NSEG;
                x = cos(phi);
                y = sin(phi);
                glVertex2d(x, y);
            }
            glEnd();

            glBegin(GL_TRIANGLE_FAN);
            glVertex2d(XMIN, YMAX);
            for (i=0; i<=NSEG; i++)
            {
                phi = (double)i*PID/(double)NSEG + PID;
                x = cos(phi);
                y = sin(phi);
                glVertex2d(x, y);
            }
            glEnd();

            glBegin(GL_TRIANGLE_FAN);
            glVertex2d(XMIN, YMIN);
            for (i=0; i<=NSEG; i++)
            {
                phi = (double)i*PID/(double)NSEG + PI;
                x = cos(phi);
                y = sin(phi);
                glVertex2d(x, y);
            }
            glEnd();

            glBegin(GL_TRIANGLE_FAN);
            glVertex2d(XMAX, YMIN);
            for (i=0; i<=NSEG; i++)
            {
                phi = (double)i*PID/(double)NSEG + 3.0*PID;
                x = cos(phi);
                y = sin(phi);
                glVertex2d(x, y);
            }
            glEnd();
            
            glBegin(GL_TRIANGLES);
            glVertex2d(XMAX, YMAX);
            glVertex2d(1.0, 0.0);
            glVertex2d(XMAX, YMIN);
            glVertex2d(XMAX, YMIN);
            glVertex2d(0.0, -1.0);
            glVertex2d(XMIN, YMIN);
            glVertex2d(XMIN, YMIN);
            glVertex2d(-1.0, 0.0);
            glVertex2d(XMIN, YMAX);
            glVertex2d(XMIN, YMAX);
            glVertex2d(0.0, 1.0);
            glVertex2d(XMAX, YMAX);
            glEnd();

            /* draw circles */
            if (BLACK) glColor3f(1.0, 1.0, 1.0);
            else glColor3f(0.0, 0.0, 0.0);
            
            draw_circle(MU, 0.0, LAMBDA, NSEG);
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
                glVertex2d(SCALING_FACTOR*x, SCALING_FACTOR*y);
            }
            glEnd ();
            break;
        }
        case (D_REULEAUX):
        {
            omega = DPI/((double)NPOLY);
            beta2 = asin(sin(omega*0.5)/LAMBDA);
            if (LAMBDA > 0.0) x2 = cos(omega*0.5) + sqrt(LAMBDA*LAMBDA - sin(omega*0.5)*sin(omega*0.5));
            else x2 = cos(omega*0.5) - sqrt(LAMBDA*LAMBDA - sin(omega*0.5)*sin(omega*0.5));
            glBegin(GL_LINE_STRIP);
            for (i=0; i<=NPOLY; i++)
            {
                for (j=0; j<NSEG; j++)
                {
                    s = 2.0*(((double)j/(double)NSEG)-0.5)*beta2;
                    x1 = x2 - LAMBDA*cos(s);
                    y1 = LAMBDA*sin(s);
                    angle = i*omega + APOLY*PID;
                    x = cos(angle)*x1 - sin(angle)*y1;
                    y = sin(angle)*x1 + cos(angle)*y1;
                    glVertex2d(SCALING_FACTOR*x, SCALING_FACTOR*y);
                }
            }
            glEnd ();
            break;
        }
        case (D_FLOWER):
        {
            compute_flower_parameters(&omega, &co, &so, &axis1, &axis2, &phimax);

            /* draw inner polygon and radial lines */
            if (DRAW_CONSTRUCTION_LINES)
            {
                glColor3f(0.5, 0.5, 0.5);
                glBegin(GL_LINE_LOOP);
                for (i=0; i<=NPOLY; i++)
                {
                    x = cos(i*omega + APOLY*PID);
                    y = sin(i*omega + APOLY*PID);
                    glVertex2d(SCALING_FACTOR*x, SCALING_FACTOR*y);
                }
                glEnd ();

                r = LAMBDA*cos(0.5*omega)/co;
                for (i=0; i<=NPOLY; i++)
                {
                    glBegin(GL_LINE_STRIP);
                    glVertex2d(0.0, 0.0);
                    x = r*cos(((double)i + 0.5)*omega + APOLY*PID);
                    y = r*sin(((double)i + 0.5)*omega + APOLY*PID);
                    glVertex2d(SCALING_FACTOR*x, SCALING_FACTOR*y);
                    glEnd ();
                }
            }
            
            /* draw billiard boundary */
            if (!PAINT_INT) 
            {
                if (BLACK) glColor3f(1.0, 1.0, 1.0);
                else glColor3f(0.0, 0.0, 0.0);
            }
            glBegin(GL_LINE_STRIP);
            for (i=0; i<=NPOLY; i++)
//             for (i=0; i<=1; i++)
            {
                for (j=0; j<NSEG; j++)
                {
//                     s = 2.0*(((double)j/(double)NSEG)-0.5)*PI;
                    s = 2.0*(((double)j/(double)NSEG)-0.5)*phimax;
                    x1 = co + axis1*cos(s);
                    y1 = axis2*sin(s);
                    angle = i*omega + APOLY*PID;
                    x = cos(angle)*x1 - sin(angle)*y1;
                    y = sin(angle)*x1 + cos(angle)*y1;
                    glVertex2d(SCALING_FACTOR*x, SCALING_FACTOR*y);
                }
            }
            glEnd ();
            
            
            break;
        }
        case (D_ALT_REU):
        {
            omega = DPI/((double)NPOLY);
            beta2 = asin(sin(omega*0.5)/LAMBDA);
            x2plus = cos(omega*0.5) + sqrt(LAMBDA*LAMBDA - sin(omega*0.5)*sin(omega*0.5));
            x2minus = cos(omega*0.5) - sqrt(LAMBDA*LAMBDA - sin(omega*0.5)*sin(omega*0.5));
            glBegin(GL_LINE_STRIP);
            for (i=0; i<=NPOLY; i++)
            {
                for (j=0; j<NSEG; j++)
                {
                    s = 2.0*(((double)j/(double)NSEG)-0.5)*beta2;
                    if (i%2==0) x1 = x2plus - LAMBDA*cos(s);
                    else x1 = x2minus + LAMBDA*cos(s);
                    y1 = LAMBDA*sin(s);
                    angle = i*omega + APOLY*PID;
                    x = cos(angle)*x1 - sin(angle)*y1;
                    y = sin(angle)*x1 + cos(angle)*y1;
                    glVertex2d(x, y);
                }
            }
            glEnd ();
            break;
        }
        case (D_ANGLE):
        {
            if (DRAW_CONSTRUCTION_LINES)
            {
                glColor3f(0.5, 0.5, 0.5);
                glBegin(GL_LINE_STRIP);
                phi = LAMBDA*PI;
                r = 2.0*YMAX;
                for (i=-(int)(PI/phi)-1; i <= (int)(PI/phi)+2; i++)
                {
//                     printf("%i \n", i);
                    glVertex2d(0.0, 0.0);
                    glVertex2d(-r*cos((double)i*phi), r*sin((double)i*phi));
                }
                glEnd();
            }
            
            init_billiard_color();
            
            glBegin(GL_LINE_STRIP);
            glVertex2d(XMIN, 0.0);
            glVertex2d(0.0, 0.0);
            glVertex2d(-YMAX/tan(LAMBDA*PI), YMAX);
            glEnd ();
            break;
        }
        case (D_LSHAPE):
        {
            glBegin(GL_LINE_LOOP);    
            glVertex2d(-1.0, -1.0);
            glVertex2d(1.0, -1.0);
            glVertex2d(1.0, 0.0);
            glVertex2d(0.0, 0.0);
            glVertex2d(0.0, 1.0);
            glVertex2d(-1.0, 1.0);
            glEnd();
            break; 
        }
        case (D_GENUSN):    /* like for polygon */
        {
            omega = DPI/((double)NPOLY);
            glBegin(GL_LINE_LOOP);
            for (i=0; i<=NPOLY; i++)
            {
                x = cos(i*omega + APOLY*PID);
                y = sin(i*omega + APOLY*PID);
                glVertex2d(SCALING_FACTOR*x, SCALING_FACTOR*y);
            }
            glEnd ();
            break;
        }
        case (D_PARABOLAS):
        {
            omega = PI/((double)NPOLY);
            a = 0.25/MU;
            b = 1.0/tan(omega);
            cc = LAMBDA + MU;
            ymax = (-b + sqrt(b*b + 4.0*a*cc))/(2.0*a);
            dy = 2.0*ymax/(double)NSEG; 
            
            if (PAINT_EXT)  /* paint billiard exterior in another color */
            {
                glColor3f(0.1, 0.1, 0.1);
                for (k=0; k<NPOLY; k++)  
                {
                    alpha = APOLY*PID + (2.0*(double)k)*omega;
                    glBegin(GL_TRIANGLE_FAN);
                    x1 = 10.0*(MU + LAMBDA);
                    x = x1*cos(alpha);
                    y = x1*sin(alpha);
                    glVertex2d(x, y);
                    for (i = 0; i < NSEG+1; i++) 
                    {
                        y1 = -ymax + dy*(double)i;
                        x1 = MU + LAMBDA - 0.25*y1*y1/MU;
                        x = x1*cos(alpha) - y1*sin(alpha);
                        y = x1*sin(alpha) + y1*cos(alpha);
                        glVertex2d(x, y);
                    }
                    glEnd();
                    
                    glBegin(GL_TRIANGLE_FAN);
                    
                    x1 = 10.0*(MU + LAMBDA);
                    x = x1*cos(alpha);
                    y = x1*sin(alpha);
                    glVertex2d(x, y);
                    
                    y1 = ymax;
                    x1 = MU + LAMBDA - 0.25*y1*y1/MU;
                    x = x1*cos(alpha) - y1*sin(alpha);
                    y = x1*sin(alpha) + y1*cos(alpha);
                    glVertex2d(x, y);
                    
                    x1 = 10.0*(MU + LAMBDA);
                    x = x1*cos(alpha + 2.0*omega);
                    y = x1*sin(alpha + 2.0*omega);
                    glVertex2d(x, y);
                    
                    glEnd();
                }
                
                
//                 glVertex2d(cc, 0.0);
//                 for (i=0; i<=NSEG; i++)
//                 {
//                     phi = (double)i*dphi;
//                     x = cc - 0.5*MU*sin(phi);
//                     y = MU*cos(phi);
//                     glVertex2d(x, y);
//                 }
                
            }

            init_billiard_color();
            glBegin(GL_LINE_LOOP);
//             glColor3f(1.0, 1.0, 1.0);
            for (k=0; k<NPOLY; k++)  
            {
//                 alpha = APOLY*PID + (2.0*(double)k+1.0)*omega;
                alpha = APOLY*PID + (2.0*(double)k)*omega;
                for (i = 0; i < NSEG+1; i++) 
                {
                    y1 = -ymax + dy*(double)i;
                    x1 = MU + LAMBDA - 0.25*y1*y1/MU;
                    x = x1*cos(alpha) - y1*sin(alpha);
                    y = x1*sin(alpha) + y1*cos(alpha);
                    glVertex2d(x, y);
                }
            }
            glEnd ();
            
            if (FOCI)
            {
                glColor3f(0.3, 0.3, 0.3);
                for (k=0; k<NPOLY; k++) 
                {
                    alpha = APOLY*PID + (2.0*(double)k)*omega;
                    draw_circle(LAMBDA*cos(alpha), LAMBDA*sin(alpha), r, NSEG);
                }
            }
            break;
        }
        case (D_PENROSE):
        {
            cc = sqrt(LAMBDA*LAMBDA - (1.0 - MU)*(1.0 - MU));
            width = 0.1*MU;
            x1 = vabs(x);
            y1 = vabs(y);
            dphi = PI/(double)NSEG;
            
            if (PAINT_EXT)  /* paint billiard exterior in another color */
            {
                rgb[0] = 0.1;   rgb[1] = 0.1;   rgb[2] = 0.1;
                erase_rectangle(LAMBDA, YMIN, XMAX, YMAX, rgb);
                erase_rectangle(-LAMBDA, YMIN, XMIN, YMAX, rgb);
                erase_rectangle(XMIN, 1.0, XMAX, YMAX, rgb);
                erase_rectangle(XMIN, -1.0, XMAX, YMIN, rgb);
                erase_rectangle(cc, -width, LAMBDA, width, rgb);
                erase_rectangle(-cc, -width, -LAMBDA, width, rgb);
                
                glColor3f(0.1, 0.1, 0.1);
                glBegin(GL_TRIANGLE_FAN);
                glVertex2d(cc, 0.0);
                for (i=0; i<=NSEG; i++)
                {
                    phi = (double)i*dphi;
                    x = cc - 0.5*PENROSE_RATIO*MU*sin(phi);
                    y = MU*cos(phi);
                    glVertex2d(x, y);
                }
                glEnd();
                
                glBegin(GL_TRIANGLE_FAN);
                glVertex2d(-cc, 0.0);
                for (i=0; i<=NSEG; i++)
                {
                    phi = (double)i*dphi;
                    x = -cc + 0.5*PENROSE_RATIO*MU*sin(phi);
                    y = MU*cos(phi);
                    glVertex2d(x, y);
                }
                glEnd();
                
                glBegin(GL_TRIANGLE_FAN);
                glVertex2d(XMAX, YMAX);
                for (i=0; i<=NSEG/2; i++)
                {
                    phi = (double)i*dphi;
                    x = LAMBDA*cos(phi);
                    y = MU + (1.0-MU)*sin(phi);
                    glVertex2d(x, y);
                }
                glEnd();
                
                glBegin(GL_TRIANGLE_FAN);
                glVertex2d(XMIN, YMAX);
                for (i=0; i<=NSEG/2; i++)
                {
                    phi = (double)i*dphi;
                    x = -LAMBDA*cos(phi);
                    y = MU + (1.0-MU)*sin(phi);
                    glVertex2d(x, y);
                }
                glEnd();
                
                glBegin(GL_TRIANGLE_FAN);
                glVertex2d(XMAX, YMIN);
                for (i=0; i<=NSEG/2; i++)
                {
                    phi = (double)i*dphi;
                    x = LAMBDA*cos(phi);
                    y = -MU - (1.0-MU)*sin(phi);
                    glVertex2d(x, y);
                }
                glEnd();
                
                glBegin(GL_TRIANGLE_FAN);
                glVertex2d(XMIN, YMIN);
                for (i=0; i<=NSEG/2; i++)
                {
                    phi = (double)i*dphi;
                    x = -LAMBDA*cos(phi);
                    y = -MU - (1.0-MU)*sin(phi);
                    glVertex2d(x, y);
                }
                glEnd();
            }
            
            init_billiard_color();
            
            glBegin(GL_LINE_LOOP);
            /* upper half ellipse */
            for (i=0; i<=NSEG; i++)
            {
                phi = (double)i*dphi;
                x = LAMBDA*cos(phi);
                y = MU + (1.0-MU)*sin(phi);
                glVertex2d(x, y);
            }
            
            /* straight parts */
            glVertex2d(-LAMBDA, width);
            glVertex2d(-cc, width);
            glVertex2d(-cc, MU);
            
            /* left half ellipse */
            for (i=0; i<=NSEG; i++)
            {
                phi = (double)i*dphi;
                x = -cc + 0.5*MU*PENROSE_RATIO*sin(phi);
                y = MU*cos(phi);
                glVertex2d(x, y);
            }
            
            /* straight parts */
            glVertex2d(-cc, -width);
            glVertex2d(-LAMBDA, -width);
            glVertex2d(-LAMBDA, -MU);
            
            /* lower half ellipse */
            for (i=0; i<=NSEG; i++)
            {
                phi = (double)i*dphi;
                x = -LAMBDA*cos(phi);
                y = -MU - (1.0-MU)*sin(phi);
                glVertex2d(x, y);
            }
            
            /* straight parts */
            glVertex2d(LAMBDA, -width);
            glVertex2d(cc, -width);
            glVertex2d(cc, -MU);
            
            /* right half ellipse */
            for (i=0; i<=NSEG; i++)
            {
                phi = (double)i*dphi;
                x = cc - 0.5*MU*PENROSE_RATIO*sin(phi);
                y = -MU*cos(phi);
                glVertex2d(x, y);
            }
            
            /* straight parts */
            glVertex2d(cc, width);
            glVertex2d(LAMBDA, width);
            glVertex2d(LAMBDA, MU);
            
            glEnd ();

            break;
        }
        case (D_CIRCLES):
        {
            rgb[0] = 0.0; rgb[1] = 0.0; rgb[2] = 0.0;
            for (k=0; k<ncircles; k++) if (circles[k].active)
            {
                if (circles[k].color == 0) draw_colored_circle(circles[k].xc, circles[k].yc, circles[k].radius, NSEG, rgb);
                else
                {
                    if (circles[k].new >= 1)
                    {
                        rgb_color_scheme_lum(circles[k].color, 0.85, rgb1);
                        circles[k].new--;
                    }
                    else rgb_color_scheme(circles[k].color, rgb1);
                    draw_colored_circle(circles[k].xc, circles[k].yc, circles[k].radius, NSEG, rgb1);
                }
            }
            init_billiard_color();
            for (k=0; k<ncircles; k++) if (circles[k].active)
                 draw_circle(circles[k].xc, circles[k].yc, circles[k].radius, NSEG);
            break; 
        }
        case (D_CIRCLES_IN_RECT):
        {
            rgb[0] = 0.0; rgb[1] = 0.0; rgb[2] = 0.0;
            for (k=0; k<ncircles; k++) if (circles[k].active)
            {
//                 printf("k = %i, color = %i\n", k, circlecolor[k]);
                if (circles[k].color == 0) draw_colored_circle(circles[k].xc, circles[k].yc, circles[k].radius, NSEG, rgb);
                else
                {
                    if (circles[k].new >= 1)
                    {
                        rgb_color_scheme_lum(circles[k].color, 0.85, rgb1);
                        circles[k].new--;
                    }
                    else rgb_color_scheme(circles[k].color, rgb1);
                    draw_colored_circle(circles[k].xc, circles[k].yc, circles[k].radius, NSEG, rgb1);
                }
            }
            init_billiard_color();
            for (k=0; k<ncircles; k++) if (circles[k].active >= 1)
            {
                if (circles[k].active == 2) 
                {
//                     hsl_to_rgb(150.0, 0.9, 0.4, rgb);
//                     glColor3f(rgb[0], rgb[1], rgb[2]);
                    glColor3f(0.0, 1.0, 0.0);
                    rgb[0] = 0.0;   rgb[1] = 0.9;   rgb[2] = 0.0;
                    draw_colored_circle(circles[k].xc, circles[k].yc, circles[k].radius, NSEG, rgb);
                }
                else draw_circle(circles[k].xc, circles[k].yc, circles[k].radius, NSEG);                
                init_billiard_color();
            }
            
            /* draw shooter position for laser pattern */
            if (CIRCLE_PATTERN == C_LASER)
            {
                hsl_to_rgb(0.0, 0.9, 0.5, rgb);
                glColor3f(rgb[0], rgb[1], rgb[2]);
                
                draw_circle(x_shooter, y_shooter, circles[ncircles-1].radius, NSEG);
            }
            
            init_billiard_color();

            glBegin(GL_LINE_LOOP);    
            glVertex2d(LAMBDA, -1.0);
            glVertex2d(LAMBDA, 1.0);
            glVertex2d(-LAMBDA, 1.0);
            glVertex2d(-LAMBDA, -1.0);
            glEnd();
            break; 
        }
        case (D_CIRCLES_IN_GENUSN):    
        {
//             for (k=0; k<ncircles; k++) if (circles[k].active >= 1)
            for (k=0; k<ncircles; k++) 
                if ((circles[k].active >= 1)&&(in_polygon(circles[k].xc, circles[k].yc, 1.0, NPOLY, APOLY)))
                {
                    if (circles[k].active == 2) 
                    {
                        hsl_to_rgb(150.0, 0.9, 0.4, rgb);
                            glColor3f(rgb[0], rgb[1], rgb[2]);
                    }
                    draw_circle(circles[k].xc, circles[k].yc, circles[k].radius, NSEG);
                    init_billiard_color();
                }
            
            /* draw shooter position for laser pattern */
            if ((CIRCLE_PATTERN == C_LASER)||(CIRCLE_PATTERN == C_LASER_GENUSN))
            {
                hsl_to_rgb(0.0, 0.9, 0.5, rgb);
                glColor3f(rgb[0], rgb[1], rgb[2]);
                
                draw_circle(x_shooter, y_shooter, circles[ncircles-1].radius, NSEG);
            }
            
            /* draw polygon */
            init_billiard_color();

            omega = DPI/((double)NPOLY);
            glBegin(GL_LINE_LOOP);
            for (i=0; i<=NPOLY; i++)
            {
                x = cos(i*omega + APOLY*PID);
                y = sin(i*omega + APOLY*PID);
                glVertex2d(SCALING_FACTOR*x, SCALING_FACTOR*y);
            }
            glEnd ();
            break;
        }
        case (D_CIRCLES_IN_TORUS):
        {
            rgb[0] = 0.0; rgb[1] = 0.0; rgb[2] = 0.0;
            for (k=0; k<ncircles; k++) if (circles[k].active)
            {
//                 printf("k = %i, color = %i\n", k, circlecolor[k]);
                if (circles[k].color == 0) draw_colored_circle(circles[k].xc, circles[k].yc, circles[k].radius, NSEG, rgb);
                else
                {
                    if (circles[k].new >= 1)
                    {
                        rgb_color_scheme_lum(circles[k].color, 0.85, rgb1);
                        circles[k].new--;
                    }
                    else rgb_color_scheme(circles[k].color, rgb1);
                    draw_colored_circle(circles[k].xc, circles[k].yc, circles[k].radius, NSEG, rgb1);
                }
            }
            init_billiard_color();
            for (k=0; k<ncircles; k++) if (circles[k].active >= 1)
            {
                if (circles[k].active == 2) 
                {
                    hsl_to_rgb(150.0, 0.9, 0.4, rgb);
                    glColor3f(rgb[0], rgb[1], rgb[2]);
                }
                draw_circle(circles[k].xc, circles[k].yc, circles[k].radius, NSEG);                
                init_billiard_color();
            }
            
            /* draw shooter position for laser pattern */
            if (CIRCLE_PATTERN == C_LASER)
            {
                hsl_to_rgb(0.0, 0.9, 0.5, rgb);
                glColor3f(rgb[0], rgb[1], rgb[2]);
                
                draw_circle(x_shooter, y_shooter, circles[ncircles-1].radius, NSEG);
            }
            
            init_billiard_color();

            glBegin(GL_LINE_LOOP);    
            glVertex2d(LAMBDA, -1.0);
            glVertex2d(LAMBDA, 1.0);
            glVertex2d(-LAMBDA, 1.0);
            glVertex2d(-LAMBDA, -1.0);
            glEnd();
            break; 
        }
        case (D_POLYLINE):
        {
            for (k=0; k<nsides; k++)
            {
                glBegin(GL_LINE_STRIP);    
                glVertex2d(polyline[k].x1, polyline[k].y1);
                glVertex2d(polyline[k].x2, polyline[k].y2);
                glEnd();
            }
            
            if (FOCI)
            {
                switch (POLYLINE_PATTERN) {
                    case (P_TOKARSKY):
                    {
                        glLineWidth(2);
                        rgb[0] = 1.0;   rgb[1] = 0.0;   rgb[2] = 0.0;
                        draw_colored_circle(-0.95, 0.0, r, NSEG, rgb);
                        rgb[0] = 0.0;   rgb[1] = 0.8;   rgb[2] = 0.2;
                        draw_colored_circle(0.95, 0.0, r, NSEG, rgb);
                        break;
                    }
                    case (P_TOKA_PRIME):
                    {
                        glLineWidth(2);
                        rgb[0] = 1.0;   rgb[1] = 0.0;   rgb[2] = 0.0;
                        draw_colored_circle(-polyline[84].x1, polyline[84].y1, r, NSEG, rgb);
                        rgb[0] = 0.0;   rgb[1] = 0.8;   rgb[2] = 0.2;
                        draw_colored_circle(polyline[84].x1, polyline[84].y1, r, NSEG, rgb);
                        break;
                    }
                    case (P_TOKA_NONSELF):
                    {
                        glLineWidth(2);
                        rgb[0] = 0.0;   rgb[1] = 0.8;   rgb[2] = 0.2;
                        draw_colored_circle(0.0, 0.0, r, NSEG, rgb);
                        break;
                    }
                }
            }
            if (ABSORBING_CIRCLES)
            {
                rgb[0] = 0.7;   rgb[1] = 0.7;   rgb[2] = 0.7;
                for (k=0; k<ncircles; k++) draw_colored_circle(circles[k].xc, circles[k].yc, circles[k].radius, NSEG, rgb);
            }
            break; 
        }        
        case (D_POLYLINE_ARCS):
        {
            for (k=0; k<nsides; k++)
            {
                glBegin(GL_LINE_STRIP);    
                glVertex2d(polyline[k].x1, polyline[k].y1);
                glVertex2d(polyline[k].x2, polyline[k].y2);
                glEnd();
            }
            
            for (k=0; k<narcs; k++)
            {
                glBegin(GL_LINE_STRIP); 
                for (i=0; i<=NSEG; i++)
                {
                    phi = arcs[k].angle1 + (double)i*(arcs[k].dangle)/(double)NSEG;
                    x = arcs[k].xc + arcs[k].radius*cos(phi);
                    y = arcs[k].yc + arcs[k].radius*sin(phi);
                    glVertex2d(x, y);
                }
                glEnd();
            }
            
            if (FOCI)
            {
                switch (POLYLINE_PATTERN) {
                    case (P_TOKARSKY):
                    {
                        glLineWidth(2);
                        rgb[0] = 1.0;   rgb[1] = 0.0;   rgb[2] = 0.0;
                        draw_colored_circle(-0.95, 0.0, r, NSEG, rgb);
                        rgb[0] = 0.0;   rgb[1] = 0.8;   rgb[2] = 0.2;
                        draw_colored_circle(0.95, 0.0, r, NSEG, rgb);
                        break;
                    }
                    case (P_TOKA_PRIME):
                    {
                        glLineWidth(2);
                        rgb[0] = 1.0;   rgb[1] = 0.0;   rgb[2] = 0.0;
                        draw_colored_circle(-polyline[84].x1, polyline[84].y1, r, NSEG, rgb);
                        rgb[0] = 0.0;   rgb[1] = 0.8;   rgb[2] = 0.2;
                        draw_colored_circle(polyline[84].x1, polyline[84].y1, r, NSEG, rgb);
                        break;
                    }
                    case (P_TOKA_NONSELF):
                    {
                        glLineWidth(2);
                        rgb[0] = 0.0;   rgb[1] = 0.8;   rgb[2] = 0.2;
                        draw_colored_circle(0.0, 0.0, r, NSEG, rgb);
                        break;
                    }
                }
            }
            if (ABSORBING_CIRCLES)
            {
                rgb[0] = 0.7;   rgb[1] = 0.7;   rgb[2] = 0.7;
                for (k=0; k<ncircles; k++) draw_colored_circle(circles[k].xc, circles[k].yc, circles[k].radius, NSEG, rgb);
            }
            break; 
        }        
        default: 
        {
            printf("Function draw_billiard not defined for this billiard \n");
        }
    }
}



/*********************************/
/* computation of the collisions */
/*********************************/

/* The variable config contains information on the state of the particle 
 * and on its next collision with the boundary: 
 * [0] position of next collision (ellipse parametrised by (LAMBDA*cos(s), sin(s))
 * [1] angle to tangent of boundary after next collision
 * [2] running time 
 * [3] initial distance to next collision
 * [4,5] initial position
 * [6,7] coordinates of next collision
 * The running time is incremented until it equals the distance to the next collision
*/


void print_config(double conf[8])  /* for debugging purposes */
{
    printf("s = %.3lg, u = %.3lg, t = %.3lg, L = %.3lg, x0 = %.3lg, y0 = %.3lg, x1 = %.3lg, y1 = %.3lg\n", conf[0], conf[1]/PI, conf[2], conf[3], conf[4], conf[5], conf[6], conf[7]);
}

void print_config_23(double conf[8])  /* for debugging purposes */
{
    printf("t = %.8f, L = %.8f\n", conf[2], conf[3]);
}

void print_colors(int color[NPARTMAX])  /* for debugging purposes */
{
    int i;
    
    for (i=0; i<NPART; i++) printf("%i ", color[i]);
    printf("\n");
}



/****************************************************************************************/
/* rectangle billiard */
/****************************************************************************************/
 
 
 int pos_rectangle(double conf[2], double pos[2], double *alpha)
 /* determine position on boundary of rectangle */
 /* the corners of the rectangle are (-LAMBDA,-1), ..., (LAMBDA,1) */
 /* returns the number of the side hit, or 0 if hitting a corner */
 {
	double s, theta;

	s = conf[0]; 
        if (s<0.0) s = 0.0;
        if (s>4.0*LAMBDA + 4.0) s = 4.0*LAMBDA + 4.0;
        
        theta = conf[1];

        /* we treat the cases of starting in one corner separately */
        /* to avoid numerical problems due to hitting a corner */
        
        /* bottom left corner */
        if ((s==0)||(s==4.0*LAMBDA + 4.0))
        {
            pos[0] = -LAMBDA;
            pos[1] = -1.0;
            if (theta > PID) theta = PID;
            *alpha = theta;
            return(0);
        }
        /* bottom right corner */
        else if (s==2.0*LAMBDA)
        {
            pos[0] = LAMBDA;
            pos[1] = -1.0;
            if (theta > PID) theta = PID;
            *alpha = theta + PID;
            return(0);            
        }
        /* top right corner */
        else if (s==2.0*LAMBDA + 2.0)
        {
            pos[0] = LAMBDA;
            pos[1] = -1.0;
            if (theta > PID) theta = PID;
            *alpha = theta + PI;
            return(0);            
        }
        /* top left corner */
        else if (s==4.0*LAMBDA + 2.0)
        {
            pos[0] = LAMBDA;
            pos[1] = -1.0;
            if (theta > PID) theta = PID;
            *alpha = theta + 3.0*PID;
            return(0);            
        }        
        /* bottom side */
        else if ((s>0)&&(s<2.0*LAMBDA))
        {
            pos[0] = s - LAMBDA;
            pos[1] = -1.0;
            *alpha = theta;
            return(1);
        }
        /* right side */
        else if (s<2.0*LAMBDA + 2.0)
        {
            pos[0] = LAMBDA;
            pos[1] = s - 2.0*LAMBDA - 1.0;
            *alpha = theta + PID;
            return(2);
        }
        /* top side */
        else if (s<4.0*LAMBDA + 2.0)
        {
            pos[0] = 3.0*LAMBDA + 2.0 - s;
            pos[1] = 1.0;
            *alpha = theta + PI;
            return(3);
        }
        /* left side */
        else
        {
            pos[0] = -LAMBDA;
            pos[1] = 4.0*LAMBDA + 3.0 - s;
            *alpha = theta + 3.0*PID;
            return(4);
        }

 }
 
 
 int vrectangle_xy(double config[8], double alpha, double pos[2])
 /* determine initial configuration for start at point pos = (x,y) */
 {
    double l, s0, c0, x1, y1, margin = 1e-12;
    int c, intb=1;

    /* initial position and velocity */

    s0 = sin(alpha);
    c0 = cos(alpha);

    /* intersection with lower part of boundary */
    if (s0<0.0)
    {
        x1 = pos[0] - c0*(1.0 + pos[1])/s0;
        y1 = -1.0;
        if ((x1>=-LAMBDA)&&(x1<=LAMBDA))
        {
            config[0] = x1 + LAMBDA;
            if ((x1 <= -LAMBDA + margin)||(x1 >= LAMBDA -margin)) config[1] = alpha + PI;   /* corners */
//             if ((x1 == -LAMBDA)||(x1 == LAMBDA)) config[1] = alpha + PI;   /* corners */
            else config[1] = -alpha;
            intb = 0;
        }
    }
    /* intersection with right-hand part of boundary */
    if (intb&&(c0>0.0))
    {
        x1 = LAMBDA;
        y1 = pos[1] + s0*(LAMBDA - pos[0])/c0;
        if ((y1>=-1.0)&&(y1<=1.0))
        {
            config[0] = 2.0*LAMBDA + 1.0 + y1;
            if ((y1 <= -1.0 + margin)||(y1 >= 1.0 -margin)) config[1] = alpha + PI;   /* corners */
//             if ((y1 == -1.0)||(y1 == 1.0)) config[1] = alpha + PI;   /* corners */
            else config[1] = PID-alpha;
            intb = 0;
        }
    }
    /* intersection with upper part of boundary */
    if (intb&&(s0>0.0))
    {
        x1 = pos[0] + c0*(1.0 - pos[1])/s0;
        y1 = 1.0;
        if ((x1>=-LAMBDA)&&(x1<=LAMBDA))
        {
            config[0] = 3.0*LAMBDA + 2.0 - x1;
            if ((x1 <= -LAMBDA + margin)||(x1 >= LAMBDA -margin)) config[1] = alpha + PI;   /* corners */
//             if ((x1 == -LAMBDA)||(x1 == LAMBDA)) config[1] = alpha + PI;   /* corners */
            else config[1] = PI-alpha;
            intb = 0;
        }
    }
    /* intersection with left-hand part of boundary */
    if (intb&&(c0<0.0))
    {
        x1 = -LAMBDA;
        y1 = pos[1] + s0*(-LAMBDA - pos[0])/c0;
        if ((y1>=-1.0)&&(y1<=1.0))
        {
            config[0] = 4.0*LAMBDA + 3.0 - y1;
            if ((y1 <= -1.0 + margin)||(y1 >= 1.0 -margin)) config[1] = alpha + PI;   /* corners */
//             if ((y1 == -1.0)||(y1 == 1.0)) config[1] = alpha + PI;   /* corners */
            else config[1] = 3.0*PID-alpha;
            intb = 0;
        }
    }
    
    if (config[1] < 0.0) config[1] += DPI;
       
    config[2] = 0.0;	/* running time */ 
    config[3] = module2(x1-pos[0], y1-pos[1]);
    config[4] = pos[0];
    config[5] = pos[1];
    config[6] = x1;
    config[7] = y1;
 
    return(c);	
 }
 
 int vrectangle(double config[8])
 {
	double pos[2], alpha;
	int c;

	/* position et vitesse de depart */

	c = pos_rectangle(config, pos, &alpha);
        
        vrectangle_xy(config, alpha, pos);

        return(c);
 }
 

/****************************************************************************************/
/* elliptic billiard */
/****************************************************************************************/

 int pos_ellipse(double conf[2], double pos[2], double *alpha)
 /* determine position on boundary of ellipse */
 {
	double theta;

        pos[0] = LAMBDA*cos(conf[0]);
        pos[1] = sin(conf[0]);
        
        theta = argument(-LAMBDA*pos[1],pos[0]/LAMBDA);
        *alpha = theta + conf[1]; 
        
        return(1);
 }
 
 
 int vellipse_xy(double config[8], double alpha, double pos[2])
 /* determine initial configuration for start at point pos = (x,y) */
 {
	double c0, s0, lam2, a, b, c, x1, y1, t, theta;
	int i;

        c0 = cos(alpha);
        s0 = sin(alpha);
        lam2 = LAMBDA*LAMBDA;
        
        /* intersection with ellipse, using parametric equation of line */
        a = c0*c0 + lam2*s0*s0;
        b = pos[0]*c0 + lam2*pos[1]*s0;
        c = pos[0]*pos[0] + lam2*pos[1]*pos[1] - lam2;
        
        t = (-b+sqrt(b*b - a*c))/a;
        x1 = pos[0] + t*c0;
        y1 = pos[1] + t*s0;
        
        /* parameter of intersection with boundary ellipse */
        config[0] = argument(x1/LAMBDA, y1);
        while (config[0] < 0.0) config[0] += DPI;
        while (config[0] > DPI) config[0] -= DPI;
        
        /* computation of outgoing angle after collision with boundary */
        theta = argument(-LAMBDA*y1,x1/LAMBDA);
        config[1] = theta - alpha; 
        while (config[1] < 0.0) config[1] += DPI;
        while (config[1] > DPI) config[1] -= DPI;
   
        config[2] = 0.0;	/* running time */ 
	config[3] = module2(x1-pos[0], y1-pos[1]);     /* distance to collision */
	config[4] = pos[0];    /* start position */
	config[5] = pos[1];
	config[6] = x1;        /* position of collision */
	config[7] = y1;
	
	return(1);
 }

 int vellipse(double config[8])
 /* determine initial configuration when starting from boundary */
 {
	double pos[2], theta, alpha;
	int i;

        pos[0] = LAMBDA*cos(config[0]);
        pos[1] = sin(config[0]);
        
        theta = argument(-LAMBDA*pos[1],pos[0]/LAMBDA);
        alpha = theta + config[1]; 
        
        vellipse_xy(config, alpha, pos);
	
	return(1);
 }
 
/****************************************************************************************/
/* stadium billiard */
/****************************************************************************************/

 int pos_stade(double conf[2], double pos[2], double *alpha)
 /* determine position on boundary of stadium */
 {
	double s, theta, l, psi0, psi;

	s = conf[0]; 
        theta = conf[1];
	l = LAMBDA/2.0;

	if (l >= 0.0)
	{
            if ((s>=0)&&(s<=LAMBDA))
            {
                pos[0] = s - l;
                pos[1] = -1.0;
                *alpha = theta;
                return(0);
            }
            else if (s<=LAMBDA+PI)
            {
                pos[0] = l + sin(s - LAMBDA);
                pos[1] = -cos(s - LAMBDA);
                *alpha = theta + s - LAMBDA;
                return(1);
            }
            else if (s<=2.0*LAMBDA+PI)
            {
                pos[0] = 3.0*l + PI - s;
                pos[1] = 1.0;
                *alpha = theta + PI;
                return(2);
            }
            else
            {
                pos[0] = -l - sin(s - 2.0*LAMBDA - PI);
                pos[1] = cos(s - 2.0*LAMBDA - PI);
                *alpha = theta + s - 2.0*LAMBDA;
                return(3);
            }
	}
	else /* for lens-shaped billiard, to be checked */
	{
            psi0 = asin(-l);
            if ((s>=0)&&(s<=PI-2.0*psi0))
            {
                psi = s + psi0;
                pos[0] = sin(psi) + l;
                pos[1] = -cos(psi);
                *alpha = theta + psi;
                return(0);
            }
            else
            {
                psi = s + 3.0*psi0 - PI;
                pos[0] = - sin(psi) - l;
                pos[1] = cos(psi);
                *alpha = theta + psi + PI;
                return(2);
            }
	}
 }
 
 
 int vstade_xy(double config[8], double alpha, double pos[2])
 /* determine initial configuration for start at point pos = (x,y) */
 {
	double l, s0, c0, t, x, y, x1, y1, a, b, res[2];
	double smin, psi, margin = 1e-12;
	int c, intb=1, intc, i;

	/* initial position and velocity */

	l = LAMBDA/2.0;
	if (l>=0.0) smin = 0.0; else smin = -l;
	s0 = sin(alpha);
	c0 = cos(alpha);

	/* intersection with lower straight part of boundary */
	if ((s0<0.0)&&(l>0))
	{
            x1 = pos[0] + c0*(-1.0 - pos[1])/s0;
            y1 = -1.0;
            if ((x1>=-l)&&(x1<=l))
            {
                config[0] = x1 + l;
                config[1] = -alpha;
                intb = 0;
            }
	}
	/* intersection with upper straight part of boundary */
	if (intb&&(s0>0.0)&&(l>0))
	{
            x1 = pos[0] + c0*(1.0 - pos[1])/s0;
            y1 = 1.0;
            if ((x1>=-l)&&(x1<=l))
            {
                config[0] = 3.0*l + PI - x1;
                config[1] = PI-alpha;
                intb = 0;
            }
	}
	/* intersection with right-hand arc of boundary */
	if (intb)
	{
            a = 2.0*pos[0]*c0 + 2.0*pos[1]*s0 - LAMBDA*c0;
            b = pos[0]*pos[0] + pos[1]*pos[1] + l*l - LAMBDA*pos[0] - 1.0;
            intc = polynome(1.0, a, b, res);
            if (intc) for(i=0; i<2; i++)
            {
                x = pos[0] + c0*res[i];
                y = pos[1] + s0*res[i];
                psi = argument(-y, x-l);
                if (intb&&(sin(psi) >= smin)&&(res[i]>margin))
                {
                    if (l>0.0) config[0] = LAMBDA + psi;
                    else config[0] = psi - asin(-l);
                    config[1] = -alpha + psi;
                    intb = 0;
                    x1 = x; y1 = y;
                }
            }
	}
	/* intersection with left-hand arc of boundary */
	if (intb)
	{
            a = 2.0*pos[0]*c0 + 2.0*pos[1]*s0 + LAMBDA*c0;
            b = pos[0]*pos[0] + pos[1]*pos[1] + l*l + LAMBDA*pos[0] - 1.0;
            intc = polynome(1.0, a, b, res);
            if (intc) for(i=0; i<2; i++)
            {
                x = pos[0] + c0*res[i];
                y = pos[1] + s0*res[i];
                psi = argument(y, -l-x);
                if (intb&&(sin(psi) >= smin)&&(res[i]>margin))
                {
                    if (l>0.0) config[0] = 2.0*LAMBDA + PI + psi;
                    else config[0] = psi - 3.0*asin(-l) + PI;
                    config[1] = -alpha + psi + PI;
                    intb = 0;
                    x1 = x; y1 = y;
                }
            }
	}

	config[2] = 0.0;	/* running time */ 
	config[3] = module2(x1-pos[0], y1-pos[1]);
	config[4] = pos[0];
	config[5] = pos[1];
	config[6] = x1;
	config[7] = y1;
 
	return(c);
 }
 
 int vstade(double config[8])
 {
	double alpha, pos[2];
	int c;

	c = pos_stade(config, pos, &alpha);
        
        vstade_xy(config, alpha, pos);

        return(c);
 }
 
/****************************************************************************************/
/* Sinai billiard */
/****************************************************************************************/
 
 int pos_sinai(double conf[2], double pos[2], double *alpha)
 /* determine position on boundary of Sinai billiard */
 /* s in [0,2 Pi) is on the circle, other s are on boundary of window */
 {
	double s, theta, psi0, psi, s1, s2, s3, s4;

	s = conf[0]; 
        theta = conf[1];
        if (conf[1] < 0.0) conf[1] += DPI;
        
        s1 = DPI + XMAX - XMIN;
        s2 = s1 + YMAX - YMIN;
        s3 = s2 + XMAX - XMIN;
        s4 = s3 + YMAX - YMIN;

        if (s < DPI)    /* circle */
        {
            pos[0] = LAMBDA*cos(s);
            pos[1] = LAMBDA*sin(s);
            theta = PID + s;
            *alpha = theta - conf[1]; 
            return(0);
        }
        else if (s < s1)    /* boundary of window */
        {
            pos[0] = XMIN + s - DPI;
            pos[1] = YMIN;
            *alpha = conf[1];
            return(-1);
        }
        else if (s < s2)
        {
            pos[0] = XMAX;
            pos[1] = YMIN + s - s1;
            *alpha = conf[1];
            return(-2);
        }
        else if (s < s3)
        {
            pos[0] = XMAX - s + s2;
            pos[1] = YMAX;
            *alpha = conf[1];
            return(-3);
        }
        else 
        {
            pos[0] = XMIN;
            pos[1] = YMAX - s + s3;
            *alpha = conf[1];
            return(-4);
        }
 }
 
 
 int vsinai_xy(double config[8], double alpha, double pos[2])
 /* determine initial configuration for start at point pos = (x,y) */
 {
	double l, s0, c0, t, t1, x, y, x1, y1, a, b, delta, res[2], s1, s2, s3, s4, s, r;
	double psi, lam2, margin = 1e-12;
	int c, intb=1, intc, i;

	/* initial position and velocity */

        c0 = cos(alpha);
        s0 = sin(alpha);
        s1 = DPI + XMAX - XMIN;
        s2 = s1 + YMAX - YMIN;
        s3 = s2 + XMAX - XMIN;
        s4 = s3 + YMAX - YMIN;
       
        /* intersection with circle, using parametric equation of line */
        b = pos[0]*c0 + pos[1]*s0;
        a = pos[0]*pos[0] + pos[1]*pos[1] - LAMBDA*LAMBDA;
        delta = b*b - a;
                
        if ((delta > margin)&&(a > margin))
        {
            t = - b - sqrt(delta);  
            x1 = pos[0] + t*c0;
            y1 = pos[1] + t*s0;
            s = argument(x1,y1);
            if (s<0.0) s += DPI;
            if (s>=DPI) s -= DPI;
            config[0] = s;
            config[1] = 3.0*PID - s + alpha;
            c = 0;
        }
        else if (c0 > 0.0)   /* intersection with boundary of window */
        {
            y1 = pos[1] + (XMAX - pos[0])*s0/c0;
            if ((y1 >= YMIN)&&(y1 <= YMAX))     /* hitting right boundary */
            {
                x1 = XMAX;
                config[0] = s3 + YMAX - y1;
                config[1] = alpha;
                c = 2;
            }
            else if (s0 > 0.0)      /* hitting upper boundary */
            {
                x1 = pos[0] + (YMAX - pos[1])*c0/s0;
                y1 = YMAX; 
                config[0] = DPI + x1 - XMIN;
                config[1] = alpha;
                c = 3;
            }
            else                    /* hitting lower boundary */
            {
                x1 = pos[0] + (YMIN - pos[1])*c0/s0;
                y1 = YMIN; 
                config[0] = s2 + XMAX - x1;
                config[1] = alpha;    
                c = 1;
            }
        }
        else if (c0 < 0.0)
        {
            y1 = pos[1] + (XMIN - pos[0])*s0/c0;
            if ((y1 >= YMIN)&&(y1 <= YMAX))     /* hitting left boundary */
            {
                x1 = XMIN;
                config[0] = s1 + y1 - YMIN;
                config[1] = alpha;
                c = 4;
            }
            else if (s0 > 0.0)      /* hitting upper boundary */
            {
                x1 = pos[0] + (YMAX - pos[1])*c0/s0;
                y1 = YMAX; 
                config[0] = DPI + x1 - XMIN;
                config[1] = alpha;
                c = 3;
            }
            else                    /* hitting lower boundary */
            {
                x1 = pos[0] + (YMIN - pos[1])*c0/s0;
                y1 = YMIN; 
                config[0] = s2 + XMAX - x1;
                config[1] = alpha;
                c = 1;
            }          
        }
        else    /* vertical motion */
        {
            if (s0 > 0.0)
            {
                x1 = pos[0];
                y1 = YMAX;
                config[0] = DPI + x1 - XMIN;
                config[1] = alpha;
                c = 3;
            }
            else
            {
                x1 = pos[0];
                y1 = YMIN;
                config[0] = s2 + XMAX - x1;
                config[1] = alpha;
                c = 1;
            }
        }
        
        if (config[1] < 0.0) config[1] += DPI;

	config[2] = 0.0;	/* running time */ 
	config[3] = module2(x1-pos[0], y1-pos[1]);
	config[4] = pos[0];
	config[5] = pos[1];
	config[6] = x1;
	config[7] = y1;
         
	return(-c);
        /* return a negative value if the disc is not hit, for color scheme */
 }
 
 int vsinai(double config[8])
 {
	double alpha, pos[2];
	int c;

	/* position et vitesse de depart */

	c = pos_sinai(config, pos, &alpha);
        
        vsinai_xy(config, alpha, pos);

        return(c);
 }
 
 
/****************************************************************************************/
/* triangle billiard */
/****************************************************************************************/
 
 
 int pos_triangle(double conf[2], double pos[2], double *alpha)
 /* determine position on boundary of triangle */
 /* the corners of the triangle are (-LAMBDA,-1), (LAMBDA,-1), (-LAMBDA,1) */
 /* we use arclength for horizontal and vertical side, x for diagonal */
 {
	double s, theta;

	s = conf[0]; 
        theta = conf[1];

        if ((s>=0)&&(s<=2.0*LAMBDA))
        {
            pos[0] = s - LAMBDA;
            pos[1] = -1.0;
            *alpha = theta;
            return(0);
        }
        else if (s<=4.0*LAMBDA)
        {
            pos[0] = 3.0*LAMBDA - s;
            pos[1] = -3.0 + s/LAMBDA;
            *alpha = theta + PI - argument(LAMBDA, 1.0);
            return(1);
        }
        else
        {
            pos[0] = -LAMBDA;
            pos[1] = 4.0*LAMBDA + 1.0 - s;
            *alpha = theta + 3.0*PID;
            return(2);
        }

 }
 
 
 int vtriangle_xy(double config[8], double alpha, double pos[2])
 /* determine initial configuration for start at point pos = (x,y) */
 /* Warning: reflection in the corners is not yet implemented correctly */
 {
	double s0, c0, t, x, y, x1, y1, psi;
	int c, intb=1, intc, i;

	/* initial position and velocity */

	s0 = sin(alpha);
	c0 = cos(alpha);

	/* intersection with lower part of boundary */
// 	if ((s0<0.0)&&(pos[1]>0.0))
	if (s0<0.0)
	{
            x1 = pos[0] - c0*(1.0 + pos[1])/s0;
            y1 = -1.0;
            if ((x1>=-LAMBDA)&&(x1<=LAMBDA))
            {
                config[0] = x1 + LAMBDA;
                config[1] = -alpha;
                intb = 0;
            }
	}
	/* intersection with left-hand part of boundary */
	if (intb&&(c0<0.0))
	{
            x1 = -LAMBDA;
            y1 = pos[1] + s0*(-LAMBDA - pos[0])/c0;
            if ((y1>=-1.0)&&(y1<=1.0))
            {
                config[0] = 4.0*LAMBDA + 1.0 - y1;
                config[1] = 3.0*PID-alpha;
                intb = 0;
            }
	}
	/* intersection with diagonal part of boundary */
	if (intb)
	{
            t = -(pos[0] + LAMBDA*pos[1])/(c0 + LAMBDA*s0);
            x1 = pos[0] + t*c0;
            y1 = pos[1] + t*s0;
            if ((x1>=-LAMBDA)&&(x1<=LAMBDA))
            {
                psi = argument(LAMBDA, 1.0);
                config[0] = 3.0*LAMBDA - x1;
                config[1] = PI - alpha - psi;
//                 config[1] = PI - alpha - atan(1.0/LAMBDA);
                intb = 0;
            }
	}

	config[2] = 0.0;	/* running time */ 
	config[3] = module2(x1-pos[0], y1-pos[1]);
	config[4] = pos[0];
	config[5] = pos[1];
	config[6] = x1;
	config[7] = y1;
 
	return(c);
 }
 
 int vtriangle(double config[8])
 {
	double alpha, pos[2];
	int c;

	/* position et vitesse de depart */

	c = pos_triangle(config, pos, &alpha);
        
        vtriangle_xy(config, alpha, pos);

        return(c);
 }
 

/****************************************************************************************/
 /* annulus billiard */
/****************************************************************************************/

 int pos_annulus(double conf[2], double pos[2], double *alpha)
 /* determine position on boundary of annulus */
 {
	double s, theta, psi0, psi, s1, s2, s3, s4;

	s = conf[0]; 
        theta = conf[1];
        if (conf[1] < 0.0) conf[1] += DPI;
 
        if (conf[0] < DPI)      /* inner circle */
        {
            pos[0] = LAMBDA*cos(conf[0]) + MU;
            pos[1] = LAMBDA*sin(conf[0]);
        
            theta = PID + conf[0];
            *alpha = theta - conf[1]; 
            return(0);
        }
        else                    /* outer circle */
        {
            pos[0] = cos(conf[0]);
            pos[1] = sin(conf[0]);
        
            theta = argument(-pos[1],pos[0]);
            *alpha = theta + conf[1]; 
            return(1);
        }
 }
 
 int vannulus_xy(double config[8], double alpha, double pos[2])
 /* determine initial configuration for start at point pos = (x,y) */
 {
	double l, s0, c0, t, t1, x, y, x1, y1, a, b, delta, res[2], s, r;
	double psi, lam2, margin = 1.0e-14, theta;
// 	double psi, lam2, margin = 1.0e-14, theta;
	int c, intb=1, intc, i;

	/* initial position and velocity */

        c0 = cos(alpha);
        s0 = sin(alpha);
       
        /* intersection with inner circle, using parametric equation of line */
        b = (pos[0]-MU)*c0 + pos[1]*s0;
        a = (pos[0]-MU)*(pos[0]-MU) + pos[1]*pos[1] - LAMBDA*LAMBDA;
        delta = b*b - a;
                
        if ((delta > margin)&&(a > margin))
        {
            t = - b - sqrt(delta);  
            x1 = pos[0] + t*c0;
            y1 = pos[1] + t*s0;
            s = argument(x1-MU,y1);
            while (s<0.0) s += DPI;
            while (s>=DPI) s -= DPI;
            config[0] = s;
            config[1] = 3.0*PID - s + alpha;
            c = 0;
        }
        else    /* intersection with outer circle, using parametric equation of line */
        {
            b = pos[0]*c0 + pos[1]*s0;
            a = pos[0]*pos[0] + pos[1]*pos[1] - 1.0;
        
            t = (-b+sqrt(b*b - a));
            x1 = pos[0] + t*c0;
            y1 = pos[1] + t*s0;
        
            /* parameter of intersection with outer circle */
            config[0] = argument(x1, y1);
            while (config[0] < DPI) config[0] += DPI;
            while (config[0] >= 2.0*DPI) config[0] -= DPI;
        
            /* computation of outgoing angle after collision with outer circle */
            theta = argument(-y1,x1);
            config[1] = theta - alpha; 
//             while (config[1] < 0.0) config[1] += DPI;
//             while (config[1] > DPI) config[1] -= DPI; 
            c = 1;
        }
   
        if (config[1] < 0.0) config[1] += DPI;

//         config[2] = 1.0e-12;	/* running time */ 
        config[2] = 0.0;	/* running time */ 
	config[3] = module2(x1-pos[0], y1-pos[1]);     /* distance to collision */
	config[4] = pos[0];    /* start position */
	config[5] = pos[1];
	config[6] = x1;        /* position of collision */
	config[7] = y1;
	
	return(c);
 }

 int vannulus(double config[8])
 /* determine initial configuration when starting from boundary */
 {
	double pos[2], alpha;
	int c;

        c = pos_annulus(config, pos, &alpha);
        
        vannulus_xy(config, alpha, pos);
	
	return(c);
 }
 
/****************************************************************************************/
/* polygonal billiard */
/****************************************************************************************/

 int pos_polygon(double conf[2], double pos[2], double *alpha)
 /* determine position on boundary of polygon */
 /* conf[0] is arclength on boundary */
 {
	double s, theta, omega, length, s1, angle, x, y;
        int c;

	s = conf[0]; 
        theta = conf[1];
        
        omega = DPI/((double)NPOLY);
        length = 2.0*sin(0.5*omega);

        c = (int)(s/length);         /* side of polygon */
        
        s1 = s - ((double)c)*length;
        
        x = 1.0 + (cos(omega) - 1.0)*s1/length;
        y = sin(omega)*s1/length;
        
        angle = (double)c*omega + PID*APOLY;
        
        pos[0] = x*cos(angle) - y*sin(angle);
        pos[1] = x*sin(angle) + y*cos(angle);
        
        *alpha = (0.5 + (double)c)*omega + theta + PID*(1.0 + APOLY);
        
        return(c);
 }
 
 int vpolygon_xy(double config[8], double alpha, double pos[2])
 /* determine initial configuration for start at point pos = (x,y) */
 {
	double s, theta, omega, length, rlength, s1, rangle, x, y, xp, yp, x1, y1, ca, sa;
	int k, c, intb=1, intc, i;

        /* dimensions/angles of polygon */
        omega = DPI/((double)NPOLY);
        length = 2.0*sin(0.5*omega);
        rlength = cos(0.5*omega);
        
        for (k=0; k<NPOLY; k++)
        {
            /* rotate position so that kth side is vertical */
            rangle = (0.5 + (double)k)*omega + APOLY*PID;
            theta = alpha - rangle;
                
            if ((cos(theta) > 0.0)&&(intb))
            {
                ca = cos(rangle);
                sa = sin(rangle);
                
                x = pos[0]*ca + pos[1]*sa;
                y = -pos[0]*sa + pos[1]*ca;
            
                xp = rlength;
                yp = y + (xp-x)*tan(theta);
                
                if (vabs(yp) < 0.5*length) 
                {
                    /* rotate back */
                    x1 = xp*ca - yp*sa;
                    y1 = xp*sa + yp*ca;
                    
                    intb = 0;
                    c = k;
                    config[0] = ((double)k + 0.5)*length + yp;
                    config[1] = PID - theta;
                }
            }
        }
       
        if (config[1] < 0.0) config[1] += DPI;

        config[2] = 0.0;	/* running time */ 
	config[3] = module2(x1-pos[0], y1-pos[1]);     /* distance to collision */
	config[4] = pos[0];    /* start position */
	config[5] = pos[1];
	config[6] = x1;        /* position of collision */
	config[7] = y1;
	
	return(c);
 }

 int vpolygon(double config[8])
 /* determine initial configuration when starting from boundary */
 {
	double pos[2], alpha;
	int c;

        c = pos_polygon(config, pos, &alpha);
        
        vpolygon_xy(config, alpha, pos);
	
	return(c);
 }
 
/****************************************************************************************/
/* Reuleaux-type and star-shaped billiard */
/****************************************************************************************/

 int pos_reuleaux(double conf[2], double pos[2], double *alpha)
 /* determine position on boundary of polygon */
 /* conf[0] is arclength on boundary */
 {
	double s, theta, omega2, beta2, beta, s1, angle, x2, x, y;
        int c;

	s = conf[0]; 
        theta = conf[1];
        
        omega2 = PI/((double)NPOLY);
        beta2 = asin(sin(omega2)/vabs(LAMBDA));
        beta = beta2*2.0;

        c = (int)(s/beta);         /* side of shape */
        
        s1 = s - ((double)c)*beta;
        
        if (LAMBDA > 0.0) x2 = cos(omega2) + sqrt(LAMBDA*LAMBDA - sin(omega2)*sin(omega2));
        else x2 = cos(omega2) - sqrt(LAMBDA*LAMBDA - sin(omega2)*sin(omega2));
        
        x = x2 - LAMBDA*cos(s1 - beta2);
        y = LAMBDA*sin(s1 - beta2);
        
        angle = 2.0*((double)c)*omega2 + PID*APOLY;
        
        pos[0] = x*cos(angle) - y*sin(angle);
        pos[1] = x*sin(angle) + y*cos(angle);
        
        *alpha = PID - s1 + beta2 + theta + 2.0*(double)c*omega2 + APOLY*PID;
        
//         printf("alpha = %.5lg\t", *alpha);
        
        return(c);
 }
 
 int vreuleaux_xy(double config[8], double alpha, double pos[2])
 /* determine initial configuration for start at point pos = (x,y) */
 {
	double s, theta, omega2, beta, s1, rangle, x, y, x1[NPOLY], y1[NPOLY], xi, yi, t, x2;
        double ca, sa, a, b, margin = 1.0e-14, tmin, tval[NPOLY], tempconf[NPOLY][2];
	int k, c, intb=1, intc, i, nt = 0, cval[NPOLY], ntmin;

        /* dimensions/angles of polygon */
        omega2 = PI/((double)NPOLY);
        beta = 2.0*asin(sin(omega2)/vabs(LAMBDA));
//         printf("beta = %.5lg\n", beta);
        
        if (LAMBDA > 0.0) x2 = cos(omega2) + sqrt(LAMBDA*LAMBDA - sin(omega2)*sin(omega2));
        else x2 = cos(omega2) - sqrt(LAMBDA*LAMBDA - sin(omega2)*sin(omega2));
//         printf("x2 = %.5lg\n", x2);
        
        for (k=0; k<NPOLY; k++)
        {
            /* rotate position so that kth side is vertical */
            rangle = 2.0*(double)k*omega2 + APOLY*PID;
            theta = alpha - rangle;
            
            ca = cos(rangle);
            sa = sin(rangle);
                                
            x = pos[0]*ca + pos[1]*sa;
            y = -pos[0]*sa + pos[1]*ca;
                
            a = (x-x2)*cos(theta) + y*sin(theta);
            b = (x-x2)*(x-x2) + y*y - LAMBDA*LAMBDA;
                
            if (a*a - b > margin)
            {
                if (LAMBDA > 0.0) t = -a - sqrt(a*a - b);
                else t = -a + sqrt(a*a - b);
                    
                xi = x + t*cos(theta);
                yi = y + t*sin(theta);
                    
                if ((t > margin)&&(vabs(yi) <= sin(omega2))) 
                {
                    cval[nt] = k;
                    tval[nt] = t;
                        
                    /* rotate back */
                    x1[nt] = xi*ca - yi*sa;
                    y1[nt] = xi*sa + yi*ca;
                    
                    tempconf[nt][0] = ((double)k + 0.5)*beta + asin(yi/LAMBDA);
                    tempconf[nt][1] = PID - asin(yi/LAMBDA) - theta;      
                    nt++;
                }
            }
        }
        
        /* find earliest intersection */
        tmin = tval[0];
        ntmin = 0;
        for (i=1; i<nt; i++) 
            if (tval[i] < tmin) 
            {
                tmin = tval[i];
                ntmin = i;
            }
            
        config[0] = tempconf[ntmin][0];
        config[1] = tempconf[ntmin][1];
        c = cval[ntmin];
 
//         printf("nt = %i\t ntmin = %i \tcmin = %i\n", nt, ntmin, c);
        

        if (config[1] < 0.0) config[1] += DPI;

        config[2] = 0.0;	/* running time */ 
	config[3] = module2(x1[ntmin]-pos[0], y1[ntmin]-pos[1]);     /* distance to collision */
	config[4] = pos[0];    /* start position */
	config[5] = pos[1];
	config[6] = x1[ntmin];        /* position of collision */
	config[7] = y1[ntmin];
	
	return(c);
 }

 int vreuleaux(double config[8])
 /* determine initial configuration when starting from boundary */
 {
	double pos[2], alpha;
	int c;

        c = pos_reuleaux(config, pos, &alpha);
        
        vreuleaux_xy(config, alpha, pos);
	
	return(c);
 }

/****************************************************************************************/
/* Bunimovich flower billiard */
/****************************************************************************************/

 int pos_flower(double conf[2], double pos[2], double *alpha)
 /* determine position on boundary of polygon */
 /* conf[0] is arclength on boundary, it belongs to [0,2*NPOLY*phimax) */
 {
	double s, theta, omega, co, so, axis1, axis2, phimax, s1, x, y, angle;
        int c;

	s = conf[0]; 
        theta = conf[1];
        
        compute_flower_parameters(&omega, &co, &so, &axis1, &axis2, &phimax);
        
        c = (int)(s/(2.0*phimax));         /* side of shape */
        
        s1 = s - (((double)c)*2.0 + 1.0)*phimax;
        
        x = co + axis1*cos(s1);
        y = axis2*sin(s1);
        
        angle = ((double)c)*omega + PID*APOLY;
//         angle = 2.0*((double)c)*omega + PID*APOLY;
        
        pos[0] = x*cos(angle) - y*sin(angle);
        pos[1] = x*sin(angle) + y*cos(angle);
        
        *alpha = argument(-axis1*sin(s1), axis2*cos(s1)) + theta + angle;
        
//         printf("alpha = %.5lg\t", *alpha);
        
        return(c);
 }
 
 int vflower_xy(double config[8], double alpha, double pos[2])
 /* determine initial configuration for start at point pos = (x,y) */
//  double config[8], alpha, pos[2];
 {
	double s, theta, omega, omega2, s1, rangle, x, y, x1, y1, xi, yi, t;
        double ca, sa, aa, bb, cc, margin = 1.0e-14, tmin;
        double co, so, co2, so2, ct, st, phimax, phi, axis1, axis2;
	int k, c, intb=1, intc, i, nt = 0, ntmin, sign;

        compute_flower_parameters(&omega, &co, &so, &axis1, &axis2, &phimax);
        
        for (k=0; k<NPOLY; k++) if (intb)
        {
            /* rotate position so that kth side is vertical */
//             rangle = (double)(2*k)*omega + APOLY*PID;
            rangle = (double)k*omega + APOLY*PID;
            theta = alpha - rangle;
            
            ca = cos(rangle);
            sa = sin(rangle);
            
            ct = cos(theta);
            st = sin(theta);
                                
            x = pos[0]*ca + pos[1]*sa;
            y = -pos[0]*sa + pos[1]*ca;
                
            /* find intersection with elliptical arc */
            aa = ct*ct/(axis1*axis1) + st*st/(axis2*axis2);
            bb = (x-co)*ct/(axis1*axis1) + y*st/(axis2*axis2);
            cc = (x-co)*(x-co)/(axis1*axis1) + y*y/(axis2*axis2) - 1.0;
                
//             if (bb*bb - aa*cc > margin) 
            if (bb*bb - aa*cc >= 0.0) 
            {
                t = (-bb + sqrt(bb*bb - aa*cc))/aa;
                    
                xi = x + t*cos(theta);
                yi = y + t*sin(theta);
                
                if (yi >= 0.0) phi = argument((xi - co)/axis1, yi/axis2);
                else phi = -argument((xi - co)/axis1, -yi/axis2);
                
                if ((t > margin)&&((vabs(phi) <= phimax)||(vabs(phi-DPI) <= phimax))) 
                {
                    intb = 0;
                    c = k;
                        
                    /* rotate back */
                    x1 = xi*ca - yi*sa;
                    y1 = xi*sa + yi*ca;
                    
                    config[0] = (double)(2*k + 1)*phimax + phi;
                    config[1] = argument(-axis1*sin(phi), axis2*cos(phi)) - theta;      
                }
            }
        }
        
//        if (nt == 0) printf("nt = %i\t ntmin = %i \tcmin = %i\n", nt, ntmin, c);

        if (config[1] < 0.0) config[1] += DPI;

        config[2] = 0.0;	/* running time */ 
	config[3] = module2(x1-pos[0], y1-pos[1]);     /* distance to collision */
	config[4] = pos[0];    /* start position */
	config[5] = pos[1];
	config[6] = x1;        /* position of collision */
	config[7] = y1;
	
	return(c);
 }

//   int old_vflower_xy(config, alpha, pos)
//  /* determine initial configuration for start at point pos = (x,y) */
//  double config[8], alpha, pos[2];
// 
//  {
// 	double s, theta, omega, omega2, s1, rangle, x, y, x1[2*NPOLY], y1[2*NPOLY], xi, yi, t;
//         double ca, sa, aa, bb, cc, margin = 1.0e-14, tmin, tval[2*NPOLY], tempconf[2*NPOLY][2];
//         double co, so, co2, so2, ct, st, phimax, phi, axis1, axis2;
// 	int k, c, intb=1, intc, i, nt = 0, cval[2*NPOLY], ntmin, sign;
// 
//         compute_flower_parameters(&omega, &co, &so, &axis1, &axis2, &phimax);
//         
//         for (k=0; k<NPOLY; k++)
//         {
//             /* rotate position so that kth side is vertical */
// //             rangle = (double)(2*k)*omega + APOLY*PID;
//             rangle = (double)k*omega + APOLY*PID;
//             theta = alpha - rangle;
//             
//             ca = cos(rangle);
//             sa = sin(rangle);
//             
//             ct = cos(theta);
//             st = sin(theta);
//                                 
//             x = pos[0]*ca + pos[1]*sa;
//             y = -pos[0]*sa + pos[1]*ca;
//                 
//             /* find intersection with elliptical arc */
//             aa = ct*ct/(axis1*axis1) + st*st/(axis2*axis2);
//             bb = (x-co)*ct/(axis1*axis1) + y*st/(axis2*axis2);
//             cc = (x-co)*(x-co)/(axis1*axis1) + y*y/(axis2*axis2) - 1.0;
//                 
// //             if (bb*bb - aa*cc > margin) 
//             if (bb*bb - aa*cc >= 0.0) 
//             {
//                 t = (-bb + sqrt(bb*bb - aa*cc))/aa;
//                     
//                 xi = x + t*cos(theta);
//                 yi = y + t*sin(theta);
//                 
//                 if (yi >= 0.0) phi = argument((xi - co)/axis1, yi/axis2);
//                 else phi = -argument((xi - co)/axis1, -yi/axis2);
//                 
// //                 phi = argument((xi - co)/axis1, yi/axis2);
// //                 if (phi > PI) phi += -DPI;
//                     
//                 if ((t > margin)&&((vabs(phi) <= phimax)||(vabs(phi-DPI) <= phimax))) 
// //                 if (((vabs(phi) <= phimax)||(vabs(phi-DPI) <= phimax))) 
// //                 if ((t > margin)) 
//                 {
//                     cval[nt] = k;
// //                     cval[nt] = 2*k;
//                     tval[nt] = t;
//                         
//                     /* rotate back */
//                     x1[nt] = xi*ca - yi*sa;
//                     y1[nt] = xi*sa + yi*ca;
//                     
//                     tempconf[nt][0] = (double)(2*k + 1)*phimax + phi;
//                     tempconf[nt][1] = argument(-axis1*sin(phi), axis2*cos(phi)) - theta;      
//                     nt++;
//                 }
//             }
//         }
//         
//         /* find earliest intersection */
//         tmin = tval[0];
//         ntmin = 0;
//         for (i=1; i<nt; i++) 
//             if (tval[i] < tmin) 
//             {
//                 tmin = tval[i];
//                 ntmin = i;
//             }
//             
//         config[0] = tempconf[ntmin][0];
//         config[1] = tempconf[ntmin][1];
//         c = cval[ntmin];
//  
//        if (nt == 0) printf("nt = %i\t ntmin = %i \tcmin = %i\n", nt, ntmin, c);
//         
// 
//         if (config[1] < 0.0) config[1] += DPI;
// 
//         config[2] = 0.0;	/* running time */ 
// 	config[3] = module2(x1[ntmin]-pos[0], y1[ntmin]-pos[1]);     /* distance to collision */
// 	config[4] = pos[0];    /* start position */
// 	config[5] = pos[1];
// 	config[6] = x1[ntmin];        /* position of collision */
// 	config[7] = y1[ntmin];
// 	
// 	return(c);
//  }
 
 int vflower(double config[8])
 /* determine initial configuration when starting from boundary */
 {
	double pos[2], alpha;
	int c;

        c = pos_flower(config, pos, &alpha);
        
        vflower_xy(config, alpha, pos);
	
	return(c);
 } 
 
 /****************************************************************************************/
/* Alternating between Reuleaux-type and star-shaped billiard */
/****************************************************************************************/

 int pos_alt_reuleaux(double conf[2], double pos[2], double *alpha)
 /* determine position on boundary of polygon */
 /* conf[0] is arclength on boundary */
 {
	double s, theta, omega2, beta2, beta, s1, angle, x2plus, x2minus, x, y;
        int c;

	s = conf[0]; 
        theta = conf[1];
        
        omega2 = PI/((double)NPOLY);
        beta2 = asin(sin(omega2)/vabs(LAMBDA));
        beta = beta2*2.0;

        c = (int)(s/beta);         /* side of shape */
        
        s1 = s - ((double)c)*beta;
        
        x2plus = cos(omega2) + sqrt(LAMBDA*LAMBDA - sin(omega2)*sin(omega2));
        x2minus = cos(omega2) - sqrt(LAMBDA*LAMBDA - sin(omega2)*sin(omega2));
        
        if (c%2 == 0) x = x2plus - LAMBDA*cos(s1 - beta2);
        else x = x2minus + LAMBDA*cos(s1 - beta2);
        
        if (c%2 == 0) y = LAMBDA*sin(s1 - beta2);
        else y = -LAMBDA*sin(s1 - beta2);
        
        /* test, to be removed */
//         x = x2plus - LAMBDA*cos(s1 - beta2);
//         y = LAMBDA*sin(s1 - beta2);
        
        angle = 2.0*((double)c)*omega2 + PID*APOLY;
        
        pos[0] = x*cos(angle) - y*sin(angle);
        pos[1] = x*sin(angle) + y*cos(angle);
        
        *alpha = PID - s1 + beta2 + theta + 2.0*(double)c*omega2 + APOLY*PID;
        
//         printf("alpha = %.5lg\t", *alpha);
        
        return(c);
 }
 
 int valt_reuleaux_xy(double config[8], double alpha, double pos[2])
 /* determine initial configuration for start at point pos = (x,y) */
 {
	double s, theta, omega2, beta, s1, rangle, x, y, x1[NPOLY], y1[NPOLY], xi, yi, t, x2plus, x2minus, arcsine;
        double ca, sa, a, b, margin = 1.0e-14, tmin, tval[NPOLY], tempconf[NPOLY][2];
	int k, c, intb=1, intc, i, nt = 0, cval[NPOLY], ntmin;

        /* dimensions/angles of polygon */
        omega2 = PI/((double)NPOLY);
        beta = 2.0*asin(sin(omega2)/vabs(LAMBDA));
//         printf("beta = %.5lg\n", beta);
        
        x2plus = cos(omega2) + sqrt(LAMBDA*LAMBDA - sin(omega2)*sin(omega2));
        x2minus = cos(omega2) - sqrt(LAMBDA*LAMBDA - sin(omega2)*sin(omega2));
//         printf("x2 = %.5lg\n", x2);
        
        for (k=0; k<NPOLY; k++)
        {
            /* rotate position so that kth side is vertical */
            rangle = 2.0*(double)k*omega2 + APOLY*PID;
            theta = alpha - rangle;
            
            ca = cos(rangle);
            sa = sin(rangle);
                                
            x = pos[0]*ca + pos[1]*sa;
            y = -pos[0]*sa + pos[1]*ca;
                
            if (k%2==0)
            {
                a = (x-x2plus)*cos(theta) + y*sin(theta);
                b = (x-x2plus)*(x-x2plus) + y*y - LAMBDA*LAMBDA;
            }
            else 
            {
                a = (x-x2minus)*cos(theta) + y*sin(theta);
                b = (x-x2minus)*(x-x2minus) + y*y - LAMBDA*LAMBDA;
            }
    
            if (a*a - b > margin)
            {
                if (k%2==0) t = -a - sqrt(a*a - b);
                else t = -a + sqrt(a*a - b);
                    
                xi = x + t*cos(theta);
                yi = y + t*sin(theta);
                    
                if ((t > margin)&&(vabs(yi) <= sin(omega2))) 
                {
                    cval[nt] = k;
                    tval[nt] = t;
                        
                    /* rotate back */
                    x1[nt] = xi*ca - yi*sa;
                    y1[nt] = xi*sa + yi*ca;
                    
                    if (k%2==0) arcsine = asin(yi/LAMBDA);
                    else arcsine = -asin(yi/LAMBDA);
                    
                    tempconf[nt][0] = ((double)k + 0.5)*beta + arcsine;
                    tempconf[nt][1] = PID - arcsine - theta;      
                    nt++;
                }
            }
        }
        
        /* find earliest intersection */
        tmin = tval[0];
        ntmin = 0;
        for (i=1; i<nt; i++) 
            if (tval[i] < tmin) 
            {
                tmin = tval[i];
                ntmin = i;
            }
            
        config[0] = tempconf[ntmin][0];
        config[1] = tempconf[ntmin][1];
        c = cval[ntmin];
 
//         printf("nt = %i\t ntmin = %i \tcmin = %i\n", nt, ntmin, c);
        

        if (config[1] < 0.0) config[1] += DPI;

        config[2] = 0.0;	/* running time */ 
	config[3] = module2(x1[ntmin]-pos[0], y1[ntmin]-pos[1]);     /* distance to collision */
	config[4] = pos[0];    /* start position */
	config[5] = pos[1];
	config[6] = x1[ntmin];        /* position of collision */
	config[7] = y1[ntmin];
	
	return(c);
 }

 int valt_reuleaux(double config[8])
 /* determine initial configuration when starting from boundary */
 {
	double pos[2], alpha;
	int c;

        c = pos_alt_reuleaux(config, pos, &alpha);
        
        valt_reuleaux_xy(config, alpha, pos);
	
	return(c);
 }

/****************************************************************************************/
/* angle billiard */
/****************************************************************************************/
 
 
 int pos_angle(double conf[2], double pos[2], double *alpha)
 /* determine position on boundary of triangle */
 /* corner of angle is at (0,0), sides are horizontal and x*sin(a) = -y*cos(a), a=LAMBDA*PI */
 /* we use arclength for horizontal and vertical side, y for diagonal */
 {
	double s, theta;

	s = conf[0]; 
        theta = conf[1];

        if (s<=0)
        {
            pos[0] = s;
            pos[1] = 0.0;
            *alpha = theta;
            return(0);
        }
        else 
        {
            pos[0] = - s/tan(LAMBDA*PI);
            pos[1] = s;
            *alpha = theta + (1.0-LAMBDA)*PI;
            return(1);
        }
 }
 
 
 int vangle_xy(double config[8], double alpha, double pos[2])
 /* determine initial configuration for start at point pos = (x,y) */
 /* Warning: reflection in the corners is not yet implemented correctly */
 {
	double s0, c0, t, x, y, x1, y1, phi;
	int c, intb=1, intc, i;

	/* initial position and velocity */

	s0 = sin(alpha);
	c0 = cos(alpha);

	/* intersection with lower part of boundary */
	if (s0<0.0)
	{
            x1 = pos[0] - pos[1]*c0/s0;
            y1 = 0.0;
            if (x1<=0.0)
            {
                config[0] = x1;
                config[1] = -alpha;
                intb = 0;
            }
	}
	/* intersection with diagonal part of boundary */
	if ((intb)&&(sin(alpha+LAMBDA*PI)>0.0))
	{
            phi = LAMBDA*PI;
            t = -(pos[0]*sin(phi) + pos[1]*cos(phi))/sin(alpha+phi);
            x1 = pos[0] + t*c0;
            y1 = pos[1] + t*s0;
            if (y1 > 0.0)
            {
                config[0] = y1;
                config[1] = PI - alpha - phi;
                intb = 0;
            }
	}
	/* other cases, assign a dummy coordinate */
	if (intb)
        {
            if (s0 > 0.0)
            {
                t = (1000.0 - pos[0])/s0;
                x1 = pos[0] + t*c0;
                y1 = 1000.0;
            }
            else
            {
                t = -(1000.0 + pos[0])/c0;
                x1 = -1000.0;
                y1 = pos[1] + t*s0;
            }
            config[0] = -1000.0;
            config[1] = PID;
        }

	config[2] = 0.0;	/* running time */ 
	config[3] = module2(x1-pos[0], y1-pos[1]);
	config[4] = pos[0];
	config[5] = pos[1];
	config[6] = x1;
	config[7] = y1;
 
	return(c);
 }
 
 int vangle(double config[8])
 {
	double alpha, pos[2];
	int c;

	/* position et vitesse de depart */

	c = pos_angle(config, pos, &alpha);
        
        vangle_xy(config, alpha, pos);

        return(c);
 }
 
/****************************************************************************************/
/* L-shaped billiard (conical singularity) */
/****************************************************************************************/
 
 
 int pos_lshape(double conf[2], double pos[2], double *alpha)
 /* determine position on boundary of L-shaped billiard */
 /* the corners of the L shape are (-1,-1), (1,-1), (1,0), (0,0), (0,1), (-1,1) */
 /* returns the number of the side hit, or 0 if hitting a corner */
 {
	double s, theta;

	s = conf[0]; 
        if (s<0.0) s = 0.0;
        if (s>8.0) s = 8.0;
        
        theta = conf[1];


        /* bottom side */
        if ((s>0)&&(s<2.0))
        {
            pos[0] = s -1.0;
            pos[1] = -1.0;
            *alpha = theta;
            return(1);
        }
        /* lower right side */
        else if (s<3.0)
        {
            pos[0] = 1.0;
            pos[1] = s - 3.0;
            *alpha = theta + PID;
            return(2);
        }
        /* lower top side */
        else if (s<4.0)
        {
            pos[0] = 4.0 - s;
            pos[1] = 0.0;
            *alpha = theta + PI;
            return(3);
        }
        /* upper right side */
        else if (s<5.0)
        {
            pos[0] = 0.0;
            pos[1] = s - 4.0;
            *alpha = theta + PID;
            return(4);
        }
        /* upper top side */
        else if (s<6.0)
        {
            pos[0] = 5.0 - s;
            pos[1] = 1.0;
            *alpha = theta + PI;
            return(5);
        }
        /* left side */
        else
        {
            pos[0] = -1.0;
            pos[1] = 7.0 - s;
            *alpha = theta - PID;
            return(6);
        }

 }
 
 
 int vlshape_xy(double config[8], double alpha, double pos[2])
 /* determine initial configuration for start at point pos = (x,y) */
 {
    double t, l, s0, c0, margin = 1e-12, tmin;
    double x1[6], y1[6], tval[6], tempconf[6][2];
    int i, c, intb=1, cval[6], nt = 0, ntmin;

    /* initial position and velocity */

    s0 = sin(alpha);
    c0 = cos(alpha);
    
//     print_config(config);

    /* intersection with lower part of boundary */
    if (s0<0.0)
    {
        t = - (1.0 + pos[1])/s0;
        x1[nt] = pos[0] + c0*t;
        y1[nt] = -1.0;
        if ((x1[nt]>=-1.0)&&(x1[nt]<=0.0))
        {
            tempconf[nt][0] = 5.0 - x1[nt];
            tempconf[nt][1] = alpha + PI;
            cval[nt] = 5;
            tval[nt] = t;
            nt++;
        }
        else if ((x1[nt]>0.0)&&(x1[nt]<=1.0))
        {
            tempconf[nt][0] = 4.0 - x1[nt];
            tempconf[nt][1] = alpha + PI;
            cval[nt] = 3;
            tval[nt] = t;
            nt++;
        }
    }
    /* intersection with lower part of right-hand boundary */
    if (c0>0.0)
    {
        t = (1.0 - pos[0])/c0;
        x1[nt] = 1.0;
        y1[nt] = pos[1] + s0*t;
        if ((y1[nt]>=-1.0)&&(y1[nt]<=0.0))
        {
            tempconf[nt][0] = 7.0 - y1[nt];
            tempconf[nt][1] = PID + alpha;
            cval[nt] = 6;
            tval[nt] = t;
            nt++;
        }
    }
    /* intersection with upper part of right-hand boundary */
    if ((pos[0] < 0.0)&&(c0>0.0))
    {
        t = - pos[0]/c0;
        x1[nt] = 0.0;
        y1[nt] = pos[1] + s0*t;
        if ((y1[nt]>=0.0)&&(y1[nt]<=1.0))
        {
            tempconf[nt][0] = 7.0 - y1[nt];
            tempconf[nt][1] = PID + alpha;
            cval[nt] = 4;
            tval[nt] = t;
            nt++;
        }
    }
    /* intersection with right-hand part of upper boundary */
    if ((pos[1] < 0.0)&&(s0>0.0))
    {
        t = - pos[1]/s0;
        x1[nt] = pos[0] + c0*t;
        y1[nt] = 0.0;
        if ((x1[nt]>=-0.0)&&(x1[nt]<=1.0))
        {
            tempconf[nt][0] = 1.0 + x1[nt];
            tempconf[nt][1] = alpha;
            cval[nt] = 1;
            tval[nt] = t;
            nt++;
        }
    }
    /* intersection with left-hand part of upper boundary */
    if (s0>0.0)
    {
        t = (1.0 - pos[1])/s0;
        x1[nt] = pos[0] + c0*t;
        y1[nt] = 1.0;
        if ((x1[nt]>=-1.0)&&(x1[nt]<=0.0))
        {
            tempconf[nt][0] = 1.0 + x1[nt];
            tempconf[nt][1] = alpha;
            cval[nt] = 1;
            tval[nt] = t;
            nt++;
        }
    }
    /* intersection with left-hand part of boundary */
    if (c0<0.0)
    {
        t = (-1.0 - pos[0])/c0;
        x1[nt] = -1.0;
        y1[nt] = pos[1] + s0*t;
        if ((y1[nt]>=0.0)&&(y1[nt]<=1.0))
        {
            tempconf[nt][0] = 4.0 + y1[nt];
            tempconf[nt][1] = alpha - PID;
            cval[nt] = 4;
            tval[nt] = t;
            nt++;
        }
        else if ((y1[nt]>=-1.0)&&(y1[nt]<0.0))
        {
            tempconf[nt][0] = 3.0 + y1[nt];
            tempconf[nt][1] = alpha - PID;
            cval[nt] = 2;
            tval[nt] = t;
            nt++;
        }
    }
    
    /* find earliest intersection */
    tmin = tval[0];
    ntmin = 0;
    for (i=1; i<nt; i++) 
        if (tval[i] < tmin) 
        {
            tmin = tval[i];
            ntmin = i;
        }
        
//         printf("nt = %i, ntmin = %i, cmin = %i\n", nt, ntmin, cval[ntmin]);
            
    config[0] = tempconf[ntmin][0];
    config[1] = tempconf[ntmin][1];
    c = cval[ntmin];

    if (config[1] < 0.0) config[1] += DPI;
    
      
    config[2] = 0.0;	/* running time */ 
    config[3] = module2(x1[ntmin]-pos[0], y1[ntmin]-pos[1]);
    config[4] = pos[0];
    config[5] = pos[1];
    config[6] = x1[ntmin];
    config[7] = y1[ntmin];
    
//     print_config(config);
 
    return(c);	
 }
 
 int vlshape(double config[8])
 {
	double pos[2], alpha;
	int c;

	/* position et vitesse de depart */

	c = pos_lshape(config, pos, &alpha);
        
        vlshape_xy(config, alpha, pos);

        return(c);
 }
 
 /****************************************************************************************/
/* genus n surface - polygon with opposite sides identified */
/****************************************************************************************/

 int pos_genusn(double conf[2], double pos[2], double *alpha)
 /* determine position on boundary of polygon */
 /* conf[0] is arclength on boundary */
 {
	double s, theta, omega, length, s1, angle, x, y;
        int c;

	s = conf[0]; 
        theta = conf[1];
        
        omega = DPI/((double)NPOLY);
        length = 2.0*sin(0.5*omega);

        c = (int)(s/length);         /* side of polygon */
        
        s1 = s - ((double)c)*length;
        
        x = 1.0 + (cos(omega) - 1.0)*s1/length;
        y = sin(omega)*s1/length;
        
        angle = (double)c*omega + PID*APOLY;
        
        pos[0] = x*cos(angle) - y*sin(angle);
        pos[1] = x*sin(angle) + y*cos(angle);
        
        *alpha = (0.5 + (double)c)*omega + theta + PID*(1.0 + APOLY);
        
        return(c);
 }
 
 int vgenusn_xy(double config[8], double alpha, double pos[2])
 /* determine initial configuration for start at point pos = (x,y) */
 {
	double s, theta, omega, length, rlength, s1, rangle, x, y, xp, yp, x1, y1, ca, sa, margin = 1.0e-14;
	int k, c, intb=1, intc, i;

        /* dimensions/angles of polygon */
        omega = DPI/((double)NPOLY);
        length = 2.0*sin(0.5*omega);
        rlength = cos(0.5*omega);
        
        for (k=0; k<NPOLY; k++)
        {
            /* rotate position so that kth side is vertical */
            rangle = (0.5 + (double)k)*omega + APOLY*PID;
            theta = alpha - rangle;
                
            if ((cos(theta) > 0.0)&&(intb))
            {
                ca = cos(rangle);
                sa = sin(rangle);
                
                x = pos[0]*ca + pos[1]*sa;
                y = -pos[0]*sa + pos[1]*ca;
            
//                 xp = -rlength;
                xp = rlength;
                yp = y + (xp-x)*tan(theta);
                
                if (vabs(yp) < 0.5*length - margin) 
                {
                    /* rotate back */
                    x1 = xp*ca - yp*sa;
                    y1 = xp*sa + yp*ca;
                    
                    intb = 0;
                    c = k + NPOLY/2;
                    if (c > NPOLY) c -= NPOLY;
                    config[0] = ((double)c + 0.5)*length - yp;
//                     if (config[0] > (double)NPOLY*length) config[0] -= (double)NPOLY*length;
                    config[1] = PID + theta;
//                     config[0] = ((double)k + 0.5)*length + yp;
//                     config[1] = PID - theta;
                }
            }
        }
        if (intb) x1 = 1.0e12;  /* to make inactive particles too close to a corner */
       
        if (config[1] < 0.0) config[1] += DPI;

        config[2] = 0.0;	/* running time */ 
	config[3] = module2(x1-pos[0], y1-pos[1]);     /* distance to collision */
	config[4] = pos[0];    /* start position */
	config[5] = pos[1];
	config[6] = x1;        /* position of collision */
	config[7] = y1;
	
	return(c);
 }

 int vgenusn(double config[8])
 /* determine initial configuration when starting from boundary */
 {
	double pos[2], alpha;
	int c;

        c = pos_genusn(config, pos, &alpha);
        
        vgenusn_xy(config, alpha, pos);
	
	return(c);
 }
 
/****************************************************************************************/
/* Polygonal billiard with parabolic sides */
/****************************************************************************************/

 int pos_parabolas(double conf[2], double pos[2], double *alpha)
 /* determine position on boundary of polygon */
 /* conf[0] is y + ymax on first side */
 {
	double s, theta, omega2, aa, bb, cc, ymax, angle, x2, x, y;
        int c;

	s = conf[0]; 
        theta = conf[1];
        
        omega2 = PI/((double)NPOLY);
        aa = 0.25/MU;
        bb = 1.0/tan(omega2);
        cc = LAMBDA + MU;
        ymax = ( - bb + sqrt(bb*bb + 4.0*aa*cc))/(2.0*aa);
        
        c = (int)(s/(2.0*ymax));         /* side of shape */
        
        y = s - (2.0*(double)c + 1.0)*ymax;
        x = LAMBDA + MU - 0.25*y*y/MU;
        
//         printf("y = %.3lg\n", y);
        
        angle = 2.0*((double)c)*omega2 + PID*APOLY;
        
        pos[0] = x*cos(angle) - y*sin(angle);
        pos[1] = x*sin(angle) + y*cos(angle);
        
        *alpha = PID + theta + atan(0.5*y/MU) + angle;
        
//         *alpha = PID + theta + atan(0.5*y/MU) + 2.0*(double)c*omega2 + APOLY*PID;
        
//         printf("alpha = %.5lg\t", *alpha);
        
        return(c);
 }
 
 int vparabolas_xy(double config[8], double alpha, double pos[2])
 /* determine initial configuration for start at point pos = (x,y) */
 {
	double s, theta, omega2, beta, s1, rangle, x, y, x1, y1, xi, yi, t, x2;
        double ca, sa, a, b, margin = 1.0e-14, tmin, aa, bb, cc, ymax;
	int k, c, intb=1, intc, i, nt = 0, ntmin;

        /* dimensions/angles of polygon */
        omega2 = PI/((double)NPOLY);
        aa = 0.25/MU;
        bb = 1.0/tan(omega2);
        cc = LAMBDA + MU;
        ymax = ( - bb + sqrt(bb*bb + 4.0*aa*cc))/(2.0*aa);
//         printf("ymax = %.3lg\n", ymax);
        
//         print_config(config);

        for (k=0; k<NPOLY; k++) if (intb)
        {
            /* rotate position so that kth side is vertical */
            rangle = 2.0*((double)k)*omega2 + APOLY*PID;
            theta = alpha - rangle;
            
//             printf("theta = %.3lg\n", theta);
            
            ca = cos(rangle);
            sa = sin(rangle);
                                
            x = pos[0]*ca + pos[1]*sa;
            y = -pos[0]*sa + pos[1]*ca;
                
            aa = sin(theta)*sin(theta);
            bb = y*sin(theta) + 2.0*MU*cos(theta);
            cc = y*y + 4.0*MU*x - 4.0*MU*(LAMBDA+MU);
            
//             printf("y = %.3lg\n", y);
            
//             if (cos(theta) == 1.0)
            if (cos(theta) >= 1.0 - margin)
            {
                
//                 t = -cc/(2.0*bb);
                xi = LAMBDA + MU - 0.25*y*y/MU;
                yi = y;
                
//                 if ((t > margin)&&(vabs(yi) <= ymax)) 
                if (vabs(yi) <= ymax)
                {
                    intb = 0;
                    
                    /* rotate back */
                    x1 = xi*ca - yi*sa;
                    y1 = xi*sa + yi*ca;
                
                    config[0] = yi + (1.0 + 2.0*(double)k)*ymax;
                    config[1] = PID - theta + atan(0.5*yi/MU); 
                }
                
//                 printf("s = %.3lg\n", config[0]);
            }
            else if (bb*bb - aa*cc > margin)
            {
                t = (-bb + sqrt(bb*bb - aa*cc))/aa;
                
                xi = x + t*cos(theta);
                yi = y + t*sin(theta);
                    
                if ((t > margin)&&(vabs(yi) <= ymax)) 
//                 if (vabs(yi) <= ymax)
                {
                    intb = 0;
                    
//                     printf("t = %.3lg\n", t);
                        
                    /* rotate back */
                    x1 = xi*ca - yi*sa;
                    y1 = xi*sa + yi*ca;
                    
                    config[0] = yi + (1.0 + 2.0*(double)k)*ymax;
                    config[1] = PID - theta + atan(0.5*yi/MU); 
                }
            }
//             if (!intb) printf("Intersection found with side %i\n", k);
        }
        
//         if (intb) printf("No intersection found!\n");

        if (config[1] < 0.0) config[1] += DPI;

        config[2] = 0.0;	/* running time */ 
	config[3] = module2(x1-pos[0], y1-pos[1]);     /* distance to collision */
	config[4] = pos[0];    /* start position */
	config[5] = pos[1];
	config[6] = x1;        /* position of collision */
	config[7] = y1;
	
//         print_config(config);
        
	return(k);
 }

 int vparabolas(double config[8])
 /* determine initial configuration when starting from boundary */
 {
	double pos[2], alpha;
	int c;

        c = pos_parabolas(config, pos, &alpha);
        
        vparabolas_xy(config, alpha, pos);
	
	return(c);
 }

/****************************************************************************************/
/* Penrose solution to illumination problem */
/****************************************************************************************/

 int pos_penrose(double conf[2], double pos[2], double *alpha)
 /* determine position on boundary of domain */
 /* conf[0] parametrization of boundary by arclength or angle */
 {
	double s, s1, theta, cc, l1, l2, width, x, y, phi; 
        int c, i;
        static int first = 1;
        static double sval[17];

	s = conf[0]; 
        theta = conf[1];
        
        cc = sqrt(LAMBDA*LAMBDA - (1.0-MU)*(1.0-MU));
        width = 0.1*MU;
        l1 = MU - width;
        l2 = LAMBDA - cc;
        
        /* s values of different boundary parts */
        if (first)
        {
            sval[0] = 0.0;  sval[1] = PI;   sval[2] = sval[1] + l1;
            sval[3] = sval[2] + l2;     sval[4] = sval[3] + l1;
            sval[5] = sval[4] + PI;     sval[6] = sval[5] + l1;
            sval[7] = sval[6] + l2;     sval[8] = sval[7] + l1;
            for (i=1; i<=8; i++) sval[8+i] = sval[8] + sval[i];
//             for (i=0; i<16; i++) printf("sval[%i] = %.3lg\n", i, sval[i]);
            first = 0;
        }
        
        if (s < sval[1])     /* upper ellipse */
        {
            pos[0] = LAMBDA*cos(s);
            pos[1] = MU + (1.0-MU)*sin(s);
            phi = argument(-LAMBDA*sin(s),(1.0-MU)*cos(s));
            *alpha = phi + theta; 
            return(0);
        }
        else if (s < sval[2])       /* upper left straight parts */
        {
            pos[0] = -LAMBDA;
            pos[1] = MU - s + PI;
            *alpha = theta - PID;
            return(1);
        }
        else if (s < sval[3])
        {
            pos[0] = -LAMBDA + s - sval[2];
            pos[1] = width;
            *alpha = theta;
            return(2);
        }
        else if (s < sval[4])
        {
            pos[0] = -cc;
            pos[1] = width + s -sval[3];
            *alpha = theta + PID;
            return(3);
        }
        else if (s < sval[5])     /* left ellipse/mushroom head */
        {
            s1 = s - sval[4];
            pos[0] = -cc + 0.5*PENROSE_RATIO*MU*sin(s1);
            pos[1] = MU*cos(s1);
            phi = argument(0.5*PENROSE_RATIO*MU*cos(s1), -MU*sin(s1));
            *alpha = phi + theta;
            return(4);
        }
        else if (s < sval[6])     /* lower left straight parts */
        {
            s1 = s - sval[5];
            pos[0] = -cc;
            pos[1] = -MU + s1;
            *alpha = theta + PID;
            return(5);
        }
        else if (s < sval[7])
        {
            s1 = s - sval[6];
            pos[0] = -cc - s1;
            pos[1] = -width;
            *alpha = theta + PI;
            return(6);
        }
        else if (s < sval[8])
        {
            s1 = s - sval[7];
            pos[0] = -LAMBDA;
            pos[1] = -width - s1;
            *alpha = theta - PID;
            return(7);
        }
        
        else if (s < sval[9])     /* lower ellipse */
        {
            s1 = s -  sval[8];
            pos[0] = -LAMBDA*cos(s1);
            pos[1] = -MU - (1.0-MU)*sin(s1);
            phi = argument(LAMBDA*sin(s1),-(1.0-MU)*cos(s1));
            *alpha = phi + theta; 
            return(8);
        }
        else if (s < sval[10])       /* lower right straight parts */
        {
            s1 = s - sval[9];
            pos[0] = LAMBDA;
            pos[1] = -MU + s1;
            *alpha = theta + PID;
            return(9);
        }
        else if (s < sval[11])
        {
            s1 = s - sval[10];
            pos[0] = LAMBDA - s1;
            pos[1] = -width;
            *alpha = theta + PI;
            return(10);
        }
        else if (s < sval[12])
        {
            s1 = s - sval[11];
            pos[0] = cc;
            pos[1] = -width -s1;
            *alpha = theta - PID;
            return(11);
        }
        else if (s < sval[13])     /* right ellipse/mushroom head */
        {
            s1 = s - sval[12];
            pos[0] = cc - 0.5*PENROSE_RATIO*MU*sin(s1);
            pos[1] = -MU*cos(s1);
            phi = argument(-0.5*PENROSE_RATIO*MU*cos(s1), MU*sin(s1));
            *alpha = phi + theta;
            return(12);
        }
        else if (s < sval[14])     /* upper right straight parts */
        {
            s1 = s - sval[13];
            pos[0] = cc;
            pos[1] = MU - s1;
            *alpha = theta - PID;
            return(13);
        }
        else if (s < sval[15])
        {
            s1 = s - sval[14];
            pos[0] = cc + s1;
            pos[1] = width;
            *alpha = theta;
            return(14);
        }
        else 
        {
            s1 = s - sval[15];
            pos[0] = LAMBDA;
            pos[1] = width + s1;
            *alpha = theta + PID;
            return(15);
        }
 }
 
 int vpenrose_xy(double config[8], double alpha, double pos[2])
 /* determine initial configuration for start at point pos = (x,y) */
 {
	double s, theta, cc, width, l1, l2, s1, rangle, x, y, x1[30], y1[30], xi, yi, t, x2;
        double ca, sa, a, b, d, margin = 1.0e-14, tmin, tval[30], tempconf[30][2], lam2, mu2; 
	int k, c, intb=1, intc, i, nt = 0, cval[30], ntmin;
        static int first = 1;
        static double sval[17];

        /* dimensions of domain */
        cc = sqrt(LAMBDA*LAMBDA - (1.0-MU)*(1.0-MU));
        width = 0.1*MU;
        l1 = MU - width;
        l2 = LAMBDA - cc;
        lam2 = LAMBDA*LAMBDA;
        mu2 = (1.0-MU)*(1.0-MU);
        
        /* s values of different boundary parts */
        /* USE STATIC DOUBLES ? */
        if (first)
        {
            sval[0] = 0.0;  sval[1] = PI;   sval[2] = sval[1] + l1;
            sval[3] = sval[2] + l2;     sval[4] = sval[3] + l1;
            sval[5] = sval[4] + PI;     sval[6] = sval[5] + l1;
            sval[7] = sval[6] + l2;     sval[8] = sval[7] + l1;
            for (i=1; i<=8; i++) sval[8+i] = sval[8] + sval[i];
//             for (i=0; i<16; i++) printf("sval[%i] = %.3lg\n", i, sval[i]);
            first = 0;
        }
        
        ca = cos(alpha);
        sa = sin(alpha);
        
        /* intersection with upper ellipse */
        a = mu2*ca*ca + lam2*sa*sa;
        b = mu2*pos[0]*ca + lam2*(pos[1]-MU)*sa;
        d = mu2*pos[0]*pos[0] + lam2*(pos[1]-MU)*(pos[1]-MU) - lam2*mu2;
        
        if (b*b > a*d)
        {
            t = (-b + sqrt(b*b - a*d))/a;
            x = pos[0] + t*ca;
            y = pos[1] + t*sa;
            
//             printf("x = %.3lg, y = %.3lg\n", x, y);
            
            if ((t > margin)&&(y >= MU))
            {
                cval[nt] = 0;
                tval[nt] = t;
                x1[nt] = x;
                y1[nt] = y;
                tempconf[nt][0] = argument(x/LAMBDA, (y-MU)/(1.0-MU));
                tempconf[nt][1] = argument(-LAMBDA*(y-MU)/(1.0-MU), (1.0-MU)*x/LAMBDA) - alpha;      
//                 tempconf[nt][1] = argument(-LAMBDA*(y-MU)/(1.0-MU), x/MU) - alpha;      
                nt++;
            }
        }
        
        /* intersection with lower ellipse */
        b = mu2*pos[0]*ca + lam2*(pos[1]+MU)*sa;
        d = mu2*pos[0]*pos[0] + lam2*(pos[1]+MU)*(pos[1]+MU) - lam2*mu2;
        
        if (b*b > a*d)
        {
            t = (-b + sqrt(b*b - a*d))/a;
            x = pos[0] + t*ca;
            y = pos[1] + t*sa;
            
            if ((t > margin)&&(y <= -MU))
            {
                cval[nt] = 8;
                tval[nt] = t;
                x1[nt] = x;
                y1[nt] = y;
                tempconf[nt][0] = sval[8] + argument(-x/LAMBDA, -(y+MU)/(1.0-MU));
                tempconf[nt][1] = argument(-LAMBDA*(y+MU)/(1.0-MU), (1.0-MU)*x/LAMBDA) - alpha;      
                nt++;
            }
        }
        
        /* intersection with right ellipse */
        a = 4.0*ca*ca + sa*sa*PENROSE_RATIO*PENROSE_RATIO;
        b = 4.0*(pos[0] - cc)*ca + pos[1]*sa*PENROSE_RATIO*PENROSE_RATIO;
        d = 4.0*(pos[0] - cc)*(pos[0] - cc) + (pos[1]*pos[1] - MU*MU)*PENROSE_RATIO*PENROSE_RATIO;
        
        if (b*b > a*d)
        {
            t = (-b + sqrt(b*b - a*d))/a;
            x = pos[0] + t*ca;
            y = pos[1] + t*sa;
            
            if ((t > margin)&&(x <= cc))
            {
                cval[nt] = 12;
                tval[nt] = t;
                x1[nt] = x;
                y1[nt] = y;
                tempconf[nt][0] = sval[12] + argument(-y/MU, 2.0*(cc-x)/(MU*PENROSE_RATIO));
                tempconf[nt][1] = argument(0.5*PENROSE_RATIO*y, 2.0*(cc-x)/PENROSE_RATIO) - alpha;      
                nt++;
            }
            
            t = (-b - sqrt(b*b - a*d))/a;
            x = pos[0] + t*ca;
            y = pos[1] + t*sa;
            
            if ((t > margin)&&(x <= cc))
            {
                cval[nt] = 12;
                tval[nt] = t;
                x1[nt] = x;
                y1[nt] = y;
                tempconf[nt][0] = sval[12] + argument(-y/MU, 2.0*(cc-x)/(MU*PENROSE_RATIO));
                tempconf[nt][1] = argument(0.5*PENROSE_RATIO*y, 2.0*(cc-x)/PENROSE_RATIO) - alpha;      
                nt++;
            }
            
        }
        
        /* intersection with left ellipse */
        b = 4.0*(pos[0] + cc)*ca + pos[1]*sa*PENROSE_RATIO*PENROSE_RATIO;
        d = 4.0*(pos[0] + cc)*(pos[0] + cc) + (pos[1]*pos[1] - MU*MU)*PENROSE_RATIO*PENROSE_RATIO;
        
        if (b*b > a*d)
        {
            t = (-b + sqrt(b*b - a*d))/a;
            x = pos[0] + t*ca;
            y = pos[1] + t*sa;
            
            if ((t > margin)&&(x >= -cc))
            {
                cval[nt] = 4;
                tval[nt] = t;
                x1[nt] = x;
                y1[nt] = y;
                tempconf[nt][0] = sval[4] + argument(y/MU, 2.0*(cc+x)/(MU*PENROSE_RATIO));
                tempconf[nt][1] = argument(0.5*PENROSE_RATIO*y, -2.0*(cc+x)/PENROSE_RATIO) - alpha;      
                nt++;
            }
            
            t = (-b - sqrt(b*b - a*d))/a;
            x = pos[0] + t*ca;
            y = pos[1] + t*sa;
            
            if ((t > margin)&&(x >= -cc))
            {
                cval[nt] = 4;
                tval[nt] = t;
                x1[nt] = x;
                y1[nt] = y;
                tempconf[nt][0] = sval[4] + argument(y/MU, 2.0*(cc+x)/(MU*PENROSE_RATIO));
                tempconf[nt][1] = argument(0.5*PENROSE_RATIO*y, -2.0*(cc+x)/PENROSE_RATIO) - alpha;      
                nt++;
            }
            
        }
        
        /* rightmost vertical segments */
        if (ca > 0.0)
        {
            t = (LAMBDA - pos[0])/ca;
            x = LAMBDA;
            y = pos[1] + t*sa;
            
            if ((y >= width)&&(y < MU))
            {
                cval[nt] = 15;
                tval[nt] = t;
                x1[nt] = x;
                y1[nt] = y;
                tempconf[nt][0] = sval[15] + y - width;
                tempconf[nt][1] = PID - alpha;      
                nt++;
            }
            
            else if ((y <= -width)&&(y > -MU))
            {
                cval[nt] = 9;
                tval[nt] = t;
                x1[nt] = x;
                y1[nt] = y;
                tempconf[nt][0] = sval[9] + y + MU;
                tempconf[nt][1] = PID - alpha;      
                nt++;
            }
        }
        
        /* leftmost vertical segments */
        if (ca < 0.0)
        {
            t = -(LAMBDA + pos[0])/ca;
            x = -LAMBDA;
            y = pos[1] + t*sa;
            
            if ((y >= width)&&(y < MU))
            {
                cval[nt] = 1;
                tval[nt] = t;
                x1[nt] = x;
                y1[nt] = y;
                tempconf[nt][0] = sval[1] + MU - y;
                tempconf[nt][1] = 3.0*PID - alpha;      
                nt++;
            }
            
            else if ((y <= -width)&&(y > -MU))
            {
                cval[nt] = 7;
                tval[nt] = t;
                x1[nt] = x;
                y1[nt] = y;
                tempconf[nt][0] = sval[7] - width - y;
                tempconf[nt][1] = 3.0*PID - alpha;      
                nt++;
            }
        }

        /* vertical segments of right mushroom head */
        if (ca < 0.0)
        {
            t = (cc - pos[0])/ca;
            x = cc;
            y = pos[1] + t*sa;
            
            if ((t > 0.0)&&(y >= width)&&(y < MU))
            {
                cval[nt] = 13;
                tval[nt] = t;
                x1[nt] = x;
                y1[nt] = y;
                tempconf[nt][0] = sval[13] - y + MU;
                tempconf[nt][1] = 3.0*PID - alpha;      
                nt++;
            }
            
            else if ((t > 0.0)&&(y <= -width)&&(y > -MU))
            {
                cval[nt] = 11;
                tval[nt] = t;
                x1[nt] = x;
                y1[nt] = y;
                tempconf[nt][0] = sval[11] - y - width;
                tempconf[nt][1] = 3.0*PID - alpha;      
                nt++;
            }
        }
        
        /* vertical segments of left mushroom head */
        if (ca > 0.0)
        {
            t = (-cc - pos[0])/ca;
            x = -cc;
            y = pos[1] + t*sa;
            
            if ((t > 0.0)&&(y >= width)&&(y < MU))
            {
                cval[nt] = 3;
                tval[nt] = t;
                x1[nt] = x;
                y1[nt] = y;
                tempconf[nt][0] = sval[3] + y - width;
                tempconf[nt][1] = PID - alpha;      
                nt++;
            }
            
            else if ((t > 0.0)&&(y <= -width)&&(y > -MU))
            {
                cval[nt] = 5;
                tval[nt] = t;
                x1[nt] = x;
                y1[nt] = y;
                tempconf[nt][0] = sval[5] + y + MU;
                tempconf[nt][1] = PID - alpha;      
                nt++;
            }
        }
        
        /* upper horizontal segments */
        if (sa < 0.0)
        {
            t = (width - pos[1])/sa;
            x = pos[0] + t*ca;
            y = width;
            
            if ((t > 0.0)&&(x >= cc)&&(x < LAMBDA))
            {
                cval[nt] = 14;
                tval[nt] = t;
                x1[nt] = x;
                y1[nt] = y;
                tempconf[nt][0] = sval[14] + x - cc;
                tempconf[nt][1] = - alpha;      
                nt++;
            }
            
            else if ((t > 0.0)&&(x <= -cc)&&(x > -LAMBDA))
            {
                cval[nt] = 2;
                tval[nt] = t;
                x1[nt] = x;
                y1[nt] = y;
                tempconf[nt][0] = sval[2] + x + LAMBDA;
                tempconf[nt][1] = - alpha;      
                nt++;
            }
        }

        /* lower horizontal segments - TO BE CORRECTED */
        if (sa > 0.0)
        {
            t = (-width - pos[1])/sa;
            x = pos[0] + t*ca;
            y = -width;
            
            if ((t > 0.0)&&(x >= cc)&&(x < LAMBDA))
            {
                cval[nt] = 10;
                tval[nt] = t;
                x1[nt] = x;
                y1[nt] = y;
                tempconf[nt][0] = sval[10] - x + LAMBDA;
                tempconf[nt][1] = PI - alpha;      
                nt++;
            }
            
            else if ((t > 0.0)&&(x <= -cc)&&(x > -LAMBDA))
            {
                cval[nt] = 6;
                tval[nt] = t;
                x1[nt] = x;
                y1[nt] = y;
                tempconf[nt][0] = sval[6] - x - cc;
                tempconf[nt][1] = PI - alpha;      
                nt++;
            }
        }
        /* find earliest intersection */
        tmin = tval[0];
        ntmin = 0;
        for (i=1; i<nt; i++) 
            if (tval[i] < tmin) 
            {
                tmin = tval[i];
                ntmin = i;
            }
            
        config[0] = tempconf[ntmin][0];
        config[1] = tempconf[ntmin][1];
        c = cval[ntmin];
 
//         printf("nt = %i\t ntmin = %i \tcmin = %i\n", nt, ntmin, c);
        

        if (config[1] < 0.0) config[1] += DPI;

        config[2] = 0.0;	/* running time */ 
	config[3] = module2(x1[ntmin]-pos[0], y1[ntmin]-pos[1]);     /* distance to collision */
	config[4] = pos[0];    /* start position */
	config[5] = pos[1];
	config[6] = x1[ntmin];        /* position of collision */
	config[7] = y1[ntmin];
        
//         print_config(config);
	
	return(c);
 }

 int vpenrose(double config[8])
 /* determine initial configuration when starting from boundary */
 {
	double pos[2], alpha;
	int c;

        c = pos_penrose(config, pos, &alpha);
        
        vpenrose_xy(config, alpha, pos);
	
	return(c);
 }


/****************************************************************************************/
/* billiard with circular scatterers */
/****************************************************************************************/

 int pos_circles(double conf[2], double pos[2], double *alpha)
 /* determine position on boundary of circle */
 /* position varies between 0 and ncircles*2Pi */
 /* returns number of hit circle */
 {
	double angle;
        int ncirc;
        
        ncirc = (int)(conf[0]/DPI);
        if (ncirc >= ncircles) ncirc = ncircles - 1;
        
        angle = conf[0] - (double)ncirc*DPI;

        pos[0] = circles[ncirc].xc + circles[ncirc].radius*cos(angle);
        pos[1] = circles[ncirc].yc + circles[ncirc].radius*sin(angle);
        
        *alpha = angle + PID + conf[1]; 
        
        return(ncirc);
 }
 
 
 int vcircles_xy(double config[8], double alpha, double pos[2])
 /* determine initial configuration for start at point pos = (x,y) */
 {
	double c0, s0, b, c, t, theta, delta, margin = 1.0e-12, tmin, rlarge = 1000.0;
        double tval[ncircles], xint[ncircles], yint[ncircles], phiint[ncircles];
	int i, nt = 0, nscat[ncircles], ntmin;

        c0 = cos(alpha);
        s0 = sin(alpha);
        
        for (i=0; i<ncircles; i++) if (circles[i].active)
        {
            b = (pos[0]-circles[i].xc)*c0 + (pos[1]-circles[i].yc)*s0;
            c = (pos[0]-circles[i].xc)*(pos[0]-circles[i].xc) + (pos[1]-circles[i].yc)*(pos[1]-circles[i].yc) -            circles[i].radius*circles[i].radius;
        
            delta = b*b - c;
            if (delta > margin)     /* there is an intersection with circle i */
            {
                t = -b - sqrt(delta);            
                if (t > margin) 
                {
                    nscat[nt] = i;
                
                    tval[nt] = t;
                    xint[nt] = pos[0] + t*c0;
                    yint[nt] = pos[1] + t*s0;
                    phiint[nt] = argument(xint[nt] - circles[i].xc, yint[nt] - circles[i].yc);

                    nt++;
                }
            }
        }
        
        if (nt > 0)     /* there is at least one intersection */
        {
            /* find earliest intersection */
            tmin = tval[0];
            ntmin = 0;
            for (i=1; i<nt; i++) 
            if (tval[i] < tmin) 
            {
                tmin = tval[i];
                ntmin = i;
            }
            while (phiint[ntmin] < 0.0) phiint[ntmin] += DPI;
            while (phiint[ntmin] >= DPI) phiint[ntmin] -= DPI;
            
            config[0] = (double)nscat[ntmin]*DPI + phiint[ntmin];
            config[1] = PID - alpha + phiint[ntmin];        /* CHECK */
            if (config[1] < 0.0) config[1] += DPI;
            if (config[1] >= PI) config[1] -= DPI;
            
            config[2] = 0.0;	/* running time */ 
            config[3] = module2(xint[ntmin]-pos[0], yint[ntmin]-pos[1]);     /* distance to collision */
            config[4] = pos[0];    /* start position */
            config[5] = pos[1];
            config[6] = xint[ntmin];        /* position of collision */
            config[7] = yint[ntmin];
            
            /* set dummy coordinates if circles are absorbing */
            if (ABSORBING_CIRCLES)
            {
                config[0] = DUMMY_ABSORBING;
                config[1] = PI;
            }

            return(nscat[ntmin]);
        }
        else    /* there is no intersection - set dummy values */
        {
            config[0] = DUMMY_ABSORBING;
            config[1] = PI;
            config[2] = 0.0;
            config[3] = rlarge;
            config[4] = pos[0];    /* start position */
            config[5] = pos[1];
            config[6] = rlarge*cos(alpha);
            config[7] = rlarge*sin(alpha);
            
            return(-1);
        }
 }

 int vcircles(double config[8])
 /* determine initial configuration when starting from boundary */
  {
	double pos[2], alpha;
	int c;

        c = pos_circles(config, pos, &alpha);
        
        c = vcircles_xy(config, alpha, pos);
	
	return(c);
 }
 
/****************************************************************************************/
/* billiard with circular scatterers in a rectangle */
/****************************************************************************************/

 int pos_circles_in_rect(double conf[2], double pos[2], double *alpha)
 /* determine position on boundary of circle */
 /* position varies between 0 and ncircles*2Pi for circles and between -BOUNDARY_SHIFT and 0 for boundary*/
 /* returns number of hit circle */
 {
	double angle;
        int ncirc, c;
        
        if (conf[0] >= 0)
        {
            ncirc = (int)(conf[0]/DPI);
            if (ncirc >= ncircles) ncirc = ncircles - 1;
        
            angle = conf[0] - (double)ncirc*DPI;

            pos[0] = circles[ncirc].xc + circles[ncirc].radius*cos(angle);
            pos[1] = circles[ncirc].yc + circles[ncirc].radius*sin(angle);
        
            *alpha = angle + PID + conf[1]; 
        
            return(ncirc);
        }
        else /* particle starts on boundary */
        {
//             conf[0] += 4.0*(LAMBDA + 1.0);
            conf[0] += BOUNDARY_SHIFT;
            c = pos_rectangle(conf, pos, alpha);
            
//             conf[0] -= 4.0*(LAMBDA + 1.0);
            conf[0] -= BOUNDARY_SHIFT;
            
            return(-c-1);
        }
 }
 
 
 int vcircles_in_rect_xy(double config[8], double alpha, double pos[2])
 /* determine initial configuration for start at point pos = (x,y) */
//  double config[8], alpha, pos[2];

 {
	double c0, s0, b, c, t, theta, delta, margin = 1.0e-12, tmin, rlarge = 1.0e10;
        double tval[ncircles], xint[ncircles], yint[ncircles], phiint[ncircles];
	int i, nt = 0, nscat[ncircles], ntmin, side;

        c0 = cos(alpha);
        s0 = sin(alpha);
        
        for (i=0; i<ncircles; i++) if (circles[i].active)
        {
            b = (pos[0]-circles[i].xc)*c0 + (pos[1]-circles[i].yc)*s0;
            c = (pos[0]-circles[i].xc)*(pos[0]-circles[i].xc) + (pos[1]-circles[i].yc)*(pos[1]-circles[i].yc) - circles[i].radius*circles[i].radius;
        
            delta = b*b - c;
            if (delta > margin)     /* there is an intersection with circle i */
            {
                t = -b - sqrt(delta);            
                if (t > margin) 
                {
                    nscat[nt] = i;
                
                    tval[nt] = t;
                    xint[nt] = pos[0] + t*c0;
                    yint[nt] = pos[1] + t*s0;
                    phiint[nt] = argument(xint[nt] - circles[i].xc, yint[nt] - circles[i].yc);

                    /* test wether intersection is in rectangle */
                    if ((vabs(xint[nt]) < LAMBDA)&&(vabs(yint[nt]) < 1.0)) nt++;
                }
            }
        }
        
        if (nt > 0)     /* there is at least one intersection */
        {
            /* find earliest intersection */
            tmin = tval[0];
            ntmin = 0;
            for (i=1; i<nt; i++) 
            if (tval[i] < tmin) 
            {
                tmin = tval[i];
                ntmin = i;
            }
            while (phiint[ntmin] < 0.0) phiint[ntmin] += DPI;
            while (phiint[ntmin] >= DPI) phiint[ntmin] -= DPI;
            
            config[0] = (double)nscat[ntmin]*DPI + phiint[ntmin];
            config[1] = PID - alpha + phiint[ntmin];        /* CHECK */
            if (config[1] < 0.0) config[1] += DPI;
            if (config[1] >= PI) config[1] -= DPI;
            
            config[2] = 0.0;	/* running time */ 
            config[3] = module2(xint[ntmin]-pos[0], yint[ntmin]-pos[1]);     /* distance to collision */
            config[4] = pos[0];    /* start position */
            config[5] = pos[1];
            config[6] = xint[ntmin];        /* position of collision */
            config[7] = yint[ntmin];
            
            
            /* set dummy coordinates if circles are absorbing */
            if (ABSORBING_CIRCLES)
            {
                config[0] = DUMMY_ABSORBING;
                config[1] = PI;
//                 return(DUMMY_SIDE_ABS);
            }
            
            if (ABSORBING_CIRCLES) return(DUMMY_SIDE_ABS);
            else return(nscat[ntmin]);
        }
        else    /* there is no intersection with the circles - compute intersection with boundary */
        {
            
            side = vrectangle_xy(config, alpha, pos);
            config[0] -= BOUNDARY_SHIFT;
//             config[0] -= 4.0*(LAMBDA+1.0);
            
//             printf("Hit side %i\n", side);
//             print_config(config);
            return(side - 5);
        }
 }

 int vcircles_in_rect(double config[8])
 /* determine initial configuration when starting from boundary */
  {
	double pos[2], alpha;
	int c;

        c = pos_circles_in_rect(config, pos, &alpha);
        
//         vcircles_in_rect_xy(config, alpha, pos);
        c = vcircles_in_rect_xy(config, alpha, pos);
	
	return(c);
 }
 

/****************************************************************************************/
/* billiard with circular scatterers in a genus n surface (polygon with identifies opposite sides) */
/****************************************************************************************/

 int pos_circles_in_genusn(double conf[2], double pos[2], double *alpha)
 /* determine position on boundary of circle */
 /* position varies between 0 and ncircles*2Pi for circles and between */
 /*  BOUNDARY_SHIFT and 0 for boundary*/
 /* returns number of hit circle */
 {
	double s, theta, omega, length, s1, angle, x, y;
        int ncirc, c;
        
        if (conf[0] >= 0)
        {
            ncirc = (int)(conf[0]/DPI);
            if (ncirc >= ncircles) ncirc = ncircles - 1;
        
            angle = conf[0] - (double)ncirc*DPI;

            pos[0] = circles[ncirc].xc + circles[ncirc].radius*cos(angle);
            pos[1] = circles[ncirc].yc + circles[ncirc].radius*sin(angle);
        
            *alpha = angle + PID + conf[1]; 
        
            return(ncirc);
        }
        else /* particle starts on boundary */
        {
            omega = DPI/((double)NPOLY);
            length = 2.0*sin(0.5*omega);
        
//             conf[0] += 2.0*length*(double)NPOLY;
            conf[0] += BOUNDARY_SHIFT;
            c = pos_genusn(conf, pos, alpha);
            
//             conf[0] -= 2.0*length*(double)NPOLY;
            conf[0] -= BOUNDARY_SHIFT;
            
            return(-c);
        }
 }
 
 
 int vcircles_in_genusn_xy(double config[8], double alpha, double pos[2])
 /* determine initial configuration for start at point pos = (x,y) */
 {
	double c0, s0, b, c, t, theta, delta, margin = 1.0e-12, tmin, rlarge = 1.0e10, omega, length, angle, cw;
        double tval[ncircles], xint[ncircles], yint[ncircles], phiint[ncircles];
	int i, k, nt = 0, nscat[ncircles], ntmin, side, condition;

        c0 = cos(alpha);
        s0 = sin(alpha);
        omega = DPI/((double)NPOLY);
        length = 2.0*sin(0.5*omega);
        cw = cos(omega*0.5);
                    
        for (i=0; i<ncircles; i++) if (circles[i].active)
        {
            b = (pos[0]-circles[i].xc)*c0 + (pos[1]-circles[i].yc)*s0;
            c = (pos[0]-circles[i].xc)*(pos[0]-circles[i].xc) + (pos[1]-circles[i].yc)*(pos[1]-circles[i].yc) - circles[i].radius*circles[i].radius;
        
            delta = b*b - c;
            if (delta > margin)     /* there is an intersection with circle i */
            {
                t = -b - sqrt(delta);            
                if (t > margin) 
                {
                    nscat[nt] = i;
                
                    tval[nt] = t;
                    xint[nt] = pos[0] + t*c0;
                    yint[nt] = pos[1] + t*s0;
                    phiint[nt] = argument(xint[nt] - circles[i].xc, yint[nt] - circles[i].yc);

                    /* test wether intersection is in polygon */
                    if (in_polygon(xint[nt], yint[nt], 1.0, NPOLY, APOLY)) nt++;
                }
            }
        }
        
        if (nt > 0)     /* there is at least one intersection */
        {
            /* find earliest intersection */
            tmin = tval[0];
            ntmin = 0;
            for (i=1; i<nt; i++) 
            if (tval[i] < tmin) 
            {
                tmin = tval[i];
                ntmin = i;
            }
            while (phiint[ntmin] < 0.0) phiint[ntmin] += DPI;
            while (phiint[ntmin] >= DPI) phiint[ntmin] -= DPI;
            
            config[0] = (double)nscat[ntmin]*DPI + phiint[ntmin];
            config[1] = PID - alpha + phiint[ntmin];        /* CHECK */
            if (config[1] < 0.0) config[1] += DPI;
            if (config[1] >= PI) config[1] -= DPI;
            
            config[2] = 0.0;	/* running time */ 
            config[3] = module2(xint[ntmin]-pos[0], yint[ntmin]-pos[1]);     /* distance to collision */
            config[4] = pos[0];    /* start position */
            config[5] = pos[1];
            config[6] = xint[ntmin];        /* position of collision */
            config[7] = yint[ntmin];
            
            
            /* set dummy coordinates if circles are absorbing */
            if (ABSORBING_CIRCLES)
            {
                config[0] = DUMMY_ABSORBING;
                config[1] = PI;
            }
            
            return(nscat[ntmin]);
        }
        else    /* there is no intersection with the circles - compute intersection with boundary */
        {
            side = vgenusn_xy(config, alpha, pos);
            
//             config[0] -= 2.0*length*(double)NPOLY;
            config[0] -= BOUNDARY_SHIFT;

            return(side);
        }
 }

 int vcircles_in_genusn(double config[8])
 /* determine initial configuration when starting from boundary */
  {
	double pos[2], alpha;
	int c;

        c = pos_circles_in_genusn(config, pos, &alpha);
        
        vcircles_in_genusn_xy(config, alpha, pos);
	
	return(c);
 }
 

/****************************************************************************************/
/* billiard with circular scatterers in a torus (rectangle with periodic boundary conditions) */
/****************************************************************************************/

 int pos_circles_in_torus(double conf[2], double pos[2], double *alpha)
 /* determine position on boundary of circle */
 /* position varies between 0 and ncircles*2Pi for circles and between -BOUNDARY_SHIFT and 0 for boundary*/
 /* returns number of hit circle */
 {
	double angle;
        int ncirc, c;
        
        if (conf[0] >= 0)
        {
            ncirc = (int)(conf[0]/DPI);
            if (ncirc >= ncircles) ncirc = ncircles - 1;
        
            angle = conf[0] - (double)ncirc*DPI;

            pos[0] = circles[ncirc].xc + circles[ncirc].radius*cos(angle);
            pos[1] = circles[ncirc].yc + circles[ncirc].radius*sin(angle);
        
            *alpha = angle + PID + conf[1]; 
        
            return(ncirc);
        }
        else /* particle starts on boundary */
        {
//             conf[0] += 4.0*(LAMBDA + 1.0);
            conf[0] += BOUNDARY_SHIFT;
            c = pos_rectangle(conf, pos, alpha);
            
//             conf[0] -= 4.0*(LAMBDA + 1.0);
            conf[0] -= BOUNDARY_SHIFT;
            
            return(-c-1);
        }
 }
 
 
 int vcircles_in_torus_xy(double config[8], double alpha, double pos[2])
 /* determine initial configuration for start at point pos = (x,y) */
//  double config[8], alpha, pos[2];

 {
	double c0, s0, b, c, t, theta, delta, margin = 1.0e-12, tmin, rlarge = 1.0e10;
        double tval[ncircles], xint[ncircles], yint[ncircles], phiint[ncircles];
	int i, nt = 0, nscat[ncircles], ntmin, side;

        c0 = cos(alpha);
        s0 = sin(alpha);
        
        for (i=0; i<ncircles; i++) if (circles[i].active)
        {
            b = (pos[0]-circles[i].xc)*c0 + (pos[1]-circles[i].yc)*s0;
            c = (pos[0]-circles[i].xc)*(pos[0]-circles[i].xc) + (pos[1]-circles[i].yc)*(pos[1]-circles[i].yc) - circles[i].radius*circles[i].radius;
        
            delta = b*b - c;
            if (delta > margin)     /* there is an intersection with circle i */
            {
                t = -b - sqrt(delta);            
                if (t > margin) 
                {
                    nscat[nt] = i;
                
                    tval[nt] = t;
                    xint[nt] = pos[0] + t*c0;
                    yint[nt] = pos[1] + t*s0;
                    phiint[nt] = argument(xint[nt] - circles[i].xc, yint[nt] - circles[i].yc);

                    /* test wether intersection is in rectangle */
                    if ((vabs(xint[nt]) < LAMBDA)&&(vabs(yint[nt]) < 1.0)) nt++;
                }
            }
        }
        
        if (nt > 0)     /* there is at least one intersection */
        {
            /* find earliest intersection */
            tmin = tval[0];
            ntmin = 0;
            for (i=1; i<nt; i++) 
            if (tval[i] < tmin) 
            {
                tmin = tval[i];
                ntmin = i;
            }
            while (phiint[ntmin] < 0.0) phiint[ntmin] += DPI;
            while (phiint[ntmin] >= DPI) phiint[ntmin] -= DPI;
            
            config[0] = (double)nscat[ntmin]*DPI + phiint[ntmin];
            config[1] = PID - alpha + phiint[ntmin];        /* CHECK */
            if (config[1] < 0.0) config[1] += DPI;
            if (config[1] >= PI) config[1] -= DPI;
            
            config[2] = 0.0;	/* running time */ 
            config[3] = module2(xint[ntmin]-pos[0], yint[ntmin]-pos[1]);     /* distance to collision */
            config[4] = pos[0];    /* start position */
            config[5] = pos[1];
            config[6] = xint[ntmin];        /* position of collision */
            config[7] = yint[ntmin];
            
            
            /* set dummy coordinates if circles are absorbing */
            if (ABSORBING_CIRCLES)
            {
                config[0] = DUMMY_ABSORBING;
                config[1] = PI;
            }
            
            return(nscat[ntmin]);
        }
        else    /* there is no intersection with the circles - compute intersection with boundary */
        {
            
            side = vrectangle_xy(config, alpha, pos);
            
            if (config[0] < 2.0*LAMBDA) config[0] = 4.0*LAMBDA + 2.0 - config[0];
            else if (config[0] < 2.0*LAMBDA + 2.0) config[0] = 6.0*LAMBDA + 4.0 - config[0];
            else if (config[0] < 4.0*LAMBDA + 2.0) config[0] = 4.0*LAMBDA + 2.0 - config[0];
            else config[0] = 6.0*LAMBDA + 4.0 - config[0];

            config[0] -= BOUNDARY_SHIFT;
            config[1] = PI - config[1]; 
//             config[0] -= 4.0*(LAMBDA+1.0);
            
//             printf("Hit side %i\n", side);
//             print_config(config);
            return(side - 5);
        }
 }

 int vcircles_in_torus(double config[8])
 /* determine initial configuration when starting from boundary */
  {
	double pos[2], alpha;
	int c;

        c = pos_circles_in_torus(config, pos, &alpha);
        
        vcircles_in_torus_xy(config, alpha, pos);
	
	return(c);
 }
 
 
/****************************************************************************************/
/* billiard in poylgonal line */
/****************************************************************************************/

 int pos_polyline(double conf[2], double pos[2], double *alpha)
 /* determine position on boundary of domain */
 /* position varies between 0 and nsides */
 /* returns number of hit setment */
 {
	double len, rlarge = 1.0e8;
        int nside;
        
        nside = (int)conf[0];
        if (nside >= nsides) nside = nsides - 1;
        
        if (nside >= 0)      /* position is on a side of the polygonal line */
        {
            len = conf[0] - (double)nside;

            pos[0] = polyline[nside].x1 + len*(polyline[nside].x2 - polyline[nside].x1);
            pos[1] = polyline[nside].y1 + len*(polyline[nside].y2 - polyline[nside].y1);
        
            *alpha = polyline[nside].angle + conf[1]; 
            return(nside);
        }
        else    /* position is on an absorbing circle */
        {
            pos[0] = rlarge;
            pos[1] = 0.0;
            
            *alpha = 0.0;
            return(-1);    
        }
 }
 
 
 int vpolyline_xy(double config[8], double alpha, double pos[2])
 /* determine initial configuration for start at point pos = (x,y) */
 {
	double c0, s0, a, b, c, t, dx, delta, s, xi, yi, margin = 1.0e-12, tmin, rlarge = 1.0e8;
        double tval[nsides + ncircles], xint[nsides + ncircles], yint[nsides + ncircles], sint[nsides + ncircles];
	int i, nt = 0, nsegment[nsides + ncircles], ntmin;

        c0 = cos(alpha);
        s0 = sin(alpha);
        
        for (i=0; i<nsides; i++) 
        {
//             printf("testing side %i\n", i);
            
            a = polyline[i].y2 - polyline[i].y1;
            b = - polyline[i].x2 + polyline[i].x1;
            c = -a*polyline[i].x1 - b*polyline[i].y1;
            
//             printf("a = %.2f, b = %.2f, c = %.2f\n", a, b, c);
            
            delta = a*c0 + b*s0;
            if (vabs(delta) > margin)     /* there is an intersection with the line containing segment i */
            {
                t = -(a*pos[0] + b*pos[1] + c)/delta; 
                if (t > margin) 
                {
                    xi = pos[0] + t*c0;
                    yi = pos[1] + t*s0;
                    dx = polyline[i].x2 - polyline[i].x1;
                    
                    if (vabs(dx) > margin) s = (xi - polyline[i].x1)/dx; 
                    else s = (yi - polyline[i].y1)/(polyline[i].y2 - polyline[i].y1); 
//                     printf("s = %.2f\n", s);
                    
                    if ((s >= 0.0)&&(s <= 1.0))     
                    /* the intersection is on the segment */
                    {
                        nsegment[nt] = i;
                        tval[nt] = t;
                        sint[nt] = s;
                        xint[nt] = pos[0] + t*c0;
                        yint[nt] = pos[1] + t*s0;
//                         printf("s = %.2f, x = %.2f, y = %.2f\n", s, xint[nt], yint[nt]);
                        nt++;                        
                    }
                }
            }
        }
        
        if (ABSORBING_CIRCLES) for (i=0; i<ncircles; i++) 
        {
            b = (pos[0]-circles[i].xc)*c0 + (pos[1]-circles[i].yc)*s0;
            c = (pos[0]-circles[i].xc)*(pos[0]-circles[i].xc) + (pos[1]-circles[i].yc)*(pos[1]-circles[i].yc) - circles[i].radius*circles[i].radius;
        
            delta = b*b - c;
            if (delta > margin)     /* there is an intersection with circle i */
            {
                t = -b - sqrt(delta);            
                if (t > margin) 
                {
                    nsegment[nt] = -1-i;
                
                    tval[nt] = t;
                    xint[nt] = pos[0] + t*c0;
                    yint[nt] = pos[1] + t*s0;
                    sint[nt] = argument(xint[nt] - circles[i].xc, yint[nt] - circles[i].yc);

                    nt++;
                }
            }
        }
        
        if (nt > 0)     /* there is at least one intersection */
        {
            /* find earliest intersection */
            tmin = tval[0];
            ntmin = 0;
            for (i=1; i<nt; i++) 
            if (tval[i] < tmin) 
            {
                tmin = tval[i];
                ntmin = i;
            }
            
//             printf("ntmin = %i\n", ntmin); 
            if (nsegment[ntmin] >= 0) 
            {
                config[0] = (double)nsegment[ntmin] + sint[ntmin];
                config[1] = polyline[nsegment[ntmin]].angle - alpha;        
                if (config[1] < 0.0) config[1] += DPI;
                if (config[1] >= PI) config[1] -= DPI;
            }
            /* set dummy coordinates if circles are absorbing */
            else if ((ABSORBING_CIRCLES)&&(nsegment[ntmin] < 0))
            {
                config[0] = DUMMY_ABSORBING;
                config[1] = PI;
            }
            config[2] = 0.0;	/* running time */ 
            config[3] = module2(xint[ntmin]-pos[0], yint[ntmin]-pos[1]);     /* distance to collision */
            config[4] = pos[0];    /* start position */
            config[5] = pos[1];
            config[6] = xint[ntmin];        /* position of collision */
            config[7] = yint[ntmin];
            
            
//             print_config(config);
            
            return(nsegment[ntmin]);
        }
        else    /* there is no intersection - set dummy values */
        {
            config[0] = DUMMY_ABSORBING;
            config[1] = PI;
            config[2] = 0.0;
            config[3] = rlarge;
            config[4] = pos[0];    /* start position */
            config[5] = pos[1];
            config[6] = rlarge*cos(alpha);
            config[7] = rlarge*sin(alpha);
            
            return(-1);
        }
 }

 int vpolyline(double config[8])
 /* determine initial configuration when starting from boundary */
  {
	double pos[2], alpha;
	int c;

        c = pos_polyline(config, pos, &alpha);
        
        vpolyline_xy(config, alpha, pos);
//         c = vpolyline_xy(config, alpha, pos);
	
	return(c);
 }
 
/****************************************************************************************/
/* billiard in poylgonal line with additional circular arcs */
/****************************************************************************************/

 int pos_polyline_arc(double conf[2], double pos[2], double *alpha)
 /* determine position on boundary of domain */
 /* position varies between 0 and nsides */
 /* returns number of hit setment */
 {
	double len, angle, rlarge = 1.0e8;
        int nside, narc;
        
        nside = (int)conf[0];
        
        if (nside < BOUNDARY_SHIFT) 
        {
            return (pos_polyline(conf, pos, alpha));    /* position is on a line segment */
        }
        else    /* position is on a circular arc */
        {
            narc = nside - BOUNDARY_SHIFT;
            if (narc > narcs) return(0);
            
            len = conf[0] - (double)nside;
            angle = arcs[narc].angle1 + len*(arcs[narc].dangle);
            
            pos[0] = arcs[narc].xc + arcs[narc].radius*cos(angle);
            pos[1] = arcs[narc].yc + arcs[narc].radius*sin(angle);
            
            *alpha = angle + PID + conf[1];
            
            return(narc + BOUNDARY_SHIFT);
        }
 }
 
 
 int vpolyline_arc_xy(double config[8], double alpha, double pos[2])
 /* determine initial configuration for start at point pos = (x,y) */
 {
	double c0, s0, x1, y1, a, b, c, t, dx, delta, s, xi, yi, angle, margin = 1.0e-12, tmin, rlarge = 1.0e8;
        double tval[nsides + ncircles], xint[nsides + ncircles], yint[nsides + ncircles], sint[nsides + ncircles], 
        aint[nsides + ncircles];
	int i, nt = 0, nsegment[nsides + ncircles], narc[nsides + ncircles], ntmin, sign, type[nsides + ncircles];

        c0 = cos(alpha);
        s0 = sin(alpha);
        
        /* look for intersections with line segments */
        for (i=0; i<nsides; i++) 
        {
//             printf("testing side %i\n", i);
            
            a = polyline[i].y2 - polyline[i].y1;
            b = - polyline[i].x2 + polyline[i].x1;
            c = -a*polyline[i].x1 - b*polyline[i].y1;
            
//             printf("a = %.2f, b = %.2f, c = %.2f\n", a, b, c);
            
            delta = a*c0 + b*s0;
            if (vabs(delta) > margin)     /* there is an intersection with the line containing segment i */
            {
                t = -(a*pos[0] + b*pos[1] + c)/delta; 
                if (t > margin) 
                {
                    xi = pos[0] + t*c0;
                    yi = pos[1] + t*s0;
                    dx = polyline[i].x2 - polyline[i].x1;
                    
                    if (vabs(dx) > margin) s = (xi - polyline[i].x1)/dx; 
                    else s = (yi - polyline[i].y1)/(polyline[i].y2 - polyline[i].y1); 
//                     printf("s = %.2f\n", s);
                    
                    if ((s >= 0.0)&&(s <= 1.0))     
                    /* the intersection is on the segment */
                    {
                        nsegment[nt] = i;
                        tval[nt] = t;
                        sint[nt] = s;
//                         xint[nt] = pos[0] + t*c0;
//                         yint[nt] = pos[1] + t*s0;
                        xint[nt] = xi;
                        yint[nt] = yi;
                        type[nt] = 0;
//                         printf("s = %.2f, x = %.2f, y = %.2f\n", s, xint[nt], yint[nt]);
                        nt++;                        
                    }
                }
            }
        }
        
        /* look for intersections with arcs */
        for (i=0; i<narcs; i++)
        {
            x1 = pos[0] - arcs[i].xc;
            y1 = pos[1] - arcs[i].yc;
            
            b = x1*c0 + y1*s0;
            c = x1*x1 + y1*y1 - arcs[i].radius*arcs[i].radius;
            
            if (b*b - c > margin) for (sign=-1; sign<=1; sign+=2)
            {
                t = -b + (double)sign*sqrt(b*b - c);
                if (t > margin)
                {
                    xi = pos[0] + t*c0;
                    yi = pos[1] + t*s0;
                    
                    angle = argument(xi - arcs[i].xc, yi - arcs[i].yc);
                    if (angle < 0.0) angle += DPI;
                    
                    if ((angle > arcs[i].angle1)&&(angle < arcs[i].angle1 + arcs[i].dangle))
                    {
                        narc[nt] = i;
                        tval[nt] = t;
                        sint[nt] = (angle - arcs[i].angle1)/(arcs[i].dangle);
                        xint[nt] = xi;
                        yint[nt] = yi;
                        aint[nt] = angle;
                        type[nt] = 1;
//                         printf("s = %.2f, x = %.2f, y = %.2f\n", sint[nt], xint[nt], yint[nt]);
                        nt++;
                    }
                }
            }
        }
        
        if (ABSORBING_CIRCLES) for (i=0; i<ncircles; i++) 
        {
            b = (pos[0]-circles[i].xc)*c0 + (pos[1]-circles[i].yc)*s0;
            c = (pos[0]-circles[i].xc)*(pos[0]-circles[i].xc) + (pos[1]-circles[i].yc)*(pos[1]-circles[i].yc) - circles[i].radius*circles[i].radius;
        
            delta = b*b - c;
            if (delta > margin)     /* there is an intersection with circle i */
            {
                t = -b - sqrt(delta);            
                if (t > margin) 
                {
                    nsegment[nt] = -1-i;
                
                    tval[nt] = t;
                    xint[nt] = pos[0] + t*c0;
                    yint[nt] = pos[1] + t*s0;
                    sint[nt] = argument(xint[nt] - circles[i].xc, yint[nt] - circles[i].yc);

                    nt++;
                }
            }
        }
        
        if (nt > 0)     /* there is at least one intersection */
        {
            /* find earliest intersection */
            tmin = tval[0];
            ntmin = 0;
            for (i=1; i<nt; i++) 
            if (tval[i] < tmin) 
            {
                tmin = tval[i];
                ntmin = i;
            }
            
//             printf("ntmin = %i\n", ntmin); 
            if ((nsegment[ntmin] >= 0)&&(type[ntmin] == 0))     /* first intersection is with a segment */
            {
                config[0] = (double)nsegment[ntmin] + sint[ntmin];
                config[1] = polyline[nsegment[ntmin]].angle - alpha;        
                if (config[1] < 0.0) config[1] += DPI;
                if (config[1] >= PI) config[1] -= DPI;
            }
//             else if ((narc[ntmin] >= 0)&&(type[ntmin] == 1))    /* first intersection is with a segment */
            else if ((type[ntmin] == 1))    /* first intersection is with an arc */
            {
                config[0] = BOUNDARY_SHIFT + (double)narc[ntmin] + sint[ntmin]; 
                config[1] = PID + aint[ntmin] - alpha;        
                if (config[1] < 0.0) config[1] += DPI;
                if (config[1] >= DPI) config[1] -= DPI;
                
//                 printf("s = %.2f, alpha = %.2f, normal = %.2f, theta = %.2f\n", config[0], alpha*180.0/PI, aint[ntmin]*180.0/PI, config[1]*180.0/PI);
            }
            /* set dummy coordinates if circles are absorbing */
            else if ((ABSORBING_CIRCLES)&&(nsegment[ntmin] < 0))
            {
                config[0] = DUMMY_ABSORBING;
                config[1] = PI;
            }
            config[2] = 0.0;	/* running time */ 
            config[3] = module2(xint[ntmin]-pos[0], yint[ntmin]-pos[1]);     /* distance to collision */
            config[4] = pos[0];    /* start position */
            config[5] = pos[1];
            config[6] = xint[ntmin];        /* position of collision */
            config[7] = yint[ntmin];
            
            
//             print_config(config);
            
            /* returned value should be positive for end of trajectory detection */
            if (nsegment[ntmin] > 0) return(nsegment[ntmin]);
            else return(BOUNDARY_SHIFT-nsegment[ntmin]);
        }
        else    /* there is no intersection - set dummy values */
        {
            config[0] = DUMMY_ABSORBING;
            config[1] = PI;
            config[2] = 0.0;
            config[3] = rlarge;
            config[4] = pos[0];    /* start position */
            config[5] = pos[1];
            config[6] = rlarge*cos(alpha);
            config[7] = rlarge*sin(alpha);
            
            return(-1);
        }
 }

 int vpolyline_arc(double config[8])
 /* determine initial configuration when starting from boundary */
  {
	double pos[2], alpha;
	int c;

        c = pos_polyline_arc(config, pos, &alpha);
        
        vpolyline_arc_xy(config, alpha, pos);
//         c = vpolyline_xy(config, alpha, pos);
	
	return(c);
 }
 
/****************************************************************************************/
/* general billiard */
/****************************************************************************************/
 
 int pos_billiard(double conf[8], double pos[2], double *alpha)
 /* determine initial configuration for start at point pos = (x,y) */
 {
    switch (B_DOMAIN) {
        case (D_RECTANGLE):
        {
            return(pos_rectangle(conf, pos, alpha));
            break;
        }
        case (D_ELLIPSE):
        {
            return(pos_ellipse(conf, pos, alpha));
            break;
        }
        case (D_STADIUM):
        {
            return(pos_stade(conf, pos, alpha));
            break;
        }
        case (D_SINAI):
        {
            return(pos_sinai(conf, pos, alpha));
            break;
        }
        case (D_TRIANGLE):
        {
            return(pos_triangle(conf, pos, alpha));
            break;
        }
        case (D_ANNULUS):
        {
            return(pos_annulus(conf, pos, alpha));
            break;
        }
        case (D_POLYGON):
        {
            return(pos_polygon(conf, pos, alpha));
            break;
        }
        case (D_REULEAUX):
        {
            return(pos_reuleaux(conf, pos, alpha));
            break;
        }
        case (D_FLOWER):
        {
            return(pos_flower(conf, pos, alpha));
            break;
        }
        case (D_ALT_REU):
        {
            return(pos_alt_reuleaux(conf, pos, alpha));
            break;
        }
        case (D_ANGLE):
        {
            return(pos_angle(conf, pos, alpha));
            break;
        }
        case (D_LSHAPE):
        {
            return(pos_lshape(conf, pos, alpha));
            break;
        }
        case (D_GENUSN):
        {
            return(pos_genusn(conf, pos, alpha));
            break;
        }
        case (D_PARABOLAS):
        {
            return(pos_parabolas(conf, pos, alpha));
            break;
        }
        case (D_PENROSE):
        {
            return(pos_penrose(conf, pos, alpha));
            break;
        }
        case (D_CIRCLES):
        {
            return(pos_circles(conf, pos, alpha));
            break;
        }
        case (D_CIRCLES_IN_RECT):
        {
            return(pos_circles_in_rect(conf, pos, alpha));
            break;
        }
        case (D_CIRCLES_IN_GENUSN):
        {
            return(pos_circles_in_genusn(conf, pos, alpha));
            break;
        }
        case (D_CIRCLES_IN_TORUS):
        {
            return(pos_circles_in_torus(conf, pos, alpha));
            break;
        }
        case (D_POLYLINE):
        {
            return(pos_polyline(conf, pos, alpha));
            break;
        }
        case (D_POLYLINE_ARCS):
        {
            return(pos_polyline_arc(conf, pos, alpha));
            break;
        }
        default: 
        {
            printf("Function pos_billiard not defined for this billiard \n");
        }
    }
 }


 
 int vbilliard_xy(double config[8], double alpha, double pos[2])
 /* determine initial configuration for start at point pos = (x,y) */
 {
    switch (B_DOMAIN) {
        case (D_RECTANGLE):
        {
            return(vrectangle_xy(config, alpha, pos));
            break;
        }
        case (D_ELLIPSE):
        {
            return(vellipse_xy(config, alpha, pos));
            break;
        }
        case (D_STADIUM):
        {
            return(vstade_xy(config, alpha, pos));
            break;
        }
        case (D_SINAI):
        {
            return(vsinai_xy(config, alpha, pos));
            break;
        }
        case (D_TRIANGLE):
        {
            return(vtriangle_xy(config, alpha, pos));
            break;
        }
        case (D_ANNULUS):
        {
            return(vannulus_xy(config, alpha, pos));
            break;
        }
        case (D_POLYGON):
        {
            return(vpolygon_xy(config, alpha, pos));
            break;
        }
        case (D_REULEAUX):
        {
            return(vreuleaux_xy(config, alpha, pos));
            break;
        }
        case (D_FLOWER):
        {
            return(vflower_xy(config, alpha, pos));
            break;
        }
        case (D_ALT_REU):
        {
            return(valt_reuleaux_xy(config, alpha, pos));
            break;
        }
        case (D_ANGLE):
        {
            return(vangle_xy(config, alpha, pos));
            break;
        }
        case (D_LSHAPE):
        {
            return(vlshape_xy(config, alpha, pos));
            break;
        }
        case (D_GENUSN):
        {
            return(vgenusn_xy(config, alpha, pos));
            break;
        }
        case (D_PARABOLAS):
        {
            return(vparabolas_xy(config, alpha, pos));
            break;
        }
        case (D_PENROSE):
        {
            return(vpenrose_xy(config, alpha, pos));
            break;
        }
        case (D_CIRCLES):
        {
            return(vcircles_xy(config, alpha, pos));
            break;
        }
        case (D_CIRCLES_IN_RECT):
        {
            return(vcircles_in_rect_xy(config, alpha, pos));
            break;
        }
        case (D_CIRCLES_IN_GENUSN):
        {
            return(vcircles_in_genusn_xy(config, alpha, pos));
            break;
        }
        case (D_CIRCLES_IN_TORUS):
        {
            return(vcircles_in_torus_xy(config, alpha, pos));
            break;
        }
        case (D_POLYLINE):
        {
            return(vpolyline_xy(config, alpha, pos));
            break;
        }
        case (D_POLYLINE_ARCS):
        {
            return(vpolyline_arc_xy(config, alpha, pos));
            break;
        }
        default: 
        {
            printf("Function vbilliard_xy not defined for this billiard \n");
        }
    }
 }

 /* TO DO: fix returned value */
 
 int vbilliard(double config[8])
 /* determine initial configuration when starting from boundary */
 {
    double pos[2], theta, alpha;
    int c;

    switch (B_DOMAIN) {
        case (D_RECTANGLE):
        {
            c = pos_rectangle(config, pos, &alpha);
        
            return(vrectangle(config));
            break;
        }
        case (D_ELLIPSE):
        {
            c = pos_ellipse(config, pos, &alpha);
         
            return(vellipse(config));
            break;
        }
        case (D_STADIUM):
        {
            c = pos_stade(config, pos, &alpha);
        
            return(vstade(config));
            break;
        }
        case (D_SINAI):
        {
            c = pos_sinai(config, pos, &alpha);
        
            return(vsinai(config));
            break;
        }
        case (D_TRIANGLE):
        {
            c = pos_triangle(config, pos, &alpha);
        
            return(vtriangle(config));
            break;
        }
        case (D_ANNULUS):
        {
            c = pos_annulus(config, pos, &alpha);
        
            return(vannulus(config));
            break;
        }
        case (D_POLYGON):
        {
            c = pos_polygon(config, pos, &alpha);
        
            return(vpolygon(config));
            break;
        }
        case (D_REULEAUX):
        {
            c = pos_reuleaux(config, pos, &alpha);
        
            return(vreuleaux(config));
            break;
        }
        case (D_FLOWER):
        {
            c = pos_flower(config, pos, &alpha);
        
            return(vflower(config));
            break;
        }
        case (D_ALT_REU):
        {
            c = pos_alt_reuleaux(config, pos, &alpha);
        
            return(valt_reuleaux(config));
            break;
        }
        case (D_ANGLE):
        {
            c = pos_angle(config, pos, &alpha);
        
            return(vangle(config));
            break;
        }
        case (D_LSHAPE):
        {
            c = pos_lshape(config, pos, &alpha);
        
            return(vlshape(config));
            break;
        }
        case (D_GENUSN):
        {
            c = pos_genusn(config, pos, &alpha);
        
            return(vgenusn(config));
            break;
        }
        case (D_PARABOLAS):
        {
            c = pos_parabolas(config, pos, &alpha);
        
            return(vparabolas(config));
            break;
        }
        case (D_PENROSE):
        {
            c = pos_penrose(config, pos, &alpha);
        
            return(vpenrose(config));
            break;
        }
        case (D_CIRCLES):
        {
            c = pos_circles(config, pos, &alpha);
        
            return(vcircles(config));
            break;
        }
        case (D_CIRCLES_IN_RECT):
        {
            c = pos_circles_in_rect(config, pos, &alpha);
        
            return(vcircles_in_rect(config));
            break;
        }
        case (D_CIRCLES_IN_GENUSN):
        {
            c = pos_circles_in_genusn(config, pos, &alpha);
        
            return(vcircles_in_genusn(config));
            break;
        }
        case (D_CIRCLES_IN_TORUS):
        {
            c = pos_circles_in_torus(config, pos, &alpha);
        
            return(vcircles_in_torus(config));
            break;
        }
        case (D_POLYLINE):
        {
            c = pos_polyline(config, pos, &alpha);
        
            return(vpolyline(config));
            break;
        }
        case (D_POLYLINE_ARCS):
        {
            c = pos_polyline_arc(config, pos, &alpha);
        
            return(vpolyline_arc(config));
            break;
        }
        default: 
        {
            printf("Function vbilliard not defined for this billiard \n");
        }
    }
 }
 
 int xy_in_billiard(double x, double y)
 /* returns 1 if (x,y) represents a point in the billiard */
 {
    double l2, r1, r2, omega, omega2, c, angle, x1, y1, x2, co, so, x2plus, x2minus, width;
    int condition, k;
 
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
        case D_ANNULUS:
        {
            l2 = LAMBDA*LAMBDA;
            r1 = x*x + y*y;
            r2 = (x-MU)*(x-MU) + y*y;
            if ((r2 > l2)&&(r1 < 1.0)) return(1);
            else return(0);
            break;
        }
        case D_POLYGON:
        {
            return(in_polygon(x, y, 1.0, NPOLY, APOLY));
            break;
        }
        case D_REULEAUX:
        {
            condition = 1;
            omega2 = PI/((double)NPOLY);
            co = cos(omega2);
            so = sin(omega2);
            if (LAMBDA > 0.0) x2 = co + sqrt(LAMBDA*LAMBDA - so*so);
            else x2 = co - sqrt(LAMBDA*LAMBDA - so*so);
                        
            for (k=0; k<NPOLY; k++)  
            {
                angle = 2.0*(double)k*omega2 + APOLY*PID;
                
                x1 = x*cos(angle) + y*sin(angle);
                y1 = -x*sin(angle) + y*cos(angle);
                if (LAMBDA > 0.0) condition = condition*((x1-x2)*(x1-x2) + y1*y1 > LAMBDA*LAMBDA);
                else condition = condition*((x1-x2)*(x1-x2) + y1*y1 < LAMBDA*LAMBDA);
            }
            return(condition);
            break;            
        }
        /* D_REULEAUX : distance to all centers of arcs should be larger than LAMBDA */
        case D_FLOWER:
        {
            /* TO DO */
            return(1);
            break;
        }
        case D_ALT_REU:
        {
            condition = 1;
            omega2 = PI/((double)NPOLY);
            co = cos(omega2);
            so = sin(omega2);
            x2plus = co + sqrt(LAMBDA*LAMBDA - so*so);
            x2minus = co - sqrt(LAMBDA*LAMBDA - so*so);
                        
            for (k=0; k<NPOLY; k++)  
            {
                angle = 2.0*(double)k*omega2 + APOLY*PID;
                
                x1 = x*cos(angle) + y*sin(angle);
                y1 = -x*sin(angle) + y*cos(angle);
                if (k%2==0) condition = condition*((x1-x2plus)*(x1-x2plus) + y1*y1 > LAMBDA*LAMBDA);
                else condition = condition*((x1-x2minus)*(x1-x2minus) + y1*y1 < LAMBDA*LAMBDA);
            }
            return(condition);
            break;            
        }
        case D_ANGLE:
        {
            if (y <= 0.0) return(0);
            else if (x*sin(LAMBDA*PI) < -y*cos(LAMBDA*PI)) return(1);
            else return(0);
            break;
        }
        case D_LSHAPE:
        {
            if ((y <= 0.0)&&(x >= -1.0)&&(x <= 1.0)) return(1);
            else if ((x >= -1.0)&&(x <= 0.0)) return(1);
            else return(0);
            break;
        }
        case D_GENUSN:      /* same as polygon */
        {
            return(in_polygon(x, y, 1.0, NPOLY, APOLY));
            break;
        }
        case D_PARABOLAS:
        {
            condition = 1;
            omega = DPI/((double)NPOLY);
            for (k=0; k<NPOLY; k++)  
            {
                angle = APOLY*PID + (double)k*omega;
                x1 = x*cos(angle) + y*sin(angle);
                y1 = -x*sin(angle) + y*cos(angle);
                condition = condition*(x1 < LAMBDA + MU - 0.25*y1*y1/MU);
            }
            return(condition);
        }
        case D_PENROSE:
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
            else if ((vabs(x) <= c)&&(4.0*(x1-c)*(x1-c)/(MU*MU*PENROSE_RATIO*PENROSE_RATIO) + y*y/(MU*MU) <= 1.0)) return(0);
            /* straight parts */
            else if ((vabs(x) >= c)&&(vabs(y) <= width)) return(0);
            else return(1);
        }
        case D_CIRCLES:      
        {
            condition = 1;
            for (k=0; k<ncircles; k++)  
                if (circles[k].active) condition = condition*out_circle(x-circles[k].xc, y-circles[k].yc, circles[k].radius);
//                 if (circleactive[k]) condition = condition*out_circle(x-circlex[k], y-circles[k].yc, circlerad[k]);
            return(condition);
            break;
        }
        case D_CIRCLES_IN_RECT:      
        {
            if ((vabs(x) >= LAMBDA)||(vabs(y) >= 1.0)) return(0);
            else 
            {
                condition = 1;
                for (k=0; k<ncircles; k++)  
                    if (circles[k].active) condition = condition*out_circle(x-circles[k].xc, y-circles[k].yc, circles[k].radius);
                return(condition);
            }
            break;
        }
        case D_CIRCLES_IN_GENUSN:      
        {
            condition = in_polygon(x, y, 1.0, NPOLY, APOLY);
            if (condition == 0) return(0);
            else    /* test whether (x,y) outside all circles */
            {
                condition = 1;
                for (k=0; k<ncircles; k++)  
                    condition = condition*circles[k].active*out_circle(x-circles[k].xc, y-circles[k].yc, circles[k].radius);
                return(condition);
            }
            break;
        }
        case D_CIRCLES_IN_TORUS:        /* same as D_CIRCLES_IN_RECT */      
        {
            if ((vabs(x) >= LAMBDA)||(vabs(y) >= 1.0)) return(0);
            else 
            {
                condition = 1;
                for (k=0; k<ncircles; k++)  
                    if (circles[k].active) condition = condition*out_circle(x-circles[k].xc, y-circles[k].yc, circles[k].radius);
                return(condition);
            }
            break;
        }
        case D_POLYLINE:
        {
            /* not easy to implement for non-convex polygons */
            if (POLYLINE_PATTERN == P_MAZE) return ((vabs(x) < 1.1*XMAX)&&(vabs(y) < 1.1*YMAX));
            else if ((POLYLINE_PATTERN == P_MAZE_DIAG)||(POLYLINE_PATTERN == P_MAZE_RANDOM)) 
                return ((vabs(x) < 1.1*XMAX)&&(vabs(y) < 1.1*YMAX));
            else if ((POLYLINE_PATTERN == P_MAZE_CIRCULAR)||(POLYLINE_PATTERN == P_MAZE_HEX)||(POLYLINE_PATTERN == P_MAZE_OCT))
                return ((vabs(x) < 1.1*XMAX)&&(vabs(y) < 1.1*YMAX));
            else return(1);
            break;
        }
        case D_POLYLINE_ARCS:
        {
            /* not easy to implement for non-convex polygons */
            if (POLYLINE_PATTERN == P_MAZE) return ((vabs(x) < 1.1*XMAX)&&(vabs(y) < 1.1*YMAX));
            else if ((POLYLINE_PATTERN == P_MAZE_DIAG)||(POLYLINE_PATTERN == P_MAZE_RANDOM)) 
                return ((vabs(x) < 1.1*XMAX)&&(vabs(y) < 1.1*YMAX));
            else if (POLYLINE_PATTERN == P_MAZE_CIRCULAR) 
                return ((vabs(x) < 1.1*XMAX)&&(vabs(y) < 1.1*YMAX));
            else return(1);
            break;
        }
        default: 
        {
            printf("Function ij_in_billiard not defined for this billiard \n");
            return(1);
        }
    }
 }
 


void init_circles(t_circle circles[NMAXCIRCLES])
{
    int i, j, k, n, ncirc0, n_p_active, ncandidates=5000, naccepted; 
    double dx, dy, xx[4], yy[4], x, y, gamma, height, phi, r0, r, dpoisson = 3.25*MU;
    short int active_poisson[NMAXCIRCLES], far;

    switch (CIRCLE_PATTERN) {
        case (C_FOUR_CIRCLES):
        {
            ncircles = 4;
            
            circles[0].xc = 1.0;
            circles[0].yc = 0.0;
            circles[0].radius = 0.8;
                        
            circles[1].xc = -1.0;
            circles[1].yc = 0.0;
            circles[1].radius = 0.8;
            
            circles[2].xc = 0.0;
            circles[2].yc = 0.8;
            circles[2].radius = 0.4;
            
            circles[3].xc = 0.0;
            circles[3].yc = -0.8;
            circles[3].radius = 0.4;
            
            for (i=0; i<4; i++) circles[i].active = 1;

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
                    circles[n].xc = ((double)(i-NCX/2) + 0.5)*dy;
                    circles[n].yc = YMIN + ((double)j + 0.5)*dy;
                    circles[n].radius = MU;
                    circles[n].active = 1;
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
                    circles[n].xc = ((double)(i-NCX/2) + 0.5)*dy;
                    circles[n].yc = YMIN + ((double)j - 0.5)*dy;
                    if ((i+NCX)%2 == 1) circles[n].yc += 0.5*dy;
                    circles[n].radius = MU;
                    circles[n].active = 1;
                }
            break;
        }
        case (C_TRI):
        {
            ncircles = NCX*(NCY+1);
            dy = (YMAX - YMIN)/((double)NCY);
            dx = dy*0.5*sqrt(3.0);
//             dx = (XMAX - XMIN)/((double)NCX);
//             dy = dx/(0.5*sqrt(3.0));
            for (i = 0; i < NCX; i++)
                for (j = 0; j < NCY+1; j++)
                {
                    n = (NCY+1)*i + j;
                    circles[n].xc = ((double)(i-NCX/2) + 0.5)*dx;
                    circles[n].yc = YMIN + ((double)j - 0.5)*dy;
                    if ((i+NCX)%2 == 1) circles[n].yc += 0.5*dy;
                    circles[n].radius = MU;
                    circles[n].active = 1;
                }
            break;
        }
        case (C_GOLDEN_MEAN):
        {
            ncircles = 200;
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
            
            for (i=0; i<NGOLDENSPIRAL; i++) 
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
            }
        break;
        }
        case (C_RAND_DISPLACED):
        {
            ncircles = NCX*NCY;
            dy = (YMAX - YMIN)/((double)NCY);
            for (i = 0; i < NCX; i++)
                for (j = 0; j < NCY; j++)
                {
                    n = NCY*i + j;
                    circles[n].xc = ((double)(i-NCX/2) + 0.5*((double)rand()/RAND_MAX - 0.5))*dy;
                    circles[n].yc = YMIN + ((double)j + 0.5 + 0.5*((double)rand()/RAND_MAX - 0.5))*dy;
                    circles[n].radius = MU*sqrt(1.0 + 0.8*((double)rand()/RAND_MAX - 0.5));
                    circles[n].active = 1;
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
        case (C_POISSON_DISC):
        {
            printf("Generating Poisson disc sample\n");
            /* generate first circle */
            circles[0].xc = LAMBDA*(2.0*(double)rand()/RAND_MAX - 1.0);
            circles[0].yc = (YMAX - YMIN)*(double)rand()/RAND_MAX + YMIN;
            active_poisson[0] = 1;
            n_p_active = 1;
            ncircles = 1;
            
            while ((n_p_active > 0)&&(ncircles < NMAXCIRCLES))
            {
                /* randomly select an active circle */
                i = rand()%(ncircles);
                while (!active_poisson[i]) i = rand()%(ncircles);                 
//                 printf("Starting from circle %i at (%.3f,%.3f)\n", i, circlex[i], circley[i]);
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
                        circles[ncircles].xc = y;
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

        case (C_LASER):
        {
            ncircles = 17;
            
            xx[0] = 0.5*(x_shooter + x_target);
            xx[1] = LAMBDA - 0.5*(x_target - x_shooter);
            if (xx[1] > LAMBDA) xx[1] = 2.0*LAMBDA - xx[1];  
            xx[2] = -xx[0];
            xx[3] = -xx[1];
            
            yy[0] = 0.5*(y_shooter + y_target);
            yy[1] = 1.0 - 0.5*(y_target - y_shooter);
            if (yy[1] > 1.0) yy[1] = 2.0 - yy[1];  
            yy[2] = -yy[0];
            yy[3] = -yy[1];

            for (i = 0; i < 4; i++)
                for (j = 0; j < 4; j++)
                {
                    circles[4*i + j].xc = xx[i];
                    circles[4*i + j].yc = yy[j];
                    
                }
                
            circles[ncircles - 1].xc = x_target;
            circles[ncircles - 1].yc = y_target;
            
            for (i=0; i<ncircles - 1; i++)
            {
                circles[i].radius = MU;
                circles[i].active = 1;
            }
            
            circles[ncircles - 1].radius = 0.5*MU;
            circles[ncircles - 1].active = 2;
            
            break;
        }
        default: 
        {
            printf("Function init_circle_config not defined for this pattern \n");
        }
    }
}

int add_rectangle_to_polyline(double xc, double yc, double width, double height, t_segment polyline[NMAXPOLY], t_circle circles[NMAXCIRCLES])
/* add a rectangle to polyline pattern */
{
    int i;
    
    polyline[nsides].x1 = xc - 0.5*width;
    polyline[nsides].y1 = yc - 0.5*height;
    polyline[nsides].length = width;
    polyline[nsides].angle = 0.0;
            
    polyline[nsides+1].x1 = xc + 0.5*width;
    polyline[nsides+1].y1 = yc - 0.5*height;
    polyline[nsides+1].length = height;
    polyline[nsides+1].angle = PID;
            
    polyline[nsides+2].x1 = xc + 0.5*width;
    polyline[nsides+2].y1 = yc + 0.5*height;
    polyline[nsides+2].length = width;
    polyline[nsides+2].angle = PI;
            
    polyline[nsides+3].x1 = xc - 0.5*width;
    polyline[nsides+3].y1 = yc + 0.5*height;
    polyline[nsides+3].length = height;
    polyline[nsides+3].angle = 3.0*PID;
            
    if (nsides+4 < NMAXPOLY) for (i=nsides; i<nsides+4; i++) 
    {
        polyline[i].color = 0;
        if (i < nsides+3) polyline[i].x2 = polyline[i+1].x1;
        else polyline[i].x2 = polyline[nsides].x1;
        if (i < nsides+3) polyline[i].y2 = polyline[i+1].y1;
        else polyline[i].y2 = polyline[nsides].y1;        
    }
    else printf("Increase NMAXPOLY\n");
            
    if (nsides+4 < NMAXCIRCLES) for (i=nsides; i<nsides+4; i++) 
    {
        circles[i].xc = polyline[i].x1;
        circles[i].yc = polyline[i].y1;
        circles[i].radius = MU;
        circles[i].active = 1;
    }
    else 
    {
        printf("Increase NMAXCIRCLES\n");
        return(0);
    }
                
    nsides += 4;
    ncircles += 4;
    return(1);
}


int axial_symmetry_tsegment(t_segment z1, t_segment z2, t_segment z, t_segment *zprime)
/* compute reflection of point z wrt axis through z1 and z2 */
{
    double r, zdotu;
    t_segment u, zparallel, zperp;
    
    /* compute unit vector parallel to z1-z2 */
    u.x1 = z2.x1 - z1.x1;
    u.y1 = z2.y1 - z1.y1;
    r = module2(u.x1, u.y1);
    if (r == 0) return(0);      /* z1 and z2 are the same */
    
    u.x1 = u.x1/r;
    u.y1 = u.y1/r;
    
    /* projection of z1z on z1z2 */
    zdotu = (z.x1 - z1.x1)*u.x1 + (z.y1 - z1.y1)*u.y1;
    zparallel.x1 = zdotu*u.x1;
    zparallel.y1 = zdotu*u.y1;
    
    /* normal vector to z1z2 */
    zperp.x1 = z.x1 - z1.x1 - zparallel.x1;
    zperp.y1 = z.y1 - z1.y1 - zparallel.y1;
    
    /* reflected point */
    zprime->x1 = z.x1 - 2.0*zperp.x1;
    zprime->y1 = z.y1 - 2.0*zperp.y1;
    
    return(1);
}


void init_polyline_circular_maze(t_segment* polyline, t_circle* circles, t_arc* arcs, t_maze* maze)
{
    int nblocks, block, i, j, n, p, q;
    double rmin, rmax, angle, r, dr, phi, dphi, ww; 
            
    /* build walls of maze */
    nblocks = NYMAZE/NXMAZE;
    rmin = 0.15;
    rmax = 1.0;
    angle = DPI/(double)nblocks;
        
    dr = (rmax - rmin)/(double)(NXMAZE);
        
    nsides = 0;
    ncircles = 0;
    
    /* add straight walls */
    for (block = 0; block < nblocks; block++)
    {
        dphi = angle;
        
        /* first circle */
        n = nmaze(0, block*NXMAZE);
        r = rmin;
        phi = (double)block*angle;
        
        if (maze[n].south)
        {
            polyline[nsides].x1 = r*cos(phi) + MAZE_XSHIFT;
            polyline[nsides].y1 = r*sin(phi);
            polyline[nsides].x2 = (r+dr)*cos(phi) + MAZE_XSHIFT;
            polyline[nsides].y2 = (r+dr)*sin(phi);
            polyline[nsides].angle = phi;
            nsides++;
        }
                
        /* second circle */
        r = rmin + dr;
        dphi *= 0.5;
        for (q=0; q<2; q++)
        {
            n = nmaze(1, block*NXMAZE + q);
            phi = (double)(block)*angle + (double)q*dphi;
            
            if (maze[n].south)
            {
                polyline[nsides].x1 = r*cos(phi) + MAZE_XSHIFT;
                polyline[nsides].y1 = r*sin(phi);
                polyline[nsides].x2 = (r+dr)*cos(phi) + MAZE_XSHIFT;
                polyline[nsides].y2 = (r+dr)*sin(phi);
                polyline[nsides].angle = phi;
                nsides++;
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
//                 printf("Segment, i = %i, dphi = %.2lg, r = %.2lg\n", i, dphi, r);
                for (q = 0; q < 2*ww; q++)
                {
                    j = block*NXMAZE + q;
                    n = nmaze(i,j);
                    phi = (double)(block)*angle + (double)q*dphi;
                    
                    if (maze[n].south)
                    {
                        polyline[nsides].x1 = r*cos(phi) + MAZE_XSHIFT;
                        polyline[nsides].y1 = r*sin(phi);
                        polyline[nsides].x2 = (r+dr)*cos(phi) + MAZE_XSHIFT;
                        polyline[nsides].y2 = (r+dr)*sin(phi);
                        polyline[nsides].angle = phi;
                        nsides++;
                    }
                }
                i++;
            }
            ww *= 2;
        }
                
    }
            
    /* add circular arcs */
    narcs = 0;
    
    for (block = 0; block < nblocks; block++)
    {
        dphi = angle;
        
        /* first circle */
        n = nmaze(0, block*NXMAZE);
        r = rmin;
        phi = (double)block*angle;
        
        if ((block > 0)&&(maze[n].west))
        {
            arcs[narcs].xc = MAZE_XSHIFT;
            arcs[narcs].yc = 0.0;
            arcs[narcs].radius = r;
            arcs[narcs].angle1 = phi;
            arcs[narcs].dangle = dphi;
            narcs++;
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
                arcs[narcs].xc = MAZE_XSHIFT;
                arcs[narcs].yc = 0.0;
                arcs[narcs].radius = r;
                arcs[narcs].angle1 = phi;
                arcs[narcs].dangle = dphi;
                narcs++;
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
                        arcs[narcs].xc = MAZE_XSHIFT;
                        arcs[narcs].yc = 0.0;
                        arcs[narcs].radius = r;
                        arcs[narcs].angle1 = phi;
                        arcs[narcs].dangle = dphi;
                        narcs++;
                    }
                }
                i++;
            }
            ww *= 2;
        }
    }
    
    /* outer boundary of maze */
    arcs[narcs].xc = MAZE_XSHIFT;
    arcs[narcs].yc = 0.0;
    arcs[narcs].radius = rmax;
    if (CLOSE_MAZE) 
    {
        arcs[narcs].angle1 = 0.0;
        arcs[narcs].dangle = DPI;
    }
    else 
    {
        arcs[narcs].angle1 = dphi;
        arcs[narcs].dangle = DPI - dphi;
    }
    narcs++;
            
}

void init_polyline_hex_maze(t_segment* polyline, t_circle* circles, t_maze* maze)
{
    int i, j, n;
    double r, r1, x1, y1, padding = 0.02, h; 
            
    r = 2.0*(YMAX - YMIN)/(3.0*((double)(NXMAZE)-0.5));
    r1 = (YMAX - YMIN)/(sqrt(3.0)*(double)(NYMAZE+1));
    if (r1 < r) r = r1;
    h = 0.5*sqrt(3.0)*r;
    
    for (i = 0; i < NXMAZE; i++)
        for (j = 0; j < NYMAZE; j++)
        {
            n = nmaze(i, j);
            x1 = YMIN + padding + 1.5*(double)i*r + MAZE_XSHIFT;
            y1 = YMIN + padding + h + 2.0*(double)j*h;
            if (i%2 == 0) y1 += h;
            
            if (((i>0)||(j!=NYMAZE/2))&&(maze[n].southwest)) 
            {
                polyline[nsides].x1 = x1 - 0.5*r;
                polyline[nsides].y1 = y1 - h;
                polyline[nsides].x2 = x1 - r;
                polyline[nsides].y2 = y1;
                polyline[nsides].angle = DPI/3.0;
                nsides++;
            }
            
            if (maze[n].south) 
            {
                polyline[nsides].x1 = x1 - 0.5*r;
                polyline[nsides].y1 = y1 - h;
                polyline[nsides].x2 = x1 + 0.5*r;
                polyline[nsides].y2 = y1 - h;
                polyline[nsides].angle = 0.0;
                nsides++;
            }
            
            if (maze[n].northwest)
            {
                polyline[nsides].x1 = x1 - r;
                polyline[nsides].y1 = y1;
                polyline[nsides].x2 = x1 - 0.5*r;
                polyline[nsides].y2 = y1 + h;
                polyline[nsides].angle = PI/3.0;
                nsides++;
            }
            
            if (((j==0)||(i==NXMAZE-1))&&(maze[n].southeast))
            {
                polyline[nsides].x1 = x1 + 0.5*r;
                polyline[nsides].y1 = y1 - h;
                polyline[nsides].x2 = x1 + r;
                polyline[nsides].y2 = y1;
                polyline[nsides].angle = PI/3.0;
                nsides++;
            }
            
            if ((j==NYMAZE-1)&&(maze[n].north))
            {
                polyline[nsides].x1 = x1 - 0.5*r;
                polyline[nsides].y1 = y1 + h;
                polyline[nsides].x2 = x1 + 0.5*r;
                polyline[nsides].y2 = y1 + h;
                polyline[nsides].angle = 0.0;
                nsides++;
            }
            
            if (((j==NYMAZE-1)||((i==NXMAZE-1)&&(j!=NYMAZE/2-1)))&&(maze[n].northeast))
            {
                polyline[nsides].x1 = x1 + r;
                polyline[nsides].y1 = y1;
                polyline[nsides].x2 = x1 + 0.5*r;
                polyline[nsides].y2 = y1 + h;
                polyline[nsides].angle = DPI/3.0;
                nsides++;
            }
        }
    /* left side of maze */
    x1 = YMIN + padding + MAZE_XSHIFT - 0.5*r;
    y1 = YMIN + padding + h;
    polyline[nsides].x1 = x1 + 1.5*r;
    polyline[nsides].y1 = YMIN - MAZE_VERTICAL_EXTENSION;
    polyline[nsides].x2 = x1 + 1.5*r;
    polyline[nsides].y2 = y1 - h;
    polyline[nsides].angle = PID;
    nsides++;
            
    y1 = YMIN + padding + 3.0*h + 2.0*(double)(NYMAZE-1)*h;
    polyline[nsides].x1 = x1;
    polyline[nsides].y1 = y1;
    polyline[nsides].x2 = x1;
    polyline[nsides].y2 = YMAX + MAZE_VERTICAL_EXTENSION;
    polyline[nsides].angle = PID;
    nsides++;
            
    /* right side of maze */
    x1 = YMIN + padding + 1.5*(double)(NXMAZE-1)*r + MAZE_XSHIFT + 0.5*r;
    y1 = YMIN + padding;
    polyline[nsides].x1 = x1;
    polyline[nsides].y1 = YMIN - MAZE_VERTICAL_EXTENSION;
    polyline[nsides].x2 = x1;
    polyline[nsides].y2 = y1;
    polyline[nsides].angle = PID;
    nsides++;
            
    y1 = YMIN + padding + 2.0*h + 2.0*(double)(NYMAZE-1)*h;
    polyline[nsides].x1 = x1 - 1.5*r;
    polyline[nsides].y1 = y1 + h;
    polyline[nsides].x2 = x1 - 1.5*r;
    polyline[nsides].y2 = YMAX + MAZE_VERTICAL_EXTENSION;
    polyline[nsides].angle = PID;
    nsides++;
    
    /* add circular arcs in corners */
    if (B_DOMAIN == D_POLYLINE_ARCS)
    {
        narcs = 0;
        
        for (i=0; i<NXMAZE; i++)
            for (j=0; j<NYMAZE; j++)
            {
                n = nmaze(i, j);
                x1 = YMIN + padding + 1.5*(double)i*r + MAZE_XSHIFT;
                y1 = YMIN + padding + h + 2.0*(double)j*h;
                if (i%2 == 0) y1 += h;
                    
                if (((i < NXMAZE-1)||(j!=NYMAZE/2-1))&&(maze[n].northeast)&&(maze[n].north)) 
                {
                    arcs[narcs].xc = x1;
                    arcs[narcs].yc = y1;
                    arcs[narcs].radius = h;
                    arcs[narcs].dangle = PI/3.0;
                    arcs[narcs].angle1 = PI/6.0;
                    narcs++;
                }
                if (((i > 0)||(j!=NYMAZE/2))&&(maze[n].north)&&(maze[n].northwest)) 
                {
                    arcs[narcs].xc = x1;
                    arcs[narcs].yc = y1;
                    arcs[narcs].radius = h;
                    arcs[narcs].dangle = PI/3.0;
                    arcs[narcs].angle1 = PID;
                    narcs++;
                }        
                if (((i > 0)||(j!=NYMAZE/2))&&(maze[n].northwest)&&(maze[n].southwest)) 
                {
                    arcs[narcs].xc = x1;
                    arcs[narcs].yc = y1;
                    arcs[narcs].radius = h;
                    arcs[narcs].dangle = PI/3.0;
                    arcs[narcs].angle1 = 5.0*PI/6.0;
                    narcs++;
                }        
                if (((i > 0)||(j!=NYMAZE/2))&&(maze[n].southwest)&&(maze[n].south)) 
                {
                    arcs[narcs].xc = x1;
                    arcs[narcs].yc = y1;
                    arcs[narcs].radius = h;
                    arcs[narcs].dangle = PI/3.0;
                    arcs[narcs].angle1 = 7.0*PI/6.0;
                    narcs++;
                }        
                if (((i > 0)||(j!=NYMAZE/2))&&(maze[n].south)&&(maze[n].southeast)) 
                {
                    arcs[narcs].xc = x1;
                    arcs[narcs].yc = y1;
                    arcs[narcs].radius = h;
                    arcs[narcs].dangle = PI/3.0;
                    arcs[narcs].angle1 = 3.0*PID;
                    narcs++;
                }        
                if (((i < NXMAZE-1)||(j!=NYMAZE/2-1))&&(maze[n].southeast)&&(maze[n].northeast)) 
                {
                    arcs[narcs].xc = x1;
                    arcs[narcs].yc = y1;
                    arcs[narcs].radius = h;
                    arcs[narcs].dangle = PI/3.0;
                    arcs[narcs].angle1 = 11.0*PI/6.0;
                    narcs++;
                }        
            }
    }
}

void init_polyline_oct_maze(t_segment* polyline, t_circle* circles, t_maze* maze)
{
    int i, j, n;
    double a, b, a2, c, dx, x1, y1, padding = 0.02; 
            
    dx = (YMAX - YMIN - 2.0*padding)/((double)NYMAZE+0.5);    
    a = dx*(2.0 - sqrt(2.0));
    b = a/sqrt(2.0); 
    a2 = 0.5*a;
    c = a2 + b;
    
    for (i = 0; i < NXMAZE; i++)
        for (j = 0; j < NYMAZE; j++)
        {
            n = nmaze(i, j);
            x1 = YMIN + padding + ((double)i + 0.5)*dx + MAZE_XSHIFT;
            y1 = YMIN + padding + ((double)j + 0.75)*dx;
            
            if ((i+j)%2 == 0)
            {
                if (((i>0)||(j!=NYMAZE/2))&&(maze[n].west))
                {
                    polyline[nsides].x1 = x1 - c;
                    polyline[nsides].y1 = y1 - a2;
                    polyline[nsides].x2 = x1 - c;
                    polyline[nsides].y2 = y1 + a2;
                    polyline[nsides].angle = PID;
                    nsides++;                    
                }
                
                if (maze[n].south)
                {
                    polyline[nsides].x1 = x1 - a2;
                    polyline[nsides].y1 = y1 - c;
                    polyline[nsides].x2 = x1 + a2;
                    polyline[nsides].y2 = y1 - c;
                    polyline[nsides].angle = 0.0;                    
                    nsides++;                    
                }
                
                if (maze[n].southwest)
                {
                    polyline[nsides].x1 = x1 - a2;
                    polyline[nsides].y1 = y1 - c;
                    polyline[nsides].x2 = x1 - c;
                    polyline[nsides].y2 = y1 - a2;
                    polyline[nsides].angle = 1.5*PID;                    
                    nsides++;                    
                }

                if (maze[n].northwest)
                {
                    polyline[nsides].x1 = x1 - c;
                    polyline[nsides].y1 = y1 + a2;
                    polyline[nsides].x2 = x1 - a2;
                    polyline[nsides].y2 = y1 + c;
                    polyline[nsides].angle = 0.5*PID;                    
                    nsides++;                    
                }
            }
            else
            {
                if (((i>0)||(j!=NYMAZE/2))&&(maze[n].west))
                {
                    polyline[nsides].x1 = x1 - a2;
                    polyline[nsides].y1 = y1 - a2;
                    polyline[nsides].x2 = x1 - a2;
                    polyline[nsides].y2 = y1 + a2;
                    polyline[nsides].angle = PID;
                    nsides++;                    
                }
                
                if (maze[n].south)
                {
                    polyline[nsides].x1 = x1 - a2;
                    polyline[nsides].y1 = y1 - a2;
                    polyline[nsides].x2 = x1 + a2;
                    polyline[nsides].y2 = y1 - a2;
                    polyline[nsides].angle = 0.0;                    
                    nsides++;                    
                }
                
                
            }
        }
        
        /* bottom of maze */
        for (i=0; i<NXMAZE; i+=2)
        {
            n = nmaze(i, 0);
            x1 = YMIN + padding + ((double)i + 0.5)*dx + MAZE_XSHIFT;
            y1 = YMIN + padding + 0.75*dx;
            
            polyline[nsides].x1 = x1 + a2;
            polyline[nsides].y1 = y1 - c;
            polyline[nsides].x2 = x1 + c;
            polyline[nsides].y2 = y1 - a2;
            polyline[nsides].angle = 0.5*PID;                    
            nsides++;                    
        }
        
        /* top of maze */
        for (i=0; i<NXMAZE; i++) 
        {
            n = nmaze(i, NYMAZE-1);
            x1 = YMIN + padding + ((double)i + 0.5)*dx + MAZE_XSHIFT;
            y1 = YMIN + padding + ((double)(NYMAZE-1) + 0.75)*dx;
            
            if ((i+NYMAZE-1)%2 == 0)
            {
                polyline[nsides].x1 = x1 + c;
                polyline[nsides].y1 = y1 + a2;
                polyline[nsides].x2 = x1 + a2;
                polyline[nsides].y2 = y1 + c;
                polyline[nsides].angle = 1.5*PID;                    
                nsides++;   
                
                polyline[nsides].x1 = x1 + a2;
                polyline[nsides].y1 = y1 + c;
                polyline[nsides].x2 = x1 - a2;
                polyline[nsides].y2 = y1 + c;
                polyline[nsides].angle = 0.0;                    
                nsides++;                    
            }
            else
            {
                polyline[nsides].x1 = x1 + a2;
                polyline[nsides].y1 = y1 + a2;
                polyline[nsides].x2 = x1 - a2;
                polyline[nsides].y2 = y1 + a2;
                polyline[nsides].angle = 0.0;                    
                nsides++;                    
            }
        }
        
        /* right side of maze */
        for (j=0; j<NYMAZE; j++)
        {
            n = nmaze(NXMAZE-1, j);
            x1 = YMIN + padding + ((double)(NXMAZE-1) + 0.5)*dx + MAZE_XSHIFT;
            y1 = YMIN + padding + ((double)j + 0.75)*dx;
            
            if ((NXMAZE-1+j)%2 == 0)
            {
                polyline[nsides].x1 = x1 + a2;
                polyline[nsides].y1 = y1 - c;
                polyline[nsides].x2 = x1 + c;
                polyline[nsides].y2 = y1 - a2;
                polyline[nsides].angle = 0.5*PID;                    
                nsides++; 
                polyline[nsides].x1 = x1 + c;
                polyline[nsides].y1 = y1 + a2;
                polyline[nsides].x2 = x1 + a2;
                polyline[nsides].y2 = y1 + c;
                polyline[nsides].angle = 1.5*PID;                    
                nsides++; 
                if ((j!=NYMAZE/2)&&(j!=NYMAZE/2+1))
                {
                    polyline[nsides].x1 = x1 + c;
                    polyline[nsides].y1 = y1 - a2;
                    polyline[nsides].x2 = x1 + c;
                    polyline[nsides].y2 = y1 + a2;
                    polyline[nsides].angle = PID;                    
                    nsides++;
                }
            }
            else
            {
                polyline[nsides].x1 = x1 + a2;
                polyline[nsides].y1 = y1 - a2;
                polyline[nsides].x2 = x1 + a2;
                polyline[nsides].y2 = y1 + a2;
                polyline[nsides].angle = PID;                    
                nsides++;
            }
        }
            
    /* add circular arcs in corners */
    if (B_DOMAIN == D_POLYLINE_ARCS)
    {
        narcs = 0;
        
        
    }
}

void compute_maze_boundaries(int type, double *xmin, double *xmax)
{
    double r, r1, h, padding = 0.02;
    
    switch (type) {
        case (P_MAZE_HEX):
        {
            r = 2.0*(YMAX - YMIN)/(3.0*((double)(NXMAZE)-0.5));
            r1 = 2.0*(YMAX - YMIN)/(sqrt(3.0)*(double)(NYMAZE));
            if (r1 < r) r = r1;
            
            *xmin = YMIN + padding + MAZE_XSHIFT - 2.0*r;
            *xmax = YMIN + padding + 1.5*(double)(NXMAZE-1)*r + MAZE_XSHIFT + 2.0*r;
            
            printf("xmin = %.3lg, xmax = %.3lg\n", *xmin, *xmax);
            break;
        }
        default:
        {
            *xmin = YMIN + MAZE_XSHIFT; 
            *xmax = YMAX + MAZE_XSHIFT;
            break;
        }
    }
    
}

void init_polyline(t_segment polyline[NMAXPOLY], t_circle circles[NMAXCIRCLES], t_arc arcs[NMAXCIRCLES])
{
    int i, j, k, l, n, z, ii, jj, terni[SDEPTH], ternj[SDEPTH], quater[SDEPTH], cond, p, q, block, nblocks, ww;
    short int vkoch[NMAXCIRCLES], turnright; 
    double ratio, omega, angle, sw, length, dist, x, y, ta, tb, a, b, ra, rb, r, 
    x1, y1, x2, y2, dx, dy, padding = 0.02, width = 0.01, xright, yright, xtop, ytop, rand_factor, rmin, rmax, dr, phi, dphi, rscat;
    double *maze_coords;
    t_maze* maze;
    
    narcs = 0;      /* default value */
    
    switch (POLYLINE_PATTERN) {
        case (P_RECTANGLE):
        {
            add_rectangle_to_polyline(0.0, 0.0, 2.0*LAMBDA, 2.0, polyline, circles);
            break;
        }
        case (P_TOKARSKY):
        {
            nsides = 26;
            ncircles = 26;
            
            polyline[0].x1 = 0.0;   polyline[0].y1 = 2.0;       polyline[0].angle = 0.0;
            polyline[1].x1 = 1.0;   polyline[1].y1 = 2.0;       polyline[1].angle = -PID;
            polyline[2].x1 = 1.0;   polyline[2].y1 = 1.0;       polyline[2].angle = 0.0;
            polyline[3].x1 = 2.0;   polyline[3].y1 = 1.0;       polyline[3].angle = -PID;
            polyline[4].x1 = 2.0;   polyline[4].y1 = 0.0;       polyline[4].angle = 0.5*PID;
            polyline[5].x1 = 3.0;   polyline[5].y1 = 1.0;       polyline[5].angle = 0.0;
            polyline[6].x1 = 4.0;   polyline[6].y1 = 1.0;       polyline[6].angle = -PID;
            polyline[7].x1 = 4.0;   polyline[7].y1 = 0.0;       polyline[7].angle = 0.5*PID;
            polyline[8].x1 = 5.0;   polyline[8].y1 = 1.0;       polyline[8].angle = 0.0;
            polyline[9].x1 = 6.0;   polyline[9].y1 = 1.0;       polyline[9].angle = -PID;
            polyline[10].x1 = 6.0;   polyline[10].y1 = 0.0;     polyline[10].angle = 0.5*PID;
            polyline[11].x1 = 7.0;   polyline[11].y1 = 1.0;     polyline[11].angle = PID;
            polyline[12].x1 = 7.0;   polyline[12].y1 = 2.0;     polyline[12].angle = 0.0;
            polyline[13].x1 = 8.0;   polyline[13].y1 = 2.0;     polyline[13].angle = PID;
            polyline[14].x1 = 8.0;   polyline[14].y1 = 3.0;     polyline[14].angle = PI;
            polyline[15].x1 = 7.0;   polyline[15].y1 = 3.0;     polyline[15].angle = 1.5*PID;
            polyline[16].x1 = 6.0;   polyline[16].y1 = 4.0;     polyline[16].angle = -PID;
            polyline[17].x1 = 6.0;   polyline[17].y1 = 3.0;     polyline[17].angle = PI;
            polyline[18].x1 = 5.0;   polyline[18].y1 = 3.0;     polyline[18].angle = -PID;
            polyline[19].x1 = 5.0;   polyline[19].y1 = 2.0;     polyline[19].angle = PI;
            polyline[20].x1 = 3.0;   polyline[20].y1 = 2.0;     polyline[20].angle = PID;
            polyline[21].x1 = 3.0;   polyline[21].y1 = 3.0;     polyline[21].angle = PI;
            polyline[22].x1 = 2.0;   polyline[22].y1 = 3.0;     polyline[22].angle = PID;
            polyline[23].x1 = 2.0;   polyline[23].y1 = 4.0;     polyline[23].angle = PI;
            polyline[24].x1 = 1.0;   polyline[24].y1 = 4.0;     polyline[24].angle = -PID;
            polyline[25].x1 = 1.0;   polyline[25].y1 = 3.0;     polyline[25].angle = 2.5*PID;
            
            ratio = (XMAX - XMIN)/8.4;
            for (i=0; i<nsides; i++)
            {
                polyline[i].x1 = ratio*(polyline[i].x1 - 4.0);
                polyline[i].y1 = ratio*(polyline[i].y1 - 2.0);
            }
                
            for (i=0; i<nsides; i++) if (i < nsides-1) 
            {   
                polyline[i].x2 = polyline[i+1].x1;
                polyline[i].y2 = polyline[i+1].y1;
            }
            
            polyline[nsides-1].x2 = polyline[0].x1;
            polyline[nsides-1].y2 = polyline[0].y1;
            
            for (i=0; i<nsides; i++) 
                polyline[i].length = module2(polyline[i].x2 - polyline[i].x1, polyline[i].y2 - polyline[i].y1);
            
            for (i=0; i<ncircles; i++)
            {
                circles[i].xc = polyline[i].x1;
                circles[i].yc = polyline[i].y1;
                circles[i].radius = MU;
                circles[i].active = 1;
            }
                        
            break;
        }
        case (P_POLYRING):
        {
            nsides = 2*NPOLY;
            ncircles = 0;
            omega = DPI/(double)NPOLY;
            sw = sin(omega/2.0);
            
            for (i=0; i<NPOLY; i++)
            {
                angle = APOLY + (double)i*omega;
                polyline[i].x1 = LAMBDA*cos(angle);
                polyline[i].y1 = LAMBDA*sin(angle);
                polyline[i].angle = angle + PID + 0.5*omega;
                polyline[i].length = 2.0*LAMBDA*sw;
            }
            for (i=0; i<NPOLY; i++)
            {
                polyline[i].x2 = polyline[(i+1)%NPOLY].x1;
                polyline[i].y2 = polyline[(i+1)%NPOLY].y1;
            }
            
            for (i=0; i<NPOLY; i++)
            {
                angle = APOLY + ((double)i+0.5)*omega;
                polyline[i+NPOLY].x1 = MU*cos(angle);
                polyline[i+NPOLY].y1 = MU*sin(angle);
                polyline[i+NPOLY].angle = angle + PID + 0.5*omega;
                polyline[i+NPOLY].length = 2.0*MU*sw;
            }
            for (i=0; i<NPOLY; i++)
            {
                polyline[i+NPOLY].x2 = polyline[(i+1)%NPOLY+NPOLY].x1;
                polyline[i+NPOLY].y2 = polyline[(i+1)%NPOLY+NPOLY].y1;
            }
            
            for (i=0; i<nsides; i++) polyline[i].color = 0;
                            
            break;
        }
        case (P_SIERPINSKI):
        {
            nsides = 0;
            ncircles = 0;
            
            add_rectangle_to_polyline(0.0, 0.0, 2.0*LAMBDA, 2.0, polyline, circles);
            
            length = 2.0/3.0;
            dist = 2.0;
            n = 1;
            
            for (k=0; k<SDEPTH; k++)
            {
                for (i=0; i<n; i++)
                    for (j=0; j<n; j++)
                    {
                        /* compute ternary expansion of i */
                        ii = i;
                        for (l=0; l<k; l++)
                        {
                            terni[l] = ii%3;
                            ii = ii - (ii%3);
                            ii = ii/3;
                        }
                        
                        /* compute ternary expansion of j */
                        jj = j;
                        for (l=0; l<k; l++)
                        {
                            ternj[l] = jj%3;
                            jj = jj - (jj%3);
                            jj = jj/3;
                        }
                        
                        /* check whether ternary expansions do not have 1 at same position */
                        cond = 1;
                        for (l=0; l<k; l++) 
                            if ((terni[l] == 1)&&(ternj[l] == 1)) cond = 0;
                            
                        if (cond)
                        {
                            x = -1.0 + dist*((double)i + 0.5);
                            y = -1.0 + dist*((double)j + 0.5);
                            add_rectangle_to_polyline(x, y, length, length, polyline, circles);
                        }
                    }
                length = length/3.0;
                dist = dist/3.0;
                n = n*3;
            }
                         
            printf("nsides = %i\n", nsides);
            
            break;
        }
        case (P_VONKOCH):
        {
            nsides = 3;
            for (k=0; k<SDEPTH; k++) nsides *= 4;
            printf("nsides = %i\n", nsides);
            ncircles = nsides;
            
            if (nsides > NMAXPOLY)
            {
                printf("NMAXPOLY has to be increased to at least %i\n", nsides);
                nsides = NMAXPOLY;
            }

            for (i=0; i<nsides/3; i++)
            {
                /* compute quaternary expansion of i */
                ii = i;
                for (l=0; l<SDEPTH; l++)
                {
                    quater[l] = ii%4;
                    ii = ii - (ii%4);
                    ii = ii/4;
                }
                
                /* find first nonzero digit */
                z = 0;
                while ((quater[z] == 0)&&(z<SDEPTH)) z++;
                
                /* compute left/right turns */
                if (i==0) vkoch[0] = 0;
                else if (z != SDEPTH)
                {   
                    if (quater[z] == 2) vkoch[i] = 0;
                    else vkoch[i] = 1;
                }
//                 printf("%i", vkoch[i]);
            }
            printf("\n");
            
            /* compute vertices */
            angle = 2.0*PI/3.0;
            x = cos(PI/6.0);
            y = -sin(PI/6.0);
            length = 2.0*sin(PI/3.0);
            
            for (k=0; k<SDEPTH; k++) length = length/3.0;
            printf("Length = %.2f\n", length);
            
            for (i=0; i<nsides; i++)
            {
                polyline[i].x1 = x;
                polyline[i].y1 = y; 
                polyline[i].angle = angle;
                
                x += length*cos(angle);
                y += length*sin(angle);
                polyline[i].length = length;
                                
                turnright = vkoch[i%(nsides/3)+1];
                if (turnright) angle -= PI/3.0;
                else angle += 2.0*PI/3.0;
                
                while (angle > DPI) angle -= DPI;
                while (angle < 0.0) angle += DPI;                
            }
            
            for (i=0; i<nsides; i++)
            {
                polyline[i].color = 0;
                if (i < nsides-1) polyline[i].x2 = polyline[i+1].x1;
                else polyline[i].x2 = polyline[0].x1;
                if (i < nsides-1) polyline[i].y2 = polyline[i+1].y1;
                else polyline[i].y2 = polyline[0].y1;  
            }
            
            for (i=0; i<ncircles; i++)
            {
                circles[i].xc = polyline[i].x1;
                circles[i].yc = polyline[i].y1;
                circles[i].radius = MU;
                circles[i].active = 1;
            }
            
            break;
        }
        case (P_POLYGON):
        {
            nsides = NPOLY;
            ncircles = 0;
            omega = DPI/(double)NPOLY;
            sw = sin(omega/2.0);
            
            for (i=0; i<NPOLY; i++)
            {
                angle = APOLY + (double)i*omega;
                polyline[i].x1 = LAMBDA*cos(angle);
                polyline[i].y1 = LAMBDA*sin(angle);
                polyline[i].angle = angle + PID + 0.5*omega;
                polyline[i].length = 2.0*LAMBDA*sw;
            }
            for (i=0; i<NPOLY; i++)
            {
                polyline[i].x2 = polyline[(i+1)%NPOLY].x1;
                polyline[i].y2 = polyline[(i+1)%NPOLY].y1;
            }
                
            for (i=0; i<nsides; i++) polyline[i].color = 0;
                            
            break;
        }
        case (P_TOKA_PRIME):
        {
            nsides = 84;
            ncircles = 84;
            
            polyline[0].x1 = 0.0;
            polyline[0].y1 = 1.0;
    
            polyline[42].x1 = 0.0;
            polyline[42].y1 = 1.0 - LAMBDA;
    
            ta = tan(0.05*PI);
            tb = tan(0.4*PI);
    
            a = LAMBDA*tb/(ta + tb);
            b = a*ta;
    
            polyline[41].x1 = b;
            polyline[41].y1 = 1.0 - a;
    
            axial_symmetry_tsegment(polyline[0], polyline[41], polyline[42], &polyline[40]);
            axial_symmetry_tsegment(polyline[0], polyline[40], polyline[41], &polyline[1]);
            axial_symmetry_tsegment(polyline[40], polyline[1], polyline[0], &polyline[84]);
            
            axial_symmetry_tsegment(polyline[1], polyline[84], polyline[40], &polyline[2]);
            for (i=2; i<39; i++)
                axial_symmetry_tsegment(polyline[i], polyline[84], polyline[i-1], &polyline[i+1]);
    
            for (i=1; i<42; i++)
            {
                polyline[84-i].x1 = -polyline[i].x1;
                polyline[84-i].y1 = polyline[i].y1;
            }
            
            printf("xc = %.5f, yc = %.5f\n", polyline[84].x1, polyline[84].y1);
    
            for (i=0; i<nsides; i++)
            {
//                 polyline[i].x1 += -MU;
                x = polyline[(i+1)%nsides].x1;
                y = polyline[(i+1)%nsides].y1;
                
                polyline[i].x2 = x;
                polyline[i].y2 = y;
                polyline[i].length = module2(x - polyline[i].x1, y - polyline[i].y1);
                polyline[i].angle = argument(x - polyline[i].x1, y - polyline[i].y1);
            }
            
            for (i=0; i<ncircles; i++)
            {
                circles[i].xc = polyline[i].x1;
                circles[i].yc = polyline[i].y1;
                circles[i].radius = MU;
                circles[i].active = 1;
            }
    
            break;
        }
        case (P_TREE):
        {
            nsides = 11;
            ncircles = 11;
            
            polyline[0].x1 = 0.85;
            polyline[0].y1 = -1.0;
            
            polyline[1].x1 = 0.4;
            polyline[1].y1 = -0.4;
            
            polyline[2].x1 = 0.65;
            polyline[2].y1 = -0.4;
            
            polyline[3].x1 = 0.25;
            polyline[3].y1 = 0.2;
            
            polyline[4].x1 = 0.5;
            polyline[4].y1 = 0.2;
            
            polyline[5].x1 = 0.0;
            polyline[5].y1 = 1.0;
            
            for (i=6; i<11; i++)
            {
                polyline[i].x1 = -polyline[10-i].x1;
                polyline[i].y1 = polyline[10-i].y1;
            }

            for (i=0; i<nsides; i++)
            {
                x = polyline[(i+1)%nsides].x1;
                y = polyline[(i+1)%nsides].y1;
                
                polyline[i].x2 = x;
                polyline[i].y2 = y;
                polyline[i].length = module2(x - polyline[i].x1, y - polyline[i].y1);
                polyline[i].angle = argument(x - polyline[i].x1, y - polyline[i].y1);
            }
            
            for (i=0; i<ncircles; i++)
            {
                circles[i].xc = polyline[i].x1;
                circles[i].yc = polyline[i].y1;
                circles[i].radius = MU;
                circles[i].active = 1;
            }
            break;
        }
        case (P_TOKA_NONSELF):
        {
            nsides = 12;
            ncircles = 12;
            
            polyline[0].x1 = 0.0;    polyline[0].y1 = 2.0;       polyline[0].angle = 3.0*PID;
            polyline[1].x1 = 0.0;    polyline[1].y1 = 1.0;       polyline[1].angle = 0.0;
            polyline[2].x1 = 1.0;    polyline[2].y1 = 1.0;       polyline[2].angle = 3.5*PID;
            polyline[3].x1 = 2.0;    polyline[3].y1 = 0.0;       polyline[3].angle = PI;
            polyline[4].x1 = 1.0;    polyline[4].y1 = 0.0;       polyline[4].angle = 3.0*PID;
            polyline[5].x1 = 1.0;    polyline[5].y1 = -1.0;      polyline[5].angle = 2.5*PID;
            polyline[6].x1 = 0.0;    polyline[6].y1 = -2.0;      polyline[6].angle = PID;
            polyline[7].x1 = 0.0;    polyline[7].y1 = -1.0;      polyline[7].angle = PI;
            polyline[8].x1 = -1.0;   polyline[8].y1 = -1.0;      polyline[8].angle = 1.5*PID;
            polyline[9].x1 = -2.0;   polyline[9].y1 = 0.0;       polyline[9].angle = 0.0;
            polyline[10].x1 = -1.0;  polyline[10].y1 = 0.0;      polyline[10].angle = PID;
            polyline[11].x1 = -1.0;  polyline[11].y1 = 1.0;     polyline[11].angle = 0.5*PID;
            
            ratio = (YMAX - YMIN)/4.5;
            for (i=0; i<nsides; i++)
            {
                polyline[i].x1 = ratio*(polyline[i].x1);
                polyline[i].y1 = ratio*(polyline[i].y1);
            }
                
            for (i=0; i<nsides; i++) if (i < nsides-1) 
            {   
                polyline[i].x2 = polyline[i+1].x1;
                polyline[i].y2 = polyline[i+1].y1;
            }
            
            polyline[nsides-1].x2 = polyline[0].x1;
            polyline[nsides-1].y2 = polyline[0].y1;
            
            for (i=0; i<nsides; i++) 
                polyline[i].length = module2(polyline[i].x2 - polyline[i].x1, polyline[i].y2 - polyline[i].y1);
            
            for (i=0; i<ncircles; i++)
            {
                circles[i].xc = polyline[i].x1;
                circles[i].yc = polyline[i].y1;
                circles[i].radius = MU;
                circles[i].active = 1;
            }
            break;
        }
        case (P_MAZE):
        {
            maze = (t_maze *)malloc(NXMAZE*NYMAZE*sizeof(t_maze));
    
            init_maze_exit(0, NYMAZE/2, maze);
            
            /* build walls of maze */
            dx = (YMAX - YMIN - 2.0*padding)/(double)(NXMAZE);
            dy = (YMAX - YMIN - 2.0*padding)/(double)(NYMAZE);
            
            nsides = 0;
            ncircles = 0;
    
            for (i=0; i<NXMAZE; i++)
                for (j=0; j<NYMAZE; j++)
                {
                    n = nmaze(i, j);
                    x1 = YMIN + padding + (double)i*dx + MAZE_XSHIFT;
                    y1 = YMIN + padding + (double)j*dy;
            
                    if (((i>0)||(j!=NYMAZE/2)||(CLOSE_MAZE))&&(maze[n].west)) 
                    {
                        polyline[nsides].x1 = x1;
                        polyline[nsides].y1 = y1;
                        polyline[nsides].x2 = x1;
                        polyline[nsides].y2 = y1 + dy;
                        polyline[nsides].angle = PID;
                        nsides++;
                    }
                    
                    if (maze[n].south) 
                    {
                        polyline[nsides].x1 = x1;
                        polyline[nsides].y1 = y1;
                        polyline[nsides].x2 = x1 + dx;
                        polyline[nsides].y2 = y1;
                        polyline[nsides].angle = 0.0;
                        nsides++;
                    }
                }
    
            /* top side of maze */
            polyline[nsides].x1 = YMIN + padding + MAZE_XSHIFT;
            polyline[nsides].y1 = YMAX - padding;
            polyline[nsides].x2 = YMAX - padding + MAZE_XSHIFT;
            polyline[nsides].y2 = YMAX - padding;
            polyline[nsides].angle = 0.0;
            nsides++;
    
            /* right side of maze */
            y1 = YMIN + padding + dy*((double)NYMAZE/2);
            x1 = YMAX - padding + MAZE_XSHIFT;
            polyline[nsides].x1 = x1;
            polyline[nsides].y1 = YMIN - 1000.0;
            polyline[nsides].x2 = x1;
            if (CLOSE_MAZE) polyline[nsides].y2 = y1;
            else polyline[nsides].y2 = y1 - dy;
            polyline[nsides].angle = PID;
            nsides++;
            
            polyline[nsides].x1 = x1;
            polyline[nsides].y1 = y1;
            polyline[nsides].x2 = x1;
            polyline[nsides].y2 = YMAX + 1000.0;
            polyline[nsides].angle = PID;
            nsides++;
    
            /* left side of maze */
            x1 = YMIN + padding + MAZE_XSHIFT;
            polyline[nsides].x1 = x1;
            polyline[nsides].y1 = YMIN - 1000.0;
            polyline[nsides].x2 = x1;
            polyline[nsides].y2 = YMIN + padding;
            polyline[nsides].angle = PID;
            nsides++;
            
            polyline[nsides].x1 = x1;
            polyline[nsides].y1 = YMAX - padding;
            polyline[nsides].x2 = x1;
            polyline[nsides].y2 = YMAX + 1000.0;
            polyline[nsides].angle = PID;
            nsides++;
            
            /* add circular arcs in corners */
            if (B_DOMAIN == D_POLYLINE_ARCS)
            {
                narcs = 0;
                
                for (i=0; i<NXMAZE; i++)
                    for (j=0; j<NYMAZE; j++)
                    {
                        n = nmaze(i, j);
                        x1 = YMIN + padding + (double)i*dx + MAZE_XSHIFT;
                        y1 = YMIN + padding + (double)j*dy;
                        
                        ra = MAZE_CORNER_RADIUS;
                        rb = 1.0 - MAZE_CORNER_RADIUS;
            
                        if (((i>0)||(j!=NYMAZE/2))&&(maze[n].south)&&(maze[n].west)) 
                        {
                            arcs[narcs].xc = x1 + ra*dx;
                            arcs[narcs].yc = y1 + ra*dy;
                            arcs[narcs].radius = ra*dx;
                            arcs[narcs].angle1 = PI;
                            arcs[narcs].dangle = PID;
                            narcs++;
                        }
                        
                        if (((i<NXMAZE-1)||(j!=NYMAZE/2-1))&&(maze[n].south)&&(maze[n].east)) 
                        {
                            arcs[narcs].xc = x1 + rb*dx;
                            arcs[narcs].yc = y1 + ra*dy;
                            arcs[narcs].radius = ra*dx;
                            arcs[narcs].angle1 = 3.0*PID;
                            arcs[narcs].dangle = PID;
                            narcs++;
                        }
                    
                        if (((i<NXMAZE-1)||(j!=NYMAZE/2-1))&&(maze[n].north)&&(maze[n].east)) 
                        {
                            arcs[narcs].xc = x1 + rb*dx;
                            arcs[narcs].yc = y1 + rb*dy;
                            arcs[narcs].radius = ra*dx;
                            arcs[narcs].angle1 = 0.0;
                            arcs[narcs].dangle = PID;
                            narcs++;
                        }
                    
                        if (((i>0)||(j!=NYMAZE/2))&&(maze[n].north)&&(maze[n].west)) 
                        {
                            arcs[narcs].xc = x1 + ra*dx;
                            arcs[narcs].yc = y1 + rb*dy;
                            arcs[narcs].radius = ra*dx;
                            arcs[narcs].angle1 = PID;
                            arcs[narcs].dangle = PID;
                            narcs++;
                        }
                    }
            }
    
            free(maze);
            break;
        }
        case (P_MAZE_DIAG):
        {
            maze = (t_maze *)malloc(NXMAZE*NYMAZE*sizeof(t_maze));
    
            init_maze_exit(0, NYMAZE/2, maze);
            
            /* build walls of maze */
            dx = (YMAX - YMIN - 2.0*padding)/(double)(NXMAZE);
            dy = (YMAX - YMIN - 2.0*padding)/(double)(NYMAZE);
            
            nsides = 0;
            ncircles = 0;
            
            /* deactivate sides for openings */
            maze[nmaze(0,NYMAZE/2)].west = 0;
    
            for (i=0; i<NXMAZE; i++)
                for (j=0; j<NYMAZE; j++)
                {
                    n = nmaze(i, j);
                    x1 = YMIN + padding + (double)i*dx + MAZE_XSHIFT;
                    y1 = YMIN + padding + (double)j*dy;
            
                    if (((i>0)||(j!=NYMAZE/2))&&(maze[n].west)) 
                    {
                        polyline[nsides].x1 = x1;
                        polyline[nsides].y1 = y1;
                        polyline[nsides].x2 = x1;
                        polyline[nsides].y2 = y1 + dy;
                        polyline[nsides].angle = PID;
                        nsides++;
                    }
                    
                    
                    if (maze[n].south) 
                    {
                        polyline[nsides].x1 = x1;
                        polyline[nsides].y1 = y1;
                        polyline[nsides].x2 = x1 + dx;
                        polyline[nsides].y2 = y1;
                        polyline[nsides].angle = 0.0;
                        nsides++;
                    }
                    
                }
                
            /* 45 degrees parts */
            for (i=0; i<NXMAZE; i++)
                for (j=0; j<NYMAZE; j++)
                {
                    n = nmaze(i, j);
                    x1 = YMIN + padding + (double)i*dx + MAZE_XSHIFT;
                    y1 = YMIN + padding + (double)j*dy;
            
                    if (((i>0)||(j!=NYMAZE/2))&&(maze[n].south)&&(maze[n].west)) 
                    {
                        polyline[nsides].x1 = x1 + 0.25*dx;
                        polyline[nsides].y1 = y1;
                        polyline[nsides].x2 = x1;
                        polyline[nsides].y2 = y1 + 0.25*dy;
                        polyline[nsides].angle = 1.5*PID;
                        nsides++;
                    }
                    
//                     if ((maze[n].south)&&(maze[n].east)) 
                    if (((i<NXMAZE-1)||(j!=NYMAZE/2-1))&&(maze[n].south)&&(maze[n].east)) {
                        polyline[nsides].x1 = x1 + 0.75*dx;
                        polyline[nsides].y1 = y1;
                        polyline[nsides].x2 = x1 + dx;
                        polyline[nsides].y2 = y1 + 0.25*dy;
                        polyline[nsides].angle = 0.5*PID;
                        nsides++;
                    }
                    
//                     if ((maze[n].north)&&(maze[n].east)) 
                    if (((i<NXMAZE-1)||(j!=NYMAZE/2-1))&&(maze[n].north)&&(maze[n].east)) 
                        {
                        polyline[nsides].x1 = x1 + dx;
                        polyline[nsides].y1 = y1 + 0.75*dy;
                        polyline[nsides].x2 = x1 + 0.75*dx;
                        polyline[nsides].y2 = y1 + dy;
                        polyline[nsides].angle = 1.5*PID;
                        nsides++;
                    }
                    
                    if (((i>0)||(j!=NYMAZE/2))&&(maze[n].north)&&(maze[n].west)) 
                    {
                        polyline[nsides].x1 = x1;
                        polyline[nsides].y1 = y1 + 0.75*dy;
                        polyline[nsides].x2 = x1 + 0.25*dx;
                        polyline[nsides].y2 = y1 + dy;
                        polyline[nsides].angle = 0.5*PID;
                        nsides++;
                    }
                }
    
            /* top side of maze */
            polyline[nsides].x1 = YMIN + padding + MAZE_XSHIFT;
            polyline[nsides].y1 = YMAX - padding;
            polyline[nsides].x2 = YMAX - padding + MAZE_XSHIFT;
            polyline[nsides].y2 = YMAX - padding;
            polyline[nsides].angle = 0.0;
            nsides++;
    
            /* right side of maze */
            y1 = YMIN + padding + dy*((double)NYMAZE/2);
            x1 = YMAX - padding + MAZE_XSHIFT;
            polyline[nsides].x1 = x1;
            polyline[nsides].y1 = YMIN - 1000.0;
            polyline[nsides].x2 = x1;
            polyline[nsides].y2 = y1 - dy;
            polyline[nsides].angle = PID;
            nsides++;
            
            polyline[nsides].x1 = x1;
            polyline[nsides].y1 = y1;
            polyline[nsides].x2 = x1;
            polyline[nsides].y2 = YMAX + 1000.0;
            polyline[nsides].angle = PID;
            nsides++;
    
            /* left side of maze */
            x1 = YMIN + padding + MAZE_XSHIFT;
            polyline[nsides].x1 = x1;
            polyline[nsides].y1 = YMIN - 1000.0;
            polyline[nsides].x2 = x1;
            polyline[nsides].y2 = YMIN + padding;
            polyline[nsides].angle = PID;
            nsides++;
            
            polyline[nsides].x1 = x1;
            polyline[nsides].y1 = YMAX - padding;
            polyline[nsides].x2 = x1;
            polyline[nsides].y2 = YMAX + 1000.0;
            polyline[nsides].angle = PID;
            nsides++;
    
            free(maze);
            break;
        }
        case (P_MAZE_RANDOM):
        {
            maze = (t_maze *)malloc(NXMAZE*NYMAZE*sizeof(t_maze));
            maze_coords = (double *)malloc(2*NXMAZE*NYMAZE*sizeof(double));
    
            init_maze_exit(0, NYMAZE/2, maze);
            
            /* build walls of maze */
            dx = (YMAX - YMIN - 2.0*padding)/(double)(NXMAZE);
            dy = (YMAX - YMIN - 2.0*padding)/(double)(NYMAZE);
            
            /* compute maze coordinates with random shift */
            for (i=0; i<NXMAZE; i++)
                for (j=0; j<NYMAZE; j++)
                {
                    n = nmaze(i, j);
                    maze_coords[2*n] = YMIN + padding + (double)i*dx + MAZE_XSHIFT;
                    maze_coords[2*n+1] = YMIN + padding + (double)j*dy;
                }  
                
            for (i=1; i<NXMAZE-1; i++)
                for (j=1; j<NYMAZE-1; j++)
                {
                    n = nmaze(i, j);
                    maze_coords[2*n] += MAZE_RANDOM_FACTOR*dx*rand()/RAND_MAX;
                    maze_coords[2*n+1] += MAZE_RANDOM_FACTOR*dy*rand()/RAND_MAX;
                }  
                
//             for (j=0; j<NYMAZE; j++) for (i=0; i<NXMAZE; i++)
//             {
//                 n = nmaze(i, j);
//                 printf("i=%i, j=%i, coord[%i]=%.2lg, coord[%i]=%.2lg\n", i, j, 2*n, maze_coords[2*n], 2*n+1, maze_coords[2*n+1]);
//             }
//             sleep(3);
            
            nsides = 0;
            ncircles = 0;
    
            for (i=0; i<NXMAZE; i++)
                for (j=0; j<NYMAZE; j++)
                {
                    n = nmaze(i, j);
                    x1 = maze_coords[2*n];
                    y1 = maze_coords[2*n+1];
                    xright = x1 + dx;
                    yright = y1;
                    xtop = x1;
                    ytop = y1 + dy;
//                     printf("i=%i, j=%i, x1=%.3lg, y1=%.3lg, x1+=%.3lg, y1+=%.3lg\n", i, j, x1, y1, xright, yright);
                    if (i<NXMAZE-1) 
                    {
                        xright = maze_coords[2*nmaze(i+1, j)];
                        yright = maze_coords[2*nmaze(i+1, j)+1];
                    }
                    if (j<NYMAZE-1) 
                    {
                        xtop = maze_coords[2*nmaze(i, j+1)];
                        ytop = maze_coords[2*nmaze(i, j+1)+1];
                    }
//                     printf("i=%i, j=%i, x1=%.3lg, y1=%.3lg, x1+=%.3lg, y1+=%.3lg\n\n", i, j, x1, y1, xright, yright);
            
                    if (((i>0)||(j!=NYMAZE/2))&&(maze[n].west)) 
                    {
                        polyline[nsides].x1 = x1;
                        polyline[nsides].y1 = y1;
                        polyline[nsides].x2 = xtop;
                        polyline[nsides].y2 = ytop;
                        polyline[nsides].angle = PID;
                        nsides++;
                    }
                    
                    if (maze[n].south) 
                    {
                        polyline[nsides].x1 = x1;
                        polyline[nsides].y1 = y1;
                        polyline[nsides].x2 = xright;
                        polyline[nsides].y2 = yright;
                        polyline[nsides].angle = 0.0;
                        nsides++;
                    }
                }
//             sleep(3);
    
            /* top side of maze */
            polyline[nsides].x1 = YMIN + padding + MAZE_XSHIFT;
            polyline[nsides].y1 = YMAX - padding;
            polyline[nsides].x2 = YMAX - padding + MAZE_XSHIFT;
            polyline[nsides].y2 = YMAX - padding;
            polyline[nsides].angle = 0.0;
            nsides++;
    
            /* right side of maze */
            y1 = YMIN + padding + dy*((double)NYMAZE/2);
            x1 = YMAX - padding + MAZE_XSHIFT;
            polyline[nsides].x1 = x1;
            polyline[nsides].y1 = YMIN - 1000.0;
            polyline[nsides].x2 = x1;
            polyline[nsides].y2 = y1 - dy;
            polyline[nsides].angle = PID;
            nsides++;
            
            polyline[nsides].x1 = x1;
            polyline[nsides].y1 = y1;
            polyline[nsides].x2 = x1;
            polyline[nsides].y2 = YMAX + 1000.0;
            polyline[nsides].angle = PID;
            nsides++;
    
            /* left side of maze */
            x1 = YMIN + padding + MAZE_XSHIFT;
            polyline[nsides].x1 = x1;
            polyline[nsides].y1 = YMIN - 1000.0;
            polyline[nsides].x2 = x1;
            polyline[nsides].y2 = YMIN + padding;
            polyline[nsides].angle = PID;
            nsides++;
            
            polyline[nsides].x1 = x1;
            polyline[nsides].y1 = YMAX - padding;
            polyline[nsides].x2 = x1;
            polyline[nsides].y2 = YMAX + 1000.0;
            polyline[nsides].angle = PID;
            nsides++;
            
//             ncircles = nsides;
//             for (i=0; i<ncircles; i++)
//             {
//                 circles[i].xc = polyline[i].x1;
//                 circles[i].yc = polyline[i].y1;
//                 circles[i].radius = MU;
//                 circles[i].active = 1;
//             }
    
            free(maze);
            free(maze_coords);
            break;
        }
        case (P_MAZE_CIRCULAR):
        {
            maze = (t_maze *)malloc(NXMAZE*NYMAZE*sizeof(t_maze));
            
            init_circular_maze(maze);
            
            init_polyline_circular_maze(polyline, circles, arcs, maze);
            
            free(maze);
            
            break;
        }
        case (P_MAZE_CIRC_SCATTERER):
        {
            
            maze = (t_maze *)malloc(NXMAZE*NYMAZE*sizeof(t_maze));
    
            init_circular_maze(maze);
            
            init_polyline_circular_maze(polyline, circles, arcs, maze);
            
            /* add scatterers */
            nblocks = NYMAZE/NXMAZE;
            rmin = 0.15;
            rmax = 1.0;
            angle = DPI/(double)nblocks;
            dphi = angle;
            
            dr = (rmax - rmin)/(double)(NXMAZE);
            rscat = 0.15*dr;
                        
            for (block = 0; block < nblocks; block++)
            {
                dphi = angle;
                
                /* first circle */
                r = rmin;
                phi = (double)block*angle;
                
                arcs[narcs].xc = rmin*cos(phi) + MAZE_XSHIFT;
                arcs[narcs].yc = rmin*sin(phi);
                arcs[narcs].radius = rscat;
                arcs[narcs].angle1 = 0.0;
                arcs[narcs].dangle = DPI;
                narcs++;
                
                /* second circle */
                dphi *= 0.5;
                for (q=0; q<2; q++)
                {
                    r = rmin + dr;
                    phi = (double)block*angle + (double)q*dphi;
                    
                    arcs[narcs].xc = r*cos(phi) + MAZE_XSHIFT;
                    arcs[narcs].yc = r*sin(phi);
                    arcs[narcs].radius = rscat;
                    arcs[narcs].angle1 = 0.0;
                    arcs[narcs].dangle = DPI;
                    narcs++;
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
                        for (q = 0; q < 2*ww; q++)
                        {
                            phi = (double)block*angle + (double)q*dphi;
                            
                            arcs[narcs].xc = r*cos(phi) + MAZE_XSHIFT;
                            arcs[narcs].yc = r*sin(phi);
                            arcs[narcs].radius = rscat;
                            arcs[narcs].angle1 = 0.0;
                            arcs[narcs].dangle = DPI;
                            narcs++;
                        }
                        i++;
                    }
                    ww *= 2;
                }
            }

            free(maze);
            break;
        }
        case (P_MAZE_HEX):
        {
            maze = (t_maze *)malloc(NXMAZE*NYMAZE*sizeof(t_maze));
            
            init_hex_maze(maze);
            
            nsides = 0;
            ncircles = 0;
    
            init_polyline_hex_maze(polyline, circles, maze);
            
            free(maze);
            
            break;
        }
        case (P_MAZE_OCT):
        {
            maze = (t_maze *)malloc(NXMAZE*NYMAZE*sizeof(t_maze));
            
            init_oct_maze(maze);
            
            nsides = 0;
            ncircles = 0;
    
            init_polyline_oct_maze(polyline, circles, maze);
            
            free(maze);
            
            break;
        }
    }
}


void init_polyline_depth(t_segment polyline[NMAXPOLY], t_circle circles[NMAXCIRCLES], int sdepth, double mu)
/* initialise polyline with variable depth parameter (for von Koch snowflake) */
{
    int i, j, k, l, n, z, ii, jj, terni[SDEPTH], ternj[SDEPTH], quater[SDEPTH], cond;
    short int vkoch[NMAXCIRCLES], turnright; 
    double ratio, omega, angle, sw, length, dist, x, y, ta, tb, a, b;
    
    switch (POLYLINE_PATTERN) {
        case (P_VONKOCH):
        {
            nsides = 3;
            for (k=0; k<sdepth; k++) nsides *= 4;
            printf("nsides = %i\n", nsides);
            ncircles = nsides;
            
            if (nsides > NMAXPOLY)
            {
                printf("NMAXPOLY has to be increased to at least %i\n", nsides);
                nsides = NMAXPOLY;
            }

            for (i=0; i<nsides/3; i++)
            {
                /* compute quaternary expansion of i */
                ii = i;
                for (l=0; l<sdepth; l++)
                {
                    quater[l] = ii%4;
                    ii = ii - (ii%4);
                    ii = ii/4;
                }
                
                /* find first nonzero digit */
                z = 0;
                while ((quater[z] == 0)&&(z<sdepth)) z++;
                
                /* compute left/right turns */
                if (i==0) vkoch[0] = 0;
                else if (z != sdepth)
                {   
                    if (quater[z] == 2) vkoch[i] = 0;
                    else vkoch[i] = 1;
                }
//                 printf("%i", vkoch[i]);
            }
            printf("\n");
            
            /* compute vertices */
            angle = 2.0*PI/3.0;
            x = cos(PI/6.0);
            y = -sin(PI/6.0);
            length = 2.0*sin(PI/3.0);
            
            for (k=0; k<sdepth; k++) length = length/3.0;
            printf("Length = %.2f\n", length);
            
            for (i=0; i<nsides; i++)
            {
                polyline[i].x1 = x;
                polyline[i].y1 = y; 
                polyline[i].angle = angle;
                
                x += length*cos(angle);
                y += length*sin(angle);
                polyline[i].length = length;
                                
                turnright = vkoch[i%(nsides/3)+1];
                if (turnright) angle -= PI/3.0;
                else angle += 2.0*PI/3.0;
                
                while (angle > DPI) angle -= DPI;
                while (angle < 0.0) angle += DPI;                
            }
            
            for (i=0; i<nsides; i++)
            {
                polyline[i].color = 0;
                if (i < nsides-1) polyline[i].x2 = polyline[i+1].x1;
                else polyline[i].x2 = polyline[0].x1;
                if (i < nsides-1) polyline[i].y2 = polyline[i+1].y1;
                else polyline[i].y2 = polyline[0].y1;  
            }
            
            for (i=0; i<ncircles; i++)
            {
                circles[i].xc = polyline[i].x1;
                circles[i].yc = polyline[i].y1;
                circles[i].radius = mu;
                circles[i].active = 1;
            }
            
            break;
        }
    }
}


int test_initial_condition(double *configs[NPARTMAX], int active[NPARTMAX], int color[NPARTMAX])
/* apply a test to initial conditions - so far, whether particle has left maze on the right */
{
    int i, j, time, nactive = 0, counter = 0, tmax[NPARTMAX], tmaxmax = 0, tmaxmin = 1000;
    double conf[8], newconf[8], cosphi, x2, pcolor;
    
    for (i=0; i<nparticles; i++)
    {
        for (j=0; j<8; j++) conf[j] = configs[i][j];
        for (j=0; j<8; j++) newconf[j] = configs[i][j];
        
        time = 0;
//         while ((time < 100000)&&(newconf[4] < 1000.0))
        while ((time < 100000)&&(newconf[0] < DUMMY_ABSORBING))
        {
            for (j=0; j<8; j++) conf[j] = newconf[j];
            vbilliard(newconf); 
            time++;
        }
        
        tmax[i] = time;
        
        printf("tmax = %i\n", time);
        
        cosphi = (conf[6] - conf[4])/conf[3];
        x2 = conf[4] + conf[2]*cosphi;
//         x2 = conf[4] + 10.0*cosphi;
        
        if ((conf[0] >= DUMMY_ABSORBING)&&(x2 > XMAX)) 
//         if (conf[0] >= DUMMY_ABSORBING)    
        {
            active[i] = 1;
            if (time > tmaxmax) tmaxmax = time;
            else if (time < tmaxmin) tmaxmin = time;
        }
        else active[i] = 0;
        
        nactive += active[i];
    }
    
    printf("%i particles, %i active particles\n", nparticles, nactive);
    printf("tmin = %i, tmax = %i\n", tmaxmin, tmaxmax);
//     sleep(5);
    
    /* reorder particles */
    for (i=0; i<nparticles; i++)
    {
        if (active[i])
        {
            for (j=0; j<8; j++) configs[counter][j] = configs[i][j];
            pcolor = sqrt((double)tmax[i] - (double)tmaxmin - 1.0)/sqrt((double)tmaxmax - (double)tmaxmin);
            color[counter] = (int)((double)NCOLORS*pcolor);
            active[counter] = 1;
            counter++;
        }
    }
    
    for (i=0; i<counter; i++) active[i] = 1;
    for (i=counter; i<nparticles; i++) active[i] = 0;

    return(counter);
}

int maze_type(int polyline_pattern)
/* group polyline pattern by type of maze */
{
    switch (POLYLINE_PATTERN) {
        case (P_MAZE): return(0);
        case (P_MAZE_DIAG): return(0);
        case (P_MAZE_CIRCULAR): return(1);
        case (P_MAZE_CIRC_SCATTERER): return(1);
        case (P_MAZE_HEX): return(2);
        case (P_MAZE_OCT): return(3);
        default: return(0);
    }
    
}

int find_maze_cell(double x, double y)
/* return maze cell number for coordinates (x, y) */
{
    int i, j, n, block, q, w, k, l;
    double padding = 0.02, x1, y1, u, v, r, phi, phi1, angle1, dphi, r1, tolerance = 0.025, x2, y2;
    static double dx, dy, rmin, rmax, angle, rr, rrtol, h, dr, dphi_table[5], x0, y0, a, b, d, a2, a2tol;
    static int first = 1, nblocks;
    
    if (first)
    {
        switch (maze_type(POLYLINE_PATTERN)) {
            case (0):   /* maze with square or rectangular cells */
            {
                dx = (YMAX - YMIN - 2.0*padding)/(double)(NXMAZE);
                dy = (YMAX - YMIN - 2.0*padding)/(double)(NYMAZE);
                break;
            }
            case (1):   /* circular maze */
            {
                nblocks = NYMAZE/NXMAZE;
                rmin = 0.15;
                rmax = 1.0;
                angle = DPI/(double)nblocks;
                dphi_table[0] = angle;
                for (i=1; i<5; i++) dphi_table[i] = 0.5*dphi_table[i-1];
        
                dr = (rmax - rmin)/(double)(NXMAZE);
                
                break;
            }
            case (2):   /* honeycomb maze */
            {
                rr = 2.0*(YMAX - YMIN)/(3.0*((double)(NXMAZE)-0.5));
                r1 = (YMAX - YMIN)/(sqrt(3.0)*(double)(NYMAZE+1));
                if (r1 < rr) rr = r1;
                h = 0.5*sqrt(3.0)*rr;
                rrtol = rr*(1.0 - tolerance);
                
                x0 = YMIN + padding + MAZE_XSHIFT;
                y0 = YMIN + padding + h;
                
                break;
            }
            case (3):   /* octahedron-square maze */
            {
                dx = (YMAX - YMIN - 2.0*padding)/((double)NYMAZE+0.5);    
                a = dx*(2.0 - sqrt(2.0));
                a2 = 0.5*a;
                b = a/sqrt(2.0);
                d = 2.0*dx;
    
                x0 = YMIN + padding + MAZE_XSHIFT - 0.5*b;
                y0 = YMIN + padding + 0.25*dx - 0.5*b;
                rr = sqrt((b+a2)*(b+a2) + a2*a2);
                rrtol = rr*(1.0 - tolerance);
                a2tol = a2*(1.0 - tolerance);
                
                break;
            }
        }
        first = 0;
    }
    
    switch (maze_type(POLYLINE_PATTERN)) {
        case (0):   /* maze with square or rectangular cells */
        {
            if (x < YMIN + MAZE_XSHIFT) return (-1);
            
            i = (int)((x - YMIN - padding - MAZE_XSHIFT)/dx);
            j = (int)((y - YMIN - padding)/dy);
            
            if (i < 0) return(-1);
            if (i > NXMAZE-1) return(-1);
            if (j < 0) return(-1);
            if (j > NYMAZE-1) return(-1);
            
            x1 = x - YMIN - padding - MAZE_XSHIFT - (double)i*dx;
            y1 = y - YMIN - padding - (double)j*dy;
            
            /* avoid finding wrong cell for particles close to wall */
            if (x1 < tolerance*dx) return(-10);
            if (x1 > dx - tolerance*dx) return(-10);
            if (y1 < tolerance*dy) return(-10);
            if (y1 > dy - tolerance*dy) return(-10);
             
            n = nmaze(i, j);
            
            if (n < NXMAZE*NYMAZE) return(n);
            else return(-1);
        }
        case (1):   /* circular maze */
        {
            r = module2(x - MAZE_XSHIFT,y);
            phi = argument(x - MAZE_XSHIFT,y);
            if (phi < 0.0) phi += DPI;
            if (phi >= DPI) phi -= DPI;
                        
            if (r < rmin*(1.0 - tolerance)) return(NXMAZE*NYMAZE);
            if (r < rmin*(1.0 + tolerance)) return(-10);
            if (r > rmax) return(-2);
            
            i = (int)((r - rmin)/dr);
            if (i > NXMAZE-1) return(-1);
            
            /* avoid finding wrong cell for particles close to wall */
            if (r - rmin - (double)i*dr < tolerance*dr) return(-10);
            if (r - rmin - (double)(i+1)*dr > -tolerance*dr) return(-10);
            
            block = (int)(phi/angle);
            phi1 = phi - (double)block*angle;
            
            /* avoid finding wrong cell for particles close to wall */
            if (phi1 < tolerance*angle)
                return(-10);
            if (phi1 - angle > -tolerance*angle)
                return(-10);
            
            if (i == 0) j = block*NXMAZE;
            else
            {
                w = 1 + (int)(log((double)i)/log(2.0));
//                 dphi = angle/ipow(2.0, w);
                dphi = dphi_table[w];
                q = (int)(phi1/dphi);
                j = block*NXMAZE + q;
                
                /* avoid finding wrong cell for particles close to wall */
                if (phi1 - (double)q*dphi < tolerance*dphi)
                    return(-1);
                if (phi1 - (double)(q+1)*dphi > - tolerance*dphi)
                    return(-1);
            }
            
            n = nmaze(i, j);
            
            if (n < NXMAZE*NYMAZE) return(n);
            else return(-1);
        }
        case (2):   /* honeycomb maze */
        {
             /* transformation to coordinate system along hexagon centers */
            u = (x-x0)/(1.5*rr);
            v = 0.5*(u + (y-y0)/h);
            
            i = (int)u;
            j = (int)v;
            
            j -= i/2;
            
            if (i < 0) return(-1);
            if (i > NXMAZE-1) return(-1);
            if (j < 0) return(-1);
            if (j > NYMAZE-1) return(-1);
            
            for (k=-1; k<2; k++) if ((i+k >= 0)&&(i+k < NXMAZE))
                for (l=-1; l<2; l++) if ((j+l >= 0)&&(j+l < NYMAZE))
                {
                    x1 = x0 + 1.5*(double)(i+k)*rr;
                    y1 = y0 + 2.0*(double)(j+l)*h;
                    if ((i+k)%2 == 0) y1 += h;
                    if (in_polygon(x - x1, y - y1, rrtol, 6, 0.0))
                    {
                        n = nmaze(i+k, j+l);
                        
                        if (i+k < 0) return(-1);
                        if (i+k > NXMAZE-1) return(-1);
                        if (j+l < 0) return(-1);
                        if (j+l > NYMAZE-1) return(-1);
                        
                        if (n < NXMAZE*NYMAZE) return(n);
                        else return(-1);
                    }
                }
           return(-10);
        }
        case (3):   /* octahedron-square maze */
        {
            i = 2*(int)((x-x0)/d);
            j = 2*(int)((y-y0)/d);
            
            x1 = x - x0 - (double)(i)*dx;
            y1 = y - y0 - (double)(j)*dx;
            
//             printf("(x,y) = (%.3lg,%.3lg), (i,j) = (%i,%i), (x1,y1) = (%.3lg,%.3lg)\n", x, y, i, j, x1/d, y1/d);
            
            if (i < 0) return(-1);
            if (i > NXMAZE-1) return(-1);
            if (j < 0) return(-1);
            if (j > NYMAZE-1) return(-1);
            
            if (x1 < b)
            {
                if (x1 + y1 < b) i--;
                else if (y1 - x1 > dx) i--; 
            }
            else if ((x1 > dx)&&(x1 < a + 2.0*b))
            {
                if (y1 - x1 < - dx) i++;
                else if (x1 + y1 > 3.0*a + 2.0*b) i++;
            }
            else if (x1 >= a + 2.0*b) i++; 
            
            if (y1 < b)
            {
                if (x1 + y1 < b) j--;
                else if (x1 - y1 > a + b) j--;
            }
            else if ((y1 > a + b)&&(y1 < a + 2.0*b))
            {
                if (x1 - y1 < - dx) j++;
                else if (x1 + y1 > 3.0*a + 2.0*b) j++;
            }
            else if (y1 >= a + 2.0*b) j++;
            
            if (i < 0) return(-1);
            if (i > NXMAZE-1) return(-1);
            if (j < 0) return(-1);
            if (j > NYMAZE-1) return(-1);

            
            /* avoid finding wrong cell for particles close to wall */
            x2 = x - x0 - (double)i*dx - b - a2;
            y2 = y - y0 - (double)j*dx - b - a2;
            if ((i+j)%2 == 0)
            {
                if (!in_polygon(x2, y2, rrtol, 8, 0.25)) return(-10);
            }
            else
            {
                if (vabs(x2) > a2tol) return(-10);
                if (vabs(y2) > a2tol) return(-10);
            }
            
            n = nmaze(i, j);
            
            if (n < NXMAZE*NYMAZE) return(n);
            else return(-1);
        }
    }
}

void draw_maze_cell(int n, int part_number, double minprop)
/* give specified color to maze cell */
{
    int i, j, block, q;
    double x, y, padding = 0.02, rgb[3], w;
    static double dx, dy, rmin, rmax, angle, r, r1, dr, phi1, dphi, h, x0, y0, a, b, a2, square_prop;
    static int first = 1, nblocks;
    
    if (first) switch (maze_type(POLYLINE_PATTERN)) {
        case (0):   /* maze with square or rectangular cells */
        {
            dx = (YMAX - YMIN - 2.0*padding)/(double)(NXMAZE);
            dy = (YMAX - YMIN - 2.0*padding)/(double)(NYMAZE);
            break;
        }
        case (1):   /* circular maze */
        {
            nblocks = NYMAZE/NXMAZE;
            rmin = 0.15;
            rmax = 1.0;
            angle = DPI/(double)nblocks;
    
            dr = (rmax - rmin)/(double)(NXMAZE);
            
            break;
        }
        case (2):   /* honeycomb maze */
        {
            r = 2.0*(YMAX - YMIN)/(3.0*((double)(NXMAZE)-0.5));
            r1 = (YMAX - YMIN)/(sqrt(3.0)*(double)(NYMAZE+1));
            if (r1 < r) r = r1;
            h = 0.5*sqrt(3.0)*r;                
            x0 = YMIN + padding + MAZE_XSHIFT;
            y0 = YMIN + padding + h;
            break;
        }
        case (3):   /* octahedron-square maze */
        {
            dx = (YMAX - YMIN - 2.0*padding)/((double)NYMAZE+0.5);    
            a = dx*(2.0 - sqrt(2.0));
            b = a/sqrt(2.0);
            a2 = 0.5*a;
            
            r = sqrt((b+a2)*(b+a2) + a2*a2);
            square_prop = 2.0*(sqrt(2.0) + 1.0); 

            x0 = YMIN + padding + MAZE_XSHIFT;
            y0 = YMIN + padding + 0.25*dx;
                
            break;
        }
        first = 0;
    }
    
    if (part_number != 0) switch (maze_type(POLYLINE_PATTERN)) {
        case (0):   /* maze with square or rectangular cells */
        {
            i = n%NXMAZE;
            j = n/NXMAZE;
            
            x = YMIN + padding + (double)i*dx + MAZE_XSHIFT;
            y = YMIN + padding + (double)j*dy;
            
            rgb_color_scheme_density(part_number, rgb, minprop); 
            erase_rectangle(x, y, x + dx, y + dy, rgb);
            break;
        }
        case (1):   /* circular maze */
        {
            if (n == NXMAZE*NYMAZE) /* inner circle */
            {
                rgb_color_scheme_density(part_number, rgb, minprop); 
                draw_filled_sector(MAZE_XSHIFT, 0.0, 0.0, rmin, 0.0, DPI, NSEG, rgb);
                break;
            }
            
            i = n%NXMAZE;
            j = n/NXMAZE;
            block = j/NXMAZE;
            
            r = rmin + (double)i*dr;
            r1 = r + dr;
                if (i == 0)    /* first circle */
            {
                dphi = angle;
                phi1 = (double)block*dphi;
            }
            else                /* other circles */
            {
                w = 1 + (int)(log((double)i)/log(2.0));
                dphi = angle*ipow(0.5, w);
                q = j - block*NXMAZE;
                phi1 = (double)block*angle + (double)q*dphi;
            }
            rgb_color_scheme_density(part_number, rgb, minprop); 
            draw_filled_sector(MAZE_XSHIFT, 0.0, r, r1, phi1, dphi, NSEG, rgb);
            break;
        }
        case (2):   /* honeycomb maze */
        {
            i = n%NXMAZE;
            j = n/NXMAZE;
            
            x = x0 + 1.5*(double)i*r;
            y = y0 + 2.0*(double)j*h;
            if (i%2 == 0) y += h;
            
            rgb_color_scheme_density(part_number, rgb, minprop); 
            draw_colored_circle(x, y, r, 6, rgb);
            break;
        }
        case (3):   /* octahedron-square maze */
        {
            i = n%NXMAZE;
            j = n/NXMAZE;
            
            x = x0 + ((double)i + 0.5)*dx;
            y = y0 + ((double)j + 0.5)*dx;
            
            /* adapt particle number in square to have constant density */
            if ((i+j)%2 == 1) part_number = (int)(square_prop*(double)part_number);
            
            rgb_color_scheme_density(part_number, rgb, minprop); 
            if ((i+j)%2 == 0) draw_colored_octahedron(x, y, r,rgb);
            else erase_rectangle(x-a2, y-a2, x+a2, y+a2, rgb);
            break;
        }
    }
        
}


/* can probably be removed */
void print_left_right_part_number_old(double *configs[NPARTMAX], int active[NPARTMAX], double xl, double yl, double xr, double yr, double xmin, double xmax)
{
    char message[50];
    int i, nleft = 0, nright = 0;
    double rgb[3], x1, y1, cosphi, sinphi;
    static int maxnleft = 0, maxnright = 0;
    
    /* count active particles, using the fact that absorbed particles have been given dummy coordinates */
    for (i=0; i<nparticles; i++) if ((configs[i][0] >= DUMMY_ABSORBING))
    {
        cosphi = (configs[i][6] - configs[i][4])/configs[i][3];
        sinphi = (configs[i][7] - configs[i][5])/configs[i][3];
        x1 = configs[i][4] + configs[i][2]*cosphi;
        y1 = configs[i][5] + configs[i][2]*cosphi;
        if (x1 < xmin) nleft++;
//         else if ((x1 > xmax)&&(y1 <= YMAX)&&(y1 >= YMIN))
        else if ((x1 > xmax))
        {
            nright++;
            printf("Detected leaving particle %i at (%.2f, %2f)\n", i, x1, y1);
        }
    }
    if (nleft > maxnleft) maxnleft = nleft;
    if (nright > maxnright) maxnright = nright;
        
    hsl_to_rgb(0.0, 0.0, 0.0, rgb);
    
    erase_area(xl, yl - 0.03, 0.45, 0.12, rgb);
    erase_area(xr, yr - 0.03, 0.35, 0.12, rgb);
    
    glColor3f(1.0, 1.0, 1.0);
    if (maxnleft > 1) sprintf(message, "%i particles", maxnleft);
    else sprintf(message, "%i particle", maxnleft);
    write_text(xl, yl, message);
    if (maxnright > 1) sprintf(message, "%i particles", maxnright);
    else sprintf(message, "%i particle", maxnright);
    write_text(xr, yr, message);
}
