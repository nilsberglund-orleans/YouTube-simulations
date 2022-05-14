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

void draw_polygon(double x, double y, double r, int nsides, double angle)
{
    int i;
    double pos[2], alpha, dalpha;
    
    dalpha = DPI/(double)nsides;
    
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
                    particles[n].xc = ((double)(i-NGRIDX/2) + 0.5)*dx;   /* is +0.5 needed? */
                    particles[n].yc = INITYMIN + ((double)j - 0.5)*dy;
                    if ((i+NGRIDX)%2 == 1) particles[n].yc += 0.5*dy;
                    particles[n].radius = MU;
                    /* activate only circles that intersect the domain */
                    if ((particles[n].yc < INITYMAX + MU)&&(particles[n].yc > INITYMIN - MU)) particles[n].active = 1;
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

void init_obstacle_config(t_obstacle obstacle[NMAXOBSTACLES])
/* initialise circular obstacle configuration */
{
    int i, j, n; 
    double x, y, dx, dy;
    
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
        default: 
        {
            printf("Function init_obstacle_config not defined for this pattern \n");
        }
    }
}
 
 
void init_segment_config(t_segment segment[NMAXSEGMENTS])
/* initialise linear obstacle configuration */
{
    int i, cycle = 0, iminus, iplus; 
    double angle, angle2, dx, width, height, a, b, length;
    
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
            
            for (i=0; i<nsegments; i++) segment[i].concave = 0;
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
            
            for (i=0; i<nsegments; i++) segment[i].concave = 0;
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
            
            for (i=0; i<nsegments; i++) segment[i].concave = 0;
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
    
    /* add one segment for S_POLY_ELLIPSE configuration */
    if (SEGMENT_PATTERN == S_POLY_ELLIPSE)
    {
        segment[nsegments].x1 = -cos(((double)i + 0.5)*angle);
        segment[nsegments].y1 = LAMBDA*sin(((double)i + 0.5)*angle);
        segment[nsegments].x2 = -cos(((double)i + 0.5)*angle);
        segment[nsegments].y2 = -LAMBDA*sin(((double)i + 0.5)*angle);
        nsegments++;
    }
    
    /* activate all segments */
    for (i=0; i<nsegments; i++) segment[i].active = 1;
    
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
    }
    
    /* deal with concave corners */
    for (i=0; i<nsegments; i++) if (segment[i].concave)
        {
            iminus = i-1;  if (iminus == 0) iminus = nsegments - 1;
            iplus = i+1;  if (iplus == nsegments - 1) iplus = 0;
            angle = argument(segment[iplus].x1 - segment[i].x1, segment[iplus].y1 - segment[i].y1) + PID;
            angle2 = argument(segment[i].x1 - segment[iminus].x1, segment[i].y1 - segment[iminus].y1) + PID;
            if (angle2 < angle) angle2 += DPI;
            segment[i].angle1 = angle;
            segment[i].angle2 = angle2;
        }
    
    /* make copy of initial values in case of rotation/translation */
    if ((ROTATE_BOUNDARY)||(MOVE_BOUNDARY)) for (i=0; i<nsegments; i++) 
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
    
    for (i=0; i<nsegments; i++) 
    {
        printf("Segment %i: (x1, y1) = (%.3lg,%.3lg), (x2, y2) = (%.3lg,%.3lg)\n (nx, ny) = (%.3lg,%.3lg), c = %.3lg, length = %.3lg\n", i, segment[i].x1, segment[i].y1, segment[i].x2, segment[i].y2, segment[i].nx, segment[i].ny, segment[i].c, segment[i].length);
        if (segment[i].concave) printf("Concave with angles %.3lg Pi, %.3lg Pi\n", segment[i].angle1/PI, segment[i].angle2/PI);
    }
}


int in_segment_region(double x, double y)
/* returns 1 if (x,y) is inside region delimited by obstacle segments */
{
    double angle, dx, height, width, theta;
    
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
 
void translate_segments(t_segment segment[NMAXSEGMENTS], double deltax, double deltay)
/* rotates the repelling segments by given angle */
{
    int i;
    
    for (i=0; i<nsegments; i++) 
    {
        segment[i].x1 = segment[i].x01 + deltax;
        segment[i].y1 = segment[i].y01 + deltay;
        segment[i].x2 = segment[i].x02 + deltax;
        segment[i].y2 = segment[i].y02 + deltay;
        segment[i].c = segment[i].nx*segment[i].x1 + segment[i].ny*segment[i].y1;
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

int compute_particle_interaction(int i, int k, double force[2], double *torque, t_particle* particle, 
                                 double distance, double krepel, double ca, double sa, double ca_rel, double sa_rel)
/* compute repelling force and torque of particle #k on particle #i */
/* returns 1 if distance between particles is smaller than NBH_DIST_FACTOR*MU */
{
    double x1, y1, x2, y2, r, f, angle, aniso, fx, fy, ff[2], dist_scaled, spin_f, ck, sk, ck_rel, sk_rel;
    static double dxhalf = 0.5*(BCXMAX - BCXMIN), dyhalf = 0.5*(BCYMAX - BCYMIN);
    int wwrapx, wwrapy;
    
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
            default: 
            {
                spin_f = particle[i].spin_freq;
                if (wwrapx||wwrapy) *torque = sin(spin_f*(-particle[k].angle - particle[i].angle))/(1.0e-8 + dist_scaled*dist_scaled);
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
    if (distance > distmin) 
    {
        printf("i = %i, hashcell %i, j = %i, hashcell %i\n", i, particle[i].hashcell, k, particle[k].hashcell);
        printf("X = (%.3lg, %.3lg)\n", particle[i].xc, particle[i].yc);
        printf("Y = (%.3lg, %.3lg) d = %.3lg\n", particle[k].xc, particle[k].yc, distance);
    }
        
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

        if (RANDOM_RADIUS) particle[i].radius = particle[i].radius*(0.75 + 0.5*((double)rand()/RAND_MAX));
        
        particle[i].mass_inv = 1.0/mass;
        if (particle[i].type == 0) particle[i].inertia_moment_inv = 1.0/PARTICLE_INERTIA_MOMENT;
        else particle[i].inertia_moment_inv = 1.0/PARTICLE_INERTIA_MOMENT;
                
        particle[i].vx = vx;
        particle[i].vy = vy;
        particle[i].energy = (particle[i].vx*particle[i].vx + particle[i].vy*particle[i].vy)*particle[i].mass_inv;
        
        particle[i].angle = DPI*(double)rand()/RAND_MAX;
        particle[i].omega = 0.0;
        
        if (particle[i].type == 1)
        {
            particle[i].interaction = INTERACTION_B;
            particle[i].eq_dist = EQUILIBRIUM_DIST_B;
            particle[i].spin_range = SPIN_RANGE_B;
            particle[i].spin_freq = SPIN_INTER_FREQUENCY_B;            
        }
    
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

void draw_one_particle(t_particle particle, double xc, double yc, double radius, double angle, int nsides, double width, double rgb[3])
/* draw one of the particles */ 
{
    double ca, sa, x1, x2, y1, y2, xc1, wangle;
    int wsign;
    
    if (CENTER_VIEW_ON_OBSTACLE) xc1 = xc - xshift;
    else xc1 = xc;
    glColor3f(rgb[0], rgb[1], rgb[2]);
    if ((particle.interaction == I_LJ_QUADRUPOLE)||(particle.interaction == I_LJ_DIPOLE)) 
        draw_colored_rhombus(xc1, yc, radius, angle + APOLY*PID, rgb);
    else draw_colored_polygon(xc1, yc, radius, nsides, angle + APOLY*PID, rgb);
        
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
            y1 = yc - MU_B*sa;
            x2 = xc1 + MU_B*ca;
            y2 = yc + MU_B*sa;
            draw_line(x1, y1, x2, y2);
            x1 = xc1 - MU_B*sa;
            y1 = yc + MU_B*ca;
            x2 = xc1 + MU_B*sa;
            y2 = yc - MU_B*ca;
            draw_line(x1, y1, x2, y2);
        }
        
    glLineWidth(width);
    glColor3f(1.0, 1.0, 1.0);
    if ((particle.interaction == I_LJ_QUADRUPOLE)||(particle.interaction == I_LJ_DIPOLE)) 
        draw_rhombus(xc1, yc, radius, angle + APOLY*PID);
    else draw_polygon(xc1, yc, radius, nsides, angle + APOLY*PID); 
    
    if (particle.interaction == I_LJ_WATER) for (wsign = -1; wsign <= 1; wsign+=2)
    {
        wangle = particle.angle + (double)wsign*DPI/3.0;
        x1 = xc1 + particle.radius*cos(wangle);
        y1 = yc + particle.radius*sin(wangle);
        draw_colored_polygon(x1, y1, 0.5*radius, nsides, angle + APOLY*PID, rgb);
        glColor3f(1.0, 1.0, 1.0);
        draw_polygon(x1, y1, 0.5*radius, nsides, angle + APOLY*PID);
    }
}


void draw_trajectory(t_tracer trajectory[TRAJECTORY_LENGTH*N_TRACER_PARTICLES], int traj_position, int traj_length)
/* draw tracer particle trajectory */
{
    int i, j, time;
    double x1, x2, y1, y2, rgb[3], lum;
    
    blank();
    glLineWidth(TRAJECTORY_WIDTH);
    printf("drawing trajectory\n");
    
    for (j=0; j<N_TRACER_PARTICLES; j++)
    {
        if (j == 0) hsl_to_rgb(HUE_TYPE1, 0.9, 0.5, rgb);
        else if (j == 1) hsl_to_rgb(HUE_TYPE2, 0.9, 0.5, rgb);
        else hsl_to_rgb(HUE_TYPE3, 0.9, 0.5, rgb);
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
    

void draw_particles(t_particle particle[NMAXCIRCLES], int plot)
{
    int i, j, k, m, width, nnbg, nsides;
    double ej, hue, huex, huey, rgb[3], rgbx[3], rgby[3], radius, x1, y1, x2, y2, angle, ca, sa, length, linkcolor, sign = 1.0, angle1, signy = 1.0, periodx, periody, x, y, lum;
    char message[100];
    
    if (!TRACER_PARTICLE) blank();
    glColor3f(1.0, 1.0, 1.0);
    
    /* draw the bonds first */
    if (plot == P_BONDS)
    {
        glLineWidth(LINK_WIDTH);
        for (j=0; j<ncircles; j++) if (particle[j].active) draw_one_particle_links(particle[j]);
    }
                
    /* determine particle color and size */
    for (j=0; j<ncircles; j++) if (particle[j].active)
    {
        switch (plot) {
            case (P_KINETIC): 
            {
                ej = particle[j].energy;
                if (ej > 0.0) 
                {
                    hue = ENERGY_HUE_MIN + (ENERGY_HUE_MAX - ENERGY_HUE_MIN)*ej/PARTICLE_EMAX;
                    if (hue > ENERGY_HUE_MIN) hue = ENERGY_HUE_MIN;
                    if (hue < ENERGY_HUE_MAX) hue = ENERGY_HUE_MAX;
                }
                radius = particle[j].radius;
                width = BOUNDARY_WIDTH;
                break;
            }
            case (P_NEIGHBOURS): 
            {
                hue = neighbour_color(particle[j].neighb);
                radius = particle[j].radius;
                width = BOUNDARY_WIDTH;
                break;
            }
            case (P_BONDS):
            {
//                 if (particle[j].type == 1) hue = 70.0;        /* to make second particle type more visible */
//                 if (particle[j].type == 1) hue = neighbour_color(7 - particle[j].neighb);
//                 else 
                hue = neighbour_color(particle[j].neighb);
                radius = particle[j].radius;
                width = 1;
                break;
            }
            case (P_ANGLE):
            {
                angle = particle[j].angle;
                hue = angle*particle[j].spin_freq/DPI;
                hue -= (double)((int)hue);
                huex = (DPI - angle)*particle[j].spin_freq/DPI;
                huex -= (double)((int)huex);
                angle = PI - angle;
                if (angle < 0.0) angle += DPI;
                huey = angle*particle[j].spin_freq/DPI;
                huey -= (double)((int)huey);
                hue = PARTICLE_HUE_MIN + (PARTICLE_HUE_MAX - PARTICLE_HUE_MIN)*hue;
                huex = PARTICLE_HUE_MIN + (PARTICLE_HUE_MAX - PARTICLE_HUE_MIN)*huex;
                huey = PARTICLE_HUE_MIN + (PARTICLE_HUE_MAX - PARTICLE_HUE_MIN)*huey;
                radius = particle[j].radius;
                width = BOUNDARY_WIDTH;
                break;
            }
            case (P_TYPE):
            {
//                 if (particle[j].type == 0) hue = 310.0;
//                 else hue = 70.0;
                if (particle[j].type <= 1) hue = HUE_TYPE0;
                else if (particle[j].type == 2) hue = HUE_TYPE1;
                else if (particle[j].type == 3) hue = HUE_TYPE2;
                else hue = HUE_TYPE3;
                radius = particle[j].radius;
                width = BOUNDARY_WIDTH;
                break;
            }
            case (P_DIRECTION): 
            {
                hue = argument(particle[j].vx, particle[j].vy);
                if (hue > DPI) hue -= DPI;
                if (hue < 0.0) hue += DPI;
                hue = PARTICLE_HUE_MIN + (PARTICLE_HUE_MAX - PARTICLE_HUE_MIN)*hue/DPI;
                radius = particle[j].radius;
                width = BOUNDARY_WIDTH;
                break;
            }
            case (P_DIRECT_ENERGY): 
            {
                hue = argument(particle[j].vx, particle[j].vy);
                if (hue > DPI) hue -= DPI;
                if (hue < 0.0) hue += DPI;
                hue = PARTICLE_HUE_MIN + (PARTICLE_HUE_MAX - PARTICLE_HUE_MIN)*hue/DPI;
                if (particle[j].energy < 0.1*PARTICLE_EMAX) lum = 10.0*particle[j].energy/PARTICLE_EMAX;
                else lum = 1.0;
                radius = particle[j].radius;
                width = BOUNDARY_WIDTH;
                break;
            }
            case (P_ANGULAR_SPEED): 
            {
                hue = 160.0*(1.0 + tanh(SLOPE*particle[j].omega));
//                 printf("omega = %.3lg, hue = %.3lg\n", particle[j].omega, hue);
                radius = particle[j].radius;
                width = BOUNDARY_WIDTH;
                break;
            }
            case (P_DIFF_NEIGHB): 
            {
                hue = (double)(particle[j].diff_neighb+1)/(double)(particle[j].neighb+1);
//                 hue = PARTICLE_HUE_MIN + (PARTICLE_HUE_MAX - PARTICLE_HUE_MIN)*hue;
                hue = 180.0*(1.0 + hue);
                radius = particle[j].radius;
                width = BOUNDARY_WIDTH;
                break;
            }
        }

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

        switch (plot) {
            case (P_KINETIC):  
            {
                hsl_to_rgb_turbo(hue, 0.9, 0.5, rgb);
                hsl_to_rgb_turbo(hue, 0.9, 0.5, rgbx);
                hsl_to_rgb_turbo(hue, 0.9, 0.5, rgby);
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
                hsl_to_rgb_twilight(hue, 0.9, 0.5, rgb);
                hsl_to_rgb_twilight(hue, 0.9, 0.5, rgbx);
                hsl_to_rgb_twilight(hue, 0.9, 0.5, rgby);
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
            default: 
            {
                hsl_to_rgb(hue, 0.9, 0.5, rgb);
                hsl_to_rgb(hue, 0.9, 0.5, rgbx);
                hsl_to_rgb(hue, 0.9, 0.5, rgby);
            }
        }
        angle = particle[j].angle + APOLY*DPI;
        
        draw_one_particle(particle[j], particle[j].xc, particle[j].yc, radius, angle, nsides, width, rgb);
                
        /* in case of periodic b.c., draw translates of particles */
        if (PERIODIC_BC)
        {
            x1 = particle[j].xc;
            y1 = particle[j].yc;
            
            for (i=-2; i<3; i++)
                for (k=-1; k<2; k++)
                    draw_one_particle(particle[j], x1 + (double)i*(BCXMAX - BCXMIN), y1 + (double)k*(BCYMAX - BCYMIN), radius,          angle, nsides, width, rgb);
        }
        else if (BOUNDARY_COND == BC_KLEIN)
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
        else if (BOUNDARY_COND == BC_BOY)
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
        else if (BOUNDARY_COND == BC_GENUS_TWO)
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
            draw_colored_circle(obstacle[i].xc, obstacle[i].yc, obstacle[i].radius, NSEG, rgb);
        glColor3f(1.0, 1.0, 1.0);
        for (i = 0; i < nobstacles; i++)
            draw_circle(obstacle[i].xc, obstacle[i].yc, obstacle[i].radius, NSEG);
    }
    if (ADD_FIXED_SEGMENTS)
    {
        glLineWidth(CONTAINER_WIDTH);
        glColor3f(1.0, 1.0, 1.0);
        for (i = 0; i < nsegments; i++) if (segment[i].active)
            draw_line(segment[i].x1, segment[i].y1, segment[i].x2, segment[i].y2);
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
                
                draw_colored_circle(x, 0.0, OBSTACLE_RADIUS, NSEG, rgb);
                glColor3f(1.0, 1.0, 1.0);
                draw_circle(x, 0.0, OBSTACLE_RADIUS, NSEG);
                
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
                
                draw_colored_circle(x, 0.0, OBSTACLE_RADIUS, NSEG, rgb);
                glColor3f(1.0, 1.0, 1.0);
                draw_circle(x, 0.0, OBSTACLE_RADIUS, NSEG);
                
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
                    draw_colored_circle(x, (double)j*(FUNNEL_WIDTH + OBSTACLE_RADIUS), OBSTACLE_RADIUS, NSEG, rgb);
                    glColor3f(1.0, 1.0, 1.0);
                    draw_circle(x, (double)j*(FUNNEL_WIDTH + OBSTACLE_RADIUS), OBSTACLE_RADIUS, NSEG);
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
        }
    }
    
}

void print_parameters(double beta, double temperature, double krepel, double lengthcontainer, double boundary_force, 
                      short int left, double pressure[N_PRESSURES], double gravity)
{
    char message[100];
    int i, j, k;
    double density, hue, rgb[3], logratio, x, y, meanpress[N_PRESSURES], phi, sphi, dphi, pprint, mean_temp;
    static double xbox, xtext, xmid, xmidtext, xxbox, xxtext, pressures[N_P_AVERAGE], meanpressure = 0.0, maxpressure = 0.0;
    static double press[N_PRESSURES][N_P_AVERAGE], temp[N_T_AVERAGE];
    static int first = 1, i_pressure, i_temp;
    
    if (first)
    {
        if (left)
        {
            xbox = XMIN + 0.4;
            xtext = XMIN + 0.08;
            xxbox = XMAX - 0.39;
            xxtext = XMAX - 0.73;
       }
        else
        {
            xbox = XMAX - 0.41;
            xtext = XMAX - 0.73;
            xxbox = XMIN + 0.4;
            xxtext = XMIN + 0.08;
        }
        xmid = 0.5*(XMIN + XMAX) - 0.1;
        xmidtext = xmid - 0.24;
        for (i=0; i<N_P_AVERAGE; i++) pressures[i] = 0.0;
        if (RECORD_PRESSURES) for (j=0; j<N_PRESSURES; j++) 
        {
            meanpress[j] = 0.0;
            for (i=0; i<N_P_AVERAGE; i++) press[j][i] = 0.0;
        }
        i_pressure = 0;
        i_temp = 0;
        for (i=0; i<N_T_AVERAGE; i++) temp[i] = 0.0;
        
        first = 0;
    }
    
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
    
    y = YMAX - 0.1;
    if (INCREASE_BETA)  /* print temperature */
    {
        logratio = log(beta/BETA)/log(0.5*BETA_FACTOR);
        if (logratio > 1.0) logratio = 1.0;
        if (BETA_FACTOR > 1.0) hue = PARTICLE_HUE_MAX - (PARTICLE_HUE_MAX - PARTICLE_HUE_MIN)*logratio;
        else hue = PARTICLE_HUE_MIN - (PARTICLE_HUE_MIN - PARTICLE_HUE_MAX)*logratio;
        erase_area_hsl_turbo(xbox, y + 0.025, 0.37, 0.05, hue, 0.9, 0.5);
//         erase_area_hsl_turbo(xmid + 0.1, y + 0.025, 0.4, 0.05, hue, 0.9, 0.5);
        if ((hue < 90)||(hue > 270)) glColor3f(1.0, 1.0, 1.0);
        else glColor3f(0.0, 0.0, 0.0);
        sprintf(message, "Temperature %.2f", 1.0/beta);
        write_text(xtext, y, message);
//         write_text(xmidtext, y, message);
//         y -= 0.1;
        
//         erase_area_hsl(xxbox, y + 0.025, 0.37, 0.05, 0.0, 0.9, 0.0);
//         glColor3f(1.0, 1.0, 1.0);
//         sprintf(message, "Pressure %.3f", meanpressure);
//         write_text(xxtext, y, message);
    }
    if (DECREASE_CONTAINER_SIZE)  /* print density */
    {
        density = (double)ncircles/((lengthcontainer)*(INITYMAX - INITYMIN));
        erase_area_hsl(xbox, y + 0.025, 0.37, 0.05, 0.0, 0.9, 0.0);
        glColor3f(1.0, 1.0, 1.0);
        sprintf(message, "Density %.3f", density);
        write_text(xtext, y, message);
        
        erase_area_hsl(xmid, y + 0.025, 0.37, 0.05, 0.0, 0.9, 0.0);
        glColor3f(1.0, 1.0, 1.0);
        sprintf(message, "Temperature %.2f", temperature);
        write_text(xmidtext, y, message);

        erase_area_hsl(xxbox, y + 0.025, 0.37, 0.05, 0.0, 0.9, 0.0);
        glColor3f(1.0, 1.0, 1.0);
        sprintf(message, "Pressure %.3f", meanpressure);
        write_text(xxtext, y, message);

    }   
    else if (INCREASE_KREPEL)  /* print force constant */
    {
        erase_area_hsl(xbox, y + 0.025, 0.22, 0.05, 0.0, 0.9, 0.0);
        glColor3f(1.0, 1.0, 1.0);
        sprintf(message, "Force %.0f", krepel);
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
    
    if ((PARTIAL_THERMO_COUPLING)&&(!INCREASE_BETA))
    {
        printf("Temperature %i in average: %.3lg\n", i_temp, temperature);
        temp[i_temp] = temperature;
        i_temp++;
        if (i_temp >= N_T_AVERAGE) i_temp = 0;
        
        mean_temp = 0.0;
        for (i=0; i<N_T_AVERAGE; i++) mean_temp += temp[i];
        mean_temp = mean_temp/N_T_AVERAGE;
        
        hue = PARTICLE_HUE_MIN + 0.5*(PARTICLE_HUE_MAX - PARTICLE_HUE_MIN)*mean_temp/PARTICLE_EMAX;
        if (hue < PARTICLE_HUE_MAX) hue = PARTICLE_HUE_MAX;
        erase_area_hsl_turbo(xbox, y + 0.025, 0.37, 0.05, hue, 0.9, 0.5);
        if ((hue < 90)||(hue > 270)) glColor3f(1.0, 1.0, 1.0);
        else glColor3f(0.0, 0.0, 0.0);
        sprintf(message, "Temperature %.2f", mean_temp);
        write_text(xtext, y, message);
    }
    
    if (INCREASE_GRAVITY)
    {
        erase_area_hsl(xmid, y + 0.025, 0.22, 0.05, 0.0, 0.9, 0.0);
        glColor3f(1.0, 1.0, 1.0);
        sprintf(message, "Gravity %.2f", gravity/GRAVITY);
        write_text(xmidtext + 0.1, y, message);
    }   
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

void print_omega(double angular_speed)
{
    char message[100];
    double rgb[3];
    static double xleftbox, xlefttext, xrightbox, xrighttext, y = YMAX - 0.1, ymin = YMIN + 0.05;
    static int first = 1;
    
    if (first)
    {
//         xleftbox = XMIN + 0.4;
//         xlefttext = xleftbox - 0.55;
        xrightbox = XMAX - 0.39;
        xrighttext = xrightbox - 0.55;
        first = 0;
    }
    

    erase_area_hsl(xrightbox, y + 0.025, 0.35, 0.05, 0.0, 0.9, 0.0);
    glColor3f(1.0, 1.0, 1.0);
    sprintf(message, "Angular speed = %.4f", DPI*angular_speed*25.0/(double)(PERIOD_ROTATE_BOUNDARY));
    write_text(xrighttext + 0.1, y, message);

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
    double xmin, xmax, ymin, ymax, padding, r, rp, r2, cphi, sphi, f, fperp = 0.0, x, y, xtube, distance, dx, dy, 
        width, ybin, angle, x1, x2, h, ytop, norm, dleft, dplus, dminus, tmp_pleft = 0.0, tmp_pright = 0.0, proj;
    
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
    if (ADD_FIXED_SEGMENTS) for (i=0; i<nsegments; i++) if (segment[i].active)
    {
        x = particle[j].xc;
        y = particle[j].yc;
        proj = (segment[i].ny*(x - segment[i].x1) - segment[i].nx*(y - segment[i].y1))/segment[i].length;
        if ((proj > 0.0)&&(proj < 1.0))
        {
            distance = segment[i].nx*x + segment[i].ny*y - segment[i].c;
            r = 1.5*particle[j].radius;
            if (vabs(distance) < r)
            {
                f = KSPRING_OBSTACLE*(r - distance);
                particle[j].fx += f*segment[i].nx;
                particle[j].fy += f*segment[i].ny;
                if (MOVE_BOUNDARY)
                {
                    segment[i].fx -= f*segment[i].nx;
                    segment[i].fy -= f*segment[i].ny;
//                 segment[i].fx += f*segment[i].nx;
//                 segment[i].fy += f*segment[i].ny;
                }
            }
        }
        
        /* compute force from concave corners */
        if (segment[i].concave)
        {
            distance = module2(x - segment[i].x1, y - segment[i].y1);
            angle = argument(x - segment[i].x1, y - segment[i].y1);
            if (angle < segment[i].angle1) angle += DPI;
            r = 1.5*particle[j].radius;
            
            if ((distance < r)&&(angle > segment[i].angle1)&&(angle < segment[i].angle2))
            {
                f = KSPRING_OBSTACLE*(r - distance);
                particle[j].fx += f*cos(angle);
                particle[j].fy += f*sin(angle);
                if (MOVE_BOUNDARY)
                {
                    segment[i].fx -= f*cos(angle);
                    segment[i].fy -= f*sin(angle);
//                     segment[i].fx += f*cos(angle);
//                     segment[i].fy += f*sin(angle);
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
    }
}


void compute_particle_force(int j, double krepel, t_particle particle[NMAXCIRCLES], t_hashgrid hashgrid[HASHX*HASHY])
/* compute force from other particles on particle j */
{
    int i0, j0, m0, k, m, q, close;
    double fx = 0.0, fy = 0.0, force[2], torque = 0.0, torque_ij, x, y;
    
    particle[j].neighb = 0;
    if (REACTION_DIFFUSION) particle[j].diff_neighb = 0;

    for (k=0; k<particle[j].hash_nneighb; k++)
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


int initialize_configuration(t_particle particle[NMAXCIRCLES], t_hashgrid hashgrid[HASHX*HASHY], 
                             t_obstacle obstacle[NMAXOBSTACLES], double px[NMAXCIRCLES], double py[NMAXCIRCLES], double pangle[NMAXCIRCLES], int tracer_n[N_TRACER_PARTICLES])
/* initialize all particles, obstacles, and the hashgrid */
{
    int i, j, k, n, nactive = 0;
    double x, y, h, xx, yy;
    
    for (i=0; i < ncircles; i++) 
    {
        /* set particle type */
        particle[i].type = 0;
        if ((TWO_TYPES)&&((double)rand()/RAND_MAX > TPYE_PROPORTION)) 
        {
            particle[i].type = 2;
            particle[i].radius = MU_B;
        }
        
        particle[i].neighb = 0;
        particle[i].diff_neighb = 0;
        particle[i].thermostat = 1;

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
        particle[i].energy = (particle[i].vx*particle[i].vx + particle[i].vy*particle[i].vy)*particle[i].mass_inv;
        
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
    }
    /* initialize dummy values in case particles are added */
    for (i=ncircles; i < NMAXCIRCLES; i++) 
    {
        particle[i].type = 0;
        particle[i].active = 0;
        particle[i].neighb = 0;
        particle[i].thermostat = 0;
        particle[i].energy = 0.0;
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
        if (((!POSITION_Y_DEPENDENCE)&&(particle[i].xc < 0))||((POSITION_Y_DEPENDENCE)&&(particle[i].yc < 0))) 
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
            if (!in_segment_region(particle[i].xc, particle[i].yc)) 
                particle[i].active = 0;
            
    /* case of reaction-diffusion equation */
    if (REACTION_DIFFUSION) for (i=0; i< ncircles; i++)
    {
        particle[i].type = 1 + (int)(RD_TYPES*(double)rand()/(double)RAND_MAX);
    }   
            
    /* count number of active particles */
    for (i=0; i< ncircles; i++) nactive += particle[i].active;
    printf("%i active particles\n", nactive);
    
    for (i=0; i<ncircles; i++) printf("particle %i of type %i\n", i, particle[i].type);
    
    return(nactive);
}


int add_particles(t_particle particle[NMAXCIRCLES], double px[NMAXCIRCLES], double py[NMAXCIRCLES], int nadd_particle)
/* add several particles to the system */
{
    static int i = 0;
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
            
    printf("Adding a particle\n\n");
    
    add_particle(XMIN - 0.5*MU, 0.0, 50.0 + 5.0*(double)i, 0.0, 2.0*PARTICLE_MASS, 0, particle);
    i++;
    
    particle[ncircles - 1].radius = 0.5*MU;
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

int partial_thermostat_coupling(t_particle particle[NMAXCIRCLES], double xmin)
/* only couple particles with x > xmin to thermostat */
{
    int condition, i, nthermo = 0;
    
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
                condition = (in_segment_region(particle[i].xc - xsegments, particle[i].yc - ysegments));
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

void update_types(t_particle particle[NMAXCIRCLES])
/* update the types in case of reaction-diffusion equation */
{
    int i, j, k;
    double distance;
    
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
}
    