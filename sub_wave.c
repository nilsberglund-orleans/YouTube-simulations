/*********************/
/* Graphics routines */
/*********************/

#include "colors_waves.c"


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


void init_circle_config()
/* initialise the arrays circlex, circley, circlerad and circleactive */
/* for billiard shape D_CIRCLES */
{
    int i, j, k, n, ncirc0, n_p_active, ncandidates=5000, naccepted; 
    double dx, dy, p, phi, r, r0, ra[5], sa[5], height, x, y = 0.0, gamma, dpoisson = 3.25*MU, xx[4], yy[4];
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
                    circlex[n] = ((double)(i-NGRIDX/2) + 0.5)*dy;
                    circley[n] = YMIN + ((double)j + 0.5)*dy;
                    circlerad[n] = MU;
                    circleactive[n] = 1;
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
                    circlex[n] = ((double)(i-NGRIDX/2) + 0.5)*dy;   /* is +0.5 needed? */
                    circley[n] = YMIN + ((double)j - 0.5)*dy;
                    if ((i+NGRIDX)%2 == 1) circley[n] += 0.5*dy;
                    circlerad[n] = MU;
                    /* activate only circles that intersect the domain */
                    if ((circley[n] < YMAX + MU)&&(circley[n] > YMIN - MU)) circleactive[n] = 1;
                    else circleactive[n] = 0;
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
                    circlex[n] = ((double)(i-NGRIDX/2) + 0.5*((double)rand()/RAND_MAX - 0.5))*dy;
                    circley[n] = YMIN + ((double)j + 0.5 + 0.5*((double)rand()/RAND_MAX - 0.5))*dy;
                    circlerad[n] = MU*sqrt(1.0 + 0.8*((double)rand()/RAND_MAX - 0.5));
                    circleactive[n] = 1;
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
                    circlex[n] = ((double)(i-NGRIDX/2) + 0.5)*dy;
                    circley[n] = YMIN + ((double)j + 0.5)*dy;
                    circlerad[n] = MU;
                    p = (double)rand()/RAND_MAX;
                    if (p < P_PERCOL) circleactive[n] = 1;
                    else circleactive[n] = 0;
                }
            break;
        }
        case (C_RAND_POISSON):
        {
            ncircles = NPOISSON;
            for (n = 0; n < NPOISSON; n++)
            {
                circlex[n] = LAMBDA*(2.0*(double)rand()/RAND_MAX - 1.0);
                circley[n] = (YMAX - YMIN)*(double)rand()/RAND_MAX + YMIN;
                circlerad[n] = MU;
                circleactive[n] = 1;
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
                    circlex[n] = r*cos(phi);
                    circley[n] = r*sin(phi);
                    circlerad[n] = MU;
                    circleactive[n] = 1;
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
                    circlex[n] = r*cos(phi);
                    circley[n] = r*sin(phi);
                    circlerad[n] = LAMBDA*ra[j];
                    circleactive[n] = 1;
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
                    circlex[4*i + j] = xx[i];
                    circley[4*i + j] = yy[j];
                    
                }
                
            circlex[ncircles - 1] = X_TARGET;
            circley[ncircles - 1] = Y_TARGET;
            
            for (i=0; i<ncircles - 1; i++)
            {
                circlerad[i] = MU;
                circleactive[i] = 1;
            }
            
            circlerad[ncircles - 1] = 0.5*MU;
            circleactive[ncircles - 1] = 2;
            
            break;
        }        
        case (C_POISSON_DISC):
        {
            printf("Generating Poisson disc sample\n");
            /* generate first circle */
            circlex[0] = LAMBDA*(2.0*(double)rand()/RAND_MAX - 1.0);
            circley[0] = (YMAX - YMIN)*(double)rand()/RAND_MAX + YMIN;
            active_poisson[0] = 1;
//             circleactive[0] = 1;
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
                    x = circlex[i] + r*cos(phi);
                    y = circley[i] + r*sin(phi);
//                        printf("Testing new circle at (%.3f,%.3f)\t", x, y);
                    far = 1;
                    for (k=0; k<ncircles; k++) if ((k!=i))
                    {
                        /* new circle is far away from circle k */
                        far = far*((x - circlex[k])*(x - circlex[k]) + (y - circley[k])*(y - circley[k]) >=     dpoisson*dpoisson);
                        /* new circle is in domain */
                        far = far*(vabs(x) < LAMBDA)*(y < YMAX)*(y > YMIN);
                    }
                    if (far)    /* accept new circle */
                    {
                        printf("New circle at (%.3f,%.3f) accepted\n", x, y);
                        circlex[ncircles] = x;
                        circley[ncircles] = y;
                        circlerad[ncircles] = MU;
                        circleactive[ncircles] = 1;
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
                circlex[n] = -LAMBDA + n*dx;
                circley[n] = y;
                y += height*gamma; 
                if (y > YMAX) y -= height;
                circlerad[n] = MU;
                circleactive[n] = 1;
            }
            
            /* test for circles that overlap top or bottom boundary */
            ncirc0 = ncircles;
            for (n=0; n < ncirc0; n++)
            {
                if (circley[n] + circlerad[n] > YMAX)
                {
                    circlex[ncircles] = circlex[n];
                    circley[ncircles] = circley[n] - height;
                    circlerad[ncircles] = MU;
                    circleactive[ncircles] = 1;
                    ncircles ++;
                }
                else if (circley[n] - circlerad[n] < YMIN)
                {
                    circlex[ncircles] = circlex[n];
                    circley[ncircles] = circley[n] + height;
                    circlerad[ncircles] = MU;
                    circleactive[ncircles] = 1;
                    ncircles ++;
                }
            }
            break;
        }
        case (C_GOLDEN_SPIRAL):
        {
            ncircles = 1;
            circlex[0] = 0.0;
            circley[0] = 0.0;
            
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
                    circlex[ncircles] = x;
                    circley[ncircles] = y;
                    ncircles++;
                }
            }
            
            for (i=0; i<ncircles; i++)
            {
                circlerad[i] = MU;
                /* inactivate circles outside the domain */
                if ((circley[i] < YMAX + MU)&&(circley[i] > YMIN - MU)) circleactive[i] = 1;
//                 printf("i = %i, circlex = %.3lg, circley = %.3lg\n", i, circlex[i], circley[i]);
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
                    circlex[n] = ((double)(i-NGRIDX/2) + 0.5)*dy;   /* is +0.5 needed? */
                    circley[n] = YMIN + ((double)j - 0.5)*dy;
                    if (((i+NGRIDX)%4 == 2)||((i+NGRIDX)%4 == 3)) circley[n] += 0.5*dy;
                    circlerad[n] = MU;
                    /* activate only circles that intersect the domain */
                    if ((circley[n] < YMAX + MU)&&(circley[n] > YMIN - MU)) circleactive[n] = 1;
                    else circleactive[n] = 0;
                }
            break;
        }
        case (C_ONE):
        {
            circlex[ncircles] = 0.0;
            circley[ncircles] = 0.0;
            circlerad[ncircles] = MU;
            circleactive[ncircles] = 1;
            ncircles += 1;
            break;
        }
        case (C_TWO):   /* used for comparison with cloak */
        {
            circlex[ncircles] = 0.0;
            circley[ncircles] = 0.0;
            circlerad[ncircles] = MU;
            circleactive[ncircles] = 2;
            ncircles += 1;

            circlex[ncircles] = 0.0;
            circley[ncircles] = 0.0;
            circlerad[ncircles] = 2.0*MU;
            circleactive[ncircles] = 1;
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


void compute_isospectral_coordinates(int type, double xshift, double yshift, double scaling, double vertex[9][2])
/* compute positions of vertices of isospectral billiards */
/* central triangle has coordinates (0,0), (1,0) and (LAMBDA,MU) fed into affine transformation */
/* defined by (xshift - 0.5), (yshift - 0.25) and scaling*/
{
    vertex[0][0] = (xshift - 0.5)*scaling;
    vertex[0][1] = (yshift - 0.25)*scaling;
    
    vertex[1][0] = (0.5 + xshift)*scaling;
    vertex[1][1] = (yshift - 0.25)*scaling;
    
    vertex[2][0] = (LAMBDA + xshift - 0.5)*scaling;
    vertex[2][1] = (MU + yshift - 0.25)*scaling; 
    
    axial_symmetry(vertex[0], vertex[2], vertex[1], vertex[3]);
    axial_symmetry(vertex[0], vertex[1], vertex[2], vertex[4]);
    axial_symmetry(vertex[1], vertex[2], vertex[0], vertex[5]);

    if (type == 0)
    {
        axial_symmetry(vertex[0], vertex[3], vertex[2], vertex[6]);
        axial_symmetry(vertex[1], vertex[4], vertex[0], vertex[7]);
        axial_symmetry(vertex[2], vertex[5], vertex[1], vertex[8]);
    }
    else
    {
        axial_symmetry(vertex[2], vertex[3], vertex[0], vertex[6]);
        axial_symmetry(vertex[0], vertex[4], vertex[1], vertex[7]);
        axial_symmetry(vertex[1], vertex[5], vertex[2], vertex[8]);
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


int xy_in_billiard(double x, double y)
/* returns 1 if (x,y) represents a point in the billiard */
// double x, y;
{
    double l2, r2, r2mu, omega, b, c, angle, z, x1, y1, x2, y2, u, v, u1, v1, dx, dy, width;
    int i, j, k, k1, k2, condition;
    static double vertex[9][2], wertex[9][2];
    static int first = 1;

    switch (B_DOMAIN) {
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
        case (D_ISOSPECTRAL):
        {
            if (first)
            {
                compute_isospectral_coordinates(0, ISO_XSHIFT_LEFT, ISO_YSHIFT_LEFT, ISO_SCALE, vertex);
                compute_isospectral_coordinates(1, ISO_XSHIFT_RIGHT, ISO_YSHIFT_RIGHT, ISO_SCALE, wertex);
                first = 0;
            }
            /* 1st triangle */
            condition = xy_in_triangle(x, y, vertex[0], vertex[1], vertex[2]);
            condition += xy_in_triangle(x, y, vertex[0], vertex[4], vertex[1]);
            condition += xy_in_triangle(x, y, vertex[1], vertex[5], vertex[2]);
            condition += xy_in_triangle(x, y, vertex[0], vertex[2], vertex[3]);
            condition += xy_in_triangle(x, y, vertex[1], vertex[4], vertex[7]);
            condition += xy_in_triangle(x, y, vertex[2], vertex[5], vertex[8]);
            condition += xy_in_triangle(x, y, vertex[0], vertex[3], vertex[6]);

            /* 2nd triangle */
            condition += xy_in_triangle(x, y, wertex[0], wertex[1], wertex[2]);
            condition += xy_in_triangle(x, y, wertex[0], wertex[4], wertex[1]);
            condition += xy_in_triangle(x, y, wertex[1], wertex[5], wertex[2]);
            condition += xy_in_triangle(x, y, wertex[0], wertex[2], wertex[3]);
            condition += xy_in_triangle(x, y, wertex[0], wertex[7], wertex[4]);
            condition += xy_in_triangle(x, y, wertex[1], wertex[8], wertex[5]);
            condition += xy_in_triangle(x, y, wertex[2], wertex[6], wertex[3]);
            return(condition >= 1);
        }
        case (D_CIRCLES):
        {
            for (i = 0; i < ncircles; i++)
                if (circleactive[i]) 
                {
                    x1 = circlex[i];
                    y1 = circley[i];
                    r2 = circlerad[i]*circlerad[i];
                    if ((x-x1)*(x-x1) + (y-y1)*(y-y1) < r2) return(0); 
                }
            return(1);
        }
        case (D_CIRCLES_IN_RECT):   /* returns 2 inside circles, 0 outside rectangle */
        {
            for (i = 0; i < ncircles; i++)
                if (circleactive[i]) 
                {
                    x1 = circlex[i];
                    y1 = circley[i];
                    r2 = circlerad[i]*circlerad[i];
                    if ((x-x1)*(x-x1) + (y-y1)*(y-y1) < r2) return(2); 
                }
            if ((vabs(x) >= LAMBDA)||(vabs(y) >= 1.0)) return(0);
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
        default:
        {
            printf("Function ij_in_billiard not defined for this billiard \n");
            return(0);
        }
    }
}

int ij_in_billiard(int i, int j)
/* returns 1 if (i,j) represents a point in the billiard */
{
    double xy[2];

    ij_to_xy(i, j, xy);

    return(xy_in_billiard(xy[0], xy[1]));
}

void toka_lineto(double x1, double y1)
/* draws boundary segments of Tokarsky billiard */
{
    double ratio, x, y, pos[2];
    
    ratio = (XMAX - XMIN)/8.4;
    x = ratio*(x1 - 4.0);
    y = ratio*(y1 - 2.0);
    xy_to_pos(x, y, pos);
    glVertex2d(pos[0], pos[1]);
}

void iso_lineto(double z[2])
/* draws boundary segments of isospectral billiard */
{
    double pos[2];
    
    xy_to_pos(z[0], z[1], pos);
    glVertex2d(pos[0], pos[1]);
}


void draw_billiard()      /* draws the billiard boundary */
{
    double x0, x, y, x1, y1, dx, dy, phi, r = 0.01, pos[2], pos1[2], alpha, dphi, omega, z, l, width, a, b, c, ymax;
    static double vertex[9][2], wertex[9][2];
    int i, j, k, k1, k2, mr2;
    static int first = 1;

    if (BLACK) glColor3f(1.0, 1.0, 1.0);
    else glColor3f(0.0, 0.0, 0.0);
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
                glColor3f(0.3, 0.3, 0.3);
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
            toka_lineto(0.0, 2.0);
            toka_lineto(1.0, 3.0);
            toka_lineto(1.0, 4.0);
            toka_lineto(2.0, 4.0);
            toka_lineto(2.0, 3.0);
            toka_lineto(3.0, 3.0);
            toka_lineto(3.0, 2.0);
            toka_lineto(5.0, 2.0);
            toka_lineto(5.0, 3.0);
            toka_lineto(6.0, 3.0);
            toka_lineto(6.0, 4.0);
            toka_lineto(7.0, 3.0);
            toka_lineto(8.0, 3.0);
            toka_lineto(8.0, 2.0);
            toka_lineto(7.0, 2.0);
            toka_lineto(7.0, 1.0);
            toka_lineto(6.0, 0.0);
            toka_lineto(6.0, 1.0);
            toka_lineto(5.0, 1.0);
            toka_lineto(4.0, 0.0);
            toka_lineto(4.0, 1.0);
            toka_lineto(3.0, 1.0);
            toka_lineto(2.0, 0.0);
            toka_lineto(2.0, 1.0);
            toka_lineto(1.0, 1.0);
            toka_lineto(1.0, 2.0);
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
        case (D_ISOSPECTRAL):
        {
            if (first) 
            {
                compute_isospectral_coordinates(0, ISO_XSHIFT_LEFT, ISO_YSHIFT_LEFT, ISO_SCALE, vertex);
                compute_isospectral_coordinates(1, ISO_XSHIFT_RIGHT, ISO_YSHIFT_RIGHT, ISO_SCALE, wertex);
//                 compute_isospectral_coordinates(0, -1.0, 0.05, 0.9, vertex);
//                 compute_isospectral_coordinates(1, 0.9, 0.05, 0.9, wertex);
//             for (i=0; i<9; i++) printf("(x%i, y%i) = (%.2f, %.2f)\n", i, i, vertex[i][0], vertex[i][1]);
                first = 0;
            }
            /* 1st triangle */
            glBegin(GL_LINE_LOOP);
            iso_lineto(vertex[0]);
            iso_lineto(vertex[4]);
            iso_lineto(vertex[7]);
            iso_lineto(vertex[1]);
            iso_lineto(vertex[5]);
            iso_lineto(vertex[8]);
            iso_lineto(vertex[2]);
            iso_lineto(vertex[3]);
            iso_lineto(vertex[6]);
            glEnd();
            
            /* inner lines */
            glBegin(GL_LINE_LOOP);
            iso_lineto(vertex[0]);
            iso_lineto(vertex[1]);
            iso_lineto(vertex[2]);
            iso_lineto(vertex[0]);
            iso_lineto(vertex[3]);
            iso_lineto(vertex[2]);
            iso_lineto(vertex[5]);
            iso_lineto(vertex[1]);
            iso_lineto(vertex[4]);
            glEnd();
            
            /* 2nd triangle */
            glBegin(GL_LINE_LOOP);
            iso_lineto(wertex[0]);
            iso_lineto(wertex[7]);
            iso_lineto(wertex[4]);
            iso_lineto(wertex[1]);
            iso_lineto(wertex[8]);
            iso_lineto(wertex[5]);
            iso_lineto(wertex[2]);
            iso_lineto(wertex[6]);
            iso_lineto(wertex[3]);
            glEnd();
            
            /* inner lines */
            glBegin(GL_LINE_LOOP);
            iso_lineto(wertex[0]);
            iso_lineto(wertex[1]);
            iso_lineto(wertex[2]);
            iso_lineto(wertex[0]);
            iso_lineto(wertex[4]);
            iso_lineto(wertex[1]);
            iso_lineto(wertex[5]);
            iso_lineto(wertex[2]);
            iso_lineto(wertex[3]);
            glEnd();
            break;
        }
        case (D_CIRCLES):
        {
            glLineWidth(BOUNDARY_WIDTH);
            for (i = 0; i < ncircles; i++) 
                if (circleactive[i]) draw_circle(circlex[i], circley[i], circlerad[i], NSEG);
            break;
        }
        case (D_CIRCLES_IN_RECT):
        {
            glLineWidth(BOUNDARY_WIDTH);
            for (i = 0; i < ncircles; i++) 
                if (circleactive[i]) draw_circle(circlex[i], circley[i], circlerad[i], NSEG);
            draw_rectangle(-LAMBDA, -1.0, LAMBDA, 1.0);
            if ((FOCI)&&(CIRCLE_PATTERN == C_LASER))
            {
                glColor3f(0.3, 0.3, 0.3);
                draw_circle(X_SHOOTER, Y_SHOOTER, r, NSEG);
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
       default:
        {
            printf("Function draw_billiard not defined for this billiard \n");
        }
    }
}
