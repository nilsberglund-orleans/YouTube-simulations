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
    double l2, r2, omega, c, angle, z;
    int i, k, k1, k2, condition;

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
        case D_YOUNG:
        {
            if ((x < -MU)||(x > MU)) return(1);
            else if ((vabs(y-LAMBDA) < MU)||(vabs(y+LAMBDA) < MU)) return (1);
            else return(0);
        }
        case D_GRATING:
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
        case D_EHRENFEST:
        {
            return(((x-1.0)*(x-1.0) + y*y < LAMBDA*LAMBDA)||((x+1.0)*(x+1.0) + y*y < LAMBDA*LAMBDA)||((vabs(x) < 1.0)&&(vabs(y) < MU)));
        }
        default:
        {
            printf("Function ij_in_billiard not defined for this billiard \n");
            return(0);
        }
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
    double x0, x, y, phi, r = 0.01, pos[2], pos1[2], alpha, dphi, omega, z;
    int i, j, k1, k2;

    if (BLACK) glColor3f(1.0, 1.0, 1.0);
    else glColor3f(0.0, 0.0, 0.0);
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
        case D_GRATING:
        {
            k1 = -(int)(-YMIN/LAMBDA);
            k2 = (int)(YMAX/LAMBDA);
            for (i=k1; i<= k2; i++)
            {
                z = (double)i*LAMBDA;
                glBegin(GL_LINE_LOOP);
                for (j=0; j<=NSEG; j++)
                {
                    phi = (double)j*DPI/(double)NSEG;
                    x = MU*cos(phi);
                    y = z + MU*sin(phi);
                    xy_to_pos(x, y, pos);
                    glVertex2d(pos[0], pos[1]);
                }
                glEnd ();
            }
            break;
        }
        case D_EHRENFEST:
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
        default:
        {
            printf("Function draw_billiard not defined for this billiard \n");
        }
    }
}
