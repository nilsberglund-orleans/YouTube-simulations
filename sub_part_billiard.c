long int global_time = 0;    /* counter to keep track of global time of simulation */
int nparticles=NPART; 

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
// 	if (alph < 0.0) alph += DPI;
	return(alph);
 }
 
 int polynome(a, b, c, r)
 double a, b, c, r[2];
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

    glOrtho(XMIN, XMAX, YMIN, YMAX , -1.0, 1.0);
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

void rgb_color_scheme(i, rgb) /* color scheme */
int i;
double rgb[3];
{
    double hue, y, r;
  
    hue = (double)(COLORSHIFT + i*360/NCOLORS);
    r = 0.9;
  
    while (hue < 0.0) hue += 360.0;
    while (hue >= 360.0) hue -= 360.0;
  
    hsl_to_rgb(hue, r, 0.5, rgb);
    /* saturation = r, luminosity = 0.5 */ 
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


void paint_billiard_interior()      /* points billiard interior, for use before draw_conf */
{
    double x0, x, y, phi, r = 0.01, alpha, dphi, omega;
    int i, k, c;
    
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
        default: 
        {
            
        }
    }
}

void draw_billiard()      /* draws the billiard boundary */
{
    double x0, x, y, phi, r = 0.01, alpha, dphi, omega, x1, y1, x2, beta2, angle, s;
    int i, j, k, c;
    
    if (PAINT_INT) glColor3f(0.5, 0.5, 0.5);
    else
    {
        if (BLACK) glColor3f(1.0, 1.0, 1.0);
        else glColor3f(0.0, 0.0, 0.0);
    }
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
                glEnable(GL_LINE_SMOOTH);
                glBegin(GL_LINE_LOOP);
                for (i=0; i<=NSEG; i++)
                {
                    phi = (double)i*DPI/(double)NSEG;
                    x = x0 + r*cos(phi);
                    y = r*sin(phi);
                    glVertex2d(x, y);
                }
                glEnd();

                glBegin(GL_LINE_LOOP);
                for (i=0; i<=NSEG; i++)
                {
                    phi = (double)i*DPI/(double)NSEG;
                    x = -x0 + r*cos(phi);
                    y = r*sin(phi);
                    glVertex2d(x, y);
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
            break; 
        }
        case D_SINAI:
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
            
            glBegin(GL_LINE_LOOP);
            for (i=0; i<=NSEG; i++)
            {
                phi = (double)i*DPI/(double)NSEG;
                x = LAMBDA*cos(phi) + MU;
                y = LAMBDA*sin(phi);
                glVertex2d(x, y);
            }
            glEnd ();
            glBegin(GL_LINE_LOOP);
            for (i=0; i<=NSEG; i++)
            {
                phi = (double)i*DPI/(double)NSEG;
                x = cos(phi);
                y = sin(phi);
                glVertex2d(x, y);
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
                glVertex2d(x, y);
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
                    glVertex2d(x, y);
                }
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


void print_config(conf)  /* for debugging purposes */
double conf[8];
{
    printf("s = %.3lg, u = %.3lg, t = %.3lg, L = %.3lg, x0 = %.3lg, y0 = %.3lg, x1 = %.3lg, y1 = %.3lg\n", conf[0], conf[1]/PI, conf[2], conf[3], conf[4], conf[5], conf[6], conf[7]);
}

void print_config_23(conf)  /* for debugging purposes */
double conf[8];
{
    printf("t = %.8f, L = %.8f\n", conf[2], conf[3]);
}

void print_colors(color)  /* for debugging purposes */
int color[NPARTMAX];
{
    int i;
    
    for (i=0; i<NPART; i++) printf("%i ", color[i]);
    printf("\n");
}


/****************************************************************************************/
/* rectangle billiard */
/****************************************************************************************/
 
 
 int pos_rectangle(conf, pos, alpha)
 /* determine position on boundary of rectangle */
 /* the corners of the rectangle are (-LAMBDA,-1), ..., (LAMBDA,1) */
 /* returns the number of the side hit, or 0 if hitting a corner */
 double conf[2], pos[2], *alpha;

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
 
 
 int vrectangle_xy(config, alpha, pos)
 /* determine initial configuration for start at point pos = (x,y) */
 double config[8], alpha, pos[2];
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
 
 int vrectangle(config)
 double config[8];

 {
	double pos[2], alpha;
	int c;

	/* position et vitesse de depart */

	c = pos_rectangle(config, pos, &alpha);
        
        vrectangle_xy(config, alpha, pos);

        return(c);
 }
 

/* elliptic billiard */

 int pos_ellipse(conf, pos, alpha)
 /* determine position on boundary of ellipse */
 double conf[2], pos[2], *alpha;

 {
	double theta;

        pos[0] = LAMBDA*cos(conf[0]);
        pos[1] = sin(conf[0]);
        
        theta = argument(-LAMBDA*pos[1],pos[0]/LAMBDA);
        *alpha = theta + conf[1]; 
        
        return(1);
 }
 
 
 int vellipse_xy(config, alpha, pos)
 /* determine initial configuration for start at point pos = (x,y) */
 double config[8], alpha, pos[2];

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

 int vellipse(config)
 /* determine initial configuration when starting from boundary */
 double config[8];

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

 int pos_stade(conf, pos, alpha)
 /* determine position on boundary of stadium */
 double conf[2], pos[2], *alpha;

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
 
 
 int vstade_xy(config, alpha, pos)
 /* determine initial configuration for start at point pos = (x,y) */
 double config[8], alpha, pos[2];

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
 
 int vstade(config)
 double config[8];
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
 
 int pos_sinai(conf, pos, alpha)
 /* determine position on boundary of Sinai billiard */
 /* s in [0,2 Pi) is on the circle, other s are on boundary of window */
 double conf[2], pos[2], *alpha;

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
 
 
 int vsinai_xy(config, alpha, pos)
 /* determine initial configuration for start at point pos = (x,y) */
 double config[8], alpha, pos[2];

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
 
 int vsinai(config)
 double config[8];

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
 
 
 int pos_triangle(conf, pos, alpha)
 /* determine position on boundary of triangle */
 /* the corners of the triangle are (-LAMBDA,-1), (LAMBDA,-1), (-LAMBDA,1) */
 /* we use arclength for horizontal and vertical side, x for diagonal */
  double conf[2], pos[2], *alpha;

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
 
 
 int vtriangle_xy(config, alpha, pos)
 /* determine initial configuration for start at point pos = (x,y) */
 /* Warning: reflection in the corners is not yet implemented correctly */
 double config[8], alpha, pos[2];

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
 
 int vtriangle(config)
 double config[8];

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

 int pos_annulus(conf, pos, alpha)
 /* determine position on boundary of annulus */
 double conf[2], pos[2], *alpha;

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
 
 int vannulus_xy(config, alpha, pos)
 /* determine initial configuration for start at point pos = (x,y) */
 double config[8], alpha, pos[2];

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

 int vannulus(config)
 /* determine initial configuration when starting from boundary */
 double config[8];

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

 int pos_polygon(conf, pos, alpha)
 /* determine position on boundary of polygon */
 /* conf[0] is arclength on boundary */
 double conf[2], pos[2], *alpha;

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
 
 int vpolygon_xy(config, alpha, pos)
 /* determine initial configuration for start at point pos = (x,y) */
 double config[8], alpha, pos[2];

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

 int vpolygon(config)
 /* determine initial configuration when starting from boundary */
 double config[8];

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

 int pos_reuleaux(conf, pos, alpha)
 /* determine position on boundary of polygon */
 /* conf[0] is arclength on boundary */
 double conf[2], pos[2], *alpha;

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
 
 int vreuleaux_xy(config, alpha, pos)
 /* determine initial configuration for start at point pos = (x,y) */
 double config[8], alpha, pos[2];

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
            
//             if ((intb))     /* check if condition is ok */
//             if ((cos(theta) > 0.0)&&(intb))     /* check if condition is ok */
            {
                ca = cos(rangle);
                sa = sin(rangle);
                
//                 printf("theta = %.5lg\n", theta);
//                 printf("rangle = %.5lg x0 = %.5lg y0 = %.5lg \n", rangle, pos[0], pos[1]);
                
                x = pos[0]*ca + pos[1]*sa;
                y = -pos[0]*sa + pos[1]*ca;
                
//                 printf("x = %.5lg\t y = %.5lg\n", x, y);
            
                a = (x-x2)*cos(theta) + y*sin(theta);
                b = (x-x2)*(x-x2) + y*y - LAMBDA*LAMBDA;
                
//                 printf("a = %.5lg\t b = %.5lg\n", a, b);
                
                if (a*a - b > margin)
                {
//                     t = vabs(a) - sqrt(a*a - b);
                    if (LAMBDA > 0.0) t = -a - sqrt(a*a - b);
                    else t = -a + sqrt(a*a - b);
                    
                    xi = x + t*cos(theta);
                    yi = y + t*sin(theta);
                    
//                     printf("t = %.5lg\t xi = %.5lg\t yi = %.5lg\n", t, xi, yi);
                
                    if ((t > margin)&&(vabs(yi) <= sin(omega2))) 
                    {
                        cval[nt] = k;
                        tval[nt] = t;
                        
                        /* rotate back */
                        x1[nt] = xi*ca - yi*sa;
                        y1[nt] = xi*sa + yi*ca;
                    
//                         intb = 0;
//                         c = k;
                        tempconf[nt][0] = ((double)k + 0.5)*beta + asin(yi/LAMBDA);
                        tempconf[nt][1] = PID - asin(yi/LAMBDA) - theta;      
//                         tempconf[nt][0] = ((double)k + 0.5)*beta + asin(yi/vabs(LAMBDA));
//                         tempconf[nt][1] = PID - asin(yi/vabs(LAMBDA)) - theta;      
                        nt++;
                    }
                }
            }
        }
//         printf("nt = %i\n", nt);
        
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

 int vreuleaux(config)
 /* determine initial configuration when starting from boundary */
 double config[8];

 {
	double pos[2], alpha;
	int c;

        c = pos_reuleaux(config, pos, &alpha);
        
        vreuleaux_xy(config, alpha, pos);
	
	return(c);
 }
 
/****************************************************************************************/
/* general billiard */
/****************************************************************************************/
 
 int pos_billiard(conf, pos, alpha)
 /* determine initial configuration for start at point pos = (x,y) */
 double conf[8], pos[2], *alpha;
 {
    switch (B_DOMAIN) {
        case (D_RECTANGLE):
        {
            return(pos_rectangle(conf, pos, &alpha));
            break;
        }
        case (D_ELLIPSE):
        {
            return(pos_ellipse(conf, pos, &alpha));
            break;
        }
        case (D_STADIUM):
        {
            return(pos_stade(conf, pos, &alpha));
            break;
        }
        case (D_SINAI):
        {
            return(pos_sinai(conf, pos, &alpha));
            break;
        }
        case (D_TRIANGLE):
        {
            return(pos_triangle(conf, pos, &alpha));
            break;
        }
        case (D_ANNULUS):
        {
            return(pos_annulus(conf, pos, &alpha));
            break;
        }
        case (D_POLYGON):
        {
            return(pos_polygon(conf, pos, &alpha));
            break;
        }
        case (D_REULEAUX):
        {
            return(pos_reuleaux(conf, pos, &alpha));
            break;
        }
        default: 
        {
            printf("Function pos_billiard not defined for this billiard \n");
        }
    }
 }


 
 int vbilliard_xy(config, alpha, pos)
 /* determine initial configuration for start at point pos = (x,y) */
 double config[8], alpha, pos[2];
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
        default: 
        {
            printf("Function vbilliard_xy not defined for this billiard \n");
        }
    }
 }

 /* TO DO: fix returned value */
 
 int vbilliard(config)
 /* determine initial configuration when starting from boundary */
 double config[8];
 {
    double pos[2], theta, alpha;
    int c;

    switch (B_DOMAIN) {
        case (D_RECTANGLE):
        {
            c = pos_rectangle(config, pos, &alpha);
        
            return(vrectangle(config, alpha, pos));
            break;
        }
        case (D_ELLIPSE):
        {
            c = pos_ellipse(config, pos, &alpha);
         
            return(vellipse(config, alpha, pos));
            break;
        }
        case (D_STADIUM):
        {
            c = pos_stade(config, pos, &alpha);
        
            return(vstade(config, alpha, pos));
            break;
        }
        case (D_SINAI):
        {
            c = pos_sinai(config, pos, &alpha);
        
            return(vsinai(config, alpha, pos));
            break;
        }
        case (D_TRIANGLE):
        {
            c = pos_triangle(config, pos, &alpha);
        
            return(vtriangle(config, alpha, pos));
            break;
        }
        case (D_ANNULUS):
        {
            c = pos_annulus(config, pos, &alpha);
        
            return(vannulus(config, alpha, pos));
            break;
        }
        case (D_POLYGON):
        {
            c = pos_polygon(config, pos, &alpha);
        
            return(vpolygon(config, alpha, pos));
            break;
        }
        case (D_REULEAUX):
        {
            c = pos_reuleaux(config, pos, &alpha);
        
            return(vreuleaux(config, alpha, pos));
            break;
        }
        default: 
        {
            printf("Function vbilliard not defined for this billiard \n");
        }
    }
 }
 
 int xy_in_billiard(x, y)
 /* returns 1 if (x,y) represents a point in the billiard */
 double x, y;
 {
    double l2, r2, omega, omega2, c, angle, x1, y1, x2, co, so;
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
            r2 = x*x + y*y;
            if ((r2 > l2)&&(r2 < 1.0)) return(1);
            else return(0);
            break;
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
            return(condition);
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
                
//                 if (!condition)
//                 {
//                     printf("x = %.5lg \t y = %.5lg \t x1 = %.5lg \t y1 = %.5lg \t angle = %.5lg \n", x, y, x1, y1, angle);
//                     printf("k = %i \t condition = %i\n", k, condition);
//                     sleep(1);
//                 }
            }
            return(condition);
            break;            
        }
        /* D_REULEAUX : distance to all centers of arcs should be larger than LAMBDA */
        default: 
        {
            printf("Function ij_in_billiard not defined for this billiard \n");
            return(0);
        }
    }
 }
 
