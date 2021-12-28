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


