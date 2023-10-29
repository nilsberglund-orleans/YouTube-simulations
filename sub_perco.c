/*********************/
/* Graphics routines */
/*********************/

#include "colors_waves.c"

#define CLUSTER_SHIFT 10    /* shift in numbering of open clusters */

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
  return 0;
}

void init()		/* initialisation of window */
{
    glLineWidth(3);

    glClearColor(0.0, 0.0, 0.0, 1.0);
    glClear(GL_COLOR_BUFFER_BIT);

    if (PLOT_3D) glOrtho(XMIN, XMAX, YMIN, YMAX , -1.0, 1.0);
    else glOrtho(0.0, NX, 0.0, NY, -1.0, 1.0);
}

void blank()
{
    if (BLACK) glClearColor(0.0, 0.0, 0.0, 1.0);
    else glClearColor(1.0, 1.0, 1.0, 1.0);
    glClear(GL_COLOR_BUFFER_BIT);
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


void save_frame_perc()
{
  static int counter = 0;
  char *name="perc.", n2[100];
  char format[6]=".%05i";

    counter++;
    strcpy(n2, name);
    sprintf(strstr(n2,"."), format, counter);
    strcat(n2, ".tif");
    printf(" saving frame %s \n",n2);
    writetiff(n2, "Wave equation in a planar domain", 0, 0,
         WINWIDTH, WINHEIGHT, COMPRESSION_LZW);

}

/*********************/
/* some basic math   */
/*********************/


double ipow(double x, int n)
 {
    double y;
    int i;
    
    y = x;
    for (i=1; i<n; i++) y *= x;
    
    return(y);
 }
 

int ipowi(int base, int n)
 {
    int p, i;
    
    if (n == 0) return(1);
    else
    {
        p = base;
        for (i=1; i<n; i++) p *= base;
    
        return(p);
    }
 }
 
 double module2(double x, double y)   /* Euclidean norm */
 {
	double m;

	m = sqrt(x*x + y*y);
	return(m);
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

    if (PLOT_3D)
    {
        pos[0] = x;
        pos[1] = y;
    }
    else
    {
        x1 = (x - XMIN)/(XMAX - XMIN);
        y1 = (y - YMIN)/(YMAX - YMIN);

        pos[0] = x1 * (double)NX;
        pos[1] = y1 * (double)NY;
    }
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

void draw_colored_rectangle(double x1, double y1, double x2, double y2, double hue)
{
    double pos[2], rgb[3];
    
    hsl_to_rgb_palette(hue, 0.9, 0.5, rgb, COLOR_PALETTE);
    glColor3f(rgb[0], rgb[1], rgb[2]);
    
    glBegin(GL_QUADS);
    xy_to_pos(x1, y1, pos);
    glVertex2d(pos[0], pos[1]);
    xy_to_pos(x2, y1, pos);
    glVertex2d(pos[0], pos[1]);
    xy_to_pos(x2, y2, pos);
    glVertex2d(pos[0], pos[1]);
    xy_to_pos(x1, y2, pos);
    glVertex2d(pos[0], pos[1]);
    glEnd();    

    glColor3f(0.0, 0.0, 0.0);
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


void draw_colored_circle(double x, double y, double r, int nseg)
{
    int i;
    double pos[2], alpha, dalpha;
    
    dalpha = DPI/(double)nseg;
    
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


int graphical_rep(int bcondition)
/* return type of drawing, depending on lattice/boundary conditions */
{
    switch (bcondition) {
        case (BC_SQUARE_DIRICHLET): return(PLOT_SQUARES);
        case (BC_SQUARE_PERIODIC): return(PLOT_SQUARES);
        case (BC_SQUARE_BOND_DIRICHLET): return(PLOT_SQUARE_BONDS);
        case (BC_HEX_SITE_DIRICHLET): return(PLOT_HEX); 
        case (BC_HEX_BOND_DIRICHLET): return(PLOT_HEX_BONDS); 
        case (BC_TRIANGLE_SITE_DIRICHLET): return(PLOT_TRIANGLE); 
        case (BC_POISSON_DISC): return(PLOT_POISSON_DISC);
//         case (BC_CUBIC_DIRICHLET): return(PLOT_SQUARES);
        case (BC_CUBIC_DIRICHLET): return(PLOT_CUBES);
        default: return(0);
    }
    
}

double pcritical(int lattice)
/* critical probability in terms of lattice and boundary condition */
{
    switch (lattice) {
        case (BC_SQUARE_DIRICHLET): return(0.59274);
        case (BC_SQUARE_PERIODIC): return(0.59274);
        case (BC_SQUARE_BOND_DIRICHLET): return(0.5);
        case (BC_HEX_BOND_DIRICHLET): return(1.0 - 2.0*sin(PI/18.0));
        case (BC_TRIANGLE_SITE_DIRICHLET): return(0.6970402);
        case (BC_CUBIC_DIRICHLET): return(0.311604);
        default: return(0.5);
        
    }
    
}

int cellnb(int i, int j, int group, int nx, int ny)
/* convert 2d coordinates to 1d */
{
    switch (LATTICE) {
        case (BC_SQUARE_DIRICHLET):
        {
            return(i*ny+j);
        }
        case (BC_SQUARE_PERIODIC):
        {
            return(i*ny+j);
        }
        case (BC_SQUARE_BOND_DIRICHLET):
        {
            if (group == 0) return(i+nx*j);
            else return (nx*(ny+1) + i*ny+j);
        }
        case (BC_HEX_SITE_DIRICHLET):
        {
            return(i+nx*j);
        }
        case (BC_HEX_BOND_DIRICHLET):
        {
            return (group*nx*(ny+1) + i+nx*j);
        }
        case (BC_TRIANGLE_SITE_DIRICHLET):
        {
            return(i+2*nx*j);
        }
        case (BC_POISSON_DISC):
        {
            return(i);
        }
    }
}

int cellnb_3d(int i, int j, int k, int group, int nx, int ny, int nz)
/* convert 3d coordinates to 1d */
{
    switch (LATTICE) {
        case (BC_CUBIC_DIRICHLET):
        {
            return(k*nx*ny + j*nx + i);
        }
    }
}

int cell_to_ij(int c, int *i, int *j, int nx, int ny)
/* convert 1d coordinates to 2d, returns group */
{
    int group;
    
    switch (LATTICE) {
        case (BC_SQUARE_DIRICHLET):
        {
            *i = c/ny;
            *j = c - *i*ny;
            return(0);
        }
        case (BC_SQUARE_PERIODIC):
        {
            *i = c/ny;
            *j = c - *i*ny;
            return(0);
        }
        case (BC_SQUARE_BOND_DIRICHLET):
        {
            if (c < nx*(ny+1)) 
            {
                *j = c/nx;
                *i = c - *j*nx;
                return(0);
            }
            else 
            {
                c -= nx*(ny+1);
                *i = c/ny;
                *j = c - *i*ny;
                return(1);
            }
        }
        case (BC_HEX_SITE_DIRICHLET):
        {
            *j = c/nx;
            *i = c - *j*nx;
            return(0);
        }
        case (BC_HEX_BOND_DIRICHLET):
        {
            group = c/(nx*(ny+1));
            c -= group*(nx*(ny+1));
            *j = c/nx;
            *i = c - *j*nx;
            return(group);
        }
        case (BC_TRIANGLE_SITE_DIRICHLET):
        {
            *j= c/(2*nx);
            *i = c - 2*(*j)*nx;
            return(0);
        }
    }
}

double p_schedule(int i)
/* percolation probability p as a function of time */
{
    double time, pstar;
    int factor;
    
    pstar = pcritical(LATTICE);
    
    factor = ipow(2, P_SCHEDULE_POWER);
    
    time = (double)i/(double)(NSTEPS-1);
    
    if (time > 0.5) return(pstar + factor*(1.0 - pstar)*ipow(time - 0.5, P_SCHEDULE_POWER));
    else return(pstar - factor*pstar*ipow(0.5 - time, P_SCHEDULE_POWER));
}

int in_plot_box(double x, double y)
{
    int pos[2];
    static double xmin, ymin;
    static int first = 1;
    
    if (first)
    {
        xy_to_ij(XMAX - 1.0, YMAX - 1.0, pos);
        xmin = (double)pos[0];
        ymin = (double)pos[1];
        first = 0;
    }
    
    return((x > xmin)&&(y > ymin));
}

int in_plot_box_screencoord(double x, double y)
{
    static double xmin, ymin;
    static int first = 1;
    
    if (first)
    {
        xmin = XMAX - 1.0;
        ymin = YMAX - 1.0;
        first = 0;
    }
    
    return((x > xmin)&&(y > ymin));
}

double size_ratio_color(int clustersize, int ncells)
/* color of cell as function of the size of its cluster */
{
    double ratio, minratio = 1.0e-2, x, p = 0.1;
    
    minratio = 1.0/(double)ncells;
    
    ratio = (double)clustersize/(double)ncells;
    if (ratio > 1.0) ratio = 1.0;
    else if (ratio < minratio) ratio = minratio;
//     x = log(ratio/minratio)/log(1.0/minratio);
//     x = log(1.0 + log(ratio/minratio))/log(1.0 - log(minratio));
//     x = pow(log(ratio/minratio)/(-log(minratio)), 0.1);
    x = (pow(ratio, p) - pow(minratio, p))/(1.0 - pow(minratio, p));
    return(CLUSTER_HUEMIN*x + CLUSTER_HUEMAX*(1.0 -x));
    
    /* other attempts that seem to bug */
//     ratio = log((double)(clustersize+1))/log((double)(ncells+1));
// //     ratio = ((double)clustersize)/((double)ncells);
// //     ratio = sqrt((double)clustersize)/sqrt((double)ncells);
//     if (ratio > 1.0) ratio = 1.0;
//     else if (ratio < 0.0) ratio = 0.0;
//     return(CLUSTER_HUEMAX*ratio + CLUSTER_HUEMIN*(1.0 -ratio));
}

void compute_cell_color(t_perco cell, int *cluster_sizes, int fade, int max_cluster_size, int kx, int nx, int kz, int nz, double rgb[3])
/* compute color of cell */
{
    int k, color, csize;
    double fade_factor = 0.15, hue;
    
    if (!cell.open) hsl_to_rgb_palette(HUE_CLOSED, 0.9, 0.5, rgb, COLOR_PALETTE);
    else
        {
            if (!cell.flooded) hsl_to_rgb_palette(HUE_OPEN, 0.9, 0.5, rgb, COLOR_PALETTE);
            else 
            {
                if (COLOR_CELLS_BY_XCOORD)
                {
                    hue = CLUSTER_HUEMIN + (CLUSTER_HUEMAX - CLUSTER_HUEMIN)*(double)kx/(double)nx;
                    hsl_to_rgb_palette(hue, 0.9, 0.5, rgb, COLOR_PALETTE);
                }
                else if (COLOR_CELLS_BY_ZCOORD)
                {
                    hue = CLUSTER_HUEMIN + (CLUSTER_HUEMAX - CLUSTER_HUEMIN)*(double)kz/(double)nz;
                    hsl_to_rgb_palette(hue, 0.9, 0.5, rgb, COLOR_PALETTE);
                }
                else hsl_to_rgb_palette(HUE_FLOODED, 0.9, 0.5, rgb, COLOR_PALETTE);
            }
            
            if ((FIND_ALL_CLUSTERS)&&(COLOR_CLUSTERS_BY_SIZE))
            {
                csize = cluster_sizes[cell.cluster]; 
                hue = size_ratio_color(csize, max_cluster_size);
                hsl_to_rgb_palette(hue, 0.9, 0.5, rgb, COLOR_PALETTE);
            }
            else if ((FIND_ALL_CLUSTERS)&&(cell.cluster > 1))
            {
                color = cell.cluster%N_CLUSTER_COLORS;
                hue = CLUSTER_HUEMIN + (CLUSTER_HUEMAX - CLUSTER_HUEMIN)*(double)color/(double)N_CLUSTER_COLORS;
                hsl_to_rgb_palette(hue, 0.9, 0.5, rgb, COLOR_PALETTE);
            }
            
//             if ((FLOOD_LEFT_BOUNDARY)&&(cell.flooded == 1)) hsl_to_rgb_palette(HUE_FLOODED, 0.9, 0.5, rgb, COLOR_PALETTE);
        }
        
    if (fade) for (k=0; k<3; k++) rgb[k] = 1.0 - fade_factor + fade_factor*rgb[k];
    
}

void set_cell_color(t_perco cell, int *cluster_sizes, int fade, int max_cluster_size, int i, int nx, int k, int nz)
/* set color of cell */
{
    double rgb[3];
    
    compute_cell_color(cell, cluster_sizes, fade, max_cluster_size, i, nx, k, nz, rgb);
    glColor3f(rgb[0], rgb[1], rgb[2]);
}

double plot_coord(double x, double xmin, double xmax)
{
    return(xmin + x*(xmax - xmin));
}

void draw_size_plot(double plot_cluster_size[NSTEPS], int i, double pcrit)
/* draw plot of cluster sizes in terms of p */
{
    int j;
    char message[100];
    static double xmin, xmax, ymin, ymax, xmid, ymid, dx, dy, plotxmin, plotxmax, plotymin, plotymax;
    double pos[2], x1, y1, x2, y2, rgb[3], x;
    static int first = 1;
    
    if (first)
    {
        xmin = XMAX - 0.95;
        xmax = XMAX - 0.05;
        ymin = YMAX - 0.95;
        ymax = YMAX - 0.05;
                
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
    
    rgb[0] = 1.0; rgb[1] = 1.0; rgb[2] = 1.0;
//     erase_area_rgb(xmid, ymid, dx, dy, rgb);
    
    glColor3f(0.0, 0.0, 0.0);
    glLineWidth(2);
    
    /* axes and labels */
    draw_line(plotxmin, plotymin, plotxmax + 0.05, plotymin);
    draw_line(plotxmin, plotymin, plotxmin, plotymax + 0.1);
    draw_line(plotxmin - 0.02, plotymax, plotxmin + 0.02, plotymax);
    x = plot_coord(pcrit, plotxmin, plotxmax);
    draw_line(x, plotymin, x, plotymax);
    draw_line(plotxmax, plotymin - 0.02, plotxmax, plotymin + 0.02);
    
    xy_to_pos(plotxmax + 0.06, plotymin - 0.03, pos);
    sprintf(message, "p");
    write_text_fixedwidth(pos[0], pos[1], message);
    
    xy_to_pos(x - 0.02, plotymin - 0.04, pos);
    sprintf(message, "pc");
    write_text_fixedwidth(pos[0], pos[1], message);
    
    hsl_to_rgb_palette(HUE_GRAPH_SIZE, 0.9, 0.5, rgb, COLOR_PALETTE);
    glColor3f(rgb[0], rgb[1], rgb[2]);
    
    xy_to_pos(plotxmax - 0.015, plotymin - 0.06, pos);
    sprintf(message, "1");
    write_text_fixedwidth(pos[0], pos[1], message);
    
    xy_to_pos(plotxmin + 0.02, plotymax + 0.1, pos);
    sprintf(message, "nflooded/nopen");
    write_text_fixedwidth(pos[0], pos[1], message);
    
    xy_to_pos(plotxmin - 0.05, plotymax - 0.01, pos);
    sprintf(message, "1");
    write_text_fixedwidth(pos[0], pos[1], message);
    
    /* plot */
    x1 = plotxmin;
    y1 = plotymin;
    for (j=0; j<i; j++)
    {
//         x2 = plot_coord((double)j/(double)NSTEPS, plotxmin, plotxmax);
        x2 = plot_coord(p_schedule(j), plotxmin, plotxmax);
        y2 = plot_coord(plot_cluster_size[j], plotymin, plotymax);
        
        draw_line(x1, y1, x2, y2);
        x1 = x2;
        y1 = y2;
    }
}

void draw_cluster_number_plot(int plot_cluster_number[NSTEPS], int max_number, int i, double pcrit)
/* draw plot of number of clusters in terms of p */
{
    int j;
    char message[100];
    static double xmin, xmax, ymin, ymax, xmid, ymid, dx, dy, plotxmin, plotxmax, plotymin, plotymax;
    double pos[2], x1, y1, x2, y2, rgb[3], x, y;
    static int first = 1;
    
    if (first)
    {
        xmin = XMAX - 0.95;
        xmax = XMAX - 0.05;
        ymin = YMAX - 0.95;
        ymax = YMAX - 0.05;
                
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
    
    rgb[0] = 1.0; rgb[1] = 1.0; rgb[2] = 1.0;
//     erase_area_rgb(xmid, ymid, dx, dy, rgb);
    
    glColor3f(0.0, 0.0, 0.0);
    glLineWidth(2);
    
    /* axes and labels */
    draw_line(plotxmin, plotymin, plotxmax + 0.05, plotymin);
    draw_line(plotxmin, plotymin, plotxmin, plotymax + 0.1);
    draw_line(plotxmin - 0.02, plotymax, plotxmin + 0.02, plotymax);
    x = plot_coord(pcrit, plotxmin, plotxmax);
    draw_line(x, plotymin, x, plotymax);
    draw_line(plotxmax, plotymin - 0.02, plotxmax, plotymin + 0.02);
    
    xy_to_pos(plotxmax + 0.06, plotymin - 0.03, pos);
    sprintf(message, "p");
    write_text_fixedwidth(pos[0], pos[1], message);
    
    xy_to_pos(x - 0.02, plotymin - 0.04, pos);
    sprintf(message, "pc");
    write_text_fixedwidth(pos[0], pos[1], message);
    
    hsl_to_rgb_palette(HUE_GRAPH_CLUSTERS, 0.9, 0.5, rgb, COLOR_PALETTE);
    glColor3f(rgb[0], rgb[1], rgb[2]);
    
    xy_to_pos(plotxmax - 0.015, plotymin - 0.06, pos);
    sprintf(message, "1");
    write_text_fixedwidth(pos[0], pos[1], message);
    
    xy_to_pos(plotxmin + 0.02, plotymax + 0.05, pos);
    sprintf(message, "clusters/cells");
    write_text_fixedwidth(pos[0], pos[1], message);
    
    xy_to_pos(plotxmin + 0.05, plotymax - 0.01, pos);
    sprintf(message, "%.2f", 1.0/(double)MAX_CLUSTER_NUMBER);
    write_text_fixedwidth(pos[0], pos[1], message);
    
    /* plot */
    x1 = plotxmin;
    y1 = plotymin;
    for (j=0; j<i; j++)
    {
//         x2 = plot_coord((double)j/(double)NSTEPS, plotxmin, plotxmax);
        x2 = plot_coord(p_schedule(j), plotxmin, plotxmax);
        y2 = plot_coord((double)plot_cluster_number[j]/(double)max_number, plotymin, plotymax);
        
        draw_line(x1, y1, x2, y2);
        x1 = x2;
        y1 = y2;
    }
}

void draw_cluster_histogram(int ncells, int *cluster_sizes, int maxclustersize, int maxclusterlabel)
/* draw histogram of cluster size distribution */
{   
    int i, bin, nbins, binwidth, *histo, maxheight = 0, csize, dh, n, max_x_axis, di;
    static int maxbins = HISTO_BINS;
    char message[100];
    static double xmin, xmax, ymin, ymax, xmid, ymid, dx, dy, plotxmin, plotxmax, plotymin, plotymax, x, y;
    double pos[2], x1, y1, x2, y2, hue, delta, logfactor = 5.0, logbinwidth;
    static int first = 1;
    
    if (first)
    {
        xmin = XMAX - 0.95;
        xmax = XMAX - 0.05;
        ymin = YMAX - 0.95;
        ymax = YMAX - 0.05;
                
        xmid = 0.5*(xmin + xmax);
        ymid = 0.5*(ymin + ymax);
        
        dx = 0.5*(xmax - xmin);
        dy = 0.5*(ymax - ymin);
        
        plotxmin = xmin + 0.15;
        plotxmax = xmax - 0.1;
        plotymin = ymin + 0.07;
        plotymax = ymax - 0.1;
        
        first = 0;
    }
        
    if (maxclustersize > 0)
    {
    
//     if (maxclustersize < maxbins) nbins = maxclustersize;
//     else 
        nbins = maxbins;
    
        histo = (int *)malloc(nbins*sizeof(int));
    
        for (bin = 0; bin < nbins; bin++) histo[bin] = 0;
    
//         binwidth = maxclustersize/nbins;
//         if (binwidth == 0) binwidth = 1;
        if (HISTO_X_LOG_SCALE) logbinwidth = log(logfactor*(double)maxclustersize)/nbins;
        else binwidth = maxclustersize/nbins + 1;
        
        if (binwidth < 1) binwidth = 1;
        if (logbinwidth < 0.2) logbinwidth = 0.2;
    
        printf("max cluster size = %i, binwidth = %i\n", maxclustersize, binwidth);
        
        /* compute histogram */
        for (i=CLUSTER_SHIFT; i<maxclusterlabel; i++) if (cluster_sizes[i] > 0)
        {
            if (HISTO_X_LOG_SCALE) 
            {
                bin = (int)(log(logfactor*(double)cluster_sizes[i])/logbinwidth) - 1;
                if (bin >= nbins) bin = nbins - 1;
                else if (bin < 0) bin = 0;
            }
            else bin = (cluster_sizes[i]-1)/binwidth;
//             printf("cluster size = %i, bin = %i\n", cluster_sizes[i], bin);
            histo[bin]++;
        }
        for (bin=0; bin<maxbins; bin++) if (histo[bin] > maxheight) maxheight = histo[bin];
        
        /* draw histogram */
        glColor3f(0.0, 0.0, 0.0);
        glLineWidth(2);
    
        x1 = plotxmin;
        y1 = plotymin;
        if (HISTO_X_LOG_SCALE) max_x_axis = log((double)maxclustersize); 
        else max_x_axis = maxclustersize;
        if (max_x_axis < HISTO_BINS) max_x_axis = HISTO_BINS;
        
        for (bin=0; bin < nbins; bin++)
        {
//             csize = bin*binwidth + binwidth/2;
            if (HISTO_X_LOG_SCALE) csize = (int)(exp((double)bin*logbinwidth)/logfactor);
            else csize = bin*binwidth;
            if (csize >= ncells) csize = ncells-1;
            hue = size_ratio_color(csize, ncells);
            
            x2 = plot_coord((double)((bin+1)*binwidth)/(double)(max_x_axis), plotxmin, plotxmax);
//             x2 = plot_coord((double)(bin+1)/(double)nbins, plotxmin, plotxmax);
            
            y = log((double)(histo[bin]+1))/log((double)(maxheight+1));
            y2 = plot_coord(y, plotymin, plotymax);
        
//             printf("x1 = %.2f, x2 %.2f\n", x1, x2);
            draw_colored_rectangle(x1, y1, x2, y2, hue);
            x1 = x2;
        }
        
        draw_line(plotxmin, plotymin, plotxmax + 0.05, plotymin);
        draw_line(plotxmin, plotymin, plotxmin, plotymax + 0.1);
        
        /* graduation of y axis */
        x = log((double)(maxheight+1))/log(10.0);
        for (i=0; i<(int)x + 1; i++)
        {
            n = ipowi(10, i);
            y = log((double)(n+1))/log((double)(maxheight+1));
            y1 = plot_coord(y, plotymin, plotymax);
            xy_to_pos(plotxmin - 0.1, y1, pos);
            draw_line(plotxmin - 0.02, y1, plotxmin + 0.02, y1);
            sprintf(message, "%i", n);
            xy_to_pos(plotxmin - 0.07 - 0.025*(double)i, y1 - 0.01, pos);
            write_text_fixedwidth(pos[0], pos[1], message);
        }
    
        /* graduation of x axis */
        if (HISTO_X_LOG_SCALE)
        {
            y = log(logfactor*(double)(maxclustersize+1))/log(10.0);
            printf("y = %.3lg\n", y);
            for (i=1; i<(int)y + 1; i++)
            {
                n = ipowi(10, i);
                x = log((double)(n+1))/(log(logfactor*(double)(maxclustersize+1)));
//                 printf("n = %i, x = %.3lg\n", n, x);
                x1 = plot_coord(x, plotxmin, plotxmax);
                xy_to_pos(x1, plotymin - 0.1, pos);
                draw_line(x1, plotymin - 0.02, x1, plotymin + 0.02);
                if (n <= 1000) sprintf(message, "%i", n); 
                else sprintf(message, "1e%i", i);
                xy_to_pos(x1 - 0.015, plotymin - 0.05, pos);
                write_text_fixedwidth(pos[0], pos[1], message);
            }
        }
        else
        {
            x = log((double)(max_x_axis+1))/log(10.0);
            n = ipowi(10, (int)x);
            y = (double)n/10.0;

            delta = plot_coord((double)n/((double)(max_x_axis+1)), plotxmin, plotxmax) - plotxmin;
            if (delta > 0.13 + 0.01*x) di = 1;
            else if (delta > 0.08 + 0.01*x) di = 2;
            else di = 5;
            for (i=di; i<10; i+=di)
            {
                x1 = plot_coord(y*(double)i*10.0/(double)(max_x_axis+1), plotxmin, plotxmax);
                if (i*n < max_x_axis*11/10)
                {
                    sprintf(message, "%i", i*n);
                    xy_to_pos(x1 + 0.005 - 0.012*x, plotymin - 0.07, pos);
                    write_text_fixedwidth(pos[0], pos[1], message);
                }
            }
        
            delta = plot_coord((double)n/(10.0*(double)(max_x_axis+1)), plotxmin, plotxmax) - plotxmin;
            for (i=0; i<100; i++) if (i*n < 11*max_x_axis)
            {
                y = (double)(i*n)/10.0;
                x1 = plot_coord(y/(double)(max_x_axis+1), plotxmin, plotxmax);
                xy_to_pos(x1, plotymin , pos);
                if (i%10 == 0) draw_line(x1, plotymin - 0.02, x1, plotymin + 0.02);
                else if (delta > 0.02) draw_line(x1, plotymin - 0.01, x1, plotymin + 0.01);
            }
        }
        
        /* for debugging */
        if (DEBUG) for (bin=0; bin < nbins; bin++) printf("Bin %i = %i\n", bin, histo[bin]);
        
        
        free(histo);
    }
}

void draw_cell(int c, int nx, int ny, int size)
/* draw a single cell, for debugging puposes */
{
    int i, j, ishift, k, group;
    double rgb[3], x, y, alpha, r, fade = 0, dsize;
    static double h;
    static int first = 1;
    
    dsize = (double)size;
    
    group = cell_to_ij(c, &i, &j, nx, ny);
    switch (graphical_rep(LATTICE)) {
        case (PLOT_SQUARES):
        {
            glBegin(GL_QUADS);
    
            glVertex2i(i*size, j*size);
            glVertex2i((i+1)*size, j*size);
            glVertex2i((i+1)*size, (j+1)*size);
            glVertex2i(i*size, (j+1)*size);

            glEnd ();
            
            break;
        }
        case (PLOT_SQUARE_BONDS):
        {
            glLineWidth(5);
            glBegin(GL_LINES);
            
            /* horizontal segments */
            if (group == 0)
            {
                glVertex2i(i*size, j*size);
                glVertex2i((i+1)*size, j*size);
            }
            
            /* vertical segments */
            else 
            {
                glVertex2i(i*size, j*size);
                glVertex2i(i*size, (j+1)*size);
            }
                
            glEnd ();
            break;
        }
        case (PLOT_HEX):
        {
            if (first)
            {
                h = 0.5*sqrt(3.0);
                first = 0;
            }
            r = (double)size*0.5/h;
    
            x = ((double)i + 0.75)*dsize;
            if (j%2 == 1) x += 0.5*dsize;
            y = h*dsize*((double)j + 1.0);
            
            glBegin(GL_TRIANGLE_FAN);
            glVertex2d(x, y);
            for (k=0; k<7; k++)
            {
                alpha = (1.0 + 2.0*(double)k)*PI/6.0;
                glVertex2d(x + r*cos(alpha), y + r*sin(alpha));
            }
            glEnd();
            break;
        }
        case (PLOT_HEX_BONDS):
        {
            if (first)
            {
                h = 0.5*sqrt(3.0);
                first = 0;
            }
            r = (double)size*0.5/h;
            
            x = ((double)i + 0.75)*dsize;
            if (j%2 == 1) x -= 0.5*dsize;
            y = h*dsize*((double)j + 1.0);
            
            glBegin(GL_LINES);
            switch (group){
                case (0):   /* vertical bonds */
                {
                    glVertex2d(x-0.5*dsize, y+0.5*r);
                    glVertex2d(x-0.5*dsize, y-0.5*r);
                    break;
                }
                case (1):   /* NE-SW bonds */
                {
                    glVertex2d(x, y-r);
                    glVertex2d(x+0.5*dsize, y-0.5*r);
                    break;
                }
                case (2):   /* NW-SE bonds */
                {
                    glVertex2d(x-0.5*dsize, y-0.5*r);
                    glVertex2d(x, y-r);
                    break;
                }
            }

            glEnd();
            break;
        }
        case (PLOT_TRIANGLE):
        {
            if (first)
            {
                h = 0.5*sqrt(3.0);
                first = 0;
            }
            r = (double)size*0.5/h;
    
            x = 0.5*((double)i + 1.25)*dsize;
            y = h*dsize*((double)j);
            
            printf("Drawing cell %i = (%i, %i) at (%.0f, %.0f)\n", c, i, j, x, y);
                    
            glBegin(GL_TRIANGLES);
            if ((i+j)%2 == 1)
            {
                glVertex2d(x-0.5*dsize, y);
                glVertex2d(x+0.5*dsize, y);
                glVertex2d(x, y+h*dsize);
            }
            else
            {
                glVertex2d(x-0.5*dsize, y+h*dsize);
                glVertex2d(x+0.5*dsize, y+h*dsize);
                glVertex2d(x, y);
            }
            glEnd();
            break;
        }
        case (PLOT_POISSON_DISC):
        {
            /* beta version, TO DO */
//             draw_circle();
        }   
    }
}

void test_neighbours(int start, t_perco *cell, int nx, int ny, int size, int ncells)
/* for debugging puposes */
{
    int i, k;
    
    for (i=start; i<ncells; i++)
    {
        printf("Testing cell %i of %i - %i neighbours\n", i, ncells, cell[i].nneighb);
        blank();
        glColor3f(1.0, 0.0, 0.0);
        draw_cell(i, nx, ny, size);
        glColor3f(0.0, 0.0, 1.0);
        for (k=0; k<cell[i].nneighb; k++) draw_cell(cell[i].nghb[k], nx, ny, size);
        glutSwapBuffers();
        sleep(1);
    }
}   


void draw_cube_ijk(int i, int j, int k, t_perco *cell, int *cluster_sizes, int nx, int ny, int nz, int size, int max_cluster_size)
/* draw one cube of 3d configuration */
{
    double dx, x, y, z, rgb[3];
    
    dx = 1.0/(double)nx;
    
    compute_cell_color(cell[k*nx*ny+j*nx+i], cluster_sizes, 0, max_cluster_size, i, nx, k, nz, rgb);
                        
    x = (double)i*dx - 0.5;
    y = (double)j*dx - 0.5;
    z = (double)k*dx - 0.5;
    draw_cube(x, y, z, dx, rgb);
}
 
int plot_cube(t_perco cell)
/* returns 1 if cube is plotted */
{
    if (PLOT_ONLY_FLOODED_CELLS) return(cell.flooded);
    else return(cell.open);
}

void draw_configuration(t_perco *cell, int *cluster_sizes, int ncells, int nx, int ny, int nz, int size, int max_cluster_size)
/* draw the configuration */
{
    int i, j, ishift, k, n, sector;
    double rgb[3], x, y, z, alpha, r, fade = 0, dsize, radius, x2, y2, dx, dy, dz, observer_angle;
    static double h, h1;
    static int first = 1;
    
    dsize = (double)size;
    
    switch (graphical_rep(LATTICE)) {
        case (PLOT_SQUARES):
        {
            blank();
            
            glBegin(GL_QUADS);
    
            for (i=0; i<nx; i++)
                for (j=0; j<ny; j++)
                {
                    if ((ADD_PLOT)&&(in_plot_box((double)(i+1)*dsize, (double)(j)*dsize))) fade = 1;
                    else fade = 0;
                    
                    set_cell_color(cell[i*ny+j], cluster_sizes, fade, max_cluster_size, 0, 1, 0, 1);
            
                    glVertex2i(i*size, j*size);
                    glVertex2i((i+1)*size, j*size);
                    glVertex2i((i+1)*size, (j+1)*size);
                    glVertex2i(i*size, (j+1)*size);
                }

            glEnd ();
            break;
        }
        
        case (PLOT_SQUARE_BONDS):
        {
            ishift = nx*(ny+1);
            
            blank();
       
            if (ADD_PLOT)
            {
                rgb[0] = 1.0; rgb[1] = 1.0; rgb[2] = 1.0;
                erase_area_rgb(XMAX - 0.5, YMAX - 0.5, 0.5, 0.5, rgb);
            }
            
            if (size < 8) glLineWidth(1 + size/4);
            else glLineWidth(3);
            
            glBegin(GL_LINES);
            
            /* horizontal segments */
            for (i=0; i<nx; i++)
                for (j=0; j<ny+1; j++)
                {
                    if ((ADD_PLOT)&&(in_plot_box((double)(i+1)*dsize, (double)(j)*dsize))) fade = 1;
                    else fade = 0;
                    
                    set_cell_color(cell[i+nx*j], cluster_sizes, fade, max_cluster_size, 0, 1, 0, 1);
                    
                    glVertex2i(i*size, j*size);
                    glVertex2i((i+1)*size, j*size);
                }
                
            /* vertical segments */
            for (i=0; i<nx+1; i++)
                for (j=0; j<ny; j++)
                {
                    if ((ADD_PLOT)&&(in_plot_box((double)(i+1)*dsize, (double)(j+1)*dsize))) fade = 1;
                    else fade = 0;
                    
                    set_cell_color(cell[ishift+ny*i+j], cluster_sizes, fade, max_cluster_size, 0, 1, 0, 1);
                    
                    glVertex2i(i*size, j*size);
                    glVertex2i(i*size, (j+1)*size);
                }
                
            glEnd ();
            break;
        }
        
        case (PLOT_HEX):
        {
            blank();
            
            if (first)
            {
                h = 0.5*sqrt(3.0);
                first = 0;
            }
            r = (double)size*0.5/h;
    
            for (j=0; j<ny; j++)
                for (i=0; i<nx; i++)
                {
                    
                    x = ((double)i + 0.75)*dsize;
                    if (j%2 == 1) x += 0.5*dsize;
                    y = h*dsize*((double)j + 1.0);
                    
                    if ((ADD_PLOT)&&(in_plot_box(x, y))) fade = 1;
                    else fade = 0;
                    
                    set_cell_color(cell[j*nx+i], cluster_sizes, fade, max_cluster_size, 0, 1, 0, 1);
                    
                    glBegin(GL_TRIANGLE_FAN);
                    
                    glVertex2d(x, y);
                    for (k=0; k<7; k++)
                    {
                        alpha = (1.0 + 2.0*(double)k)*PI/6.0;
                        glVertex2d(x + r*cos(alpha), y + r*sin(alpha));
                    }
                    glEnd();
                }
            break;
        }
        case (PLOT_HEX_BONDS):
        {
            blank();
            
            ishift = nx*(ny+1);
            
            if (first)
            {
                h = 0.5*sqrt(3.0);
                first = 0;
            }
            r = (double)size*0.5/h;
    
            if (ADD_PLOT)
            {
                rgb[0] = 1.0; rgb[1] = 1.0; rgb[2] = 1.0;
                erase_area_rgb(XMAX - 0.5, YMAX - 0.5, 0.5, 0.5, rgb);
            }
            
            if (size < 8) glLineWidth(1 + size/4);
            else glLineWidth(3);
            
            glBegin(GL_LINES);
            for (i=0; i<nx; i++)
                for (j=0; j<ny; j++)
                {
                    x = ((double)i + 0.5)*dsize;
                    if (j%2 == 1) x -= 0.5*dsize;
                    y = h*dsize*((double)j + 0.75);
                    
                    if ((ADD_PLOT)&&(in_plot_box(x, y))) fade = 1;
                    else fade = 0;
                    
                    /* vertical bonds */
                    if (j<ny-1)
                    {
                        set_cell_color(cell[i+nx*j], cluster_sizes, fade, max_cluster_size, 0, 1, 0, 1);
                    
                        glVertex2d(x-0.5*dsize, y+0.5*r);
                        glVertex2d(x-0.5*dsize, y-0.5*r);
                    }
                    
                    if ((ADD_PLOT)&&(in_plot_box(x, y-0.5*h*dsize))) fade = 1;
                    else fade = 0;
                    
                    /* NW-SE bonds */
                    set_cell_color(cell[2*ishift + i+nx*j], cluster_sizes, fade, max_cluster_size, 0, 1, 0, 1);
                    
                    glVertex2d(x-0.5*dsize, y-0.5*r);
                    glVertex2d(x, y-r);
                    
                    if ((ADD_PLOT)&&(in_plot_box(x+0.5*dsize, y-0.5*h*dsize))) fade = 1;
                    else fade = 0;
                    
                    /* NE-SW bonds */
                    set_cell_color(cell[ishift + i+nx*j], cluster_sizes, fade, max_cluster_size, 0, 1, 0, 1);
                    
                    glVertex2d(x, y-r);
                    glVertex2d(x+0.5*dsize, y-0.5*r);

                    
                }
           
            glEnd();
            break;
        }
        case (PLOT_TRIANGLE):
        {
            blank();
            
            if (first)
            {
                h = 0.5*sqrt(3.0);
                h1 = 0.5/sqrt(3.0);
                first = 0;
            }
            r = (double)size*0.5/h;
    
            for (j=0; j<ny; j++)
                for (i=0; i<2*nx; i++)
                {
                    
                    x = 0.5*((double)i + 1.25)*dsize;
                    y = h*dsize*((double)j);
                    
                    if ((ADD_PLOT)&&(in_plot_box(x, y+0.5*h*dsize))) fade = 1;
                    else fade = 0;
                    
                    set_cell_color(cell[2*j*nx+i], cluster_sizes, fade, max_cluster_size, 0, 1, 0, 1);
                    
                    glBegin(GL_TRIANGLES);
                    if ((i+j)%2 == 1)
                    {
                        glVertex2d(x-0.5*dsize, y);
                        glVertex2d(x+0.5*dsize, y);
                        glVertex2d(x, y+h*dsize);
                    }
                    else
                    {
                        glVertex2d(x-0.5*dsize, y+h*dsize);
                        glVertex2d(x+0.5*dsize, y+h*dsize);
                        glVertex2d(x, y);
                    }
                    glEnd();
                }
            break;
        }
        case (PLOT_POISSON_DISC):
        {
            radius = sqrt((XMAX - XMIN)*(YMAX - YMIN)/(PI*(double)nx));; 
            
            blank();
            
            if (ADD_PLOT)
            {
                rgb[0] = 1.0; rgb[1] = 1.0; rgb[2] = 1.0;
                erase_area_rgb(XMAX - 0.5, YMAX - 0.5, 0.5, 0.5, rgb);
            }
            
            for (i=0; i<ncells; i++)
            {
                x = cell[i].x;
                y = cell[i].y;
                
                if ((ADD_PLOT)&&(in_plot_box_screencoord(x, y))) fade = 1;
                else fade = 0;
                    
                set_cell_color(cell[i], cluster_sizes, fade, max_cluster_size, 0, 1, 0, 1);
                
                for (k=0; k<cell[i].nneighb; k++)
                {
                    n = cell[i].nghb[k];
//                     printf("Cell %i, %ith neighbour cell %i\n", i, k, n);
                    x2 = cell[n].x;
                    y2 = cell[n].y;
                    draw_line(x, y, 0.5*(x + x2), 0.5*(y + y2));
                }
                    
                draw_colored_circle(x, y, radius, NSEG);
            }
            break;
        }
        
        case (PLOT_CUBES):  /* beta version */
        {
            blank();
            
//             glBegin(GL_QUADS);
            
            if (ADD_PLOT)
            {
                rgb[0] = 1.0; rgb[1] = 1.0; rgb[2] = 1.0;
                erase_area_rgb(XMAX - 0.5, YMAX - 0.5, 0.5, 0.5, rgb);
            }
            
            for (k=0; k<nz; k++) 
            {
//                 if (ROTATE_VIEW)
                {
                    observer_angle = argument(observer[0], observer[1]);
//                     observer_angle += 0.1*PID;
                    if (observer_angle < 0.0) observer_angle += DPI;
                    sector = (int)(observer_angle*2.0/PID);
//                     printf("Observer_angle = %.3lg\n", observer_angle*360.0/DPI); 
                    
                    switch (sector) {
                        case (0):
                        {
                            for (i=0; i<nx; i++) 
                                for (j=0; j<ny; j++) if (plot_cube(cell[k*nx*ny+j*nx+i]))
                                    draw_cube_ijk(i, j, k, cell, cluster_sizes, nx, ny, nz, size, max_cluster_size);
                            break;
                        }
                        case (1):
                        {
                            for (j=0; j<ny; j++) 
                                for (i=0; i<nx; i++) if (plot_cube(cell[k*nx*ny+j*nx+i]))
                                    draw_cube_ijk(i, j, k, cell, cluster_sizes, nx, ny, nz, size, max_cluster_size);
                            break;
                        }
                        case (2):
                        {
                            for (j=0; j<ny; j++) 
                                for (i=nx-1; i>=0; i--) if (plot_cube(cell[k*nx*ny+j*nx+i]))
                                    draw_cube_ijk(i, j, k, cell, cluster_sizes, nx, ny, nz, size, max_cluster_size);
                            break;
                        }
                        case (3):
                        {
                            for (i=nx-1; i>= 0; i--) 
                                for (j=0; j<ny; j++) if (plot_cube(cell[k*nx*ny+j*nx+i]))
                                    draw_cube_ijk(i, j, k, cell, cluster_sizes, nx, ny, nz, size, max_cluster_size);
                            break;
                        }
                        case (4):
                        {
                            for (i=nx-1; i>= 0; i--) 
                                for (j=ny-1; j>=0; j--) if (plot_cube(cell[k*nx*ny+j*nx+i]))
                                    draw_cube_ijk(i, j, k, cell, cluster_sizes, nx, ny, nz, size, max_cluster_size);
                            break;
                        }
                        case (5):
                        {
                            for (j=ny-1; j>=0; j--) 
                                for (i=nx-1; i>=0; i--) if (plot_cube(cell[k*nx*ny+j*nx+i]))
                                    draw_cube_ijk(i, j, k, cell, cluster_sizes, nx, ny, nz, size, max_cluster_size);
                            break;
                        }
                        case (6):
                        {
                            for (j=ny-1; j>=0; j--) 
                                for (i=0; i<nx; i++) if (plot_cube(cell[k*nx*ny+j*nx+i]))
                                    draw_cube_ijk(i, j, k, cell, cluster_sizes, nx, ny, nz, size, max_cluster_size);
                            break;
                        }
                        case (7):
                        {
                            for (i=0; i<nx; i++)  
                                for (j=ny-1; j>=0; j--) if (plot_cube(cell[k*nx*ny+j*nx+i]))
                                    draw_cube_ijk(i, j, k, cell, cluster_sizes, nx, ny, nz, size, max_cluster_size);
                            break;
                        }
                    }
                }
//                 else
//                 {
//                     for (i=0; i < nx; i++)
//                         for (j=0; j<ny; j++) if (plot_cube(cell[k*nx*ny+j*nx+i]))
//                             draw_cube_ijk(i, j, k, cell, cluster_sizes, nx, ny, nz, size, max_cluster_size);
//                 }
                
            }
            break;
        }
        
        
    }
}


void print_p(double p)
{
    char message[100];
    double y = YMAX - 0.13, pos[2], rgb[3] = {1.0, 1.0, 1.0};
    static double xleftbox, xlefttext, xrightbox, xrighttext;
    static int first = 1;
    
    if (first)
    {
        xleftbox = XMIN + 0.2;
        xlefttext = xleftbox - 0.45;
        if (PLOT_CLUSTER_HISTOGRAM) xrightbox = XMAX - 0.41;
        else xrightbox = XMAX - 0.27;
        xrighttext = xrightbox - 0.45;
        first = 0;
    }
    
    if (!PLOT_CLUSTER_SIZE)
    {
        if (NX > 1280) erase_area_rgb(xrightbox, y + 0.025, 0.15, 0.05, rgb);
        else erase_area_rgb(xrightbox, y + 0.025, 0.22, 0.05, rgb);
    }
    glColor3f(0.0, 0.0, 0.0);
    if (NX > 1280) xy_to_pos(xrighttext + 0.35, y, pos);
    else xy_to_pos(xrighttext + 0.3, y, pos);
    sprintf(message, "p = %.4f", p);
    write_text(pos[0], pos[1], message);
}

void print_nclusters(int nclusters)
{
    char message[100];
    double y = YMAX - 0.25, pos[2], rgb[3] = {1.0, 1.0, 1.0};
    static double xleftbox, xlefttext, xrightbox, xrighttext;
    static int first = 1;
    
    if (first)
    {
        xleftbox = XMIN + 0.2;
        xlefttext = xleftbox - 0.45;
        xrightbox = XMAX - 0.31;
        xrighttext = xrightbox - 0.48;
        first = 0;
    }
    
    if (!PLOT_CLUSTER_SIZE)
    {
        if (NX > 1280) erase_area_rgb(xrightbox, y + 0.025, 0.18, 0.05, rgb);
        else erase_area_rgb(xrightbox, y + 0.025, 0.25, 0.05, rgb);
    }
    glColor3f(0.0, 0.0, 0.0);
    if (NX > 1280) xy_to_pos(xrighttext + 0.35, y, pos);
    else xy_to_pos(xrighttext + 0.3, y, pos);
    if (nclusters == 1) sprintf(message, "%i cluster", nclusters);
    else sprintf(message, "%i clusters", nclusters);
    write_text_fixedwidth(pos[0], pos[1], message);
}

void print_largest_cluster_size(int max_cluster_size)
{
    char message[100];
    double y = YMAX - 0.25, pos[2], rgb[3] = {1.0, 1.0, 1.0};
    static double xleftbox, xlefttext, xrightbox, xrighttext;
    static int first = 1;
    
    if (first)
    {
        xleftbox = XMIN + 0.2;
        xlefttext = xleftbox - 0.45;
        xrightbox = XMAX - 0.41;
        xrighttext = xrightbox - 0.48;
        first = 0;
    }
    
    if (!PLOT_CLUSTER_SIZE)
    {
        if (NX > 1280) erase_area_rgb(xrightbox, y + 0.025, 0.18, 0.05, rgb);
        else erase_area_rgb(xrightbox, y + 0.025, 0.25, 0.05, rgb);
    }
    glColor3f(0.0, 0.0, 0.0);
    if (NX > 1280) xy_to_pos(xrighttext + 0.35, y, pos);
    else xy_to_pos(xrighttext + 0.3, y, pos);
    sprintf(message, "max size %i", max_cluster_size);
    write_text_fixedwidth(pos[0], pos[1], message);
}


/*********************/
/* animation part    */
/*********************/

int generate_poisson_discs(t_perco *cell, double dpoisson, int nmaxcells)
/* generate Poisson disc configuration */
{
    double x, y, r, phi;
    int i, j, k, n_p_active, ncircles, ncandidates=NPOISSON_CANDIDATES, naccepted; 
    short int *active_poisson, far;
    
    active_poisson = (short int *)malloc(2*(nmaxcells)*sizeof(short int));
    
    printf("Generating Poisson disc sample\n");
    /* generate first circle */
    cell[0].x = (XMAX - XMIN)*(double)rand()/RAND_MAX + XMIN;
    cell[0].y = (YMAX - YMIN)*(double)rand()/RAND_MAX + YMIN;
    active_poisson[0] = 1;
    n_p_active = 1;
    ncircles = 1;
        
    while ((n_p_active > 0)&&(ncircles < nmaxcells))
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
            x = cell[i].x + r*cos(phi);
            y = cell[i].y + r*sin(phi);
//                       printf("Testing new circle at (%.3f,%.3f)\t", x, y);
            far = 1;
            for (k=0; k<ncircles; k++) if ((k!=i))
            {
                /* new circle is far away from circle k */
                far = far*((x - cell[k].x)*(x - cell[k].x) + (y - cell[k].y)*(y - cell[k].y) >= dpoisson*dpoisson);
                /* new circle is in domain */
                far = far*(x < XMAX)*(x > XMIN)*(y < YMAX)*(y > YMIN);
            }
            if (far)    /* accept new circle */
            {
                printf("New circle at (%.3f,%.3f) accepted\n", x, y);
                cell[ncircles].x = x;
                cell[ncircles].y = y;
                cell[ncircles].active = 1;
                active_poisson[ncircles] = 1;
                ncircles++;
                n_p_active++;
                naccepted++;
            }
//                        else printf("Rejected\n");
        }
        if (naccepted == 0)    /* inactivate circle i */ 
        {
            active_poisson[i] = 0;
            n_p_active--;
        }
        printf("%i active circles\n", n_p_active);
    }
        
    printf("Generated %i circles\n", ncircles);
    
    free(active_poisson);
    
    return(ncircles);
}


int cell_number(int nx, int ny, int nz)
/* total number of cells in graph */
{
    switch (LATTICE) {
        case (BC_SQUARE_DIRICHLET): return(nx*ny);
        case (BC_SQUARE_PERIODIC): return(nx*ny);
        case (BC_SQUARE_BOND_DIRICHLET): return(2*nx*ny + nx + ny);
        case (BC_HEX_SITE_DIRICHLET): return((int)((double)(nx*ny)*2.0/sqrt(3.0))); /* hex lattice requires more vertical space! */
        case (BC_HEX_BOND_DIRICHLET): return(3*(int)((double)((nx+2)*(ny+2))*2.0/sqrt(3.0))); 
        /* hex lattice requires more vertical space! */
        case (BC_TRIANGLE_SITE_DIRICHLET): return((int)((double)(2*nx*ny)*2.0/sqrt(3.0))); /* triangle lattice requires more vertical space! */
        case (BC_POISSON_DISC): return(nx*ny);   /* TO IMPROVE */
        case (BC_CUBIC_DIRICHLET): return(nx*ny*nz);
    }
}


void compute_nxnynz(int size, int *nx, int *ny, int *nz)
/* compute the number of rows and columns depending on lattice */
{
    switch (LATTICE) {
        case (BC_HEX_SITE_DIRICHLET):
        {
            *nx = NX/size - 1;
            *ny = (int)((double)NY*2.0/((double)size*sqrt(3.0)));
            *nz = 1;
            break;
        }
        case (BC_HEX_BOND_DIRICHLET):
        {
            *nx = NX/size + 1;
            *ny = (int)((double)NY*2.0/((double)size*sqrt(3.0))) + 1;
            *nz = 1;
            break;
        }
        case (BC_TRIANGLE_SITE_DIRICHLET):
        {
            *nx = NX/size - 1;
            *ny = (int)((double)NY*2.0/((double)size*sqrt(3.0))) + 1;
            *nz = 1;
            break;
        }
        case (BC_POISSON_DISC):
        {
            /* for Poisson disc configuration, use a 1d labelling */
            *nx = NX*NY/(size*size);
            *ny = 1;
            *nz = 1;
            break;
        }
        case (BC_CUBIC_DIRICHLET):
        {
            *nx = NX/size;
            *ny = NY/size;
            *nz = NZ/size;
            break;
        }
        default:
        {
            *nx = NX/size;
            *ny = NY/size;
            *nz = 1;
        }
    }
    
}

int init_cell_lattice(t_perco *cell, int nx, int ny, int nz)
/* initialize the graph of connected cells - returns the number of cells */
{
    int i, j, k, iplus, iminus, ishift, n;
    int ncells;
    double dpoisson, radius;
    
    printf("Initializing cell lattice ...");
    
    switch (LATTICE) {
        case (BC_SQUARE_DIRICHLET):
        {
            /* neighbours in the bulk */
            #pragma omp parallel for private(i,j)
            for (i=1; i<nx-1; i++){
                for (j=1; j<ny-1; j++){
                    n = cellnb(i, j, 0, nx, ny);
                    cell[n].nneighb = 4;
                    cell[n].nghb[0] = cellnb(i+1, j, 0, nx, ny);
                    cell[n].nghb[1] = cellnb(i-1, j, 0, nx, ny);
                    cell[n].nghb[2] = cellnb(i, j+1, 0, nx, ny);
                    cell[n].nghb[3] = cellnb(i, j-1, 0, nx, ny);
                }
            }
    
            /* left boundary */
            #pragma omp parallel for private(j)
            for (j=1; j<ny-1; j++)
            {
                n = cellnb(0, j, 0, nx, ny);
                cell[n].nneighb = 3;
                cell[n].nghb[0] = cellnb(1, j, 0, nx, ny);
                cell[n].nghb[1] = cellnb(0, j+1, 0, nx, ny);
                cell[n].nghb[2] = cellnb(0, j-1, 0, nx, ny);
            }
            /* right boundary */
            #pragma omp parallel for private(j)
            for (j=1; j<ny-1; j++)
            {
                n = cellnb(nx-1, j, 0, nx, ny);
                cell[n].nneighb = 3;
                cell[n].nghb[0] = cellnb(nx-2, j, 0, nx, ny);
                cell[n].nghb[1] = cellnb(nx-1, j+1, 0, nx, ny);
                cell[n].nghb[2] = cellnb(nx-1, j-1, 0, nx, ny);
            }
            
            /* top boundary */
            #pragma omp parallel for private(i)
            for (i=1; i<nx-1; i++)
            {
                n = cellnb(i, ny-1, 0, nx, ny);
                cell[n].nneighb = 3;
                cell[n].nghb[0] = cellnb(i+1, ny-1, 0, nx, ny);
                cell[n].nghb[1] = cellnb(i-1, ny-1, 0, nx, ny);
                cell[n].nghb[2] = cellnb(i, ny-2, 0, nx, ny);
            }
            
            /* bottom boundary */
            #pragma omp parallel for private(i)
            for (i=1; i<nx-1; i++)
            {
                n = cellnb(i, 0, 0, nx, ny);
                cell[n].nneighb = 3;
                cell[n].nghb[0] = cellnb(i+1, 0, 0, nx, ny);
                cell[n].nghb[1] = cellnb(i-1, 0, 0, nx, ny);
                cell[n].nghb[2] = cellnb(i, 1, 0, nx, ny);
            }
            
            /* corners */
            n = cellnb(0, 0, 0, nx, ny);
            cell[n].nneighb = 2;
            cell[n].nghb[0] = cellnb(1, 0, 0, nx, ny);
            cell[n].nghb[1] = cellnb(0, 1, 0, nx, ny);
            
            n = cellnb(0, ny-1, 0, nx, ny);
            cell[n].nneighb = 2;
            cell[n].nghb[0] = cellnb(0, ny-2, 0, nx, ny);
            cell[n].nghb[1] = cellnb(1, ny-1, 0, nx, ny);            
            
            n = cellnb(nx-1, 0, 0, nx, ny);
            cell[n].nneighb = 2;
            cell[n].nghb[0] = cellnb(nx-1, 1, 0, nx, ny);
            cell[n].nghb[1] = cellnb(nx-2, 0, 0, nx, ny);
            
            n = cellnb(nx-1, ny-1, 0, nx, ny);
            cell[n].nneighb = 2;
            cell[n].nghb[0] = cellnb(nx-1, ny-2, 0, nx, ny);
            cell[n].nghb[1] = cellnb(nx-2, ny-1, 0, nx, ny);
            
            return(nx*ny);
            break;
        }
        case (BC_SQUARE_PERIODIC):
        {
            /* neighbours in the bulk */
            #pragma omp parallel for private(i,j)
            for (i=1; i<nx-1; i++){
                for (j=1; j<ny-1; j++){
                    n = cellnb(i, j, 0, nx, ny);
                    cell[n].nneighb = 4;
                    cell[n].nghb[0] = cellnb(i+1, j, 0, nx, ny);
                    cell[n].nghb[1] = cellnb(i-1, j, 0, nx, ny);
                    cell[n].nghb[2] = cellnb(i, j+1, 0, nx, ny);
                    cell[n].nghb[3] = cellnb(i, j-1, 0, nx, ny);
                }
            }
    
            /* left boundary */
            #pragma omp parallel for private(j)
            for (j=1; j<ny-1; j++)
            {
                n = cellnb(0, j, 0, nx, ny);
                cell[n].nneighb = 4;
                cell[n].nghb[0] = cellnb(1, j, 0, nx, ny);
                cell[n].nghb[1] = cellnb(nx-1, j, 0, nx, ny);
                cell[n].nghb[2] = cellnb(0, j+1, 0, nx, ny);
                cell[n].nghb[3] = cellnb(0, j+1, 0, nx, ny);
            }
            
            /* right boundary */
            #pragma omp parallel for private(j)
            for (j=1; j<ny-1; j++)
            {
                n = cellnb(nx-1, j, 0, nx, ny);
                cell[n].nneighb = 4;
                cell[n].nghb[0] = cellnb(nx-2, j, 0, nx, ny);
                cell[n].nghb[1] = cellnb(0, j, 0, nx, ny);
                cell[n].nghb[2] = cellnb(nx-1, j+1, 0, nx, ny);
                cell[n].nghb[3] = cellnb(nx-1, j-1, 0, nx, ny);
            }
            
            /* top boundary */
            #pragma omp parallel for private(i,iplus,iminus)
            for (i=0; i<nx; i++)
            {
                n = cellnb(i, ny-1, 0, nx, ny);
                iplus = (i+1) % nx;
                iminus = (i-1) % nx;
                if (iminus < 0) iminus += nx;
                    
                cell[n].nneighb = 4;
                cell[n].nghb[0] = cellnb(iplus, ny-1, 0, nx, ny);
                cell[n].nghb[1] = cellnb(iminus, ny-1, 0, nx, ny);
                cell[n].nghb[2] = cellnb(i, ny-2, 0, nx, ny);
                cell[n].nghb[3] = cellnb(i, 0, 0, nx, ny);
            }
            
            /* bottom boundary */
            #pragma omp parallel for private(i,iplus,iminus)
            for (i=0; i<nx; i++)
            {
                n = cellnb(i, 0, 0, nx, ny);
                iplus = (i+1) % nx;
                iminus = (i-1) % nx;
                if (iminus < 0) iminus += nx;
                    
                cell[n].nneighb = 4;
                cell[n].nghb[0] = cellnb(iplus, 0, 0, nx, ny);
                cell[n].nghb[1] = cellnb(iminus, 0, 0, nx, ny);
                cell[n].nghb[2] = cellnb(i, 1, 0, nx, ny);
                cell[n].nghb[3] = cellnb(i, ny-1, 0, nx, ny);
            }
            
            return(nx*ny);
            break;
        }
        case (BC_SQUARE_BOND_DIRICHLET): 
        {
            ishift = nx*(ny+1); 
            
            /* horizontal bonds */
            /* neighbours in the bulk  */
            #pragma omp parallel for private(i,j)
            for (i=1; i<nx-1; i++){
                for (j=1; j<ny; j++){
                    n = cellnb(i, j, 0, nx, ny);
                    cell[n].nneighb = 6;
                    cell[n].nghb[0] = cellnb(i+1, j, 0, nx, ny);
                    cell[n].nghb[1] = cellnb(i-1, j, 0, nx, ny);
                    cell[n].nghb[2] = cellnb(i, j-1, 1, nx, ny);
                    cell[n].nghb[3] = cellnb(i, j, 1, nx, ny);
                    cell[n].nghb[4] = cellnb(i+1, j-1, 1, nx, ny);
                    cell[n].nghb[5] = cellnb(i+1, j, 1, nx, ny);
                }
            }
            
            /* bottom boundary */
            #pragma omp parallel for private(i)
            for (i=1; i<nx-1; i++)
            {
                n = cellnb(i, 0, 0, nx, ny);
                cell[n].nneighb = 4;
                cell[n].nghb[0] = cellnb(i+1, 0, 0, nx, ny);
                cell[n].nghb[1] = cellnb(i-1, 0, 0, nx, ny);
                cell[n].nghb[2] = cellnb(i, 0, 1, nx, ny);
                cell[n].nghb[3] = cellnb(i+1, 0, 1, nx, ny);
            }

            /* top boundary */
            #pragma omp parallel for private(i)
            for (i=1; i<nx-1; i++)
            {
                n = cellnb(i, ny, 0, nx, ny);
                cell[n].nneighb = 4;
                cell[n].nghb[0] = cellnb(i+1, ny, 0, nx, ny);
                cell[n].nghb[1] = cellnb(i-1, ny, 0, nx, ny);
                cell[n].nghb[2] = cellnb(i, ny-1, 1, nx, ny);
                cell[n].nghb[3] = cellnb(i+1, ny-1, 1, nx, ny);
            }
            
            /* left boundary */
            #pragma omp parallel for private(j)
            for (j=1; j<ny; j++){
                n = cellnb(0, j, 0, nx, ny);
                cell[n].nneighb = 5;
                cell[n].nghb[0] = cellnb(1, j, 0, nx, ny);
                cell[n].nghb[1] = cellnb(0, j-1, 1, nx, ny);
                cell[n].nghb[2] = cellnb(0, j, 1, nx, ny);
                cell[n].nghb[3] = cellnb(1, j-1, 1, nx, ny);
                cell[n].nghb[4] = cellnb(1, j, 1, nx, ny);
            }
            
            /* right boundary */
            #pragma omp parallel for private(j)
            for (j=1; j<ny; j++){
                n = cellnb(nx-1, j, 0, nx, ny);
                cell[n].nneighb = 5;
                cell[n].nghb[0] = cellnb(nx-2, j, 0, nx, ny);
                cell[n].nghb[1] = cellnb(nx-1, j-1, 1, nx, ny);
                cell[n].nghb[2] = cellnb(nx-1, j, 1, nx, ny);
                cell[n].nghb[3] = cellnb(nx, j-1, 1, nx, ny);
                cell[n].nghb[4] = cellnb(nx, j, 1, nx, ny);
            }
            
            /* corners */
            n = cellnb(0, 0, 0, nx, ny);
            cell[n].nneighb = 3;
            cell[n].nghb[0] = cellnb(1, 0, 0, nx, ny);
            cell[n].nghb[1] = cellnb(0, 0, 1, nx, ny);
            cell[n].nghb[2] = cellnb(1, 0, 1, nx, ny);
            
            n = cellnb(nx-1, 0, 0, nx, ny);
            cell[n].nneighb = 3;
            cell[n].nghb[0] = cellnb(nx-2, 0, 0, nx, ny);
            cell[n].nghb[1] = cellnb(nx-1, 0, 1, nx, ny);
            cell[n].nghb[2] = cellnb(nx, 0, 1, nx, ny);
            
            n = cellnb(0, ny, 0, nx, ny);
            cell[n].nneighb = 3;
            cell[n].nghb[0] = cellnb(1, ny, 0, nx, ny);
            cell[n].nghb[1] = cellnb(0, ny-1, 1, nx, ny);
            cell[n].nghb[2] = cellnb(1, ny-1, 1, nx, ny);
            
            n = cellnb(nx-1, ny, 0, nx, ny);
            cell[n].nneighb = 3;
            cell[n].nghb[0] = cellnb(nx-2, ny, 0, nx, ny);
            cell[n].nghb[1] = cellnb(nx-1, ny-1, 1, nx, ny);
            cell[n].nghb[2] = cellnb(nx, ny-1, 1, nx, ny);

            /* vertical bonds */
            /* neighbours in the bulk */
            #pragma omp parallel for private(i,j)
            for (i=1; i<nx; i++){
                for (j=1; j<ny-1; j++){
                    n = cellnb(i, j, 1, nx, ny);
                    cell[n].nneighb = 6;
                    cell[n].nghb[0] = cellnb(i, j-1, 1, nx, ny);
                    cell[n].nghb[1] = cellnb(i, j+1, 1, nx, ny);
                    cell[n].nghb[2] = cellnb(i, j, 0, nx, ny);
                    cell[n].nghb[3] = cellnb(i-1, j, 0, nx, ny);
                    cell[n].nghb[4] = cellnb(i, j+1, 0, nx, ny);
                    cell[n].nghb[5] = cellnb(i-1, j+1, 0, nx, ny);
                }
            }
            
            /* left boundary */
            #pragma omp parallel for private(j)
            for (j=1; j<ny-1; j++){
                n = cellnb(0, j, 1, nx, ny);
                cell[n].nneighb = 4;
                cell[n].nghb[0] = cellnb(0, j-1, 1, nx, ny);
                cell[n].nghb[1] = cellnb(0, j+1, 1, nx, ny);
                cell[n].nghb[2] = cellnb(0, j, 0, nx, ny);
                cell[n].nghb[3] = cellnb(0, j+1, 0, nx, ny);
            }
            
            /* right boundary */
            #pragma omp parallel for private(j)
            for (j=1; j<ny-1; j++){
                n = cellnb(nx, j, 1, nx, ny);
                cell[n].nneighb = 4;
                cell[n].nghb[0] = cellnb(nx, j-1, 1, nx, ny);
                cell[n].nghb[1] = cellnb(nx, j+1, 1, nx, ny);
                cell[n].nghb[2] = cellnb(nx-1, j, 0, nx, ny);
                cell[n].nghb[3] = cellnb(nx-1, j+1, 0, nx, ny);
            }
            
            /* bottom boundary */
            #pragma omp parallel for private(i)
            for (i=1; i<nx; i++){
                n = cellnb(i, 0, 1, nx, ny);
                cell[n].nneighb = 5;
                cell[n].nghb[0] = cellnb(i, 1, 1, nx, ny);
                cell[n].nghb[1] = cellnb(i-1, 0, 0, nx, ny);
                cell[n].nghb[2] = cellnb(i, 0, 0, nx, ny);
                cell[n].nghb[3] = cellnb(i-1, 1, 0, nx, ny);
                cell[n].nghb[4] = cellnb(i, 1, 0, nx, ny);
            }
            
            /* top boundary */
            #pragma omp parallel for private(i)
            for (i=1; i<nx; i++){
                n = cellnb(i, ny-1, 1, nx, ny);
                cell[n].nneighb = 5;
                cell[n].nghb[0] = cellnb(i, ny-2, 1, nx, ny);
                cell[n].nghb[1] = cellnb(i-1, ny-1, 0, nx, ny);
                cell[n].nghb[2] = cellnb(i, ny-1, 0, nx, ny);
                cell[n].nghb[3] = cellnb(i-1, ny, 0, nx, ny);
                cell[n].nghb[4] = cellnb(i, ny, 0, nx, ny);
            }
            
            /* corners */
            n = cellnb(0, 0, 1, nx, ny);
            cell[n].nneighb = 3;
            cell[n].nghb[0] = cellnb(0, 1, 1, nx, ny);
            cell[n].nghb[1] = cellnb(0, 0, 0, nx, ny);
            cell[n].nghb[2] = cellnb(0, 1, 0, nx, ny);
            
            n = cellnb(0, ny-1, 1, nx, ny);
            cell[n].nneighb = 3;
            cell[n].nghb[0] = cellnb(0, ny-2, 1, nx, ny);
            cell[n].nghb[1] = cellnb(0, ny-1, 0, nx, ny);
            cell[n].nghb[2] = cellnb(0, ny, 0, nx, ny);

            n = cellnb(nx, 0, 1, nx, ny);
            cell[n].nneighb = 3;
            cell[n].nghb[0] = cellnb(nx, 1, 1, nx, ny);
            cell[n].nghb[1] = cellnb(nx-1, 0, 0, nx, ny);
            cell[n].nghb[2] = cellnb(nx-1, 1, 0, nx, ny);
                    
            n = cellnb(nx, ny-1, 1, nx, ny);
            cell[n].nneighb = 3;
            cell[n].nghb[0] = cellnb(nx, ny-2, 1, nx, ny);
            cell[n].nghb[1] = cellnb(nx-1, ny, 0, nx, ny);
            cell[n].nghb[2] = cellnb(nx-1, ny-1, 0, nx, ny);
            
            return(nx + ny + 2*nx*ny);
            break;
        }
        case (BC_HEX_SITE_DIRICHLET):
        {
            /* neighbours in the bulk  */
            #pragma omp parallel for private(i,j)
            for (i=1; i<nx-1; i++){
                for (j=1; j<ny-1; j+=2){
                    n = cellnb(i, j, 0, nx, ny);
                    cell[n].nneighb = 6;
                    cell[n].nghb[0] = cellnb(i+1, j, 0, nx, ny);
                    cell[n].nghb[1] = cellnb(i-1, j, 0, nx, ny);
                    cell[n].nghb[2] = cellnb(i+1, j+1, 0, nx, ny);
                    cell[n].nghb[3] = cellnb(i, j+1, 0, nx, ny);
                    cell[n].nghb[4] = cellnb(i+1, j-1,0, nx, ny);
                    cell[n].nghb[5] = cellnb(i, j-1, 0, nx, ny);
                }
            }
            #pragma omp parallel for private(i,j)
            for (i=1; i<nx-1; i++){
                for (j=2; j<ny-1; j+=2){
                    n = cellnb(i, j, 0, nx, ny);
                    cell[n].nneighb = 6;
                    cell[n].nghb[0] = cellnb(i+1, j, 0, nx, ny);
                    cell[n].nghb[1] = cellnb(i-1, j, 0, nx, ny);
                    cell[n].nghb[2] = cellnb(i-1, j+1, 0, nx, ny);
                    cell[n].nghb[3] = cellnb(i, j+1, 0, nx, ny);
                    cell[n].nghb[4] = cellnb(i-1, j-1,0, nx, ny);
                    cell[n].nghb[5] = cellnb(i, j-1, 0, nx, ny);
                }
            }
            
            /* left boundary */
            #pragma omp parallel for private(j)
            for (j=1; j<ny-1; j+=2){
                n = cellnb(0, j, 0, nx, ny);
                cell[n].nneighb = 5;
                cell[n].nghb[0] = cellnb(1, j, 0, nx, ny);
                cell[n].nghb[1] = cellnb(0, j+1, 0, nx, ny);
                cell[n].nghb[2] = cellnb(1, j+1, 0, nx, ny);
                cell[n].nghb[3] = cellnb(0, j-1, 0, nx, ny);
                cell[n].nghb[4] = cellnb(1, j-1, 0, nx, ny);
            }
            #pragma omp parallel for private(j)
            for (j=2; j<ny-1; j+=2){
                n = cellnb(0, j, 0, nx, ny);
                cell[n].nneighb = 3;
                cell[n].nghb[0] = cellnb(1, j, 0, nx, ny);
                cell[n].nghb[1] = cellnb(0, j+1, 0, nx, ny);
                cell[n].nghb[2] = cellnb(0, j-1, 0, nx, ny);
            }
            
            /* right boundary */
            #pragma omp parallel for private(j)
            for (j=1; j<ny-1; j+=2){
                n = cellnb(nx-1, j, 0, nx, ny);
                cell[n].nneighb = 3;
                cell[n].nghb[0] = cellnb(nx-2, j, 0, nx, ny);
                cell[n].nghb[1] = cellnb(nx-1, j+1, 0, nx, ny);
                cell[n].nghb[2] = cellnb(nx-1, j-1, 0, nx, ny);
            }
            #pragma omp parallel for private(j)
            for (j=2; j<ny-1; j+=2){
                n = cellnb(nx-1, j, 0, nx, ny);
                cell[n].nneighb = 5;
                cell[n].nghb[0] = cellnb(nx-2, j, 0, nx, ny);
                cell[n].nghb[1] = cellnb(nx-2, j+1, 0, nx, ny);
                cell[n].nghb[2] = cellnb(nx-1, j+1, 0, nx, ny);
                cell[n].nghb[3] = cellnb(nx-2, j-1, 0, nx, ny);
                cell[n].nghb[4] = cellnb(nx-1, j-1, 0, nx, ny);
            }
            
            /* bottom boundary */
            #pragma omp parallel for private(i)
            for (i=1; i<nx-1; i++){
                n = cellnb(i, 0, 0, nx, ny);
                cell[n].nneighb = 4;
                cell[n].nghb[0] = cellnb(i+1, 0, 0, nx, ny);
                cell[n].nghb[1] = cellnb(i-1, 0, 0, nx, ny);
                cell[n].nghb[2] = cellnb(i-1, 1, 0, nx, ny);
                cell[n].nghb[3] = cellnb(i, 1, 0, nx, ny);
            }
            
            /* top boundary */
            #pragma omp parallel for private(i)
            for (i=1; i<nx-1; i++){
                n = cellnb(i, ny-1, 0, nx, ny);
                cell[n].nneighb = 4;
                cell[n].nghb[0] = cellnb(i+1, ny-1, 0, nx, ny);
                cell[n].nghb[1] = cellnb(i-1, ny-1, 0, nx, ny);
                cell[n].nghb[2] = cellnb(i-1, ny-2, 0, nx, ny);
                cell[n].nghb[3] = cellnb(i, ny-2, 0, nx, ny);
            }
            
            /* corners */
            n = cellnb(0, 0, 0, nx, ny);
            cell[n].nneighb = 2;
            cell[n].nghb[0] = cellnb(1, 0, 0, nx, ny);
            cell[n].nghb[1] = cellnb(0, 1, 0, nx, ny);
            
            n = cellnb(nx-1, 0, 0, nx, ny);
            cell[n].nneighb = 3;
            cell[n].nghb[0] = cellnb(nx-2, 0, 0, nx, ny);
            cell[n].nghb[1] = cellnb(nx-1, 1, 0, nx, ny);
            cell[n].nghb[2] = cellnb(nx-2, 1, 0, nx, ny);
            
            n = cellnb(0, ny-1, 0, nx, ny);
            cell[n].nneighb = 2;
            cell[n].nghb[0] = cellnb(1, ny-1, 0, nx, ny);
            cell[n].nghb[1] = cellnb(0, ny-2, 0, nx, ny);
            
            n = cellnb(nx-1, ny-1, 0, nx, ny);
            cell[n].nneighb = 3;
            cell[n].nghb[0] = cellnb(nx-2, ny-1, 0, nx, ny);
            cell[n].nghb[1] = cellnb(nx-1, ny-2, 0, nx, ny);
            cell[n].nghb[2] = cellnb(nx-2, ny-2, 0, nx, ny);
            
            return(nx*ny);
            break;
        }
        case (BC_HEX_BOND_DIRICHLET):
        {
            /* neighbours in the bulk */
            /* vertical bonds */
            #pragma omp parallel for private(i,j)
            for (i=1; i<nx; i++){
                for (j=1; j<ny; j+=2){
                    n = cellnb(i, j, 0, nx, ny);
                    cell[n].nneighb = 4;
                    cell[n].nghb[0] = cellnb(i-1, j, 1, nx, ny);
                    cell[n].nghb[1] = cellnb(i, j, 2, nx, ny);
                    cell[n].nghb[2] = cellnb(i-1, j+1, 1, nx, ny);
                    cell[n].nghb[3] = cellnb(i-1, j+1, 2, nx, ny);
                }
            }
            #pragma omp parallel for private(i,j)
            for (i=1; i<nx; i++){
                for (j=0; j<ny; j+=2){
                    n = cellnb(i, j, 0, nx, ny);
                    cell[n].nneighb = 4;
                    cell[n].nghb[0] = cellnb(i-1, j, 1, nx, ny);
                    cell[n].nghb[1] = cellnb(i, j, 2, nx, ny);
                    cell[n].nghb[2] = cellnb(i, j+1, 1, nx, ny);
                    cell[n].nghb[3] = cellnb(i, j+1, 2, nx, ny);
                }
            }
            /* NE-SW bonds */
            #pragma omp parallel for private(i,j)
            for (i=0; i<nx-1; i++){
                for (j=1; j<ny-1; j+=2){
                    n = cellnb(i, j, 1, nx, ny);
                    cell[n].nneighb = 4;
                    cell[n].nghb[0] = cellnb(i+1, j, 0, nx, ny);
                    cell[n].nghb[1] = cellnb(i+1, j, 2, nx, ny);
                    cell[n].nghb[2] = cellnb(i, j-1, 0, nx, ny);
                    cell[n].nghb[3] = cellnb(i, j, 2, nx, ny);
                }
            }
            #pragma omp parallel for private(i,j)
            for (i=0; i<nx-1; i++){
                for (j=2; j<ny-1; j+=2){
                    n = cellnb(i, j, 1, nx, ny);
                    cell[n].nneighb = 4;
                    cell[n].nghb[0] = cellnb(i+1, j, 0, nx, ny);
                    cell[n].nghb[1] = cellnb(i+1, j, 2, nx, ny);
                    cell[n].nghb[2] = cellnb(i+1, j-1, 0, nx, ny);
                    cell[n].nghb[3] = cellnb(i, j, 2, nx, ny);
                }
            }
            
            /* NW-SE bonds */
            #pragma omp parallel for private(i,j)
            for (i=1; i<nx-1; i++){
                for (j=1; j<ny-1; j+=2){
                    n = cellnb(i, j, 2, nx, ny);
                    cell[n].nneighb = 4;
                    cell[n].nghb[0] = cellnb(i, j, 0, nx, ny);
                    cell[n].nghb[1] = cellnb(i, j, 1, nx, ny);
                    cell[n].nghb[2] = cellnb(i, j-1, 0, nx, ny);
                    cell[n].nghb[3] = cellnb(i-1, j, 1, nx, ny);
                }
            }
            #pragma omp parallel for private(i,j)
            for (i=0; i<nx; i++){
                for (j=2; j<ny-1; j+=2){
                    n = cellnb(i, j, 2, nx, ny);
                    cell[n].nneighb = 4;
                    cell[n].nghb[0] = cellnb(i, j, 0, nx, ny);
                    cell[n].nghb[1] = cellnb(i, j, 1, nx, ny);
                    cell[n].nghb[2] = cellnb(i+1, j-1, 0, nx, ny);
                    cell[n].nghb[3] = cellnb(i-1, j, 1, nx, ny);
                }
            }
            
            /* left boundary */
            /* vertical bonds */
            #pragma omp parallel for private(j)
            for (j=0; j<ny; j+=2)
            {
                n = cellnb(0, j, 0, nx, ny);
                cell[n].nneighb = 2;
                cell[n].nghb[0] = cellnb(0, j, 2, nx, ny);
                cell[n].nghb[1] = cellnb(0, j+1, 1, nx, ny);
            }
            /* NE-SW bonds */
            #pragma omp parallel for private(j)
            for (j=1; j<ny-1; j+=2)
            {
                n = cellnb(0, j, 1, nx, ny);
                cell[n].nneighb = 3;
                cell[n].nghb[0] = cellnb(1, j, 0, nx, ny);
                cell[n].nghb[1] = cellnb(0, j-1, 0, nx, ny);
                cell[n].nghb[2] = cellnb(1, j, 2, nx, ny);
            }
            /* NW-SE bonds */
            #pragma omp parallel for private(j)
            for (j=2; j<ny-1; j+=2)
            {
                n = cellnb(0, j, 2, nx, ny);
                cell[n].nneighb = 3;
                cell[n].nghb[0] = cellnb(0, j, 0, nx, ny);
                cell[n].nghb[1] = cellnb(1, j-1, 0, nx, ny);
                cell[n].nghb[2] = cellnb(0, j, 1, nx, ny);
            }
            
            /* right boundary */
            /* vertical bonds */
            #pragma omp parallel for private(j)
            for (j=0; j<ny; j+=2)
            {
                n = cellnb(nx-1, j, 0, nx, ny);
                cell[n].nneighb = 2;
                cell[n].nghb[0] = cellnb(nx-1, j+1, 2, nx, ny);
                cell[n].nghb[1] = cellnb(nx-2, j, 1, nx, ny);
            }
            /* NE-SW bonds */
            #pragma omp parallel for private(j)
            for (j=1; j<ny-1; j+=2)
            {
                n = cellnb(nx-1, j, 1, nx, ny);
                cell[n].nneighb = 3;
                cell[n].nghb[0] = cellnb(nx, j, 0, nx, ny);
                cell[n].nghb[1] = cellnb(nx-1, j-1, 0, nx, ny);
                cell[n].nghb[2] = cellnb(nx-1, j, 2, nx, ny);
            }
            #pragma omp parallel for private(j)
            for (j=2; j<ny-1; j+=2)
            {
                n = cellnb(nx-1, j, 1, nx, ny);
                cell[n].nneighb = 3;
                cell[n].nghb[0] = cellnb(nx, j, 0, nx, ny);
                cell[n].nghb[1] = cellnb(nx, j, 2, nx, ny);
                cell[n].nghb[2] = cellnb(nx, j-1, 0, nx, ny);
            }
            /* NW-SE bonds */
            #pragma omp parallel for private(j)
            for (j=1; j<ny-1; j+=2)
            {
                n = cellnb(nx-1, j, 2, nx, ny);
                cell[n].nneighb = 3;
                cell[n].nghb[0] = cellnb(nx-1, j-1, 0, nx, ny);
                cell[n].nghb[1] = cellnb(nx-1, j, 0, nx, ny);
                cell[n].nghb[2] = cellnb(nx-2, j, 1, nx, ny);
            }
            #pragma omp parallel for private(j)
            for (j=2; j<ny-1; j+=2)
            {
                n = cellnb(nx-1, j, 2, nx, ny);
                cell[n].nneighb = 3;
                cell[n].nghb[0] = cellnb(nx-1, j-1, 0, nx, ny);
                cell[n].nghb[1] = cellnb(nx-1, j, 0, nx, ny);
                cell[n].nghb[2] = cellnb(nx-2, j, 1, nx, ny);
            }
            
            /* bottom boundary */
            /* NE-SW bonds */
            #pragma omp parallel for private(i)
            for (i=1; i<nx-1; i++)
            {
                n = cellnb(i, 0, 1, nx, ny);
                cell[n].nneighb = 3;
                cell[n].nghb[0] = cellnb(i+1, 0, 0, nx, ny);
                cell[n].nghb[1] = cellnb(i, 0, 2, nx, ny);
                cell[n].nghb[2] = cellnb(i+1, 0, 2, nx, ny);
            }
            /* NW-SE bonds */
            #pragma omp parallel for private(i)
            for (i=1; i<nx-1; i++)
            {
                n = cellnb(i, 0, 2, nx, ny);
                cell[n].nneighb = 3;
                cell[n].nghb[0] = cellnb(i, 0, 0, nx, ny);
                cell[n].nghb[1] = cellnb(i-1, 0, 1, nx, ny);
                cell[n].nghb[2] = cellnb(i, 0, 1, nx, ny);
            }
            
            /* top boundary */
            /* vertical bonds */
            #pragma omp parallel for private(i)
            for (i=1; i<nx-1; i++)
            {
                n = cellnb(i, ny-1, 0, nx, ny);
                cell[n].nneighb = 2;
                cell[n].nghb[0] = cellnb(i-1, ny-1, 1, nx, ny);
                cell[n].nghb[1] = cellnb(i, ny-1, 2, nx, ny);
            }
            /* NE-SW bonds */
            #pragma omp parallel for private(i)
            for (i=1; i<nx-1; i++)
            {
                n = cellnb(i, ny-1, 1, nx, ny);
                cell[n].nneighb = 3;
                cell[n].nghb[0] = cellnb(i+1, ny-1, 2, nx, ny);
                cell[n].nghb[1] = cellnb(i, ny-1, 2, nx, ny);
                cell[n].nghb[2] = cellnb(i+1, ny-2, 0, nx, ny);
            }
            /* NW-SE bonds */
            #pragma omp parallel for private(i)
            for (i=1; i<nx-1; i++)
            {
                n = cellnb(i, ny-1, 2, nx, ny);
                cell[n].nneighb = 3;
                cell[n].nghb[0] = cellnb(i-1, ny-1, 1, nx, ny);
                cell[n].nghb[1] = cellnb(i, ny-1, 1, nx, ny);
                cell[n].nghb[2] = cellnb(i+1, ny-2, 0, nx, ny);
            }
            
            /* corners */
            n = cellnb(0, 0, 2, nx, ny);
            cell[n].nneighb = 2;
            cell[n].nghb[0] = cellnb(0, 0, 0, nx, ny);
            cell[n].nghb[1] = cellnb(0, 0, 1, nx, ny);
            
            return(3*(nx+1)*(ny+1));    /* TO BE CHECKED */
            break;
        }
        case (BC_TRIANGLE_SITE_DIRICHLET):
        {
            /* neighbours in the bulk */
            #pragma omp parallel for private(i,j)
            for (i=1; i<2*nx-1; i++){
                for (j=1; j<ny-1; j++){
                    n = cellnb(i, j, 0, nx, ny);
                    cell[n].nneighb = 3;
                    cell[n].nghb[0] = cellnb(i+1,j,0,nx,ny);
                    cell[n].nghb[1] = cellnb(i-1,j,0,nx,ny);
                    if ((i+j)%2 == 0) cell[n].nghb[2] = cellnb(i,j+1,0,nx,ny);
                    else cell[n].nghb[2] = cellnb(i,j-1,0,nx,ny);
                }
            }
            
            /* left boundary */
            #pragma omp parallel for private(j)
            for (j=0; j<ny; j++){
                n = cellnb(0, j, 0, nx, ny);
                cell[n].nneighb = 2;
                cell[n].nghb[0] = cellnb(1, j, 0, nx, ny);
                if (j%2 == 0) cell[n].nghb[1] = cellnb(0, j+1, 0, nx, ny);
                else cell[n].nghb[1] = cellnb(0, j-1, 0, nx, ny);
            }
            
            /* right boundary */
            #pragma omp parallel for private(j)
            for (j=1; j<ny-1; j++){
                n = cellnb(2*nx-1, j, 0, nx, ny);
                cell[n].nneighb = 2;
                cell[n].nghb[0] = cellnb(2*nx-2, j, 0, nx, ny);
                if (j%2 == 1) cell[n].nghb[1] = cellnb(2*nx-1, j+1, 0, nx, ny);
                else cell[n].nghb[1] = cellnb(2*nx-1, j-1, 0, nx, ny);
            }
            
            /* bottom boundary */
            #pragma omp parallel for private(i,j)
            for (i=1; i<2*nx-1; i++){
                n = cellnb(i, 0, 0, nx, ny);
                if (i%2 == 0)
                {
                    cell[n].nneighb = 3;
                    cell[n].nghb[0] = cellnb(i-1, 0, 0, nx, ny);
                    cell[n].nghb[1] = cellnb(i+1, 0, 0, nx, ny);
                    cell[n].nghb[2] = cellnb(i, 1, 0, nx, ny);
                }
                else 
                {
                    cell[n].nneighb = 2;
                    cell[n].nghb[0] = cellnb(i+1, 0, 0, nx, ny);
                    cell[n].nghb[1] = cellnb(i-1, 0, 0, nx, ny);
                }
            }
            
            /* top boundary */
            #pragma omp parallel for private(i,j)
            for (i=1; i<2*nx-1; i++){
                n = cellnb(i, ny-1, 0, nx, ny);
                if (i%2 == 1)
                {
                    cell[n].nneighb = 3;
                    cell[n].nghb[0] = cellnb(i-1, ny-1, 0, nx, ny);
                    cell[n].nghb[1] = cellnb(i+1, ny-1, 0, nx, ny);
                    cell[n].nghb[2] = cellnb(i, ny-2, 0, nx, ny);
                }
                else 
                {
                    cell[n].nneighb = 2;
                    cell[n].nghb[0] = cellnb(i-1, ny-1, 0, nx, ny);
                    cell[n].nghb[1] = cellnb(i+1, ny-1, 0, nx, ny);
                }
            }
            
            
            /* corners (at the right) */
            n = cellnb(2*nx-1, 0, 0, nx, ny);
            cell[n].nneighb = 1;
            cell[n].nghb[0] = cellnb(2*nx-2, 0, 0, nx, ny);
            
            n = cellnb(2*nx-1, ny-1, 0, nx, ny);
            cell[n].nneighb = 1;
            cell[n].nghb[0] = cellnb(2*nx-2, ny-1, 0, nx, ny);
            
            
            return(2*nx*ny);
            break;
        }
        case (BC_POISSON_DISC):
        {
            dpoisson = 3.0*sqrt((XMAX - XMIN)*(YMAX - YMIN)/(PI*(double)nx));
            
//             printf("nx = %i, dpoisson = %.5lg\n", nx, dpoisson);
            ncells = generate_poisson_discs(cell, dpoisson, nx);
            radius = 1.8*dpoisson; 
            
            for (i=0; i<ncells; i++)
            {
                k = 0;
                
                for (j=0; j<ncells; j++) if ((j != i)&&(k < MAX_NEIGHB))
                {
                    if (module2(cell[j].x - cell[i].x, cell[j].y - cell[i].y) < radius)
                    {
                        cell[i].nghb[k] = j;
//                         printf("%i ", j);
                        k++;
                    }
                }
                cell[i].nneighb = k;
            }
            
            return(ncells);
            break;
        }
        case (BC_CUBIC_DIRICHLET):
        {
            /* neighbours in the bulk */
            #pragma omp parallel for private(i,j,k)
            for (i=1; i<nx-1; i++){
                for (j=1; j<ny-1; j++){
                    for (k=1; k<nz-1; k++){
                        n = cellnb_3d(i, j, k, 0, nx, ny, nz);
                        cell[n].nneighb = 6;
                        cell[n].nghb[0] = cellnb_3d(i+1, j, k, 0, nx, ny, nz);
                        cell[n].nghb[1] = cellnb_3d(i-1, j, k, 0, nx, ny, nz);
                        cell[n].nghb[2] = cellnb_3d(i, j+1, k, 0, nx, ny, nz);
                        cell[n].nghb[3] = cellnb_3d(i, j-1, k, 0, nx, ny, nz);
                        cell[n].nghb[4] = cellnb_3d(i, j, k+1, 0, nx, ny, nz);
                        cell[n].nghb[5] = cellnb_3d(i, j, k-1, 0, nx, ny, nz);
                    }
                }
            }
            
            /* back and front face */
            #pragma omp parallel for private(j,k)
            for (j=1; j<ny-1; j++){
                for (k=1; k<nz-1; k++){
                    n = cellnb_3d(0, j, k, 0, nx, ny, nz);
                    cell[n].nneighb = 5;
                    cell[n].nghb[0] = cellnb_3d(0, j+1, k, 0, nx, ny, nz);
                    cell[n].nghb[1] = cellnb_3d(0, j-1, k, 0, nx, ny, nz);
                    cell[n].nghb[2] = cellnb_3d(0, j, k+1, 0, nx, ny, nz);
                    cell[n].nghb[3] = cellnb_3d(0, j, k-1, 0, nx, ny, nz);
                    cell[n].nghb[4] = cellnb_3d(1, j, k, 0, nx, ny, nz);
                    
                    n = cellnb_3d(nx-1, j, k, 0, nx, ny, nz);
                    cell[n].nneighb = 5;
                    cell[n].nghb[0] = cellnb_3d(nx-1, j+1, k, 0, nx, ny, nz);
                    cell[n].nghb[1] = cellnb_3d(nx-1, j-1, k, 0, nx, ny, nz);
                    cell[n].nghb[2] = cellnb_3d(nx-1, j, k+1, 0, nx, ny, nz);
                    cell[n].nghb[3] = cellnb_3d(nx-1, j, k-1, 0, nx, ny, nz);
                    cell[n].nghb[4] = cellnb_3d(nx-2, j, k, 0, nx, ny, nz);
                }
            }
            
            /* left and right face */
            #pragma omp parallel for private(i,k)
            for (i=1; i<nx-1; i++){
                for (k=1; k<nz-1; k++){
                    n = cellnb_3d(i, 0, k, 0, nx, ny, nz);
                    cell[n].nneighb = 5;
                    cell[n].nghb[0] = cellnb_3d(i+1, 0, k, 0, nx, ny, nz);
                    cell[n].nghb[1] = cellnb_3d(i-1, 0, k, 0, nx, ny, nz);
                    cell[n].nghb[2] = cellnb_3d(i, 0, k+1, 0, nx, ny, nz);
                    cell[n].nghb[3] = cellnb_3d(i, 0, k-1, 0, nx, ny, nz);
                    cell[n].nghb[4] = cellnb_3d(i, 1, k, 0, nx, ny, nz);
                    
                    n = cellnb_3d(i, ny-1, k, 0, nx, ny, nz);
                    cell[n].nneighb = 5;
                    cell[n].nghb[0] = cellnb_3d(i+1, ny-1, k, 0, nx, ny, nz);
                    cell[n].nghb[1] = cellnb_3d(i-1, ny-1, k, 0, nx, ny, nz);
                    cell[n].nghb[2] = cellnb_3d(i, ny-1, k+1, 0, nx, ny, nz);
                    cell[n].nghb[3] = cellnb_3d(i, ny-1, k-1, 0, nx, ny, nz);
                    cell[n].nghb[4] = cellnb_3d(i, ny-2, k, 0, nx, ny, nz);
                }
            }
            
            /* lower and upper face */
            #pragma omp parallel for private(i,j)
            for (i=1; i<nx-1; i++){
                for (j=1; j<ny-1; j++){
                    n = cellnb_3d(i, j, 0, 0, nx, ny, nz);
                    cell[n].nneighb = 5;
                    cell[n].nghb[0] = cellnb_3d(i+1, j, 0, 0, nx, ny, nz);
                    cell[n].nghb[1] = cellnb_3d(i-1, j, 0, 0, nx, ny, nz);
                    cell[n].nghb[2] = cellnb_3d(i, j+1, 0, 0, nx, ny, nz);
                    cell[n].nghb[3] = cellnb_3d(i, j-1, 0, 0, nx, ny, nz);
                    cell[n].nghb[4] = cellnb_3d(i, j, 1, 0, nx, ny, nz);
                    
                    n = cellnb_3d(i, j, nz-1, 0, nx, ny, nz);
                    cell[n].nneighb = 5;
                    cell[n].nghb[0] = cellnb_3d(i+1, j, nz-1, 0, nx, ny, nz);
                    cell[n].nghb[1] = cellnb_3d(i-1, j, nz-1, 0, nx, ny, nz);
                    cell[n].nghb[2] = cellnb_3d(i, j+1, nz-1, 0, nx, ny, nz);
                    cell[n].nghb[3] = cellnb_3d(i, j-1, nz-1, 0, nx, ny, nz);
                    cell[n].nghb[4] = cellnb_3d(i, j, nz-2, 0, nx, ny, nz);
                }
            }
            
            /* back to front edges */
            #pragma omp parallel for private(i)
            for (i=1; i<nx-1; i++){
                n = cellnb_3d(i, 0, 0, 0, nx, ny, nz);
                cell[n].nneighb = 4;
                cell[n].nghb[0] = cellnb_3d(i+1, 0, 0, 0, nx, ny, nz);
                cell[n].nghb[1] = cellnb_3d(i-1, 0, 0, 0, nx, ny, nz);
                cell[n].nghb[2] = cellnb_3d(i, 1, 0, 0, nx, ny, nz);
                cell[n].nghb[3] = cellnb_3d(i, 0, 1, 0, nx, ny, nz);
                    
                n = cellnb_3d(i, ny-1, 0, 0, nx, ny, nz);
                cell[n].nneighb = 4;
                cell[n].nghb[0] = cellnb_3d(i+1, ny-1, 0, 0, nx, ny, nz);
                cell[n].nghb[1] = cellnb_3d(i-1, ny-1, 0, 0, nx, ny, nz);
                cell[n].nghb[2] = cellnb_3d(i, ny-2, 0, 0, nx, ny, nz);
                cell[n].nghb[3] = cellnb_3d(i, ny-1, 1, 0, nx, ny, nz);
                    
                n = cellnb_3d(i, 0, nz-1, 0, nx, ny, nz);
                cell[n].nneighb = 4;
                cell[n].nghb[0] = cellnb_3d(i+1, 0, nz-1, 0, nx, ny, nz);
                cell[n].nghb[1] = cellnb_3d(i-1, 0, nz-1, 0, nx, ny, nz);
                cell[n].nghb[2] = cellnb_3d(i, 1, nz-1, 0, nx, ny, nz);
                cell[n].nghb[3] = cellnb_3d(i, 0, nz-2, 0, nx, ny, nz);
                    
                n = cellnb_3d(i, ny-1, nz-1, 0, nx, ny, nz);
                cell[n].nneighb = 4;
                cell[n].nghb[0] = cellnb_3d(i+1, ny-1, nz-1, 0, nx, ny, nz);
                cell[n].nghb[1] = cellnb_3d(i-1, ny-1, nz-1, 0, nx, ny, nz);
                cell[n].nghb[2] = cellnb_3d(i, ny-2, nz-1, 0, nx, ny, nz);
                cell[n].nghb[3] = cellnb_3d(i, ny-1, nz-2, 0, nx, ny, nz);
            }
            
            /* left to right edges */
            #pragma omp parallel for private(j)
            for (j=1; j<ny-1; j++){
                n = cellnb_3d(0, j, 0, 0, nx, ny, nz);
                cell[n].nneighb = 4;
                cell[n].nghb[0] = cellnb_3d(0, j+1, 0, 0, nx, ny, nz);
                cell[n].nghb[1] = cellnb_3d(0, j-1, 0, 0, nx, ny, nz);
                cell[n].nghb[2] = cellnb_3d(1, j, 0, 0, nx, ny, nz);
                cell[n].nghb[3] = cellnb_3d(0, j, 1, 0, nx, ny, nz);
                    
                n = cellnb_3d(nx-1, j, 0, 0, nx, ny, nz);
                cell[n].nneighb = 4;
                cell[n].nghb[0] = cellnb_3d(nx-1, j+1, 0, 0, nx, ny, nz);
                cell[n].nghb[1] = cellnb_3d(nx-1, j-1, 0, 0, nx, ny, nz);
                cell[n].nghb[2] = cellnb_3d(nx-2, j, 0, 0, nx, ny, nz);
                cell[n].nghb[3] = cellnb_3d(nx-1, j, 1, 0, nx, ny, nz);
                    
                n = cellnb_3d(0, j, nz-1, 0, nx, ny, nz);
                cell[n].nneighb = 4;
                cell[n].nghb[0] = cellnb_3d(0, j+1, nz-1, 0, nx, ny, nz);
                cell[n].nghb[1] = cellnb_3d(0, j-1, nz-1, 0, nx, ny, nz);
                cell[n].nghb[2] = cellnb_3d(1, j, nz-1, 0, nx, ny, nz);
                cell[n].nghb[3] = cellnb_3d(0, j, nz-2, 0, nx, ny, nz);
                    
                n = cellnb_3d(nx-1, j, nz-1, 0, nx, ny, nz);
                cell[n].nneighb = 4;
                cell[n].nghb[0] = cellnb_3d(nx-1, j+1, nz-1, 0, nx, ny, nz);
                cell[n].nghb[1] = cellnb_3d(nx-1, j-1, nz-1, 0, nx, ny, nz);
                cell[n].nghb[2] = cellnb_3d(nx-2, j, nz-1, 0, nx, ny, nz);
                cell[n].nghb[3] = cellnb_3d(nx-1, j, nz-2, 0, nx, ny, nz);
            }
            
            /* top to bottom edges */
            #pragma omp parallel for private(k)
            for (k=1; k<nz-1; k++){
                n = cellnb_3d(0, 0, k, 0, nx, ny, nz);
                cell[n].nneighb = 4;
                cell[n].nghb[0] = cellnb_3d(0, 0, k+1, 0, nx, ny, nz);
                cell[n].nghb[1] = cellnb_3d(0, 0, k-1, 0, nx, ny, nz);
                cell[n].nghb[2] = cellnb_3d(1, 0, k, 0, nx, ny, nz);
                cell[n].nghb[3] = cellnb_3d(0, 1, k, 0, nx, ny, nz);
                    
                n = cellnb_3d(nx-1, 0, k, 0, nx, ny, nz);
                cell[n].nneighb = 4;
                cell[n].nghb[0] = cellnb_3d(nx-1, 0, k+1, 0, nx, ny, nz);
                cell[n].nghb[1] = cellnb_3d(nx-1, 0, k-1, 0, nx, ny, nz);
                cell[n].nghb[2] = cellnb_3d(nx-2, 0, k, 0, nx, ny, nz);
                cell[n].nghb[3] = cellnb_3d(nx-1, 1, k, 0, nx, ny, nz);
                    
                n = cellnb_3d(0, ny-1, k, 0, nx, ny, nz);
                cell[n].nneighb = 4;
                cell[n].nghb[0] = cellnb_3d(0, ny-1, k+1, 0, nx, ny, nz);
                cell[n].nghb[1] = cellnb_3d(0, ny-1, k-1, 0, nx, ny, nz);
                cell[n].nghb[2] = cellnb_3d(1, ny-1, k, 0, nx, ny, nz);
                cell[n].nghb[3] = cellnb_3d(0, ny-2, k, 0, nx, ny, nz);
                   
                n = cellnb_3d(nx-1, ny-1, k, 0, nx, ny, nz);
                cell[n].nneighb = 4;
                cell[n].nghb[0] = cellnb_3d(nx-1, ny-1, k+1, 0, nx, ny, nz);
                cell[n].nghb[1] = cellnb_3d(nx-1, ny-1, k-1, 0, nx, ny, nz);
                cell[n].nghb[2] = cellnb_3d(nx-2, ny-1, k, 0, nx, ny, nz);
                cell[n].nghb[3] = cellnb_3d(nx-1, ny-2, k, 0, nx, ny, nz);
            }
            
            /* corners */
            n = cellnb_3d(0, 0, 0, 0, nx, ny, nz);
            cell[n].nneighb = 3;
            cell[n].nghb[0] = cellnb_3d(1, 0, 0, 0, nx, ny, nz);
            cell[n].nghb[1] = cellnb_3d(0, 1, 0, 0, nx, ny, nz);
            cell[n].nghb[2] = cellnb_3d(0, 0, 1, 0, nx, ny, nz);
            
            n = cellnb_3d(nx-1, 0, 0, 0, nx, ny, nz);
            cell[n].nneighb = 3;
            cell[n].nghb[0] = cellnb_3d(nx-2, 0, 0, 0, nx, ny, nz);
            cell[n].nghb[1] = cellnb_3d(nx-1, 1, 0, 0, nx, ny, nz);
            cell[n].nghb[2] = cellnb_3d(nx-1, 0, 1, 0, nx, ny, nz);
            
            n = cellnb_3d(0, ny-1, 0, 0, nx, ny, nz);
            cell[n].nneighb = 3;
            cell[n].nghb[0] = cellnb_3d(1, ny-1, 0, 0, nx, ny, nz);
            cell[n].nghb[1] = cellnb_3d(0, ny-2, 0, 0, nx, ny, nz);
            cell[n].nghb[2] = cellnb_3d(0, ny-1, 1, 0, nx, ny, nz);
            
            n = cellnb_3d(nx-1, ny-1, 0, 0, nx, ny, nz);
            cell[n].nneighb = 3;
            cell[n].nghb[0] = cellnb_3d(nx-2, ny-1, 0, 0, nx, ny, nz);
            cell[n].nghb[1] = cellnb_3d(nx-1, ny-2, 0, 0, nx, ny, nz);
            cell[n].nghb[2] = cellnb_3d(nx-1, ny-1, 1, 0, nx, ny, nz);
            
            n = cellnb_3d(0, 0, nz-1, 0, nx, ny, nz);
            cell[n].nneighb = 3;
            cell[n].nghb[0] = cellnb_3d(1, 0, nz-1, 0, nx, ny, nz);
            cell[n].nghb[1] = cellnb_3d(0, 1, nz-1, 0, nx, ny, nz);
            cell[n].nghb[2] = cellnb_3d(0, 0, nz-2, 0, nx, ny, nz);
            
            n = cellnb_3d(nx-1, 0, nz-1, 0, nx, ny, nz);
            cell[n].nneighb = 3;
            cell[n].nghb[0] = cellnb_3d(nx-2, 0, nz-1, 0, nx, ny, nz);
            cell[n].nghb[1] = cellnb_3d(nx-1, 1, nz-1, 0, nx, ny, nz);
            cell[n].nghb[2] = cellnb_3d(nx-1, 0, nz-2, 0, nx, ny, nz);
            
            n = cellnb_3d(0, ny-1, nz-1, 0, nx, ny, nz);
            cell[n].nneighb = 3;
            cell[n].nghb[0] = cellnb_3d(1, ny-1, nz-1, 0, nx, ny, nz);
            cell[n].nghb[1] = cellnb_3d(0, ny-2, nz-1, 0, nx, ny, nz);
            cell[n].nghb[2] = cellnb_3d(0, ny-1, nz-2, 0, nx, ny, nz);
            
            n = cellnb_3d(nx-1, ny-1, nz-1, 0, nx, ny, nz);
            cell[n].nneighb = 3;
            cell[n].nghb[0] = cellnb_3d(nx-2, ny-1, nz-1, 0, nx, ny, nz);
            cell[n].nghb[1] = cellnb_3d(nx-1, ny-2, nz-1, 0, nx, ny, nz);
            cell[n].nghb[2] = cellnb_3d(nx-1, ny-1, nz-2, 0, nx, ny, nz);
                        
            return(nx*ny*nz);
        }
    }
    printf("Done\n");
}


void init_cell_probabilities(t_perco *cell, int ncells)
/* initialize the probabilities of cells being open */
{
    int i;
    
    printf("Initializing cell probabilities ...");
    #pragma omp parallel for private(i)
    for (i=0; i<ncells; i++)
    {
        cell[i].proba = (double)rand()/RAND_MAX;
    }
    printf("Done\n");
}

void init_cell_state(t_perco *cell, double p, int ncells, int first)
/* initialize the probabilities of cells being open */
{
    int i, delta_i;
    
    delta_i = ncells/100;
    
    printf("Initializing cell state ...");
    #pragma omp parallel for private(i)
    for (i=0; i<ncells; i++)
    {
        if ((VERBOSE)&&(i%delta_i == 0)) printf("%i ", i/delta_i);
        if (cell[i].proba < p) 
        {
            cell[i].open = 1;
            cell[i].active = 1;
        }
        else 
        {
            cell[i].open = 0;
            cell[i].active = 0;
        }
        cell[i].flooded = 0;
        if (first) cell[i].previous_cluster = 0;
        else cell[i].previous_cluster = cell[i].cluster; 
        cell[i].cluster = 0;
        cell[i].tested = 0;
    }
    printf("Done\n");
}

int init_flooded_cells(t_perco *cell, int ncells, int nx, int ny, int nz, int bottom, t_perco* *pstack)
/* make the left row of cells flooded, returns number of flooded cells */
{
    int i, j, k, n, ishift, c;
    double pdisc_prop;
    
//     printf("Initializing flooded cells ...");
    switch (graphical_rep(LATTICE)) {
        case (PLOT_SQUARES): 
        {
            #pragma omp parallel for private(j)
            for (j=0; j<ny; j++)
            {
                cell[j].open = 1;
                cell[j].flooded = 1;
                cell[j].cluster = 1;
                cell[j].previous_cluster = 1;
                cell[j].active = 1;
                pstack[j] = &cell[j];
            }
            return(ny);
        }
        case (PLOT_SQUARE_BONDS): 
        {
            ishift = nx*(ny+1);
            
            #pragma omp parallel for private(j)
            for (j=0; j<ny; j++)
            {
                cell[ishift + j].open = 1;
                cell[ishift + j].flooded = 1;
                cell[ishift + j].cluster = 1;
                cell[ishift + j].previous_cluster = 1;
                cell[ishift + j].active = 1;
                pstack[j] = &cell[ishift + j];
            }
            return(ny);
        }
        case (PLOT_HEX): 
        {
            #pragma omp parallel for private(j)
            for (j=0; j<ny; j++)
            {
                cell[j*nx].open = 1;
                cell[j*nx].flooded = 1;
                cell[j*nx].cluster = 1;
                cell[j*nx].previous_cluster = 1;
                cell[j*nx].active = 1;
                pstack[j] = &cell[j*nx];
            }
            return(ny);
        }
        case (PLOT_HEX_BONDS): 
        {
            #pragma omp parallel for private(j)
            for (j=0; j<ny; j+=2)
            {
                cell[j*nx].open = 1;
                cell[j*nx].flooded = 1;
                cell[j*nx].cluster = 1;
                cell[j*nx].previous_cluster = 1;
                cell[j*nx].active = 1;
                pstack[j] = &cell[j*nx];
            }
            return(ny);
        }
        case (PLOT_TRIANGLE): 
        {
            #pragma omp parallel for private(j)
            for (j=0; j<ny; j++)
            {
                cell[2*j*nx].open = 1;
                cell[2*j*nx].flooded = 1;
                cell[2*j*nx].cluster = 1;
                cell[2*j*nx].previous_cluster = 1;
                cell[2*j*nx].active = 1;
                pstack[j] = &cell[2*j*nx];
            }
            return(ny);
        }
        case (PLOT_POISSON_DISC):
        {
            n = 0;
            pdisc_prop = 0.5/sqrt((double)ncells);
            #pragma omp parallel for private(j)
            for (j=0; j<ncells; j++) if (cell[j].x - XMIN < (XMAX - XMIN)*pdisc_prop)
            {
//                 printf("Flooding cell %i\n", j);
                cell[j].open = 1;
                cell[j].flooded = 1;
                cell[j].cluster = 1;
                cell[j].previous_cluster = 1;
                cell[j].active = 1;
                pstack[n] = &cell[j];
                n++;
            }
            printf("Flooded %i cells\n", n);
            return(n);
        }
        case (PLOT_CUBES): 
        {
            if (bottom)
            {
                #pragma omp parallel for private(i, j)
                for (i=0; i<nx; i++)
                    for (j=0; j<ny; j++)
                    {
                        c = cellnb_3d(i, j, 0, 0, nx, ny, nz);
                        cell[c].open = 1;
                        cell[c].flooded = 1;
                        cell[c].cluster = 1;
                        cell[c].previous_cluster = 1;
                        cell[c].active = 1;
                        pstack[j*nx+i] = &cell[c];
                    }
                return(nx*ny);                
            }
            else
            {
                #pragma omp parallel for private(j,k)
                for (j=0; j<ny; j++)
                    for (k=0; k<nz; k++)
                    {
                        c = cellnb_3d(0, j, k, 0, nx, ny, nz);
                        cell[c].open = 1;
                        cell[c].flooded = 1;
                        cell[c].cluster = 1;
                        cell[c].previous_cluster = 1;
                        cell[c].active = 1;
                        pstack[j*nz+k] = &cell[c];
                    }
                return(ny*nz);
            }
        }
        
    }
//     printf("Done\n");
}


int count_open_cells(t_perco *cell, int ncells)
/* count the number of open cells */
{
    int n = 0, i, delta_i;
    
    delta_i = ncells/100;
    for (i=0; i<ncells; i++) 
    {
        if (cell[i].open) n++;
        if ((VERBOSE)&&(i%delta_i == 0)) printf("%i ", i/delta_i);
    }
    return(n);
}

int count_flooded_cells(t_perco *cell, int ncells)
/* count the number of flooded cells */
{
    int n = 0, i;
    
    for (i=0; i<ncells; i++) if ((cell[i].open)&&(cell[i].flooded)) n++;
    return(n);
}

int count_active_cells(t_perco *cell, int ncells)
/* count the number of active cells */
{
    int n = 0, i;
    
    for (i=0; i<ncells; i++) 
        if ((cell[i].active)&&(cell[i].flooded)) n++;
        
    return(n);
}


int find_open_cells(t_perco *cell, int ncells, int position, t_perco* *pstack, int *stacksize, int cluster)
/* look for open neighbours of cell[position], returns difference in open cells */
{
    int k, nopen = 0, n;
    
//     printf("Position %i, open %i\n", position, pstack[position]->open);
    
    if (!pstack[position]->open) return(-1);
    
    pstack[position]->cluster = cluster;
    
    for (k=0; k<pstack[position]->nneighb; k++)
    {
        n = pstack[position]->nghb[k];
        if ((cell[n].open)&&(cell[n].cluster != cluster))
        {
            cell[n].flooded = 1;
            cell[n].tested = 1;
            cell[n].cluster = cluster;
            cell[n].active = 1;
            nopen++;
            if (*stacksize < ncells)
            {
                (*stacksize)++;
                pstack[*stacksize-1] = &cell[n];
            }
        }
    }

    if (nopen == 0) 
    {
        pstack[position]->active = 0;
        return(-1);
    }
    else return(nopen);
}


int find_percolation_cluster(t_perco *cell, int ncells, t_perco* *pstack, int nstack)
/* new version of cluster search algorithm, using stacks; returns number of flooded cells */
{
    int position = 0, i, stacksize = nstack, nactive = nstack; 
    
    while (nactive > 0) 
    {
        /* find an active cell in stack */
        while (!pstack[position]->active) position++;
        if (position == stacksize) position = 0;
        
        nactive += find_open_cells(cell, ncells, position, pstack, &stacksize, 1);
    }
    printf("Stack size %i\n", stacksize);
    
    return(stacksize);
}


int find_all_clusters(t_perco *cell, int ncells, int reset, int *max_cluster_label)
/* find all clusters of open cells; returns number of clusters */
{
    int nclusters = 0, i, j, nactive, stacksize, position, newcluster = 1, maxlabel = 0, delta_i;
    t_perco **cstack;
    static int cluster;
    
    cstack = (t_perco* *)malloc(ncells*sizeof(struct t_perco *));
    
    if (reset) cluster = CLUSTER_SHIFT;
    else cluster = (*max_cluster_label) + CLUSTER_SHIFT;   /* give higher labels than before to new clusters */

    delta_i = ncells/1000;
    
    for (i=0; i<ncells; i++) if ((cell[i].open)&&(cell[i].cluster != 1)&&(!cell[i].tested))
    {
        position = 0;
        stacksize = 1;
        nactive = 1;
        cstack[0] = &cell[i];
        
        /* choice of new color */ 
        if (cell[i].previous_cluster >= CLUSTER_SHIFT) newcluster = cell[i].previous_cluster;
        else 
        {
            newcluster = cluster;
            cluster++;
        }
        cell[i].cluster = newcluster;
        cell[i].active = 1;
        cell[i].tested = 1;
        
        while (nactive > 0) 
        {
            /* find an active cell in stack */
            while (!cstack[position]->active) position++;
            if (position == stacksize) position = 0;
        
            nactive += find_open_cells(cell, ncells, position, cstack, &stacksize, newcluster);
            
//             if ((VERBOSE)&&(nactive%100 == 99)) printf ("%i active cells\n", nactive + 1);
        }
        
        if ((VERBOSE)&&(i%delta_i == 0)) printf("Found cluster %i with %i active cells\n", cluster, stacksize);

        nclusters++;
    }
    
    free(cstack);
    
    /* determine maximal cluster label */
//     #pragma omp parallel for private(i)
    for (i=0; i<ncells; i++) if ((cell[i].open)&&(cell[i].cluster > maxlabel)) maxlabel = cell[i].cluster;
    
    printf("Cluster = %i, maxlabel = %i\n", cluster, maxlabel);
    
    *max_cluster_label = maxlabel;
    return(nclusters);
}

void print_cluster_sizes(t_perco *cell, int ncells, int *cluster_sizes)
/* for debugging purposes */
{
    int j;
    
    for (j=0; j<ncells; j++) if (cell[j].open == 1) 
            printf("(cell %i: cluster %i: size %i)\n", j, cell[j].cluster, cluster_sizes[cell[j].cluster]);
    sleep(1);
    printf("\n\n");
}

int find_cluster_sizes(t_perco *cell, int ncells, int *cluster_sizes, int *maxclusterlabel)
/* determine the sizes of all clusters; returns the size of largest cluster */
{
    int i, max = 0, maxlabel = 0;
    
//     if (maxclusterlabel > ncells) maxclusterlabel = ncells;
    
    /* determine label of larges cluster */
    for (i=0; i<ncells; i++) if (cell[i].cluster > maxlabel) maxlabel = cell[i].cluster;
    *maxclusterlabel = maxlabel + 1;

    #pragma omp parallel for private(i)
    for (i=0; i<*maxclusterlabel; i++) cluster_sizes[i] = 0;
    
    for (i=0; i<ncells; i++) if (cell[i].open) cluster_sizes[cell[i].cluster]++;
    
    if (DEBUG)
    {
        printf("Computed %i cluster sizes\n", *maxclusterlabel);
        print_cluster_sizes(cell, ncells, cluster_sizes);
    }
    
    for (i=0; i<*maxclusterlabel; i++) if (cluster_sizes[i] > max) max = cluster_sizes[i];
    
//     for (i=0; i<ncells; i++) if (cell[i].open) printf("%i ", cluster_sizes[cell[i].cluster]);
    
    printf("Max cluster label %i\n", *maxclusterlabel);
    printf("Largest cluster has %i cells\n", max); 
    return(max);
}
