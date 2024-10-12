void init_sphere_2D()		/* initialisation of window */
{
    glLineWidth(3);

    glClearColor(0.0, 0.0, 0.0, 1.0);
    glClear(GL_COLOR_BUFFER_BIT);

    glOrtho(0.0, NX, DPOLE, NY-DPOLE, -1.0, 1.0);
}

void draw_vertex_sphere(double xyz[3])
{
    double xy_screen[2];
    
    xyz_to_xy(xyz[0], xyz[1], xyz[2], xy_screen);
    glVertex2d(xy_screen[0], xy_screen[1]);
}

void init_wave_sphere(t_wave_sphere wsphere[NX*NY])
/* initialize sphere data */
{
    int i, j;
    double dphi, dtheta, theta0, xy[2], phishift;
    
    printf("Initializing wsphere\n");
    
    dphi = DPI/(double)NX;
//     dtheta = PI/(double)NY;
    dtheta = PI/(double)(NY-2*(DPOLE));
    theta0 = (double)(DPOLE)*dtheta;
    phishift = PHISHIFT*(XMAX-XMIN)/360.0;
    
    #pragma omp parallel for private(i,j)
    for (i=0; i<NX; i++)
    {
        for (j=DPOLE; j<NY-DPOLE; j++)
            wsphere[i*NY+j].theta = (double)j*dtheta - theta0;
        for (j=0; j<DPOLE; j++) wsphere[i*NY+j].theta = 0.0;
        for (j=NY-DPOLE; j<NY; j++) wsphere[i*NY+j].theta = PI;
        
        for (j=0; j<NY; j++)
        {
            wsphere[i*NY+j].phi = (double)i*dphi;
//             wsphere[i*NY+j].theta = (double)j*dtheta;
            
            wsphere[i*NY+j].cphi = cos(wsphere[i*NY+j].phi);
            wsphere[i*NY+j].sphi = sin(wsphere[i*NY+j].phi);
            
            wsphere[i*NY+j].ctheta = cos(wsphere[i*NY+j].theta);
            wsphere[i*NY+j].stheta = sin(wsphere[i*NY+j].theta);
            
            wsphere[i*NY+j].x = wsphere[i*NY+j].cphi*wsphere[i*NY+j].stheta;
            wsphere[i*NY+j].y = wsphere[i*NY+j].sphi*wsphere[i*NY+j].stheta;
            wsphere[i*NY+j].z = -wsphere[i*NY+j].ctheta;
            
            wsphere[i*NY+j].radius = 1.0;
            wsphere[i*NY+j].radius_dem = 1.0;
            
            ij_to_xy(NX-1-i,j,xy);
//             xy[0] = XMIN + ((double)(NX-i-1))*(XMAX-XMIN)/((double)NX);
//             xy[1] = YMIN + ((double)(j-DPOLE))*(YMAX-YMIN)/((double)(NY-2*DPOLE));
            
            xy[0] += phishift;
            if (xy[0] > XMAX) xy[0] += XMIN - XMAX;
            
            xy[1] *= (double)NY/(double)(NY-2*DPOLE);
            
            wsphere[i*NY+j].x2d = xy[0]; 
            wsphere[i*NY+j].y2d = xy[1]; 
            
            wsphere[i*NY+j].cos_angle_sphere = wsphere[i*NY+j].x*light[0] + wsphere[i*NY+j].y*light[1] + wsphere[i*NY+j].z*light[2];
        }
        
        /* cotangent, taking care of not dividing by zero */
        /* TODO clean up cottheta range ? */
        for (j=DPOLE; j<NY-DPOLE; j++) wsphere[i*NY+j].cottheta = wsphere[i*NY+j].ctheta/wsphere[i*NY+j].stheta;
        for (j=0; j<DPOLE; j++) wsphere[i*NY+j].cottheta = wsphere[i*NY+DPOLE].cottheta;
        for (j=NY-DPOLE; j<NY; j++) wsphere[i*NY+j].cottheta = wsphere[i*NY+DPOLE-1].cottheta;        
        
        /* old version */
//         for (j=1; j<NY-1; j++) wsphere[i*NY+j].cottheta = wsphere[i*NY+j].ctheta/wsphere[i*NY+j].stheta;
//         wsphere[i*NY].cottheta = wsphere[i*NY+1].cottheta;
//         wsphere[i*NY + NY-1].cottheta = wsphere[i*NY+1].cottheta;
    }
}


void read_negative_dem_values(double *height_values, t_wave_sphere wsphere[NX*NY])
/* init bathymetric data */
{
    int i, j, k, ii, jj, nx, ny, maxrgb, nmaxpixels = 6480000, hmin, hmax, ishift, nshift, sshift, rgbval, scan, rgbtot;
    int *rgb_values, *int_height_values;
    int hcont = 50, rgbdiff;
    double cratio, rx, ry, height;
    double *height_values_tmp, *height_values_tmp2;
    FILE *image_file;
    
    printf("Reading bathymetric data\n");
    
    rgb_values = (int *)malloc(3*nmaxpixels*sizeof(int));
    int_height_values = (int *)malloc(3*nmaxpixels*sizeof(int));
    height_values_tmp = (double *)malloc(NX*NY*sizeof(double));
    height_values_tmp2 = (double *)malloc(NX*NY*sizeof(double));
    
//     image_file = fopen("bathymetry_gebco_3600x1800_color.ppm", "r");
    image_file = fopen("bathymetry_gebco_2560_1280_mod2_color.ppm", "r");
    scan = fscanf(image_file,"%i %i\n", &nx, &ny);
    scan = fscanf(image_file,"%i\n", &maxrgb);
    
    for (i=0; i<NX*NY; i++) 
    {
        height_values_tmp[i] = 0.0;
        height_values_tmp2[i] = 0.0;
    }
    
    hmin = maxrgb;
    hmax = 0;
    
    if (nx*ny > nmaxpixels)
    {
        printf("bathymetric data file too large, increase nmaxpixels in read_negative_dem_values()\n");
        exit(0);
    }
    
    /* shift due to min/max latitudes of image */
    sshift = 0 + DPOLE;
    nshift = 0 + DPOLE; 
    
    /* choice of zero meridian */
    ishift = (int)(nx*ZERO_MERIDIAN/360.0);
    
    /* read rgb values */
    for (j=0; j<ny; j++)
        for (i=0; i<nx; i++)
        {
            rgbtot = 0;
            for (k=0; k<3; k++)
            {
                scan = fscanf(image_file,"%i\n", &rgbval);
                rgb_values[3*(j*nx+i)+k] = rgbval;
                rgbtot += rgbval;
            }
            if ((rgbtot < hmin)&&(rgbtot > hcont)) hmin = rgbtot;
            if (rgbtot > hmax) hmax = rgbtot;
            int_height_values[3*(j*nx+i)] = rgbtot;
        }
    
    printf("hmin = %i, hmax = %i\n", hmin, hmax);
    
    /* remove remaining black continents */
    for (i=0; i<3*nx*ny; i++) if (int_height_values[i] < hcont) int_height_values[i] = hmax;
    
    cratio = 1.0/(double)(hmax-hmin);
    rx = (double)nx/(double)NX;
    ry = (double)ny/(double)(NY - sshift - nshift);
    
    /* option: set continents color to white */
//     for (i=0; i<NX; i++)
//         for (j=0; j<NY; j++) 
//         {
//             wsphere[i*NY+j].r = 0.9;
//             wsphere[i*NY+j].g = 0.9;
//             wsphere[i*NY+j].b = 0.9;
//         }
    
    /* build underwater height table */
    for (i=0; i<NX; i++)
        for (j=0; j<NY; j++) 
        {
            ii = (int)(rx*(double)(NX-1 - i)) + nx/2 + ishift;
            if (ii > nx-1) ii -= nx;
            jj = (int)(ry*(double)(NY-nshift - j));
            if (jj > ny-1) jj = ny-1;
            if (jj < 0) jj = 0;
                        
            if (wsphere[i*NY+j].indomain)
            {
                /* set height to zero if color is black (due to black patches in bathymetric map) */
                if (int_height_values[3*(jj*nx+ii)] < hcont) height = 0.0;
                else height = -1.0 + (double)(int_height_values[3*(jj*nx+ii)])*cratio;
                if (height > 0.0) height = 0.0;
                height_values_tmp[i*NY+j] = height;
                wsphere[i*NY+j].altitude = height;

                if (int_height_values[3*(jj*nx+ii)] > hcont)
                {
                    wsphere[i*NY+j].r = 0.9*(double)rgb_values[3*(jj*nx+ii)]*cratio;
                    wsphere[i*NY+j].g = 0.9*(double)rgb_values[3*(jj*nx+ii)+1]*cratio;
                    wsphere[i*NY+j].b = 0.9*(double)rgb_values[3*(jj*nx+ii)+2]*cratio;
                }
                else 
                {
                    wsphere[i*NY+j].r = 0.29;
                    wsphere[i*NY+j].g = 0.29;
                    wsphere[i*NY+j].b = 0.29;
                }
            }
            else 
            {
                height_values_tmp[i*NY+j] = 0.0;
                height_values_tmp2[i*NY+j] = 0.0;
            }
        }
    
    /* smoothen values at low depth */
    if (SMOOTH_DEM) for (k=1; k<DEM_SMOOTH_STEPS; k++)
    {
        printf("Smoothing step %i\n", k);
        for (i=1; i<NX-1; i++)
            for (j=1; j<NY-1; j++) 
                if ((wsphere[i*NY+j].indomain)&&(height_values[i*NY+j] >= -0.25))
                {
                    height_values_tmp2[i*NY+j] = height_values_tmp[i*NY+j] + 0.1*(height_values_tmp[(i+1)*NY+j] + height_values_tmp[(i-1)*NY+j] + height_values_tmp[i*NY+j+1] + height_values_tmp[i*NY+j-1] - 4.0*height_values_tmp[i*NY+j]);
            
                    height_values_tmp[i*NY+j] = height_values_tmp2[i*NY+j] + 0.1*(height_values_tmp2[(i+1)*NY+j] + height_values_tmp2[(i-1)*NY+j] + height_values_tmp2[i*NY+j+1] + height_values_tmp2[i*NY+j-1] - 4.0*height_values_tmp2[i*NY+j]);
                }
    }
            
    if (SMOOTH_DEM) for (i=1; i<NX-1; i++)
        for (j=1; j<NY-1; j++) 
            if ((wsphere[i*NY+j].indomain)&&(height_values[i*NY+j] >= -0.25))
            {
                wsphere[i*NY+j].altitude = height_values_tmp[i*NY+j];
            }
    
    for (i=0; i<NX; i++)
        for (j=0; j<NY; j++) 
            height_values[i*NY+j] = height_values_tmp[i*NY+j];
    
    free(rgb_values);
    free(int_height_values);
    free(height_values_tmp);
    free(height_values_tmp2);
}

void init_dem(t_wave_sphere wsphere[NX*NY], int dem_number)
/* init heights from digital elevation map */
{
    int i, j, ii, jj, k, nx, ny, maxrgb, nmaxpixels = 4915200, scan, rgbval, diff, sshift, nshift, hmin, hmax, ishift, hsum;
    int *rgb_values;
    double cratio, rx, ry, cy, dx, dy, pscal, norm, vscale1, vscale2, gradx, grady, deltar, deltai[3], deltaj[3], dphi, dtheta, n[3], hsea, hmean, vshift;
    double *height_values, *height_values_tmp;
    FILE *image_file;
    
    printf("Reading digital elevation model\n");
    
    switch (dem_number) {
        case (DEM_EARTH): 
        {
            nmaxpixels = 4915200;
            image_file = fopen("digital_elevation_model_large.ppm", "r");
            hsea = 12.0;  /* sea level #0c0c0c */
            break;
        }
        case (DEM_MARS):
        {
            nmaxpixels = 8388608;
            image_file = fopen("marscyl2.ppm", "r");
            hsea = 6.0 + 255.0*PLANET_SEALEVEL/DEM_MAXHEIGHT;
            break;
        }
        case (DEM_MOON):
        {
            nmaxpixels = 2097152;
            image_file = fopen("Moon_LRO_LOLA_global_LDEM_1024.ppm", "r");
            hsea = 255.0*PLANET_SEALEVEL/(DEM_MAXHEIGHT - DEM_MAXDEPTH);
            break;
        }
        case (DEM_VENUS):
        {
            nmaxpixels = 4096*2048;
            image_file = fopen("Venus_Magellan_Topography_Global_4641m_v02_scaled2.ppm", "r");
            hsea = 255.0*PLANET_SEALEVEL/DEM_MAXHEIGHT;
            break;
        }
        case (DEM_MERCURY):
        {
            nmaxpixels = 1280*380;
            image_file = fopen("Mercury_Messenger_DEM_Global_665m_1024_1_cropped.ppm", "r");
            hsea = 255.0*PLANET_SEALEVEL/DEM_MAXHEIGHT;
            break;
        }
    }
    
    rgb_values = (int *)malloc(3*nmaxpixels*sizeof(int));
    height_values = (double *)malloc(NX*NY*sizeof(double));
    height_values_tmp = (double *)malloc(NX*NY*sizeof(double));
    
//     image_file = fopen("digital_elevation_model_large.ppm", "r");
    scan = fscanf(image_file,"%i %i\n", &nx, &ny);
    scan = fscanf(image_file,"%i\n", &maxrgb);
    
    printf("nx = %i, ny = %i, maxrgb = %i\n", nx, ny, maxrgb);
    sleep(1);
    
    hmin = maxrgb;
    hmax = 0;
    
    if (nx*ny > nmaxpixels)
    {
        printf("DEM too large, increase nmaxpixels in init_dem()\n");
        exit(0);
    }
    
    /* shift due to min/max latitudes of image */
    sshift = 0 + DPOLE;
    nshift = 0 + DPOLE; 
    
    /* choice of zero meridian */
    ishift = (int)(nx*ZERO_MERIDIAN/360.0);
    
    printf("Reading RGB values\n");
    
    /* read rgb values */
    for (j=0; j<ny; j++)
        for (i=0; i<nx; i++)
            for (k=0; k<3; k++)
                {
                    scan = fscanf(image_file,"%i\n", &rgbval);
                    rgb_values[3*(j*nx+i)+k] = rgbval;
                    if (rgbval < hmin) 
                    {
                        if (B_DOMAIN == D_SPHERE_VENUS)
                        {
                            if (rgbval > 0) hmin = rgbval;
                        }
                        else hmin = rgbval;
                    }
                    if (rgbval > hmax) hmax = rgbval;
                }
    
    printf("hmin = %i, hmax = %i, hsea = %i\n", hmin, hmax, (int)hsea);
    
    if (B_DOMAIN == D_SPHERE_VENUS)
    {
        hsum = 0;
        for (j=0; j<ny; j++)
            for (i=0; i<nx; i++)
                hsum += rgb_values[3*(j*nx+i)];
        hmean = (double)hsum/(double)(nx*ny);
        printf("hmean = %.2f\n", hmean);
    }
    
    cratio = 1.0/(double)(hmax-hmin);
    rx = (double)nx/(double)NX;
    ry = (double)ny/(double)(NY - sshift - nshift);
    
    /* build height table */
    vshift = PLANET_SEALEVEL/(DEM_MAXHEIGHT - DEM_MAXDEPTH);
    for (i=0; i<NX; i++)
        for (j=0; j<NY; j++)
        {
            ii = (int)(rx*(double)(NX-1 - i)) + nx/2 + ishift;
            if (ii > nx-1) ii -= nx;
            if (ii < 0) ii = 0;
            jj = (int)(ry*(double)(NY-nshift - j));
            if (jj > ny-1) jj = ny-1;
            if (jj < 0) jj = 0;
            height_values[i*NY+j] = ((double)rgb_values[3*(jj*nx+ii)]-hsea)*cratio;
            wsphere[i*NY+j].altitude = ((double)rgb_values[3*(jj*nx+ii)]-hsea)*cratio;
            
            /* take care of black areas (missing data) on venus */
            if ((B_DOMAIN == D_SPHERE_VENUS)&&(rgb_values[3*(jj*nx+ii)] == 0))
            {
                height_values[i*NY+j] = VENUS_NODATA_FACTOR*hmean*cratio;
                wsphere[i*NY+j].altitude = VENUS_NODATA_FACTOR*hmean*cratio;
            }
            
            if (OTHER_PLANET)
                wsphere[i*NY+j].indomain = (wsphere[i*NY+j].altitude < vshift);
            
//             if (wsphere[i*NY+j].indomain) printf("rgb = %i, altitude = %.3lg\n", rgb_values[3*(jj*nx+ii)], height_values[i*NY+j]);
        }
        
    /* smooth values in case of high resolution */
    if ((NX > nx)||(NY > ny))
    {
        for (i=1; i<NX-1; i++)
            for (j=1; j<NY-1; j++)
            {
                height_values[i*NY+j] *= 0.2; 
                height_values[i*NY+j] += 0.2*height_values[(i+1)*NY+j];
                height_values[i*NY+j] += 0.2*height_values[(i-1)*NY+j];
                height_values[i*NY+j] += 0.2*height_values[i*NY+j-1];
                height_values[i*NY+j] += 0.2*height_values[i*NY+j+1];
                
                wsphere[i*NY+j].altitude *= 0.2; 
                wsphere[i*NY+j].altitude += 0.2*wsphere[(i+1)*NY+j].altitude;
                wsphere[i*NY+j].altitude += 0.2*wsphere[(i-1)*NY+j].altitude;
                wsphere[i*NY+j].altitude += 0.2*wsphere[i*NY+j-1].altitude;
                wsphere[i*NY+j].altitude += 0.2*wsphere[i*NY+j+1].altitude;
            }
            
        /* i = 0 */
        for (j=1; j<NY-1; j++)
        {
            height_values[j] *= 0.2; 
            height_values[j] += 0.2*height_values[NY+j];
            height_values[j] += 0.2*height_values[(NX-1)*NY+j];
            height_values[j] += 0.2*height_values[j-1];
            height_values[j] += 0.2*height_values[j+1];
                
            wsphere[j].altitude *= 0.2; 
            wsphere[j].altitude += 0.2*wsphere[NY+j].altitude;
            wsphere[j].altitude += 0.2*wsphere[(NX-1)*NY+j].altitude;
            wsphere[j].altitude += 0.2*wsphere[j-1].altitude;
            wsphere[j].altitude += 0.2*wsphere[j+1].altitude;
        }
        
        /* i = NY-1 */
        for (j=1; j<NY-1; j++)
        {
            height_values[(NY-1)*NY+j] *= 0.2; 
            height_values[(NY-1)*NY+j] += 0.2*height_values[j];
            height_values[(NY-1)*NY+j] += 0.2*height_values[(NY-2)*NY+j];
            height_values[(NY-1)*NY+j] += 0.2*height_values[(NY-1)*NY+j-1];
            height_values[(NY-1)*NY+j] += 0.2*height_values[(NY-1)*NY+j+1];
                
            wsphere[(NY-1)*NY+j].altitude *= 0.2; 
            wsphere[(NY-1)*NY+j].altitude += 0.2*wsphere[j].altitude;
            wsphere[(NY-1)*NY+j].altitude += 0.2*wsphere[(NY-2)*NY+j].altitude;
            wsphere[(NY-1)*NY+j].altitude += 0.2*wsphere[(NY-1)*NY+j-1].altitude;
            wsphere[(NY-1)*NY+j].altitude += 0.2*wsphere[(NY-1)*NY+j+1].altitude;
        }
    }
        
    printf("Closing rgb_values\n");
    fclose(image_file);
    free(rgb_values);
    
    /* smoothen values at low altitude */
    if (SMOOTH_DEM) for (k=1; k<DEM_SMOOTH_STEPS; k++)
    {
        printf("Smoothing step %i\n", k);
        for (i=1; i<NX-1; i++)
            for (j=1; j<NY-1; j++) 
                if ((!wsphere[i*NY+j].indomain)&&(height_values[i*NY+j] <= DEM_SMOOTH_HEIGHT))
                {
                    height_values_tmp[i*NY+j] = height_values[i*NY+j] + 0.1*(height_values[(i+1)*NY+j] + height_values[(i-1)*NY+j] + height_values[i*NY+j+1] + height_values[i*NY+j-1] - 4.0*height_values[i*NY+j]);
            
                    height_values[i*NY+j] = height_values_tmp[i*NY+j] + 0.1*(height_values_tmp[(i+1)*NY+j] + height_values_tmp[(i-1)*NY+j] + height_values_tmp[i*NY+j+1] + height_values_tmp[i*NY+j-1] - 4.0*height_values_tmp[i*NY+j]);
                }
        /* i = 0 */
        for (j=1; j<NY-1; j++) 
            if ((!wsphere[j].indomain)&&(height_values[j] <= DEM_SMOOTH_HEIGHT))
            {
                height_values_tmp[j] = height_values[j] + 0.1*(height_values[NY+j] + height_values[(NX-1)*NY+j] + height_values[j+1] + height_values[j-1] - 4.0*height_values[j]);
            
                height_values[j] = height_values_tmp[j] + 0.1*(height_values_tmp[NY+j] + height_values_tmp[(NX-1)*NY+j] + height_values_tmp[j+1] + height_values_tmp[j-1] - 4.0*height_values_tmp[j]);
            }
        /* i = NY-1 */
        for (j=1; j<NY-1; j++) 
            if ((!wsphere[(NX-1)*NY+j].indomain)&&(height_values[(NX-1)*NY+j] <= DEM_SMOOTH_HEIGHT))
            {
                height_values_tmp[(NX-1)*NY+j] = height_values[(NX-1)*NY+j] + 0.1*(height_values[j] + height_values[(NX-2)*NY+j] + height_values[(NX-1)*NY+j+1] + height_values[(NX-1)*NY+j-1] - 4.0*height_values[(NX-1)*NY+j]);
            
                height_values[(NX-1)*NY+j] = height_values_tmp[(NX-1)*NY+j] + 0.1*(height_values_tmp[j] + height_values_tmp[(NX-2)*NY+j] + height_values_tmp[(NX-1)*NY+j+1] + height_values_tmp[(NX-1)*NY+j-1] - 4.0*height_values_tmp[(NX-1)*NY+j]);
            }
    }
            
    if (SMOOTH_DEM) for (i=0; i<NX; i++)
        for (j=1; j<NY-1; j++) 
            if ((!wsphere[i*NY+j].indomain)&&(wsphere[i*NY+j].altitude <= DEM_SMOOTH_HEIGHT))
            {
                wsphere[i*NY+j].altitude = height_values[i*NY+j];
            }
        
//     for (j=0; j<NY; j+= 100)
//         for (i=0; i<NX; i++) 
//             if ((wsphere[i*NY+j].altitude > 0.0)&&(wsphere[i*NY+j].altitude <= 0.15))
//                 printf("Smoothed altitude at (%i, %i): %.3lg\n", i, j, wsphere[i*NY+j].altitude);
        
    if ((ADD_NEGATIVE_DEM)&&(dem_number == DEM_EARTH)) 
        read_negative_dem_values(height_values, wsphere);
            
    /* set radius */
    for (i=0; i<NX; i++)
        for (j=0; j<NY; j++)
//             if (!wsphere[i*NY+j].indomain) wsphere[i*NY+j].radius = 1.0 + RSCALE_DEM*wsphere[i*NY+j].altitude;
//             else wsphere[i*NY+j].radius = 1.0;
                wsphere[i*NY+j].radius_dem = 1.0 + RSCALE_DEM*(wsphere[i*NY+j].altitude - vshift);
    
    /* set domain of evolution for simulations with flooding */
    for (i=0; i<NX; i++)
        for (j=0; j<NY; j++)
            wsphere[i*NY+j].evolve_wave = (wsphere[i*NY+j].altitude < vshift + 0.1); 
    
    /* compute light angle */  
    dx = 2.0*(XMAX - XMIN)/(double)NX;
    dy = 2.0*(YMAX - YMIN)/(double)NY;
    vscale1 = 0.1*SHADE_SCALE_2D;
    vscale2 = vscale1*vscale1;
    
    if (SHADE_2D)
    {
        for (i=1; i<NX-1; i++)
            for (j=1; j<NY-1; j++)
            {
                gradx = (wsphere[(i+1)*NY+j].radius_dem - wsphere[(i-1)*NY+j].radius_dem)/dx;
                grady = (wsphere[i*NY+j+1].radius_dem - wsphere[i*NY+j-1].radius_dem)/dy;
                
                norm = sqrt(vscale2 + gradx*gradx + grady*grady);
                pscal = -gradx*light[0] - grady*light[1] + vscale1;
                
                wsphere[i*NY+j].cos_angle = pscal/norm;
            }
        /* i = 0 */
        for (j=1; j<NY-1; j++)
        {
            gradx = (wsphere[NY+j].radius_dem - wsphere[(NX-1)*NY+j].radius_dem)/dx;
            grady = (wsphere[j+1].radius_dem - wsphere[j-1].radius_dem)/dy;
            
            norm = sqrt(vscale2 + gradx*gradx + grady*grady);
            pscal = -gradx*light[0] - grady*light[1] + vscale1;
            
            wsphere[j].cos_angle = pscal/norm;
        }
        /* i = N-1 */
        for (j=1; j<NY-1; j++)
        {
            gradx = (wsphere[j].radius_dem - wsphere[(NX-2)*NY+j].radius_dem)/dx;
            grady = (wsphere[(NX-1)*NY+j+1].radius_dem - wsphere[(NX-1)*NY+j-1].radius_dem)/dy;
            
            norm = sqrt(vscale2 + gradx*gradx + grady*grady);
            pscal = -gradx*light[0] - grady*light[1] + vscale1;
            
            wsphere[(NX-1)*NY+j].cos_angle = pscal/norm;
        }
    }
    else if (SHADE_3D)
    {
        dphi = DPI/(double)NX;
        dtheta = PI/(double)NY;
        
        for (i=1; i<NX-1; i++)
            for (j=1; j<NY-1; j++)
            {                
                /* computation of tangent vectors */
                deltar = (wsphere[(i+1)*NY+j].radius_dem - wsphere[i*NY+j].radius_dem)/dphi;
                    
                deltai[0] = -wsphere[i*NY+j].radius_dem*wsphere[i*NY+j].sphi;
                deltai[0] += deltar*wsphere[i*NY+j].cphi;
                    
                deltai[1] = wsphere[i*NY+j].radius_dem*wsphere[i*NY+j].cphi;
                deltai[1] += deltar*wsphere[i*NY+j].sphi;
                    
                deltai[2] = -deltar*wsphere[i*NY+j].cottheta;
                    
                deltar = (wsphere[i*NY+j+1].radius_dem - wsphere[i*NY+j].radius_dem)/dtheta;
                    
                deltaj[0] = wsphere[i*NY+j].radius_dem*wsphere[i*NY+j].cphi*wsphere[i*NY+j].ctheta;
                deltaj[0] += deltar*wsphere[i*NY+j].cphi*wsphere[i*NY+j].stheta;
                    
                deltaj[1] = wsphere[i*NY+j].radius_dem*wsphere[i*NY+j].sphi*wsphere[i*NY+j].ctheta;
                deltaj[1] += deltar*wsphere[i*NY+j].sphi*wsphere[i*NY+j].stheta;
                    
                deltaj[2] = wsphere[i*NY+j].radius_dem*wsphere[i*NY+j].stheta;
                deltaj[2] += -deltar*wsphere[i*NY+j].ctheta;
                    
                /* computation of normal vector */
                n[0] = deltai[1]*deltaj[2] - deltai[2]*deltaj[1];
                n[1] = deltai[2]*deltaj[0] - deltai[0]*deltaj[2];
                n[2] = deltai[0]*deltaj[1] - deltai[1]*deltaj[0];
                                        
                norm = sqrt(n[0]*n[0] + n[1]*n[1] + n[2]*n[2]);
                    
                pscal = n[0]*light[0] + n[1]*light[1] + n[2]*light[2];
                
                wsphere[i*NY+j].cos_angle = pscal/norm;
            }
        /* i = 0 */
        for (j=1; j<NY-1; j++)
        {                
            /* computation of tangent vectors */
            deltar = (wsphere[NY+j].radius_dem - wsphere[j].radius_dem)/dphi;
                    
            deltai[0] = -wsphere[j].radius_dem*wsphere[j].sphi;
            deltai[0] += deltar*wsphere[j].cphi;
                    
            deltai[1] = wsphere[j].radius_dem*wsphere[j].cphi;
            deltai[1] += deltar*wsphere[j].sphi;
                    
            deltai[2] = -deltar*wsphere[j].cottheta;
                    
            deltar = (wsphere[j+1].radius_dem - wsphere[j].radius_dem)/dtheta;
                    
            deltaj[0] = wsphere[j].radius_dem*wsphere[j].cphi*wsphere[j].ctheta;
            deltaj[0] += deltar*wsphere[j].cphi*wsphere[j].stheta;
                    
            deltaj[1] = wsphere[j].radius_dem*wsphere[j].sphi*wsphere[j].ctheta;
            deltaj[1] += deltar*wsphere[j].sphi*wsphere[j].stheta;
                    
            deltaj[2] = wsphere[j].radius_dem*wsphere[j].stheta;
            deltaj[2] += -deltar*wsphere[j].ctheta;
                    
            /* computation of normal vector */
            n[0] = deltai[1]*deltaj[2] - deltai[2]*deltaj[1];
            n[1] = deltai[2]*deltaj[0] - deltai[0]*deltaj[2];
            n[2] = deltai[0]*deltaj[1] - deltai[1]*deltaj[0];
                                        
            norm = sqrt(n[0]*n[0] + n[1]*n[1] + n[2]*n[2]);
                    
            pscal = n[0]*light[0] + n[1]*light[1] + n[2]*light[2];
                
            wsphere[j].cos_angle = pscal/norm;
        }
        /* i = NX-1 */
        for (j=1; j<NY-1; j++)
        {                
            /* computation of tangent vectors */
            deltar = (wsphere[j].radius_dem - wsphere[(NX-1)*NY+j].radius_dem)/dphi;
                    
            deltai[0] = -wsphere[(NX-1)*NY+j].radius_dem*wsphere[(NX-1)*NY+j].sphi;
            deltai[0] += deltar*wsphere[(NX-1)*NY+j].cphi;
                    
            deltai[1] = wsphere[(NX-1)*NY+j].radius_dem*wsphere[(NX-1)*NY+j].cphi;
            deltai[1] += deltar*wsphere[(NX-1)*NY+j].sphi;
                    
            deltai[2] = -deltar*wsphere[(NX-1)*NY+j].cottheta;
                    
            deltar = (wsphere[(NX-1)*NY+j+1].radius_dem - wsphere[(NX-1)*NY+j].radius_dem)/dtheta;
                    
            deltaj[0] = wsphere[(NX-1)*NY+j].radius_dem*wsphere[(NX-1)*NY+j].cphi*wsphere[(NX-1)*NY+j].ctheta;
            deltaj[0] += deltar*wsphere[(NX-1)*NY+j].cphi*wsphere[(NX-1)*NY+j].stheta;
                    
            deltaj[1] = wsphere[(NX-1)*NY+j].radius_dem*wsphere[(NX-1)*NY+j].sphi*wsphere[(NX-1)*NY+j].ctheta;
            deltaj[1] += deltar*wsphere[(NX-1)*NY+j].sphi*wsphere[(NX-1)*NY+j].stheta;
                    
            deltaj[2] = wsphere[(NX-1)*NY+j].radius_dem*wsphere[(NX-1)*NY+j].stheta;
            deltaj[2] += -deltar*wsphere[(NX-1)*NY+j].ctheta;
                    
            /* computation of normal vector */
            n[0] = deltai[1]*deltaj[2] - deltai[2]*deltaj[1];
            n[1] = deltai[2]*deltaj[0] - deltai[0]*deltaj[2];
            n[2] = deltai[0]*deltaj[1] - deltai[1]*deltaj[0];
                                        
            norm = sqrt(n[0]*n[0] + n[1]*n[1] + n[2]*n[2]);
                    
            pscal = n[0]*light[0] + n[1]*light[1] + n[2]*light[2];
                
            wsphere[(NX-1)*NY+j].cos_angle = pscal/norm;
        }
    }
    
    free(height_values);
    free(height_values_tmp);
}

void init_earth_map(t_wave_sphere wsphere[NX*NY])
/* init file from earth map */
{
    int i, j, ii, jj, k, nx, ny, maxrgb, nmaxpixels = 4915200, scan, rgbval, diff, sshift, nshift, ishift;
    int *rgb_values;
    double cratio, rx, ry, cy;
    FILE *image_file;
    
    printf("Reading Earth map\n");
    
    rgb_values = (int *)malloc(3*nmaxpixels*sizeof(int));
    
    image_file = fopen("Earth_Map_Blue_Marble_2002_large.ppm", "r");
    scan = fscanf(image_file,"%i %i\n", &nx, &ny);
    scan = fscanf(image_file,"%i\n", &maxrgb);

    if (nx*ny > nmaxpixels)
    {
        printf("Image too large, increase nmaxpixels in init_earth_map()\n");
        exit(0);
    }
    
    /* shift due to min/max latitudes of image */
    sshift = 0 + DPOLE;
    nshift = 0 + DPOLE; 
    
    /* choice of zero meridian */
    ishift = (int)(nx*ZERO_MERIDIAN/360.0);
    
    /* read rgb values */
    for (j=0; j<ny; j++)
        for (i=0; i<nx; i++)
            for (k=0; k<3; k++)
                {
                    scan = fscanf(image_file,"%i\n", &rgbval);
                    rgb_values[3*(j*nx+i)+k] = rgbval;
                }
                    
    cratio = 1.0/(double)maxrgb;
    rx = (double)nx/(double)NX;
    ry = (double)ny/(double)(NY - sshift - nshift);
//     cy = rx*(double)(NY - nshift);
    
    /* build wave table */
    for (i=0; i<NX; i++)
        for (j=0; j<NY; j++)
        {
            ii = (int)(rx*(double)(NX-1 - i)) + nx/2 + ishift;
            if (ii > nx-1) ii -= nx;
//             jj = (int)(-ry*(double)j + cy);
//             jj = (int)(ry*(double)(NY-nshift - j)) + sshift;
            jj = (int)(ry*(double)(NY-nshift - j));
            if (jj > ny-1) jj = ny-1;
            if (jj < 0) jj = 0;
            wsphere[i*NY+j].r = (double)rgb_values[3*(jj*nx+ii)]*cratio;
            wsphere[i*NY+j].g = (double)rgb_values[3*(jj*nx+ii)+1]*cratio;
            wsphere[i*NY+j].b = (double)rgb_values[3*(jj*nx+ii)+2]*cratio;
            
            /* decide which points are in the Sea */
            diff = iabs(rgb_values[3*(jj*nx+ii)] - 10);
            diff += iabs(rgb_values[3*(jj*nx+ii)+1] - 10);
            diff += iabs(rgb_values[3*(jj*nx+ii)+2] - 51);
            wsphere[i*NY+j].indomain = (diff < 15);
        }
    
    /* smooth colors in case of high resolution */
    if ((NX > nx)||(NY > ny))
        for (i=1; i<NX-1; i++)
            for (j=1; j<NY-1; j++)
            {
                wsphere[i*NY+j].r *= 0.2; 
                wsphere[i*NY+j].r += 0.2*wsphere[(i+1)*NY+j].r;
                wsphere[i*NY+j].r += 0.2*wsphere[(i-1)*NY+j].r;
                wsphere[i*NY+j].r += 0.2*wsphere[i*NY+j-1].r;
                wsphere[i*NY+j].r += 0.2*wsphere[i*NY+j+1].r;

                wsphere[i*NY+j].g *= 0.2; 
                wsphere[i*NY+j].g += 0.2*wsphere[(i+1)*NY+j].g;
                wsphere[i*NY+j].g += 0.2*wsphere[(i-1)*NY+j].g;
                wsphere[i*NY+j].g += 0.2*wsphere[i*NY+j-1].g;
                wsphere[i*NY+j].g += 0.2*wsphere[i*NY+j+1].g;

                wsphere[i*NY+j].b *= 0.2; 
                wsphere[i*NY+j].b += 0.2*wsphere[(i+1)*NY+j].b;
                wsphere[i*NY+j].b += 0.2*wsphere[(i-1)*NY+j].b;
                wsphere[i*NY+j].b += 0.2*wsphere[i*NY+j-1].b;
                wsphere[i*NY+j].b += 0.2*wsphere[i*NY+j+1].b;
            }
    
    free(rgb_values);
    fclose(image_file);
    
    if (ADD_DEM) init_dem(wsphere, DEM_EARTH);
}

void init_planet_map(t_wave_sphere wsphere[NX*NY], int planet)
/* init file from planetary map */
{
    int i, j, ii, jj, k, nx, ny, maxrgb, nmaxpixels, scan, rgbval, diff, sshift, nshift, ishift, dem_number;
    int *rgb_values;
    double cratio, rx, ry, cy, vshift;
    FILE *image_file;
    
    switch (planet){
        case (D_SPHERE_MARS): 
        {
            printf("Reading Mars map\n");
            nmaxpixels = 8388608;
            image_file = fopen("Mars_Viking_ClrMosaic_global_925m_scaled.ppm", "r");
            dem_number = DEM_MARS;
            break;
        }
        case (D_SPHERE_MOON):
        {
            printf("Reading Moon map\n");
            nmaxpixels = 2048*1024;
            image_file = fopen("Moon_photo_map.ppm", "r");
            dem_number = DEM_MOON;
            break;
        }
        case (D_SPHERE_VENUS):
        {
            printf("Reading Venus map\n");
            nmaxpixels = 1440*720;
            image_file = fopen("Venus_map_NASA_JPL_Magellan-Venera-Pioneer.ppm", "r");
            dem_number = DEM_VENUS;
            break;
        }
        case (D_SPHERE_MERCURY):
        {
            printf("Reading Mercury map\n");
            nmaxpixels = 2304*1152;
            image_file = fopen("Mercury_color_photo.ppm", "r");
            dem_number = DEM_MERCURY;
            break;
        }
    }
    
    scan = fscanf(image_file,"%i %i\n", &nx, &ny);
    scan = fscanf(image_file,"%i\n", &maxrgb);
    printf("nx*ny = %i\n", nx*ny);
    rgb_values = (int *)malloc(3*nmaxpixels*sizeof(int));
    
    if (nx*ny > nmaxpixels)
    {
        printf("Image too large, increase nmaxpixels in init_planet_map()\n");
        exit(0);
    }
    
    /* shift due to min/max latitudes of image */
    sshift = 0 + DPOLE;
    nshift = 0 + DPOLE; 
    
    /* choice of zero meridian */
    ishift = (int)(nx*ZERO_MERIDIAN/360.0);
    
    /* read rgb values */
    for (j=0; j<ny; j++)
        for (i=0; i<nx; i++)
            for (k=0; k<3; k++)
                {
                    scan = fscanf(image_file,"%i\n", &rgbval);
                    rgb_values[3*(j*nx+i)+k] = rgbval;
                }
                    
    cratio = 1.0/(double)maxrgb;
    rx = (double)nx/(double)NX;
    ry = (double)ny/(double)(NY - sshift - nshift);
    
    printf("cratio = %.3lg, rx = %.3lg, ry = %.3lg\n", cratio, rx, ry);
    
//     cy = rx*(double)(NY - nshift);
    
    /* build wave table */
    for (i=0; i<NX; i++)
        for (j=0; j<NY; j++)
        {
            ii = (int)(rx*(double)(NX-1 - i)) + nx/2 + ishift;
            if (ii > nx-1) ii -= nx;
//             jj = (int)(-ry*(double)j + cy);
//             jj = (int)(ry*(double)(NY-nshift - j)) + sshift;
            jj = (int)(ry*(double)(NY-nshift - j));
            if (jj > ny-1) jj = ny-1;
            if (jj < 0) jj = 0;
            wsphere[i*NY+j].r = (double)rgb_values[3*(jj*nx+ii)]*cratio;
            wsphere[i*NY+j].g = (double)rgb_values[3*(jj*nx+ii)+1]*cratio;
            wsphere[i*NY+j].b = (double)rgb_values[3*(jj*nx+ii)+2]*cratio;
            
//             printf("RGB = (%.2f, %.2f, %.2f)\n", wsphere[i*NY+j].r, wsphere[i*NY+j].g, wsphere[i*NY+j].b);
            
            /* decide which points are in the Sea */
            wsphere[i*NY+j].indomain = 1;
            wsphere[i*NY+j].draw_wave = 1;
//             wsphere[i*NY+j].indomain = 0;
//             wsphere[i*NY+j].indomain = (diff < 15);
        }
    
    /* smooth colors in case of high resolution */
    if ((NX > nx)||(NY > ny))
        for (i=1; i<NX-1; i++)
            for (j=1; j<NY-1; j++)
            {
                wsphere[i*NY+j].r *= 0.2; 
                wsphere[i*NY+j].r += 0.2*wsphere[(i+1)*NY+j].r;
                wsphere[i*NY+j].r += 0.2*wsphere[(i-1)*NY+j].r;
                wsphere[i*NY+j].r += 0.2*wsphere[i*NY+j-1].r;
                wsphere[i*NY+j].r += 0.2*wsphere[i*NY+j+1].r;

                wsphere[i*NY+j].g *= 0.2; 
                wsphere[i*NY+j].g += 0.2*wsphere[(i+1)*NY+j].g;
                wsphere[i*NY+j].g += 0.2*wsphere[(i-1)*NY+j].g;
                wsphere[i*NY+j].g += 0.2*wsphere[i*NY+j-1].g;
                wsphere[i*NY+j].g += 0.2*wsphere[i*NY+j+1].g;

                wsphere[i*NY+j].b *= 0.2; 
                wsphere[i*NY+j].b += 0.2*wsphere[(i+1)*NY+j].b;
                wsphere[i*NY+j].b += 0.2*wsphere[(i-1)*NY+j].b;
                wsphere[i*NY+j].b += 0.2*wsphere[i*NY+j-1].b;
                wsphere[i*NY+j].b += 0.2*wsphere[i*NY+j+1].b;
            }
    
    free(rgb_values);
    fclose(image_file);
    
    if (ADD_DEM) init_dem(wsphere, dem_number);
    
    vshift = PLANET_SEALEVEL/(DEM_MAXHEIGHT - DEM_MAXDEPTH);
    for (i=0; i<NX; i++)
        for (j=0; j<NY; j++)
            wsphere[i*NY+j].indomain = (wsphere[i*NY+j].altitude < vshift);
}

int ij_to_sphere(int i, int j, double r, t_wave_sphere wsphere[NX*NY], double xyz[3], int use_wave_radius)
/* convert spherical to rectangular coordinates */
{
    double pscal, newr;
    static double norm_observer;
    static int first = 1;
    
    if (first)
    {
        norm_observer = sqrt(observer[0]*observer[0] + observer[1]*observer[1] + observer[2]*observer[2]);
        first = 0;
    }
    
    xyz[0] = wsphere[i*NY+j].x;
    xyz[1] = wsphere[i*NY+j].y;
    xyz[2] = wsphere[i*NY+j].z;
    
    pscal = xyz[0]*observer[0] + xyz[1]*observer[1] + xyz[2]*observer[2];
    
    if (use_wave_radius)
    {
        newr = wsphere[i*NY+j].radius;
        xyz[0] *= newr;
        xyz[1] *= newr;
        xyz[2] *= newr;
    }
    else
    {
        newr = wsphere[i*NY+j].radius_dem;
        xyz[0] *= newr;
        xyz[1] *= newr;
        xyz[2] *= newr;
    }
    
    return(pscal/norm_observer > COS_VISIBLE);
}


int ij_to_sphere_r(int i, int j, double r, t_wave_sphere wsphere[NX*NY], double xyz[3])
/* convert spherical to rectangular coordinates */
{
    double pscal, newr;
    static double norm_observer;
    static int first = 1;
    
    if (first)
    {
        norm_observer = sqrt(observer[0]*observer[0] + observer[1]*observer[1] + observer[2]*observer[2]);
        first = 0;
    }
    
    xyz[0] = wsphere[i*NY+j].x;
    xyz[1] = wsphere[i*NY+j].y;
    xyz[2] = wsphere[i*NY+j].z;
    
    pscal = xyz[0]*observer[0] + xyz[1]*observer[1] + xyz[2]*observer[2];
    
    xyz[0] *= r;
    xyz[1] *= r;
    xyz[2] *= r;
    
    return(pscal/norm_observer > COS_VISIBLE);
}


int init_circle_sphere(t_circles_sphere *circles, int circle_pattern)
/* initialise circles on sphere */
/* for billiard shape D_SPHERE_CIRCLES */
{
    int ncircles, k;
    double alpha, beta, gamma;
    
    switch (circle_pattern) {
        case (C_SPH_DODECA):
        {
            alpha = acos(sqrt(5.0)/3.0);
            beta = acos(1.0/3.0);
            gamma = asin(sqrt(3.0/8.0));
            
            circles[0].theta = 0.0;
            circles[0].phi = 0.0;
            
            for (k=0; k<3; k++)
            {
                circles[1+k].theta = alpha;
                circles[1+k].phi = (double)k*DPI/3.0;
            }

            for (k=0; k<3; k++)
            {
                circles[4+k].theta = beta;
                circles[4+k].phi = (double)k*DPI/3.0 + gamma;
                circles[7+k].theta = beta;
                circles[7+k].phi = (double)k*DPI/3.0 - gamma;
            }
            
            for (k=0; k<3; k++)
            {
                circles[10+k].theta = PI - beta;
                circles[10+k].phi = (double)k*DPI/3.0 + PI/3.0 + gamma;
                circles[13+k].theta = PI - beta;
                circles[13+k].phi = (double)k*DPI/3.0 + PI/3.0 - gamma;
            }

            for (k=0; k<3; k++)
            {
                circles[16+k].theta = PI - alpha;
                circles[16+k].phi = (double)k*DPI/3.0 + PI/3.0;
            }
            
            circles[19].theta = PI;
            circles[19].phi = 0.0;

            for (k=0; k<20; k++)
            {
                circles[k].radius = MU;
                circles[k].x = sin(circles[k].theta)*cos(circles[k].phi);
                circles[k].y = sin(circles[k].theta)*sin(circles[k].phi);
                circles[k].z = cos(circles[k].theta);
            }
            
            ncircles = 20;
            break;
        }
        case (C_SPH_ICOSA):
        {
            alpha = acos(1.0/sqrt(5.0));
            
            circles[0].theta = 0.0;
            circles[0].phi = 0.0;
       
            for (k=0; k<5; k++)
            {
                circles[1+k].theta = alpha;
                circles[1+k].phi = (double)k*DPI/5.0;
            }
            
            for (k=0; k<5; k++)
            {
                circles[6+k].theta = PI - alpha;
                circles[6+k].phi = (double)k*DPI/5.0 + PI/5.0;
            }

            circles[11].theta = PI;
            circles[11].phi = 0.0;

            for (k=0; k<12; k++) 
            {
                circles[k].radius = MU;
                circles[k].x = sin(circles[k].theta)*cos(circles[k].phi);
                circles[k].y = sin(circles[k].theta)*sin(circles[k].phi);
                circles[k].z = cos(circles[k].theta);
            }
            
            ncircles = 12;
            break;
        }
        case (C_SPH_OCTA):
        {
            circles[0].theta = 0.0;
            circles[0].phi = 0.0;
       
            for (k=0; k<4; k++)
            {
                circles[1+k].theta = PID;
                circles[1+k].phi = (double)k*PID;
            }
            
            circles[5].theta = PI;
            circles[5].phi = 0.0;

            for (k=0; k<6; k++) 
            {
                circles[k].radius = MU;
                circles[k].x = sin(circles[k].theta)*cos(circles[k].phi);
                circles[k].y = sin(circles[k].theta)*sin(circles[k].phi);
                circles[k].z = cos(circles[k].theta);
            }
            
            ncircles = 6;
            break;
        }
        case (C_SPH_CUBE):
        {
            alpha = acos(1.0/3.0);
            
            circles[0].theta = 0.0;
            circles[0].phi = 0.0;
       
            for (k=0; k<3; k++)
            {
                circles[1+k].theta = alpha;
                circles[1+k].phi = (double)k*DPI/3.0;
                circles[4+k].theta = PI - alpha;
                circles[4+k].phi = (double)k*DPI/3.0 + PI/3.0;
            }
            
            circles[7].theta = PI;
            circles[7].phi = 0.0;

            for (k=0; k<8; k++) 
            {
                circles[k].radius = MU;
                circles[k].x = sin(circles[k].theta)*cos(circles[k].phi);
                circles[k].y = sin(circles[k].theta)*sin(circles[k].phi);
                circles[k].z = cos(circles[k].theta);
            }
            
            ncircles = 8;
            break;
        }
        default: 
        {
            printf("Function init_circle_sphere not defined for this pattern \n");
        }
    }
    return(ncircles);
}

void init_sphere_maze(t_wave_sphere wsphere[NX*NY], int spiral, int wave, int npole)
/* initialize maze on sphere */
{
    int nblocks, block, i, j, k, n, p, q, np, na, inrect;
    double rmin, rmax, angle, r, dr, phi, dphi, phi1, ww, width = 0.02, theta, theta1, phaseshift; 
    t_maze* maze;
    t_rectangle* polyrect;
    
    maze = (t_maze *)malloc(NXMAZE*NYMAZE*sizeof(t_maze));
    polyrect = (t_rectangle *)malloc(NMAXPOLY*sizeof(t_rectangle));
    
    printf("Initializing maze\n");
    init_circular_maze(maze);
    printf("Maze initialized\n");
    
    printf("Building maze walls\n");
    /* build walls of maze */
    
    np = 0;
    na = 0;
    
    /* build walls of maze */
    nblocks = NYMAZE/NXMAZE;
    rmin = 0.3;
    rmax = PID + 0.2;
    angle = DPI/(double)nblocks;
        
    dr = (rmax - rmin)/(double)(NXMAZE);
    
    /* add straight walls */
    for (block = 0; block < nblocks; block++)
    {
        printf("Block %i\n", block);
        dphi = angle;
        
        /* first circle */
        n = nmaze(0, block*NXMAZE);
        r = rmin - 0.5*width;
        phi = (double)block*angle;
            
        if (maze[n].south)
        {
            polyrect[np].x1 = phi - 0.5*width;
            polyrect[np].y1 = r;
            polyrect[np].x2 = phi + 0.5*width;
            polyrect[np].y2 = r+dr+width;
//             printf("Adding rectangle at (%.3f, %.3f) - (%.3f, %.3f)\n", polyrect[np].x1, polyrect[np].y1, polyrect[np].x2, polyrect[np].y2);
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
                polyrect[np].x1 = phi - 0.5*width;
                polyrect[np].y1 = r;
                polyrect[np].x2 = phi + 0.5*width;
                polyrect[np].y2 = r+dr+width;
//                 printf("Adding rectangle at (%.3f, %.3f) - (%.3f, %.3f)\n", polyrect[np].x1, polyrect[np].y1, polyrect[np].x2, polyrect[np].y2);
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
                        polyrect[np].x1 = phi - 0.5*width;
                        polyrect[np].y1 = r;
                        polyrect[np].x2 = phi + 0.5*width;
                        polyrect[np].y2 = r+dr+width;
//                         printf("Adding rectangle at (%.3f, %.3f) - (%.3f, %.3f)\n", polyrect[np].x1, polyrect[np].y1, polyrect[np].x2, polyrect[np].y2);
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
            polyrect[np].x1 = phi;
            polyrect[np].y1 = r - 0.5*width;
            polyrect[np].x2 = phi + dphi + 0.05*width;
            polyrect[np].y2 = r + 0.5*width;
            np++;
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
                polyrect[np].x1 = phi;
                polyrect[np].y1 = r - 0.5*width;
                polyrect[np].x2 = phi + dphi + 0.05*width;
                polyrect[np].y2 = r + 0.5*width;
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
                r = rmin + (double)i*dr;
                printf("Circle, i = %i, dphi = %.2lg, r = %.2lg\n", i, dphi, r);
                for (q = 0; q < 2*ww; q++)
                {
                    j = block*NXMAZE + q;
                    n = nmaze(i,j);
                    phi = (double)(block)*angle + (double)q*dphi;
                    
                    if (maze[n].west)
                    {
                        polyrect[np].x1 = phi;
                        polyrect[np].y1 = r - 0.5*width;
                        polyrect[np].x2 = phi + dphi + 0.05*width;
                        polyrect[np].y2 = r + 0.5*width;
                        np++;
                    }
                }
                i++;
            }
            ww *= 2;
        }
    }
    
    /* outer boundary of maze */
    polyrect[np].x1 = dphi;
    polyrect[np].y1 = rmax - 0.5*width;
    polyrect[np].x2 = DPI;
    polyrect[np].y2 = rmax + 0.5*width;
    np++;
    

    printf("Maze walls built\n");
    
    phaseshift = ZERO_MERIDIAN*PI/180.0;
    for (i=0; i<NX; i++)
    {
        phi = DPI*(double)i/(double)NX + phaseshift;
        if (phi > DPI) phi -= DPI;
        for (j=0; j<NY; j++)
        {
            theta = PI*(1.0 - (double)j/(double)NY);
            if (spiral) 
            {
                phi -= PID/(double)NY;
                if (phi < 0.0) phi += DPI;
            }
            if (wave)
            {
                theta1 = theta; 
                theta += 0.03*cos(6.0*phi)*theta;
                if (theta > PI) theta = PI; 
                phi += 0.00075*cos(20.0*theta1)*theta1;
//                 if (phi > DPI) phi -= DPI;
                if (phi < 0.0) phi += DPI;
            }
            inrect = 0;
            for (n=0; n<np; n++) for (k=-1; k<2; k++)
            {
                phi1 = phi + k*DPI;
                if (((phi1 > polyrect[n].x1)&&(phi1 < polyrect[n].x2)&&(theta > polyrect[n].y1)&&(theta < polyrect[n].y2))) inrect = 1;
            }
            wsphere[i*NY+j].indomain = 1-inrect;
        }
        if (npole) for (j=NY-DPOLE*3/2; j<NY; j++) wsphere[i*NY+j].indomain = 0;
    }
    
    free(maze);
    free(polyrect);
    
}

int xy_in_billiard_sphere(int i, int j, t_wave_sphere wsphere[NX*NY])
/* returns 1 if (x,y) represents a point in the billiard */
{
    int k;
    double pscal, dist, r, u, v, u1, v1;
    static double cos_rot, sin_rot;
    static int first = 1;
    
    if (first)
    {
        if (B_DOMAIN == D_SPHERE_EARTH) init_earth_map(wsphere);
        else if (OTHER_PLANET) init_planet_map(wsphere, B_DOMAIN);
        else if ((B_DOMAIN == D_SPHERE_JULIA)||(B_DOMAIN == D_SPHERE_JULIA_INV)||(B_DOMAIN == D_SPHERE_JULIA_CUBIC))
        {
            cos_rot = cos(JULIA_ROT*DPI/360.0);
            sin_rot = sin(JULIA_ROT*DPI/360.0);
        }
        else if (B_DOMAIN == D_SPHERE_MAZE) init_sphere_maze(wsphere, 0, 0, 0);
        else if (B_DOMAIN == D_SPHERE_MAZE_SPIRAL) init_sphere_maze(wsphere, 1, 0, 1);
        else if (B_DOMAIN == D_SPHERE_MAZE_WAVE) init_sphere_maze(wsphere, 0, 1, 1);
        first = 0;
    }
    
    switch (B_DOMAIN) {
        case (D_NOTHING):
        {
            return(1);
        }
        case (D_LATITUDE):
        {
            return(vabs(wsphere[i*NY+j].theta - PID) < LAMBDA*PID);
        }
        case (D_SPHERE_CIRCLES):
        {
            for (k=0; k<ncircles; k++)
            {
                pscal = cos(wsphere[i*NY+j].phi - circ_sphere[k].phi)*wsphere[i*NY+j].stheta*sin(circ_sphere[k].theta);
                pscal += wsphere[i*NY+j].ctheta*cos(circ_sphere[k].theta);
            
                dist = acos(pscal);
                if (dist < circ_sphere[k].radius) return(0);
            }
            return(1);
        }
        case (D_SPHERE_JULIA):
        {
            if (wsphere[i*NY+j].z == 1.0) return(1);
            if (wsphere[i*NY+j].z == -1.0) return(0);
            r = (1.0 + wsphere[i*NY+j].ctheta)/wsphere[i*NY+j].stheta;
            u1 = r*wsphere[i*NY+j].cphi*JULIA_SCALE;
            v1 = r*wsphere[i*NY+j].sphi*JULIA_SCALE;
            u = u1*cos_rot + v1*sin_rot;
            v = -u1*sin_rot + v1*cos_rot;
            i = 0;
            while ((i<MANDELLEVEL)&&(u*u+v*v < 1000.0*MANDELLIMIT))
            {
                u1 = u*u - v*v + JULIA_RE;
                v = 2.0*u*v + JULIA_IM;
                u = u1;
                i++;
            }
            if (u*u + v*v < MANDELLIMIT) return(1);
            return(0);
        }
        case (D_SPHERE_JULIA_INV):
        {
            if (wsphere[i*NY+j].z == 1.0) return(1);
            if (wsphere[i*NY+j].z == -1.0) return(0);
            r = (1.0 - wsphere[i*NY+j].ctheta)/wsphere[i*NY+j].stheta;
            u1 = r*wsphere[i*NY+j].cphi*JULIA_SCALE;
            v1 = r*wsphere[i*NY+j].sphi*JULIA_SCALE;
            u = u1*cos_rot + v1*sin_rot;
            v = -u1*sin_rot + v1*cos_rot;
            i = 0;
            while ((i<MANDELLEVEL)&&(u*u+v*v < 1000.0*MANDELLIMIT))
            {
                u1 = u*u - v*v + JULIA_RE;
                v = 2.0*u*v + JULIA_IM;
                u = u1;
                i++;
            }
            if (u*u + v*v < MANDELLIMIT) return(0);
            return(1);
        }
        case (D_SPHERE_JULIA_CUBIC):
        {
            if (wsphere[i*NY+j].z == 1.0) return(1);
            if (wsphere[i*NY+j].z == -1.0) return(0);
            r = (1.0 + wsphere[i*NY+j].ctheta)/wsphere[i*NY+j].stheta;
            u1 = r*wsphere[i*NY+j].cphi*JULIA_SCALE;
            v1 = r*wsphere[i*NY+j].sphi*JULIA_SCALE;
            u = u1*cos_rot + v1*sin_rot;
            v = -u1*sin_rot + v1*cos_rot;
            i = 0;
            while ((i<MANDELLEVEL)&&(u*u+v*v < 1000.0*MANDELLIMIT))
            {
                u1 = u*u*u - 3.0*u*v*v + JULIA_RE;
                v = 3.0*u*u*v - v*v*v + JULIA_IM;
                u = u1;
                i++;
            }
            if (u*u + v*v < MANDELLIMIT) return(1);
            return(0);
        }
        default:
        {
            return(wsphere[i*NY+j].indomain);
        }
    }
}

void init_speed_dissipation_sphere(short int xy_in[NX*NY], double tc[NX*NY], double tcc[NX*NY], double tgamma[NX*NY], t_wave_sphere wsphere[NX*NY])
/* initialise fields for wave speed and dissipation */
{
    int i, j, k, n, inlens;
    double courant2 = COURANT*COURANT, courantb2 = COURANTB*COURANTB, lambda1, mu1;
    double u, v, u1, x, y, xy[2], norm2, speed, r2, c, salpha, h, ll, ca, sa, x1, y1, dx, dy, sum, sigma, x0, y0, rgb[3];
    double xc[NGRIDX*NGRIDY], yc[NGRIDX*NGRIDY], height[NGRIDX*NGRIDY];
    t_circle circles[NMAXCIRCLES];
    
    if (VARIABLE_IOR)
    {
        switch (IOR) {
            case (IOR_EARTH_DEM):
            {
//                 #pragma omp parallel for private(i,j,courant2,courantb2,speed)
                for (i=0; i<NX; i++){
                    for (j=DPOLE; j<NY-DPOLE; j++){
                        if (wsphere[i*NY+j].indomain) 
                        {
                            tcc[i*NY+j] = courant2;
                            tc[i*NY+j] = COURANT;
                            tgamma[i*NY+j] = GAMMA;
                        }
                        else 
                        {
                            speed = 1.0 - 10.0*wsphere[i*NY+j].altitude;
                            if (speed < 0.0) speed = 0.0;
                            tcc[i*NY+j] = courant2*speed + courantb2*(1.0 - speed); 
                            tc[i*NY+j] = COURANT*speed + COURANTB*(1.0 - speed); 
                            tgamma[i*NY+j] = GAMMA*speed + GAMMAB*(1.0 - speed);
                        }
                    }
                    for (j=0; j<DPOLE; j++) 
                    {
                        if (wsphere[i*NY+j].indomain)
                        {
                            tcc[i*NY+j] = courant2;
                            tc[i*NY+j] = COURANT;
                            tgamma[i*NY+j] = GAMMA;
                        }
                        else
                        {
                            tcc[i*NY+j] = courantb2;
                            tc[i*NY+j] = COURANTB;
                            tgamma[i*NY+j] = GAMMAB;
                        }    
                    }
                    for (j=NY-DPOLE; j<NY; j++) 
                    {
                        if (wsphere[i*NY+j].indomain)
                        {
                            tcc[i*NY+j] = courant2;
                            tc[i*NY+j] = COURANT;
                            tgamma[i*NY+j] = GAMMA;
                        }
                        else
                        {
                            tcc[i*NY+j] = courantb2;
                            tc[i*NY+j] = COURANTB;
                            tgamma[i*NY+j] = GAMMAB;
                        }    
                    }
                }
                break;
            }
            default:
            {
                for (i=0; i<NX; i++){
                    for (j=0; j<NY; j++){
                        tc[i*NY+j] = COURANT;
                        tcc[i*NY+j] = courant2;
                        tgamma[i*NY+j] = GAMMA;
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
                if (xy_in[i*NY+j] != 0)
                {
                    tc[i*NY+j] = COURANT;
                    tcc[i*NY+j] = courant2;
                    if (xy_in[i*NY+j] == 1) tgamma[i*NY+j] = GAMMA;
                    else tgamma[i*NY+j] = GAMMAB;
                }
                else if (TWOSPEEDS)
                {
                    tc[i*NY+j] = COURANTB;
                    tcc[i*NY+j] = courantb2;
                    tgamma[i*NY+j] = GAMMAB;
                }
            }
        }
    }
}

void init_wave_flat_sphere(double phi[NX*NY], double psi[NX*NY], short int xy_in[NX*NY], t_wave_sphere wsphere[NX*NY])
/* initialise field with drop - phi is wave height, psi is phi at time t-1 */
/* phi0 is longidude, theta0 is latitude */
{
    int i, j;
    
    printf("Initializing wave\n"); 
//     #pragma omp parallel for private(i,j,xy,dist2)
    for (i=0; i<NX; i++)
    {
        if (i%100 == 0) printf("Initializing column %i of %i\n", i, NX);
        for (j=0; j<NY; j++)
        {
            xy_in[i*NY+j] = xy_in_billiard_sphere(i, j, wsphere);
            
	    phi[i*NY+j] = 0.0;
            psi[i*NY+j] = 0.0;
        }
    }
}


void init_circular_wave_sphere(double phi0, double theta0, double phi[NX*NY], double psi[NX*NY], short int xy_in[NX*NY], t_wave_sphere wsphere[NX*NY])
/* initialise field with drop - phi is wave height, psi is phi at time t-1 */
/* phi0 is longidude, theta0 is latitude */
{
    int i, j, jmin, jmax;
    double xy[2], pscal, dist, dist2, stheta, ctheta, phishift;

    ctheta = cos(theta0 + PID);
    stheta = sin(theta0 + PID);
    phishift = ZERO_MERIDIAN*PI/180.0;
    
    /* safety distance to avoid singularity at poles */
    jmin = DPOLE;
    jmax = NY-DPOLE;
        
    printf("Initializing wave\n"); 
//     #pragma omp parallel for private(i,j,xy,dist2)
    for (i=0; i<NX; i++)
    {
        if (i%100 == 0) printf("Initializing column %i of %i\n", i, NX);
        for (j=jmin; j<jmax; j++)
        {
            pscal = cos(wsphere[i*NY+j].phi - phi0 - phishift)*wsphere[i*NY+j].stheta*stheta;
            pscal += wsphere[i*NY+j].ctheta*ctheta;
            
            dist = acos(pscal);
            dist2 = dist*dist;
            
	    xy_in[i*NY+j] = xy_in_billiard_sphere(i, j, wsphere);
            
	    if ((xy_in[i*NY+j])||(TWOSPEEDS)) 
                phi[i*NY+j] = INITIAL_AMP*exp(-dist2/INITIAL_VARIANCE)*cos(-sqrt(dist2)/INITIAL_WAVELENGTH);
            else phi[i*NY+j] = 0.0;
            psi[i*NY+j] = 0.0;
        }
        for (j=0; j<jmin; j++)
        {
            phi[i*NY+j] = 0.0;
            psi[i*NY+j] = 0.0;
        }
        for (j=jmax; j<NY; j++)
        {
            phi[i*NY+j] = 0.0;
            psi[i*NY+j] = 0.0;
        }
    }
}

void add_circular_wave_sphere(double amp, double phi0, double theta0, double phi[NX*NY], double psi[NX*NY], short int xy_in[NX*NY], t_wave_sphere wsphere[NX*NY])
/* initialise field with drop - phi is wave height, psi is phi at time t-1 */
/* phi0 is longidude, theta0 is latitude */
{
    int i, j, jmin, jmax;
    double xy[2], pscal, dist, dist2, stheta, ctheta, phishift;

    ctheta = cos(theta0 + PID);
    stheta = sin(theta0 + PID);
    phishift = ZERO_MERIDIAN*PI/180.0;
    
    /* safety distance to avoid singularity at poles */
    jmin = DPOLE;
    jmax = NY-DPOLE;
    
    #pragma omp parallel for private(i,j,xy,dist2)
    for (i=0; i<NX; i++)
    {
        for (j=jmin; j<jmax; j++)
        {
            pscal = cos(wsphere[i*NY+j].phi - phi0 - phishift)*wsphere[i*NY+j].stheta*stheta;
            pscal += wsphere[i*NY+j].ctheta*ctheta;
            
            dist = acos(pscal);
            dist2 = dist*dist;
            
	    if ((xy_in[i*NY+j])||(TWOSPEEDS)) 
                phi[i*NY+j] += amp*INITIAL_AMP*exp(-dist2/INITIAL_VARIANCE)*cos(-sqrt(dist2)/INITIAL_WAVELENGTH);
        }
    }
}


void init_tidal_wave_sphere(double phi0, double land_factor, double phi[NX*NY], double psi[NX*NY], short int xy_in[NX*NY], t_wave_sphere wsphere[NX*NY])
/* initialise field for "tidal" simulation - phi is wave height, psi is phi at time t-1 */
/* phi0 is longidude, theta0 is latitude */
{
    int i, j, jmin, jmax;
    double xy[2], pscal, dist, dist2, stheta, ctheta, phishift, phaseshift;

    phishift = ZERO_MERIDIAN*PI/180.0;
    phaseshift = 0.0;
    
    /* safety distance to avoid singularity at poles */
    jmin = DPOLE;
    jmax = NY-DPOLE;
        
    printf("Initializing wave\n"); 
//     #pragma omp parallel for private(i,j,xy,dist2)
    for (i=0; i<NX; i++)
    {
        if (i%100 == 0) printf("Initializing column %i of %i\n", i, NX);
        for (j=jmin; j<jmax; j++)
        {
	    xy_in[i*NY+j] = xy_in_billiard_sphere(i, j, wsphere);
            
	    if ((xy_in[i*NY+j])||(TWOSPEEDS)) 
//             if (xy_in[i*NY+j]) 
            {
                phi[i*NY+j] = INITIAL_AMP*wsphere[i*NY+j].stheta*sin(2.0*(wsphere[i*NY+j].phi - phi0 - phishift));
                if (xy_in[i*NY+j] == 0) phi[i*NY+j] *= land_factor;
            }
            else phi[i*NY+j] = 0.0;
            psi[i*NY+j] = 0.0;
        }
        for (j=0; j<jmin; j++)
        {
            phi[i*NY+j] = 0.0;
            psi[i*NY+j] = 0.0;
        }
        for (j=jmax; j<NY; j++)
        {
            phi[i*NY+j] = 0.0;
            psi[i*NY+j] = 0.0;
        }
    }
}


void init_moving_tidal_wave_sphere(double phi0, double voverc, double land_factor, double phi[NX*NY], double psi[NX*NY], short int xy_in[NX*NY], t_wave_sphere wsphere[NX*NY])
/* initialise field for "tidal" simulation - phi is wave height, psi is phi at time t-1 */
/* phi0 is longidude, theta0 is latitude */
{
    int i, j, jmin, jmax, inew;
    double xy[2], pscal, dist, dist2, stheta, ctheta, phishift, deltaphi, factor;

    phishift = ZERO_MERIDIAN*PI/180.0;
    deltaphi = voverc*COURANT/(double)NX;
    
//     phi1 = phi0 - vphi*DPI/(double)NX;
    
    /* safety distance to avoid singularity at poles */
    jmin = DPOLE;
    jmax = NY-DPOLE;
    
    printf("Initializing wave\n"); 
//     #pragma omp parallel for private(i,j,xy,dist2)
    for (i=0; i<NX; i++)
    {
        if (i%100 == 0) printf("Initializing column %i of %i\n", i, NX);
        for (j=jmin; j<jmax; j++)
        {
	    xy_in[i*NY+j] = xy_in_billiard_sphere(i, j, wsphere);
            
	    if ((xy_in[i*NY+j])||(TWOSPEEDS)) 
            {
                phi[i*NY+j] = INITIAL_AMP*wsphere[i*NY+j].stheta*sin(2.0*(wsphere[i*NY+j].phi - phi0 - phishift));
                psi[i*NY+j] = INITIAL_AMP*wsphere[i*NY+j].stheta*sin(2.0*(wsphere[i*NY+j].phi - phi0 - phishift - deltaphi));
//                 phi[i*NY+j] += 0.1*INITIAL_AMP*sin(2.0*wsphere[i*NY+j].theta)*sin(7.0*(wsphere[i*NY+j].phi - phi0 - phishift));
//                 psi[i*NY+j] += 0.1*INITIAL_AMP*sin(2.0*wsphere[i*NY+j].theta)*sin(7.0*(wsphere[i*NY+j].phi - phi0 - phishift - deltaphi));
            }
            else 
            {
                phi[i*NY+j] = 0.0;
                psi[i*NY+j] = 0.0;
            }
        }
        for (j=0; j<jmin; j++)
        {
            phi[i*NY+j] = 0.0;
            psi[i*NY+j] = 0.0;
        }
        for (j=jmax; j<NY; j++)
        {
            phi[i*NY+j] = 0.0;
            psi[i*NY+j] = 0.0;
        }
    }
        
    for (i=0; i<NX; i++)
    {
//         if (i%100 == 0) printf("Initializing column %i of %i\n", i, NX);
        for (j=jmin; j<jmax; j++) if (xy_in[i*NY+j] == 0) 
        {
            factor = 1.0 - land_factor*wsphere[i*NY+j].altitude;
            if (factor < 0.0) factor = 0.0;
            phi[i*NY+j] *= factor;
            psi[i*NY+j] *= factor;
        }
    }
}

void compute_light_angle_sphere(double phi[NX*NY], short int xy_in[NX*NY], t_wave wave[NX*NY], t_wave_sphere wsphere[NX*NY], int movie, int transparent)
/* computes cosine of angle between normal vector and vector light */
{
    int i, j;
    double x, y, z, norm, pscal, deltai[3], deltaj[3], deltar, n[3], r;
    static double dphi, dtheta, vshift;
    static int first = 1;
    
    if (first)
    {
        dphi = DPI/(double)NX;
        dtheta = PI/(double)NY;
        vshift = PLANET_SEALEVEL/(DEM_MAXHEIGHT - DEM_MAXDEPTH);
        first = 0;
    }
    
    #pragma omp parallel for private(i,j)
    for (i=0; i<NX; i++)
    {
        for (j=0; j<NY - DPOLE; j++) if ((TWOSPEEDS)||(xy_in[i*NY+j]))
        {
//             r = 1.01 + RSCALE*(*wave[i*NY+j].p_zfield[movie]);
            r = 1.0 + RSCALE*(*wave[i*NY+j].p_zfield[movie]);
            if (r > RMAX) r = RMAX;
            if (r < RMIN) r = RMIN;
            wsphere[i*NY+j].radius = r;
            
            if (FLOODING) 
                wsphere[i*NY+j].draw_wave = (phi[i*NY+j] >= wsphere[i*NY+j].altitude - vshift);
//                 wsphere[i*NY+j].draw_wave = ((wsphere[i*NY+j].indomain)||(phi[i*NY+j] >= wsphere[i*NY+j].altitude - vshift));
        }
        /* TODO ? Avoid artifacts due to singularity at north pole */
        for (j=NY - DPOLE; j<NY; j++) if ((TWOSPEEDS)||(xy_in[i*NY+j]))
        {
//             r = 1.0 + RSCALE*(*wave[i*NY+j].p_zfield[movie]);
            r = 1.01 + RSCALE*(*wave[i*NY+j].p_zfield[movie]);
            if (r > RMAX) r = RMAX;
            if (r < 1.0) r = 1.0;
            wsphere[i*NY+j].radius = r;
            
            if (FLOODING) 
                wsphere[i*NY+j].draw_wave = (phi[i*NY+j] >= wsphere[i*NY+j].altitude - vshift);
//                 wsphere[i*NY+j].draw_wave = ((wsphere[i*NY+j].indomain)||(phi[i*NY+j] >= wsphere[i*NY+j].altitude - vshift));
        }
    }
    
    if (SHADE_WAVE)
    {
        #pragma omp parallel for private(i,j,norm,pscal,deltar,deltai,deltaj,n)
        for (i=0; i<NX-1; i++)
            for (j=0; j<NY-1; j++)
            {
                if ((TWOSPEEDS)||(xy_in[i*NY+j]))
                {
                    /* computation of tangent vectors */
                    if (transparent)
                    {
                        deltar = (wsphere[(i+1)*NY+j].radius_dem - wsphere[i*NY+j].radius_dem)/dphi;
                    
                        deltai[0] = -wsphere[i*NY+j].radius_dem*wsphere[i*NY+j].sphi;
                        deltai[0] += deltar*wsphere[i*NY+j].cphi;
                    
                        deltai[1] = wsphere[i*NY+j].radius_dem*wsphere[i*NY+j].cphi;
                        deltai[1] += deltar*wsphere[i*NY+j].sphi;
                    
                        deltai[2] = -deltar*wsphere[i*NY+j].cottheta;
                    
                        deltar = (wsphere[i*NY+j+1].radius_dem - wsphere[i*NY+j].radius_dem)/dtheta;
                    
                        deltaj[0] = wsphere[i*NY+j].radius_dem*wsphere[i*NY+j].cphi*wsphere[i*NY+j].ctheta;
                        deltaj[0] += deltar*wsphere[i*NY+j].cphi*wsphere[i*NY+j].stheta;
                    
                        deltaj[1] = wsphere[i*NY+j].radius_dem*wsphere[i*NY+j].sphi*wsphere[i*NY+j].ctheta;
                        deltaj[1] += deltar*wsphere[i*NY+j].sphi*wsphere[i*NY+j].stheta;
                    
                        deltaj[2] = wsphere[i*NY+j].radius_dem*wsphere[i*NY+j].stheta;
                        deltaj[2] += -deltar*wsphere[i*NY+j].ctheta;
                    }
                    else
                    {
                        deltar = (wsphere[(i+1)*NY+j].radius - wsphere[i*NY+j].radius)/dphi;
                    
                        deltai[0] = -wsphere[i*NY+j].radius*wsphere[i*NY+j].sphi;
                        deltai[0] += deltar*wsphere[i*NY+j].cphi;
                    
                        deltai[1] = wsphere[i*NY+j].radius*wsphere[i*NY+j].cphi;
                        deltai[1] += deltar*wsphere[i*NY+j].sphi;
                    
                        deltai[2] = -deltar*wsphere[i*NY+j].cottheta;
                    
                        deltar = (wsphere[i*NY+j+1].radius - wsphere[i*NY+j].radius)/dtheta;
                    
                        deltaj[0] = wsphere[i*NY+j].radius*wsphere[i*NY+j].cphi*wsphere[i*NY+j].ctheta;
                        deltaj[0] += deltar*wsphere[i*NY+j].cphi*wsphere[i*NY+j].stheta;
                    
                        deltaj[1] = wsphere[i*NY+j].radius*wsphere[i*NY+j].sphi*wsphere[i*NY+j].ctheta;
                        deltaj[1] += deltar*wsphere[i*NY+j].sphi*wsphere[i*NY+j].stheta;
                    
                        deltaj[2] = wsphere[i*NY+j].radius*wsphere[i*NY+j].stheta;
                        deltaj[2] += -deltar*wsphere[i*NY+j].ctheta;
                    }

                    /* computation of normal vector */
                    n[0] = deltai[1]*deltaj[2] - deltai[2]*deltaj[1];
                    n[1] = deltai[2]*deltaj[0] - deltai[0]*deltaj[2];
                    n[2] = deltai[0]*deltaj[1] - deltai[1]*deltaj[0];
                                        
                    norm = sqrt(n[0]*n[0] + n[1]*n[1] + n[2]*n[2]);
                    
                    pscal = n[0]*light[0] + n[1]*light[1] + n[2]*light[2];
                
                    wave[i*NY+j].cos_angle = pscal/norm;
                }
//                 else wave[i*NY+j].cos_angle = wsphere[i*NY+j].cos_angle;
                else 
                {
                    pscal = wsphere[i*NY+j].x*light[0] + wsphere[i*NY+j].y*light[1] + wsphere[i*NY+j].z*light[2];
                
                    wave[i*NY+j].cos_angle = pscal;
                }
            }
            
        for (i=0; i<NX-1; i++) wave[i*NY+NY-1].cos_angle = wave[i*NY+NY-2].cos_angle;
//         for (j=0; j<NY-1; j++) wave[(NX-1)*NY+j].cos_angle = wave[(NX-2)*NY+j].cos_angle;
        wave[(NX-1)*NY+NY-1].cos_angle = wave[(NX-1)*NY+NY-2].cos_angle;
        
        /* i = NX-1 */
        for (j=0; j<NY-1; j++)
            {
                if ((TWOSPEEDS)||(xy_in[(NX-1)*NY+j]))
                {
                    /* computation of tangent vectors */
                    if (transparent)
                    {
                        deltar = (wsphere[j].radius_dem - wsphere[(NX-1)*NY+j].radius_dem)/dphi;
                    
                        deltai[0] = -wsphere[(NX-1)*NY+j].radius_dem*wsphere[(NX-1)*NY+j].sphi;
                        deltai[0] += deltar*wsphere[(NX-1)*NY+j].cphi;
                    
                        deltai[1] = wsphere[(NX-1)*NY+j].radius_dem*wsphere[(NX-1)*NY+j].cphi;
                        deltai[1] += deltar*wsphere[(NX-1)*NY+j].sphi;
                    
                        deltai[2] = -deltar*wsphere[(NX-1)*NY+j].cottheta;
                    
                        deltar = (wsphere[(NX-1)*NY+j+1].radius_dem - wsphere[(NX-1)*NY+j].radius_dem)/dtheta;
                    
                        deltaj[0] = wsphere[(NX-1)*NY+j].radius_dem*wsphere[(NX-1)*NY+j].cphi*wsphere[(NX-1)*NY+j].ctheta;
                        deltaj[0] += deltar*wsphere[(NX-1)*NY+j].cphi*wsphere[(NX-1)*NY+j].stheta;
                    
                        deltaj[1] = wsphere[(NX-1)*NY+j].radius_dem*wsphere[(NX-1)*NY+j].sphi*wsphere[(NX-1)*NY+j].ctheta;
                        deltaj[1] += deltar*wsphere[(NX-1)*NY+j].sphi*wsphere[(NX-1)*NY+j].stheta;
                    
                        deltaj[2] = wsphere[(NX-1)*NY+j].radius_dem*wsphere[(NX-1)*NY+j].stheta;
                        deltaj[2] += -deltar*wsphere[(NX-1)*NY+j].ctheta;
                    }
                    else
                    {
                        deltar = (wsphere[j].radius - wsphere[(NX-1)*NY+j].radius)/dphi;
                    
                        deltai[0] = -wsphere[(NX-1)*NY+j].radius*wsphere[(NX-1)*NY+j].sphi;
                        deltai[0] += deltar*wsphere[(NX-1)*NY+j].cphi;
                    
                        deltai[1] = wsphere[(NX-1)*NY+j].radius*wsphere[(NX-1)*NY+j].cphi;
                        deltai[1] += deltar*wsphere[(NX-1)*NY+j].sphi;
                    
                        deltai[2] = -deltar*wsphere[(NX-1)*NY+j].cottheta;
                    
                        deltar = (wsphere[(NX-1)*NY+j+1].radius - wsphere[(NX-1)*NY+j].radius)/dtheta;
                    
                        deltaj[0] = wsphere[(NX-1)*NY+j].radius*wsphere[(NX-1)*NY+j].cphi*wsphere[(NX-1)*NY+j].ctheta;
                        deltaj[0] += deltar*wsphere[(NX-1)*NY+j].cphi*wsphere[(NX-1)*NY+j].stheta;
                    
                        deltaj[1] = wsphere[(NX-1)*NY+j].radius*wsphere[(NX-1)*NY+j].sphi*wsphere[(NX-1)*NY+j].ctheta;
                        deltaj[1] += deltar*wsphere[(NX-1)*NY+j].sphi*wsphere[(NX-1)*NY+j].stheta;
                    
                        deltaj[2] = wsphere[(NX-1)*NY+j].radius*wsphere[(NX-1)*NY+j].stheta;
                        deltaj[2] += -deltar*wsphere[(NX-1)*NY+j].ctheta;
                    }

                    /* computation of normal vector */
                    n[0] = deltai[1]*deltaj[2] - deltai[2]*deltaj[1];
                    n[1] = deltai[2]*deltaj[0] - deltai[0]*deltaj[2];
                    n[2] = deltai[0]*deltaj[1] - deltai[1]*deltaj[0];
                                        
                    norm = sqrt(n[0]*n[0] + n[1]*n[1] + n[2]*n[2]);
                    
                    pscal = n[0]*light[0] + n[1]*light[1] + n[2]*light[2];
                
                    wave[(NX-1)*NY+j].cos_angle = pscal/norm;
                }
//                 else wave[i*NY+j].cos_angle = wsphere[i*NY+j].cos_angle;
                else 
                {
                    pscal = wsphere[(NX-1)*NY+j].x*light[0] + wsphere[(NX-1)*NY+j].y*light[1] + wsphere[(NX-1)*NY+j].z*light[2];
                
                    wave[(NX-1)*NY+j].cos_angle = pscal;
                }
            }
    }
    else
    {
        #pragma omp parallel for private(i,j,pscal)
        for (i=0; i<NX; i++)
            for (j=0; j<NY; j++)
            {
                if ((TWOSPEEDS)||(xy_in[i*NY+j]))
                {
                    pscal = wsphere[i*NY+j].x*light[0] + wsphere[i*NY+j].y*light[1] + wsphere[i*NY+j].z*light[2];
                
                    wave[i*NY+j].cos_angle = pscal;
                }
            }
    }
}

void compute_light_angle_sphere_2d(short int xy_in[NX*NY], double phi[NX*NY], t_wave wave[NX*NY], t_wave_sphere wsphere[NX*NY], int movie, int transparent)
/* computes cosine of angle between normal vector and vector light */
{
    int i, j;
    short int draw;
    double gradx, grady, norm, pscal;
    static double dx, dy, vscale2, vshift;
    static int first = 1;
    
    if (first)
    {
        dx = 2.0*(XMAX - XMIN)/(double)NX;
        dy = 2.0*(YMAX - YMIN)/(double)NY;
        vscale2 = SHADE_SCALE_2D*SHADE_SCALE_2D;
        vshift = PLANET_SEALEVEL/(DEM_MAXHEIGHT-DEM_MAXDEPTH);
        first = 0;
    }
    
    #pragma omp parallel for private(i,j,gradx, grady, norm, pscal)
    for (i=1; i<NX-1; i++)
        for (j=1; j<NY-1; j++)
        {
//             if ((TWOSPEEDS)||(xy_in[i*NY+j]))
//             if ((wsphere[i*NY+j].indomain)||(phi[i*NY+j] >= wsphere[i*NY+j].altitude + 0.001))
            if (FLOODING) 
            {
                draw = (phi[i*NY+j] >= wsphere[i*NY+j].altitude - vshift);
                wsphere[i*NY+j].draw_wave = draw;
            }
//             else draw = ((TWOSPEEDS)||(xy_in[i*NY+j]));
            else draw = ((xy_in[i*NY+j]));
            if (draw)
            {
                if (transparent)
                {
                    gradx = (wsphere[(i+1)*NY+j].radius_dem - wsphere[(i-1)*NY+j].radius_dem)/dx;
                    grady = (wsphere[i*NY+j+1].radius_dem - wsphere[i*NY+j-1].radius_dem)/dy;
                }
                else
                {
                    gradx = (*wave[(i+1)*NY+j].p_zfield[movie] - *wave[(i-1)*NY+j].p_zfield[movie])/dx;
                    grady = (*wave[i*NY+j+1].p_zfield[movie] - *wave[i*NY+j-1].p_zfield[movie])/dy;
                }
                
                norm = sqrt(vscale2 + gradx*gradx + grady*grady);
                pscal = -gradx*light[0] - grady*light[1] + SHADE_SCALE_2D;
                
                wave[i*NY+j].cos_angle = pscal/norm;
            }
            else wave[i*NY+j].cos_angle = wsphere[i*NY+j].cos_angle;
        }
    
    /* i=0 */
    for (j=1; j<NY-1; j++)
    {
//         if ((TWOSPEEDS)||(xy_in[j]))
        if (FLOODING) 
        {
            draw = (phi[j] >= wsphere[j].altitude - vshift);
            wsphere[j].draw_wave = draw;
        }
        else draw = ((xy_in[j]));
        
        if (draw)
        {
            if (transparent)
            {
                gradx = (wsphere[NY+j].radius_dem - wsphere[(NX-1)*NY+j].radius_dem)/dx;
                grady = (wsphere[j+1].radius_dem - wsphere[j-1].radius_dem)/dy;
            }
            else
            {
                gradx = (*wave[NY+j].p_zfield[movie] - *wave[(NY-1)*NY+j].p_zfield[movie])/dx;
                grady = (*wave[j+1].p_zfield[movie] - *wave[j-1].p_zfield[movie])/dy;
            }
            
            norm = sqrt(vscale2 + gradx*gradx + grady*grady);
            pscal = -gradx*light[0] - grady*light[1] + SHADE_SCALE_2D;
                
            wave[j].cos_angle = pscal/norm;
        }
        else wave[j].cos_angle = wsphere[j].cos_angle;
    }
    
    /* i=NX-1 */
    for (j=1; j<NY-1; j++)
    {
//         if ((TWOSPEEDS)||(xy_in[(NX-1)*NY+j]))
        if (FLOODING) 
        {
            draw = (phi[(NX-1)*NY+j] >= wsphere[(NX-1)*NY+j].altitude - vshift);
            wsphere[(NX-1)*NY+j].draw_wave = draw;
        }
        else draw = ((xy_in[(NX-1)*NY+j]));
        
        if (draw)
        {
            if (transparent)
            {
                gradx = (wsphere[j].radius_dem - wsphere[(NX-2)*NY+j].radius_dem)/dx;
                grady = (wsphere[(NX-1)*NY+j+1].radius_dem - wsphere[(NX-1)*NY+j-1].radius_dem)/dy;
            }
            else
            {
                gradx = (*wave[NY+j].p_zfield[movie] - *wave[(NY-1)*NY+j].p_zfield[movie])/dx;
                grady = (*wave[j+1].p_zfield[movie] - *wave[j-1].p_zfield[movie])/dy;
            }
            
            norm = sqrt(vscale2 + gradx*gradx + grady*grady);
            pscal = -gradx*light[0] - grady*light[1] + SHADE_SCALE_2D;
                
            wave[(NX-1)*NY+j].cos_angle = pscal/norm;
        }
        else wave[(NX-1)*NY+j].cos_angle = wsphere[(NX-1)*NY+j].cos_angle;
    }
}


void compute_cfield_sphere(short int xy_in[NX*NY], int cplot, int palette, t_wave wave[NX*NY], int fade, double fade_value, int movie)
/* compute the colors */
{
    int i, j, k;
    double ca;
    static double fact_max;
    static int first = 1;
    
    if (first)
    {
        fact_max = 0.5*COS_LIGHT_MAX;
        first = 0;
    }
       
    #pragma omp parallel for private(i,j,ca)
    for (i=0; i<NX; i++) for (j=0; j<NY; j++)
    {
        compute_field_color(*wave[i*NY+j].p_cfield[movie], *wave[i*NY+j].p_cfield[movie+2], cplot, palette, wave[i*NY+j].rgb);
        if ((SHADE_3D)||(SHADE_2D))
        {
            ca = wave[i*NY+j].cos_angle;
//             ca = (ca + 1.0)*0.4 + 0.2;
            ca = (ca + 1.0)*fact_max + COS_LIGHT_MIN;
            if ((FADE_IN_OBSTACLE)&&(!xy_in[i*NY+j])) ca = (ca + 0.1)*1.6;
            for (k=0; k<3; k++) wave[i*NY+j].rgb[k] *= ca;
        }
        if (fade) for (k=0; k<3; k++) wave[i*NY+j].rgb[k] *= fade_value;
    }
}


void draw_wave_sphere_2D(int movie, double phi[NX*NY], double psi[NX*NY], short int xy_in[NX*NY], t_wave wave[NX*NY], t_wave_sphere wsphere[NX*NY], int zplot, int cplot, int palette, int fade, double fade_value, int refresh)
{
    int i, j, ii;
    short int draw;
    double x, y, ca;
    static int ishift;
    static double dx, dy;
    
    if (1)
    {
        compute_wave_fields(phi, psi, xy_in, zplot, cplot, wave);
        if (SHADE_3D) compute_light_angle_sphere(phi, xy_in, wave, wsphere, movie, TRANSPARENT_WAVE);
        else if (SHADE_2D) 
            compute_light_angle_sphere_2d(xy_in, phi, wave, wsphere, movie, TRANSPARENT_WAVE);
        compute_cfield_sphere(xy_in, cplot, palette, wave, fade, fade_value, movie);
        dx = (XMAX - XMIN)/(double)NX;
        dy = (YMAX - YMIN)/(double)(NY-2*DPOLE);
        ishift = (int)((double)NX*PHISHIFT/360.0);
    }
    
    glBegin(GL_QUADS);
    
    for (i=0; i<NX; i++)
        for (j=0; j<NY; j++)
        {
//             if (FLOODING) draw = (wsphere[i*NY+j].indomain)||(phi[i*NY+j] >= wsphere[i*NY+j].altitude + 0.001);
//             else draw = (TWOSPEEDS)||(xy_in[i*NY+j]);
            if (FLOODING) draw = wsphere[i*NY+j].draw_wave;
            else draw = (xy_in[i*NY+j]);
            if (draw)
                glColor3f(wave[i*NY+j].rgb[0], wave[i*NY+j].rgb[1], wave[i*NY+j].rgb[2]);
            else
            {
                ca = wave[i*NY+j].cos_angle;
                ca = (ca + 1.0)*0.4 + 0.2;
                if (fade) ca *= fade_value;
                if (PLANET)
                    glColor3f(wsphere[i*NY+j].r*ca, wsphere[i*NY+j].g*ca, wsphere[i*NY+j].b*ca);
                else glColor3f(COLOR_OUT_R*ca, COLOR_OUT_G*ca, COLOR_OUT_B*ca);
            }

            x = wsphere[i*NY+j].x2d;
            y = wsphere[i*NY+j].y2d;
            
            ii = NX-i-1+ishift;
            if (ii > NX) ii -= NX;
            
            glVertex2i(ii, j);
            glVertex2i(ii+1, j);
            glVertex2i(ii+1, j+1);
            glVertex2i(ii, j+1);
        }
    glEnd ();
    
    if (DRAW_MOON_POSITION)
    {
        ii = NX-moon_position-1+ishift;
        if (ii > NX) ii -= NX;
        else if (ii < 0) ii += NX;
        glColor3f(fade_value, fade_value, fade_value);
        glBegin(GL_LINE_STRIP);
        glVertex2i(ii, 0);
        glVertex2i(ii, NY);
        glEnd();
    }
}

void draw_wave_sphere_ij(int i, int iplus, int j, int jplus, int jcolor, int movie, double phi[NX*NY], double psi[NX*NY], short int xy_in[NX*NY], t_wave wave[NX*NY], t_wave_sphere wsphere[NX*NY], int zplot, int cplot, int palette, int fade, double fade_value)
/* draw wave at simulation grid point (i,j) */
{
    int k, l;
    short int draw, notdraw, draw_bc=1;
    double xyz[3], ca;
    
    if (NON_DIRICHLET_BC) 
        draw_bc = (xy_in[i*NY+j])&&(xy_in[iplus*NY+j])&&(xy_in[i*NY+jplus])&&(xy_in[iplus*NY+jplus]);
    
//     if ((TWOSPEEDS)||(xy_in[i*NY+j]))
//     if (wsphere[i*NY+j].indomain)
//     if ((wsphere[i*NY+j].indomain)||(phi[i*NY+j] >= wsphere[i*NY+j].altitude + 0.001))
//     if (FLOODING) draw = (wsphere[i*NY+j].indomain)||(phi[i*NY+j] >= wsphere[i*NY+j].altitude + 0.001);
//     if (FLOODING) draw = (wsphere[i*NY+j].indomain)||(phi[i*NY+j] >= wsphere[i*NY+j].altitude);
    if (FLOODING) draw = wsphere[i*NY+j].draw_wave;
//     else draw = ((TWOSPEEDS)||(xy_in[i*NY+j]));
    else draw = (xy_in[i*NY+j]);
    if (draw) glColor3f(wave[i*NY+jcolor].rgb[0], wave[i*NY+jcolor].rgb[1], wave[i*NY+jcolor].rgb[2]);
    else
    {
//         ca = wave[i*NY+j].cos_angle;
        ca = wsphere[i*NY+j].cos_angle;
        ca = (ca + 1.0)*0.4 + 0.2;
        if (fade) ca *= fade_value;
        if (PLANET)
            glColor3f(wsphere[i*NY+j].r*ca, wsphere[i*NY+j].g*ca, wsphere[i*NY+j].b*ca);
        else glColor3f(COLOR_OUT_R*ca, COLOR_OUT_G*ca, COLOR_OUT_B*ca);
    }
    if (FLOODING)
    {
        glBegin(GL_TRIANGLE_FAN);
        if (ij_to_sphere(i, j, *wave[i*NY+j].p_zfield[movie], wsphere, xyz, wsphere[i*NY+j].draw_wave))
            draw_vertex_sphere(xyz);
        if (ij_to_sphere(iplus, j, *wave[iplus*NY+j].p_zfield[movie], wsphere, xyz, wsphere[iplus*NY+j].draw_wave))
            draw_vertex_sphere(xyz);
        if (ij_to_sphere(iplus, jplus, *wave[iplus*NY+j+1].p_zfield[movie], wsphere, xyz, wsphere[iplus*NY+j+1].draw_wave))
            draw_vertex_sphere(xyz);
        if (ij_to_sphere(i, jplus, *wave[i*NY+j+1].p_zfield[movie], wsphere, xyz, wsphere[i*NY+j+1].draw_wave))
            draw_vertex_sphere(xyz);
        glEnd ();        
    }
    else if (draw_bc)
    {
        notdraw = (!draw);
        glBegin(GL_TRIANGLE_FAN);
        if (ij_to_sphere(i, j, *wave[i*NY+j].p_zfield[movie], wsphere, xyz, notdraw))
            draw_vertex_sphere(xyz);
        if (ij_to_sphere(iplus, j, *wave[iplus*NY+j].p_zfield[movie], wsphere, xyz, notdraw))
            draw_vertex_sphere(xyz);
        if (ij_to_sphere(iplus, jplus, *wave[iplus*NY+j+1].p_zfield[movie], wsphere, xyz, notdraw))
            draw_vertex_sphere(xyz);
        if (ij_to_sphere(i, jplus, *wave[i*NY+j+1].p_zfield[movie], wsphere, xyz, notdraw))
            draw_vertex_sphere(xyz);
        glEnd ();
    }
}


void draw_wave_sphere_3D(int movie, double phi[NX*NY], double psi[NX*NY], short int xy_in[NX*NY], t_wave wave[NX*NY], t_wave_sphere wsphere[NX*NY], int zplot, int cplot, int palette, int fade, double fade_value, int refresh)
{
    int i, j, imax, imin, jmax;
    double observer_angle, angle2;
    
    blank();
            
    if (refresh)
    {
        compute_wave_fields(phi, psi, xy_in, zplot, cplot, wave);
        if (SHADE_3D) compute_light_angle_sphere(phi, xy_in, wave, wsphere, movie, TRANSPARENT_WAVE);
        else if (SHADE_2D) compute_light_angle_sphere_2d(xy_in, phi, wave, wsphere, movie, TRANSPARENT_WAVE);
        compute_cfield_sphere(xy_in, cplot, palette, wave, fade, fade_value, movie);
    }
    
    observer_angle = argument(observer[0], observer[1]);
    if (observer_angle < 0.0) observer_angle += DPI;
    
    angle2 = observer_angle + PI;
    if (angle2 > DPI) angle2 -= DPI;
    
    imin = (int)(observer_angle*(double)NX/DPI);
    imax = (int)(angle2*(double)NX/DPI);
    if (imin >= NX-1) imin = NX-2;
    if (imax >= NX-1) imax = NX-2;
    
//     jmax = NY-25;
    jmax = NY-3;
    
//     printf("Angle = %.5lg, angle2 = %.5lg, imin = %i, imax = %i\n", observer_angle, angle2, imin, imax);
    
    if (observer[2] > 0.0)
    {
        if (imin < imax)
        {
            for (i=imax; i>imin; i--)
                for (j=0; j<=jmax; j++)
                    draw_wave_sphere_ij(i, i+1, j, j+1, j+1, movie, phi, psi, xy_in, wave, wsphere, zplot, cplot, palette, fade, fade_value);
        
            for (i=imax+1; i<NX-1; i++)
                for (j=0; j<=jmax; j++)
                    draw_wave_sphere_ij(i, i+1, j, j+1, j+1, movie, phi, psi, xy_in, wave, wsphere, zplot, cplot, palette, fade, fade_value);
        
            for (j=0; j<=jmax; j++)
                draw_wave_sphere_ij(NX-1, 0, j, j+1, j+1, movie, phi, psi, xy_in, wave, wsphere, zplot, cplot, palette, fade, fade_value);
        
            for (i=0; i<=imin; i++)
                for (j=0; j<=jmax; j++)
                    draw_wave_sphere_ij(i, i+1, j, j+1, j+1, movie, phi, psi, xy_in, wave, wsphere, zplot, cplot, palette, fade, fade_value);
        }
        else
        {
            for (i=imax; i<imin; i++)
                for (j=0; j<=jmax; j++)
                    draw_wave_sphere_ij(i, i+1, j, j+1, j+1, movie, phi, psi, xy_in, wave, wsphere, zplot, cplot, palette, fade, fade_value);
        
            for (i=imax-1; i>=0; i--)
                for (j=0; j<=jmax; j++)
                    draw_wave_sphere_ij(i, i+1, j, j+1, j+1, movie, phi, psi, xy_in, wave, wsphere, zplot, cplot, palette, fade, fade_value);
        
            for (j=0; j<=jmax; j++)
                draw_wave_sphere_ij(NX-1, 0, j, j+1, j+1, movie, phi, psi, xy_in, wave, wsphere, zplot, cplot, palette, fade, fade_value);
        
            for (i=NX-2; i>=imin; i--)
                for (j=0; j<=jmax; j++)
                    draw_wave_sphere_ij(i, i+1, j, j+1, j+1, movie, phi, psi, xy_in, wave, wsphere, zplot, cplot, palette, fade, fade_value);
        }
    
        /* North pole */
        for (i=0; i<NX-1; i++) 
            draw_wave_sphere_ij(i, i+1, NY-3, NY-2, NY-2, movie, phi, psi, xy_in, wave, wsphere, zplot, cplot, palette, fade, fade_value);
        
        draw_wave_sphere_ij(NX-1, 0, NY-3, NY-2, NY-DPOLE, movie, phi, psi, xy_in, wave, wsphere, zplot, cplot, palette, fade, fade_value);
    }
    else
    {
        if (imin < imax)
        {
            for (i=imax; i>imin; i--)
                for (j=jmax; j>=0; j--)
                    draw_wave_sphere_ij(i, i+1, j, j+1, j+1, movie, phi, psi, xy_in, wave, wsphere, zplot, cplot, palette, fade, fade_value);
        
            for (i=imax+1; i<NX-1; i++)
                for (j=jmax; j>=0; j--)
                    draw_wave_sphere_ij(i, i+1, j, j+1, j+1, movie, phi, psi, xy_in, wave, wsphere, zplot, cplot, palette, fade, fade_value);
        
            for (j=jmax; j>=0; j--)
                draw_wave_sphere_ij(NX-1, 0, j, j+1, j+1, movie, phi, psi, xy_in, wave, wsphere, zplot, cplot, palette, fade, fade_value);
        
            for (i=0; i<=imin; i++)
                for (j=jmax; j>=0; j--)
                    draw_wave_sphere_ij(i, i+1, j, j+1, j+1, movie, phi, psi, xy_in, wave, wsphere, zplot, cplot, palette, fade, fade_value);
        }
        else
        {
            for (i=imax; i<imin; i++)
                for (j=jmax; j>=0; j--)
                    draw_wave_sphere_ij(i, i+1, j, j+1, j+1, movie, phi, psi, xy_in, wave, wsphere, zplot, cplot, palette, fade, fade_value);
        
            for (i=imax-1; i>=0; i--)
                for (j=jmax; j>=0; j--)
                    draw_wave_sphere_ij(i, i+1, j, j+1, j+1, movie, phi, psi, xy_in, wave, wsphere, zplot, cplot, palette, fade, fade_value);
        
            for (j=jmax; j>=0; j--)
                draw_wave_sphere_ij(NX-1, 0, j, j+1, j+1, movie, phi, psi, xy_in, wave, wsphere, zplot, cplot, palette, fade, fade_value);
        
//             for (i=NX-2; i>=imin; i--)
            for (i=NX-2; i>=imin-1; i--)
                for (j=jmax; j>=0; j--)
                    draw_wave_sphere_ij(i, i+1, j, j+1, j+1, movie, phi, psi, xy_in, wave, wsphere, zplot, cplot, palette, fade, fade_value);
        }
    
        /* South pole */
        for (i=0; i<NX-1; i++) for (j=2; j>0; j--)
            draw_wave_sphere_ij(i, i+1, j-1, j, j, movie, phi, psi, xy_in, wave, wsphere, zplot, cplot, palette, fade, fade_value);
        
        for (j=2; j>0; j--)
            draw_wave_sphere_ij(NX-1, 0, j-1, j, DPOLE, movie, phi, psi, xy_in, wave, wsphere, zplot, cplot, palette, fade, fade_value);
    }
}

void draw_wave_sphere(int movie, double phi[NX*NY], double psi[NX*NY], short int xy_in[NX*NY], t_wave wave[NX*NY], t_wave_sphere wsphere[NX*NY], int zplot, int cplot, int palette, int fade, double fade_value, int refresh)
{
    if (PLOT_2D) draw_wave_sphere_2D(movie, phi, psi, xy_in, wave, wsphere, zplot, cplot, palette, fade,fade_value, refresh);
    else draw_wave_sphere_3D(movie, phi, psi, xy_in, wave, wsphere, zplot, cplot, palette, fade,fade_value, refresh);
}
