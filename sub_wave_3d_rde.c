/*********************/
/* animation part    */
/*********************/



void init_3d()		/* initialisation of window */
{
    glLineWidth(3);

    glClearColor(0.0, 0.0, 0.0, 1.0);
    glClear(GL_COLOR_BUFFER_BIT);

    glOrtho(XMIN, XMAX, YMIN, YMAX , -1.0, 1.0);
}


void init_wave_sphere_rde(t_wave_sphere *wsphere, int res)
/* initialize sphere data, taken from sub_sphere.c */
/* wsphere is assumed to have size res*res*NX*NY */
{
    int i, j;
    double dphi, dtheta, theta0, xy[2], phishift, reg_cot;
    
    printf("Initializing wsphere\n");
    
    dphi = DPI/(double)(res*NX);
    dtheta = PI/(double)(res*NY);
//     dtheta = PI/(double)(NY-2*(DPOLE));
//     theta0 = (double)(DPOLE)*dtheta;
    theta0 = 0;
    phishift = PHISHIFT*(XMAX-XMIN)/360.0;
    
    #pragma omp parallel for private(i,j)
    for (i=0; i<res*NX; i++)
    {
//         for (j=DPOLE; j<NY-DPOLE; j++)
//             wsphere[i*NY+j].theta = (double)j*dtheta - theta0;
//         for (j=0; j<DPOLE; j++) wsphere[i*NY+j].theta = 0.0;
//         for (j=NY-DPOLE; j<NY; j++) wsphere[i*NY+j].theta = PI;
                
        for (j=0; j<res*NY; j++)
        {
            wsphere[i*res*NY+j].phi = (double)i*dphi;
            wsphere[i*res*NY+j].theta = (double)j*dtheta;
            
            wsphere[i*res*NY+j].cphi = cos(wsphere[i*res*NY+j].phi);
            wsphere[i*res*NY+j].sphi = sin(wsphere[i*res*NY+j].phi);
            
            wsphere[i*res*NY+j].ctheta = cos(wsphere[i*res*NY+j].theta);
            wsphere[i*res*NY+j].stheta = sin(wsphere[i*res*NY+j].theta);
            
            wsphere[i*res*NY+j].x = wsphere[i*res*NY+j].cphi*wsphere[i*res*NY+j].stheta;
            wsphere[i*res*NY+j].y = wsphere[i*res*NY+j].sphi*wsphere[i*res*NY+j].stheta;
            wsphere[i*res*NY+j].z = -wsphere[i*res*NY+j].ctheta;
            
            wsphere[i*res*NY+j].radius = 1.0;
            wsphere[i*res*NY+j].radius_dem = 1.0;
            
            ij_to_xy(NX-1-i,j,xy);
//             xy[0] = XMIN + ((double)(NX-i-1))*(XMAX-XMIN)/((double)NX);
//             xy[1] = YMIN + ((double)(j-DPOLE))*(YMAX-YMIN)/((double)(NY-2*DPOLE));
            
            xy[0] += phishift;
            if (xy[0] > XMAX) xy[0] += XMIN - XMAX;
            
            xy[1] *= (double)NY/(double)(NY-2*DPOLE);
            
            wsphere[i*res*NY+j].x2d = xy[0]; 
            wsphere[i*res*NY+j].y2d = xy[1]; 
            
            wsphere[i*res*NY+j].cos_angle_sphere = wsphere[i*res*NY+j].x*light[0] + wsphere[i*res*NY+j].y*light[1] + wsphere[i*res*NY+j].z*light[2];
            
            /* default value, to be changed by init_dem */
            wsphere[i*res*NY+j].evolve_wave = 1;
        }
        
        /* cotangent, taking care of not dividing by zero */
        /* TODO clean up cottheta range ? */
        for (j=DPOLE; j<res*NY-DPOLE; j++) wsphere[i*res*NY+j].cottheta = wsphere[i*res*NY+j].ctheta/wsphere[i*res*NY+j].stheta;
        for (j=0; j<DPOLE; j++) wsphere[i*res*NY+j].cottheta = wsphere[i*res*NY+DPOLE].cottheta;
        for (j=res*NY-res*DPOLE; j<NY; j++) wsphere[i*res*NY+j].cottheta = wsphere[i*res*NY+DPOLE-1].cottheta;       
        
               
        
    }
    
    /* regularized cotangent */
    
    for (j=0; j<res*NY; j++) 
    {
        reg_cot = wsphere[j].ctheta/sqrt(1.0 + SMOOTHCOTPOLE*wsphere[j].stheta*wsphere[j].stheta);
        for (i=0; i<res*NX; i++)
        {
            wsphere[i*res*NY+j].reg_cottheta = reg_cot;
        }   
    }
}


void read_negative_dem_values_rde(double *height_values, t_wave_sphere *wsphere)
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

//                 if (int_height_values[3*(jj*nx+ii)] > hcont)
//                 {
//                     wsphere[i*NY+j].r = 0.9*(double)rgb_values[3*(jj*nx+ii)]*cratio;
//                     wsphere[i*NY+j].g = 0.9*(double)rgb_values[3*(jj*nx+ii)+1]*cratio;
//                     wsphere[i*NY+j].b = 0.9*(double)rgb_values[3*(jj*nx+ii)+2]*cratio;
//                 }
//                 else 
//                 {
//                     wsphere[i*NY+j].r = 0.29;
//                     wsphere[i*NY+j].g = 0.29;
//                     wsphere[i*NY+j].b = 0.29;
//                 }
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

void init_dem_rde(t_wave_sphere *wsphere, int dem_number, int res)
/* init heights from digital elevation map */
{
    int i, j, ii, jj, k, nx, ny, maxrgb, nmaxpixels = 4915200, scan, rgbval, diff, sshift, nshift, hmin, hmax, ishift, hsum, rnx, rny;
    int *rgb_values;
    double cratio, rx, ry, cy, dx, dy, pscal, norm, vscale1, vscale2, gradx, grady, deltar, deltai[3], deltaj[3], dphi, dtheta, n[3], hsea, hmean, vshift, altitude;
    double *height_values, *height_values_tmp;
    FILE *image_file;
    
    printf("Reading digital elevation model\n");
    rnx = res*NX;
    rny = res*NY;
    
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
    height_values = (double *)malloc(rnx*rny*sizeof(double));
    height_values_tmp = (double *)malloc(rnx*rny*sizeof(double));
    
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
    rx = (double)nx/(double)(rnx);
    ry = (double)ny/(double)(rny - sshift - nshift);
    
    /* build height table */
    vshift = PLANET_SEALEVEL/(DEM_MAXHEIGHT - DEM_MAXDEPTH);
    for (i=0; i<rnx; i++)
        for (j=0; j<rny; j++)
        {
            ii = (int)(rx*(double)(rnx-1 - i)) + nx/2 + ishift;
            if (ii > nx-1) ii -= nx;
            if (ii < 0) ii = 0;
            jj = (int)(ry*(double)(rny-nshift - j));
            if (jj > ny-1) jj = ny-1;
            if (jj < 0) jj = 0;
            height_values[i*rny+j] = ((double)rgb_values[3*(jj*nx+ii)]-hsea)*cratio;
            wsphere[i*rny+j].altitude = ((double)rgb_values[3*(jj*nx+ii)]-hsea)*cratio;
            
            /* take care of black areas (missing data) on venus */
            if ((B_DOMAIN == D_SPHERE_VENUS)&&(rgb_values[3*(jj*nx+ii)] == 0))
            {
                height_values[i*rny+j] = VENUS_NODATA_FACTOR*hmean*cratio;
                wsphere[i*rny+j].altitude = VENUS_NODATA_FACTOR*hmean*cratio;
            }
            
            if (OTHER_PLANET)
                wsphere[i*rny+j].indomain = (wsphere[i*rny+j].altitude < vshift);
            
//             if (wsphere[i*rny+j].indomain) printf("rgb = %i, altitude = %.3lg\n", rgb_values[3*(jj*nx+ii)], height_values[i*rny+j]);
        }
        
    /* smooth values in case of high resolution */
    if ((rnx > nx)||(rny > ny))
    {
        for (i=1; i<rnx-1; i++)
            for (j=1; j<rny-1; j++)
            {
                height_values[i*rny+j] *= 0.2; 
                height_values[i*rny+j] += 0.2*height_values[(i+1)*rny+j];
                height_values[i*rny+j] += 0.2*height_values[(i-1)*rny+j];
                height_values[i*rny+j] += 0.2*height_values[i*rny+j-1];
                height_values[i*rny+j] += 0.2*height_values[i*rny+j+1];
                
                wsphere[i*rny+j].altitude *= 0.2; 
                wsphere[i*rny+j].altitude += 0.2*wsphere[(i+1)*rny+j].altitude;
                wsphere[i*rny+j].altitude += 0.2*wsphere[(i-1)*rny+j].altitude;
                wsphere[i*rny+j].altitude += 0.2*wsphere[i*rny+j-1].altitude;
                wsphere[i*rny+j].altitude += 0.2*wsphere[i*rny+j+1].altitude;
            }
            
        /* i = 0 */
        for (j=1; j<rny-1; j++)
        {
            height_values[j] *= 0.2; 
            height_values[j] += 0.2*height_values[rny+j];
            height_values[j] += 0.2*height_values[(rnx-1)*rny+j];
            height_values[j] += 0.2*height_values[j-1];
            height_values[j] += 0.2*height_values[j+1];
                
            wsphere[j].altitude *= 0.2; 
            wsphere[j].altitude += 0.2*wsphere[rny+j].altitude;
            wsphere[j].altitude += 0.2*wsphere[(rnx-1)*rny+j].altitude;
            wsphere[j].altitude += 0.2*wsphere[j-1].altitude;
            wsphere[j].altitude += 0.2*wsphere[j+1].altitude;
        }
        
        /* i = rny-1 */
        for (j=1; j<rny-1; j++)
        {
            height_values[(rny-1)*rny+j] *= 0.2; 
            height_values[(rny-1)*rny+j] += 0.2*height_values[j];
            height_values[(rny-1)*rny+j] += 0.2*height_values[(rny-2)*rny+j];
            height_values[(rny-1)*rny+j] += 0.2*height_values[(rny-1)*rny+j-1];
            height_values[(rny-1)*rny+j] += 0.2*height_values[(rny-1)*rny+j+1];
                
            wsphere[(rny-1)*rny+j].altitude *= 0.2; 
            wsphere[(rny-1)*rny+j].altitude += 0.2*wsphere[j].altitude;
            wsphere[(rny-1)*rny+j].altitude += 0.2*wsphere[(rny-2)*rny+j].altitude;
            wsphere[(rny-1)*rny+j].altitude += 0.2*wsphere[(rny-1)*rny+j-1].altitude;
            wsphere[(rny-1)*rny+j].altitude += 0.2*wsphere[(rny-1)*rny+j+1].altitude;
        }
    }
        
    printf("Closing rgb_values\n");
    fclose(image_file);
    free(rgb_values);
    
    /* smoothen values at low altitude */
    if (SMOOTH_DEM) for (k=1; k<DEM_SMOOTH_STEPS; k++)
    {
        printf("Smoothing step %i\n", k);
        for (i=1; i<rnx-1; i++)
            for (j=1; j<rny-1; j++) 
                if ((!wsphere[i*rny+j].indomain)&&(height_values[i*rny+j] <= DEM_SMOOTH_HEIGHT))
                {
                    height_values_tmp[i*rny+j] = height_values[i*rny+j] + 0.1*(height_values[(i+1)*rny+j] + height_values[(i-1)*rny+j] + height_values[i*rny+j+1] + height_values[i*rny+j-1] - 4.0*height_values[i*rny+j]);
            
                    height_values[i*rny+j] = height_values_tmp[i*rny+j] + 0.1*(height_values_tmp[(i+1)*rny+j] + height_values_tmp[(i-1)*rny+j] + height_values_tmp[i*rny+j+1] + height_values_tmp[i*rny+j-1] - 4.0*height_values_tmp[i*rny+j]);
                }
        /* i = 0 */
        for (j=1; j<rny-1; j++) 
            if ((!wsphere[j].indomain)&&(height_values[j] <= DEM_SMOOTH_HEIGHT))
            {
                height_values_tmp[j] = height_values[j] + 0.1*(height_values[rny+j] + height_values[(rnx-1)*rny+j] + height_values[j+1] + height_values[j-1] - 4.0*height_values[j]);
            
                height_values[j] = height_values_tmp[j] + 0.1*(height_values_tmp[rny+j] + height_values_tmp[(rnx-1)*rny+j] + height_values_tmp[j+1] + height_values_tmp[j-1] - 4.0*height_values_tmp[j]);
            }
        /* i = rny-1 */
        for (j=1; j<rny-1; j++) 
            if ((!wsphere[(rnx-1)*rny+j].indomain)&&(height_values[(rnx-1)*rny+j] <= DEM_SMOOTH_HEIGHT))
            {
                height_values_tmp[(rnx-1)*rny+j] = height_values[(rnx-1)*rny+j] + 0.1*(height_values[j] + height_values[(rnx-2)*rny+j] + height_values[(rnx-1)*rny+j+1] + height_values[(rnx-1)*rny+j-1] - 4.0*height_values[(rnx-1)*rny+j]);
            
                height_values[(rnx-1)*rny+j] = height_values_tmp[(rnx-1)*rny+j] + 0.1*(height_values_tmp[j] + height_values_tmp[(rnx-2)*rny+j] + height_values_tmp[(rnx-1)*rny+j+1] + height_values_tmp[(rnx-1)*rny+j-1] - 4.0*height_values_tmp[(rnx-1)*rny+j]);
            }
    }
            
    if (SMOOTH_DEM) for (i=0; i<rnx; i++)
        for (j=1; j<rny-1; j++) 
            if ((!wsphere[i*rny+j].indomain)&&(wsphere[i*rny+j].altitude <= DEM_SMOOTH_HEIGHT))
            {
                wsphere[i*rny+j].altitude = height_values[i*rny+j];
            }
                
    if ((ADD_NEGATIVE_DEM)&&(dem_number == DEM_EARTH)) 
        read_negative_dem_values_rde(height_values, wsphere);
            
    /* set radius */
//     for (i=0; i<rnx; i++)
//         for (j=0; j<rny; j++)
//         {
//             altitude = wsphere[i*rny+j].altitude - vshift;
//             if (!wsphere[i*rny+j].indomain) wsphere[i*rny+j].radius_dem = 1.0 + RSCALE_DEM*altitude;
// //             if (altitude >= 0.0) wsphere[i*rny+j].radius_dem = 1.0 + RSCALE_DEM*altitude;
//             else wsphere[i*rny+j].radius_dem = 1.0;
//         }

    /* set domain in which wave is evolved */
    for (i=0; i<rnx; i++)
        for (j=0; j<rny; j++)
            wsphere[i*rny+j].evolve_wave = (wsphere[i*rny+j].altitude < vshift + 0.01); 
    
    /* compute light angle */  
    dx = 2.0*(XMAX - XMIN)/(double)rnx;
    dy = 2.0*(YMAX - YMIN)/(double)rny;
    vscale1 = 0.1*SHADE_SCALE_2D;
    vscale2 = vscale1*vscale1;
    
    if (SHADE_2D)
    {
        for (i=1; i<rnx-1; i++)
            for (j=1; j<rny-1; j++)
            {
                gradx = (wsphere[(i+1)*rny+j].radius_dem - wsphere[(i-1)*rny+j].radius_dem)/dx;
                grady = (wsphere[i*rny+j+1].radius_dem - wsphere[i*rny+j-1].radius_dem)/dy;
                
                norm = sqrt(vscale2 + gradx*gradx + grady*grady);
                pscal = -gradx*light[0] - grady*light[1] + vscale1;
                
                wsphere[i*rny+j].cos_angle = pscal/norm;
            }
        /* i = 0 */
        for (j=1; j<rny-1; j++)
        {
            gradx = (wsphere[rny+j].radius_dem - wsphere[(rnx-1)*rny+j].radius_dem)/dx;
            grady = (wsphere[j+1].radius_dem - wsphere[j-1].radius_dem)/dy;
            
            norm = sqrt(vscale2 + gradx*gradx + grady*grady);
            pscal = -gradx*light[0] - grady*light[1] + vscale1;
            
            wsphere[j].cos_angle = pscal/norm;
        }
        /* i = N-1 */
        for (j=1; j<rny-1; j++)
        {
            gradx = (wsphere[j].radius_dem - wsphere[(rnx-2)*rny+j].radius_dem)/dx;
            grady = (wsphere[(rnx-1)*rny+j+1].radius_dem - wsphere[(rnx-1)*rny+j-1].radius_dem)/dy;
            
            norm = sqrt(vscale2 + gradx*gradx + grady*grady);
            pscal = -gradx*light[0] - grady*light[1] + vscale1;
            
            wsphere[(rnx-1)*rny+j].cos_angle = pscal/norm;
        }
    }
    else if (SHADE_3D)
    {
        dphi = DPI/(double)rnx;
        dtheta = PI/(double)rny;
        
        for (i=1; i<rnx-1; i++)
            for (j=1; j<rny-1; j++)
            {                
                /* computation of tangent vectors */
                deltar = (wsphere[(i+1)*rny+j].radius_dem - wsphere[i*rny+j].radius_dem)/dphi;
                    
                deltai[0] = -wsphere[i*rny+j].radius_dem*wsphere[i*rny+j].sphi;
                deltai[0] += deltar*wsphere[i*rny+j].cphi;
                    
                deltai[1] = wsphere[i*rny+j].radius_dem*wsphere[i*rny+j].cphi;
                deltai[1] += deltar*wsphere[i*rny+j].sphi;
                    
                deltai[2] = -deltar*wsphere[i*rny+j].cottheta;
                    
                deltar = (wsphere[i*rny+j+1].radius_dem - wsphere[i*rny+j].radius_dem)/dtheta;
                    
                deltaj[0] = wsphere[i*rny+j].radius_dem*wsphere[i*rny+j].cphi*wsphere[i*rny+j].ctheta;
                deltaj[0] += deltar*wsphere[i*rny+j].cphi*wsphere[i*rny+j].stheta;
                    
                deltaj[1] = wsphere[i*rny+j].radius_dem*wsphere[i*rny+j].sphi*wsphere[i*rny+j].ctheta;
                deltaj[1] += deltar*wsphere[i*rny+j].sphi*wsphere[i*rny+j].stheta;
                    
                deltaj[2] = wsphere[i*rny+j].radius_dem*wsphere[i*rny+j].stheta;
                deltaj[2] += -deltar*wsphere[i*rny+j].ctheta;
                    
                /* computation of normal vector */
                n[0] = deltai[1]*deltaj[2] - deltai[2]*deltaj[1];
                n[1] = deltai[2]*deltaj[0] - deltai[0]*deltaj[2];
                n[2] = deltai[0]*deltaj[1] - deltai[1]*deltaj[0];
                                        
                norm = sqrt(n[0]*n[0] + n[1]*n[1] + n[2]*n[2]);
                    
                pscal = n[0]*light[0] + n[1]*light[1] + n[2]*light[2];
                
                wsphere[i*rny+j].cos_angle = pscal/norm;
            }
        /* i = 0 */
        for (j=1; j<rny-1; j++)
        {                
            /* computation of tangent vectors */
            deltar = (wsphere[rny+j].radius_dem - wsphere[j].radius_dem)/dphi;
                    
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
        /* i = rnx-1 */
        for (j=1; j<rny-1; j++)
        {                
            /* computation of tangent vectors */
            deltar = (wsphere[j].radius_dem - wsphere[(rnx-1)*rny+j].radius_dem)/dphi;
                    
            deltai[0] = -wsphere[(rnx-1)*rny+j].radius_dem*wsphere[(rnx-1)*rny+j].sphi;
            deltai[0] += deltar*wsphere[(rnx-1)*rny+j].cphi;
                    
            deltai[1] = wsphere[(rnx-1)*rny+j].radius_dem*wsphere[(rnx-1)*rny+j].cphi;
            deltai[1] += deltar*wsphere[(rnx-1)*rny+j].sphi;
                    
            deltai[2] = -deltar*wsphere[(rnx-1)*rny+j].cottheta;
                    
            deltar = (wsphere[(rnx-1)*rny+j+1].radius_dem - wsphere[(rnx-1)*rny+j].radius_dem)/dtheta;
                    
            deltaj[0] = wsphere[(rnx-1)*rny+j].radius_dem*wsphere[(rnx-1)*rny+j].cphi*wsphere[(rnx-1)*rny+j].ctheta;
            deltaj[0] += deltar*wsphere[(rnx-1)*rny+j].cphi*wsphere[(rnx-1)*rny+j].stheta;
                    
            deltaj[1] = wsphere[(rnx-1)*rny+j].radius_dem*wsphere[(rnx-1)*rny+j].sphi*wsphere[(rnx-1)*rny+j].ctheta;
            deltaj[1] += deltar*wsphere[(rnx-1)*rny+j].sphi*wsphere[(rnx-1)*rny+j].stheta;
                    
            deltaj[2] = wsphere[(rnx-1)*rny+j].radius_dem*wsphere[(rnx-1)*rny+j].stheta;
            deltaj[2] += -deltar*wsphere[(rnx-1)*rny+j].ctheta;
                    
            /* computation of normal vector */
            n[0] = deltai[1]*deltaj[2] - deltai[2]*deltaj[1];
            n[1] = deltai[2]*deltaj[0] - deltai[0]*deltaj[2];
            n[2] = deltai[0]*deltaj[1] - deltai[1]*deltaj[0];
                                        
            norm = sqrt(n[0]*n[0] + n[1]*n[1] + n[2]*n[2]);
                    
            pscal = n[0]*light[0] + n[1]*light[1] + n[2]*light[2];
                
            wsphere[(rnx-1)*rny+j].cos_angle = pscal/norm;
        }
    }
    
    free(height_values);
    free(height_values_tmp);
}

void init_earth_map_rde(t_wave_sphere *wsphere, int res)
/* init file from earth map */
{
    int i, j, ii, jj, k, nx, ny, maxrgb, nmaxpixels = 4915200, scan, rgbval, diff, sshift, nshift, ishift;
    int *rgb_values;
    double cratio, rx, ry, cy, vshift;
    FILE *image_file;
    
    printf("Reading Earth map at resolution %i\n", res);
    
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
    rx = (double)nx/(double)(res*NX);
    ry = (double)ny/(double)(res*NY - sshift - nshift);
//     cy = rx*(double)(NY - nshift);
    
    /* build wave table */
    for (i=0; i<res*NX; i++)
        for (j=0; j<res*NY; j++)
        {
            ii = (int)(rx*(double)(res*NX-1 - i)) + nx/2 + ishift;
            if (ii > nx-1) ii -= nx;
//             jj = (int)(-ry*(double)j + cy);
//             jj = (int)(ry*(double)(NY-nshift - j)) + sshift;
            jj = (int)(ry*(double)(res*NY-nshift - j));
            if (jj > ny-1) jj = ny-1;
            if (jj < 0) jj = 0;
            wsphere[i*res*NY+j].r = (double)rgb_values[3*(jj*nx+ii)]*cratio;
            wsphere[i*res*NY+j].g = (double)rgb_values[3*(jj*nx+ii)+1]*cratio;
            wsphere[i*res*NY+j].b = (double)rgb_values[3*(jj*nx+ii)+2]*cratio;
            
//             printf("RGB at (%i, %i) = (%.3lg, %3.lg, %.3lg)\n", i, j, wsphere[i*NY+j].r, wsphere[i*NY+j].g, wsphere[i*NY+j].b);
            
            /* decide which points are in the Sea */
            diff = iabs(rgb_values[3*(jj*nx+ii)] - 10);
            diff += iabs(rgb_values[3*(jj*nx+ii)+1] - 10);
            diff += iabs(rgb_values[3*(jj*nx+ii)+2] - 51);
            wsphere[i*res*NY+j].indomain = (diff < 15);
        }
    
    /* smooth colors in case of high resolution */
    if ((res*NX > nx)||(res*NY > ny))
        for (i=1; i<res*NX-1; i++)
            for (j=1; j<res*NY-1; j++)
            {
                wsphere[i*res*NY+j].r *= 0.2; 
                wsphere[i*res*NY+j].r += 0.2*wsphere[(i+1)*res*NY+j].r;
                wsphere[i*res*NY+j].r += 0.2*wsphere[(i-1)*res*NY+j].r;
                wsphere[i*res*NY+j].r += 0.2*wsphere[i*res*NY+j-1].r;
                wsphere[i*res*NY+j].r += 0.2*wsphere[i*res*NY+j+1].r;

                wsphere[i*res*NY+j].g *= 0.2; 
                wsphere[i*res*NY+j].g += 0.2*wsphere[(i+1)*res*NY+j].g;
                wsphere[i*res*NY+j].g += 0.2*wsphere[(i-1)*res*NY+j].g;
                wsphere[i*res*NY+j].g += 0.2*wsphere[i*res*NY+j-1].g;
                wsphere[i*res*NY+j].g += 0.2*wsphere[i*res*NY+j+1].g;

                wsphere[i*res*NY+j].b *= 0.2; 
                wsphere[i*res*NY+j].b += 0.2*wsphere[(i+1)*res*NY+j].b;
                wsphere[i*res*NY+j].b += 0.2*wsphere[(i-1)*res*NY+j].b;
                wsphere[i*res*NY+j].b += 0.2*wsphere[i*res*NY+j-1].b;
                wsphere[i*res*NY+j].b += 0.2*wsphere[i*res*NY+j+1].b;
            }
    
    free(rgb_values);
    fclose(image_file);
    
//     if (ADD_DEM) 
    init_dem_rde(wsphere, DEM_EARTH, res);
    
    /* set radius */
    vshift = PLANET_SEALEVEL/(DEM_MAXHEIGHT - DEM_MAXDEPTH);
    for (i=0; i<res*NX; i++)
        for (j=0; j<res*NY; j++)
        {
//             wsphere[i*res*NY+j].indomain = (wsphere[i*res*NY+j].altitude < vshift + 1.0e-7);
            if (!wsphere[i*res*NY+j].indomain)
                wsphere[i*res*NY+j].radius_dem = 1.0 + RSCALE_DEM*(wsphere[i*res*NY+j].altitude - vshift);
            else wsphere[i*res*NY+j].radius_dem = 1.0;
//             printf("Radius_dem at (%i,%i) = %.3lg\n", i, j, wsphere[i*NY+j].radius_dem);
        }
}


void init_planet_map_rde(t_wave_sphere wsphere[NX*NY], int planet, int res)
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
    rx = (double)nx/(double)(res*NX);
    ry = (double)ny/(double)(res*NY - sshift - nshift);
    
    printf("cratio = %.3lg, rx = %.3lg, ry = %.3lg\n", cratio, rx, ry);
    
//     cy = rx*(double)(NY - nshift);
    
    /* build wave table */
    for (i=0; i<res*NX; i++)
        for (j=0; j<res*NY; j++)
        {
            ii = (int)(rx*(double)(res*NX-1 - i)) + nx/2 + ishift;
            if (ii > nx-1) ii -= nx;
//             jj = (int)(-ry*(double)j + cy);
//             jj = (int)(ry*(double)(NY-nshift - j)) + sshift;
            jj = (int)(ry*(double)(res*NY-nshift - j));
            if (jj > ny-1) jj = ny-1;
            if (jj < 0) jj = 0;
            wsphere[i*res*NY+j].r = (double)rgb_values[3*(jj*nx+ii)]*cratio;
            wsphere[i*res*NY+j].g = (double)rgb_values[3*(jj*nx+ii)+1]*cratio;
            wsphere[i*res*NY+j].b = (double)rgb_values[3*(jj*nx+ii)+2]*cratio;
            
//             printf("RGB at (%i, %i) = (%.3lg, %3.lg, %.3lg)\n", i, j, wsphere[i*NY+j].r, wsphere[i*NY+j].g, wsphere[i*NY+j].b);
            
            /* decide which points are in the Sea */
            wsphere[i*NY+j].indomain = 1;
            wsphere[i*NY+j].draw_wave = 1;
        }

    
    /* smooth colors in case of high resolution */
    if ((res*NX > nx)||(res*NY > ny))
        for (i=1; i<res*NX-1; i++)
            for (j=1; j<res*NY-1; j++)
            {
                wsphere[i*res*NY+j].r *= 0.2; 
                wsphere[i*res*NY+j].r += 0.2*wsphere[(i+1)*res*NY+j].r;
                wsphere[i*res*NY+j].r += 0.2*wsphere[(i-1)*res*NY+j].r;
                wsphere[i*res*NY+j].r += 0.2*wsphere[i*res*NY+j-1].r;
                wsphere[i*res*NY+j].r += 0.2*wsphere[i*res*NY+j+1].r;

                wsphere[i*res*NY+j].g *= 0.2; 
                wsphere[i*res*NY+j].g += 0.2*wsphere[(i+1)*res*NY+j].g;
                wsphere[i*res*NY+j].g += 0.2*wsphere[(i-1)*res*NY+j].g;
                wsphere[i*res*NY+j].g += 0.2*wsphere[i*res*NY+j-1].g;
                wsphere[i*res*NY+j].g += 0.2*wsphere[i*res*NY+j+1].g;

                wsphere[i*res*NY+j].b *= 0.2; 
                wsphere[i*res*NY+j].b += 0.2*wsphere[(i+1)*res*NY+j].b;
                wsphere[i*res*NY+j].b += 0.2*wsphere[(i-1)*res*NY+j].b;
                wsphere[i*res*NY+j].b += 0.2*wsphere[i*res*NY+j-1].b;
                wsphere[i*res*NY+j].b += 0.2*wsphere[i*res*NY+j+1].b;
            }

    
    free(rgb_values);
    fclose(image_file);
    
    if (ADD_DEM) init_dem_rde(wsphere, dem_number, res);
    
    vshift = PLANET_SEALEVEL/(DEM_MAXHEIGHT - DEM_MAXDEPTH);
    for (i=0; i<res*NX; i++)
        for (j=0; j<res*NY; j++)
            wsphere[i*res*NY+j].indomain = (wsphere[i*res*NY+j].altitude < vshift);
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

int ij_to_sphere_hres(int i, int j, double r, t_wave_sphere wsphere[HRES*HRES*NX*NY], double xyz[3], int use_wave_radius)
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
    
    xyz[0] = wsphere[i*HRES*NY+j].x;
    xyz[1] = wsphere[i*HRES*NY+j].y;
    xyz[2] = wsphere[i*HRES*NY+j].z;
    
    pscal = xyz[0]*observer[0] + xyz[1]*observer[1] + xyz[2]*observer[2];
    
    if (use_wave_radius)
    {
        newr = wsphere[i*HRES*NY+j].radius;
        xyz[0] *= newr;
        xyz[1] *= newr;
        xyz[2] *= newr;
    }
    else
    {
        newr = wsphere[i*HRES*NY+j].radius_dem;
        xyz[0] *= newr;
        xyz[1] *= newr;
        xyz[2] *= newr;
    }
    
    return(pscal/norm_observer > COS_VISIBLE);
}

void xyz_to_xy(double x, double y, double z, double xy_out[2])
{
    int i;
    double s, t, xinter[3];
    static double n2, m2, d, sm2, sn2, v[3], h[2], plane_ratio = 0.5;
    static int first = 1;
    
    if (((first)&&(REPRESENTATION_3D == REP_PROJ_3D))||(reset_view))
    {
        m2 = observer[0]*observer[0] + observer[1]*observer[1];
        n2 = m2 + observer[2]*observer[2];
        d = plane_ratio*n2;
        sm2 = sqrt(m2);
        sn2 = sqrt(n2);
        h[0] = observer[1]/sm2;
        h[1] = -observer[0]/sm2;
        v[0] = -observer[0]*observer[2]/(sn2*sm2);
        v[1] = -observer[1]*observer[2]/(sn2*sm2);
        v[2] = m2/(sn2*sm2);
        first = 0;
        reset_view = 0;
//         printf("h = (%.3lg, %.3lg)\n", h[0], h[1]);
//         printf("v = (%.3lg, %.3lg, %.3lg)\n", v[0], v[1], v[2]);
    }
    
    switch (REPRESENTATION_3D) {
        case (REP_AXO_3D):
        {
            for (i=0; i<2; i++)
                xy_out[i] = x*u_3d[i] + y*v_3d[i] + z*w_3d[i];
            break;
        }
        case (REP_PROJ_3D):
        {
            if (z > ZMAX_FACTOR*n2) z = ZMAX_FACTOR*n2;
            z *= Z_SCALING_FACTOR;
            s = observer[0]*x + observer[1]*y + observer[2]*z;
            t = (d - s)/(n2 - s);
            xinter[0] = t*observer[0] + (1.0-t)*x;
            xinter[1] = t*observer[1] + (1.0-t)*y;
            xinter[2] = t*observer[2] + (1.0-t)*z;
            
            xy_out[0] = XSHIFT_3D + XY_SCALING_FACTOR*(xinter[0]*h[0] + xinter[1]*h[1]);
            xy_out[1] = YSHIFT_3D + XY_SCALING_FACTOR*(xinter[0]*v[0] + xinter[1]*v[1] + xinter[2]*v[2]);
            break;
        }
    }
}

void draw_vertex_in_spherical_coords(double x, double y, double r, int i)
{
    double phi, theta, x1, y1, z1, xy_screen[2];
    static double phi_ratio, theta_ratio, phi_offset, theta_offset;
    static int first = 1;
    
    if (first)
    {
        phi_ratio = DPI/(XMAX - XMIN);
        theta_ratio = PI/(YMAX - YMIN); 
        phi_offset = phi_ratio*XMIN;
        theta_offset = theta_ratio*YMIN;
        first = 0;
    }
    
    phi = phi_ratio*x - phi_offset;
    theta = theta_ratio*y - theta_offset;
//     phi = DPI*(x - XMIN)/(XMAX - XMIN);
//     theta = PI*(y - YMIN)/(YMAX - YMIN); 
    x1 = r*cos(phi)*sin(theta);
    y1 = r*sin(phi)*sin(theta);
    z1 = -r*cos(theta);
    
//     printf("(phi, theta) = (%.5lg, %.5lg)\n", phi, theta);
//     printf("(x1, y1, z1) = (%.5lg, %.5lg, %.5lg)\n", x1, y1, z1);
                
    xyz_to_xy(x1, y1, z1, xy_screen);
    glVertex2d(xy_screen[0], xy_screen[1]);
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
//         if (B_DOMAIN == D_SPHERE_EARTH) init_earth_map(wsphere);
//         else if (OTHER_PLANET) init_planet_map(wsphere, B_DOMAIN);
//         else 
        if ((B_DOMAIN == D_SPHERE_JULIA)||(B_DOMAIN == D_SPHERE_JULIA_INV)||(B_DOMAIN == D_SPHERE_JULIA_CUBIC))
        {
            cos_rot = cos(JULIA_ROT*DPI/360.0);
            sin_rot = sin(JULIA_ROT*DPI/360.0);
        }
//         else if (B_DOMAIN == D_SPHERE_MAZE) init_sphere_maze(wsphere, 0, 0, 0);
//         else if (B_DOMAIN == D_SPHERE_MAZE_SPIRAL) init_sphere_maze(wsphere, 1, 0, 1);
//         else if (B_DOMAIN == D_SPHERE_MAZE_WAVE) init_sphere_maze(wsphere, 0, 1, 1);
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

void draw_vertex_sphere(double xyz[3])
{
    double xy_screen[2];
    
    xyz_to_xy(xyz[0], xyz[1], xyz[2], xy_screen);
    glVertex2d(xy_screen[0], xy_screen[1]);
}



void draw_vertex_ij(int i, int j)
{
    double xy[2];
    
    ij_to_xy(i, j, xy);
//     if (xy[1] > 0.0) printf("(i,j) = (%i,%i), (x,y) = (%.2lg,%.2lg)\n", i, j, xy[0], xy[1]);
    glVertex2d(xy[0], xy[1]);
}


void draw_vertex_xyz(double xy[2], double z)
{
    double xy_screen[2];
    
    xyz_to_xy(xy[0], xy[1], z, xy_screen);
    glVertex2d(xy_screen[0], xy_screen[1]);
}
    
void draw_vertex_xyz_shift(double xy[2], double z, int shiftx, int shifty)
{
    double xy_screen[2];
    
    xyz_to_xy(xy[0] + (double)shiftx*(XMAX - XMIN), xy[1] + (double)shifty*(YMAX - YMIN), z, xy_screen);
    glVertex2d(xy_screen[0], xy_screen[1]);
}
    
void draw_vertex_x_y_z(double x, double y, double z)
{
    double xy_screen[2];
    
    xyz_to_xy(x, y, z, xy_screen);
    glVertex2d(xy_screen[0], xy_screen[1]);
}

void draw_rectangle_3d(double x1, double y1, double x2, double y2)
{
    glBegin(GL_LINE_LOOP);
    draw_vertex_x_y_z(x1, y1, 0.0);
    draw_vertex_x_y_z(x2, y1, 0.0);
    draw_vertex_x_y_z(x2, y2, 0.0);
    draw_vertex_x_y_z(x1, y2, 0.0);    
    glEnd();    
}


void draw_rectangle_noscale(double x1, double y1, double x2, double y2)
{
    glBegin(GL_LINE_LOOP);
    glVertex2d(x1, y1);
    glVertex2d(x2, y1);
    glVertex2d(x2, y2);
    glVertex2d(x1, y2);
    glEnd();    
}

void draw_circle_3d(double x, double y, double r, int nseg)
{
    int i;
    double pos[2], alpha, dalpha;
    
    dalpha = DPI/(double)nseg;
    
    glBegin(GL_LINE_LOOP);
    for (i=0; i<=nseg; i++)
    {
        alpha = (double)i*dalpha;
        draw_vertex_x_y_z(x + r*cos(alpha), y + r*sin(alpha), 0.0);
    }
    glEnd();
}

void tvertex_lineto_3d(t_vertex z)
/* draws boundary segments of isospectral billiard */
{
    draw_vertex_x_y_z(z.x, z.y, 0.0);
}


void draw_billiard_3d(int fade, double fade_value)      /* draws the billiard boundary */
{
    double x0, x, y, x1, y1, dx, dy, phi, r = 0.01, pos[2], pos1[2], alpha, dphi, omega, z, l, width, a, b, c, ymax;
    int i, j, k, k1, k2, mr2;
    static int first = 1, nsides;

    if (BLACK) 
    {
        if (fade) glColor3f(fade_value, fade_value, fade_value);
        else glColor3f(1.0, 1.0, 1.0);
    }
    else 
    {
        if (fade) glColor3f(1.0 - fade_value, 1.0 - fade_value, 1.0 - fade_value);
        else glColor3f(0.0, 0.0, 0.0);
    }
    glLineWidth(BOUNDARY_WIDTH);

    glEnable(GL_LINE_SMOOTH);

    switch (B_DOMAIN) {
        case (D_RECTANGLE):
        {
            glBegin(GL_LINE_LOOP);
            draw_vertex_x_y_z(LAMBDA, -1.0, 0.0);
            draw_vertex_x_y_z(LAMBDA, 1.0, 0.0);
            draw_vertex_x_y_z(-LAMBDA, 1.0, 0.0);
            draw_vertex_x_y_z(-LAMBDA, -1.0, 0.0);
            glEnd();
            break;
        }
        case (D_ELLIPSE):
        {
            glBegin(GL_LINE_LOOP);
            for (i=0; i<=NSEG; i++)
            {
                phi = (double)i*DPI/(double)NSEG;
                x = 1.05*LAMBDA*cos(phi);
                y = 1.05*sin(phi);
                draw_vertex_x_y_z(x, y, 0.0);
            }
            glEnd ();

            /* draw foci */
//             if (FOCI)
//             {
//                 glColor3f(0.3, 0.3, 0.3);
//                 x0 = sqrt(LAMBDA*LAMBDA-1.0);
// 
//                 glLineWidth(2);
//                 glEnable(GL_LINE_SMOOTH);
//                 
//                 draw_circle(x0, 0.0, r, NSEG);
//                 draw_circle(-x0, 0.0, r, NSEG);
//             }
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
                draw_vertex_x_y_z(x, y, 0.0);
            }
            for (i=0; i<=NSEG; i++)
            {
                phi = PID + (double)i*PI/(double)NSEG;
                x = -0.5*LAMBDA + cos(phi);
                y = sin(phi);
                draw_vertex_x_y_z(x, y, 0.0);
            }
            glEnd();
            break;
        }
        case (D_SINAI):
        {
            draw_circle_3d(0.0, 0.0, LAMBDA, NSEG);
            break;
        }
        case (D_POLYGON):
        {
            omega = DPI/((double)NPOLY);
            glBegin(GL_LINE_LOOP);
            for (i=0; i<=NPOLY; i++)
            {
                x = 1.0075*cos(i*omega + APOLY*PID);
                y = 1.0075*sin(i*omega + APOLY*PID);
                draw_vertex_x_y_z(x, y, 0.0);
            }
            glEnd ();
            break;
        }
        case (D_YOUNG):
        {
            if (FILL_BILLIARD_COMPLEMENT)
            {
                if (fade) glColor3f(0.75*fade_value, 0.75*fade_value, 0.75*fade_value);
                else glColor3f(0.75, 0.75, 0.75);
                
                glBegin(GL_TRIANGLE_FAN);
                draw_vertex_x_y_z(-MU, YMIN, 0.0);
                draw_vertex_x_y_z(-MU, -LAMBDA-MU, 0.0);
                draw_vertex_x_y_z(MU, -LAMBDA-MU, 0.0);
                draw_vertex_x_y_z(MU, YMIN, 0.0);
                glEnd();
            
                glBegin(GL_TRIANGLE_FAN);
                draw_vertex_x_y_z(-MU, YMAX, 0.0);
                draw_vertex_x_y_z(-MU, LAMBDA+MU, 0.0);
                draw_vertex_x_y_z(MU, LAMBDA+MU, 0.0);
                draw_vertex_x_y_z(MU, YMAX, 0.0);
                glEnd();

                glBegin(GL_TRIANGLE_FAN);
                draw_vertex_x_y_z(-MU, -LAMBDA+MU, 0.0);
                draw_vertex_x_y_z(-MU, LAMBDA-MU, 0.0);
                draw_vertex_x_y_z(MU, LAMBDA-MU, 0.0);
                draw_vertex_x_y_z(MU, -LAMBDA+MU, 0.0);
                glEnd();
            }
            else
            {
                glBegin(GL_LINE_STRIP);
                draw_vertex_x_y_z(-MU, YMIN, 0.0);
                draw_vertex_x_y_z(-MU, -LAMBDA-MU, 0.0);
                draw_vertex_x_y_z(MU, -LAMBDA-MU, 0.0);
                draw_vertex_x_y_z(MU, YMIN, 0.0);
                glEnd();
            
                glBegin(GL_LINE_STRIP);
                draw_vertex_x_y_z(-MU, YMAX, 0.0);
                draw_vertex_x_y_z(-MU, LAMBDA+MU, 0.0);
                draw_vertex_x_y_z(MU, LAMBDA+MU, 0.0);
                draw_vertex_x_y_z(MU, YMAX, 0.0);
                glEnd();

                glBegin(GL_LINE_LOOP);
                draw_vertex_x_y_z(-MU, -LAMBDA+MU, 0.0);
                draw_vertex_x_y_z(-MU, LAMBDA-MU, 0.0);
                draw_vertex_x_y_z(MU, LAMBDA-MU, 0.0);
                draw_vertex_x_y_z(MU, -LAMBDA+MU, 0.0);
                glEnd();
            }
        }
        case (D_GRATING):
        {
            k1 = -(int)(-YMIN/LAMBDA);
            k2 = (int)(YMAX/LAMBDA);
            for (i=k1; i<= k2; i++)
            {
                z = (double)i*LAMBDA;
                draw_circle_3d(0.0, z, MU, NSEG);
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
                x = 1.0 + (LAMBDA + 0.01)*cos(phi);
                y = (LAMBDA + 0.01)*sin(phi);
                draw_vertex_x_y_z(x, y, 0.0);
            }
            phi = PI - alpha;
            x = 1.0 + (LAMBDA + 0.01)*cos(phi);
            y = 0.01 + (LAMBDA + 0.01)*sin(phi);
            draw_vertex_x_y_z(x, y, 0.0);
            phi = alpha;
            x = -1.0 + (LAMBDA + 0.01)*cos(phi);
            y = 0.01 + (LAMBDA + 0.01)*sin(phi);
            draw_vertex_x_y_z(x, y, 0.0);
            for (i=0; i<=NSEG; i++)
            {
                phi = alpha + (double)i*dphi;
                x = -1.0 + (LAMBDA + 0.01)*cos(phi);
                y = (LAMBDA + 0.01)*sin(phi);
                draw_vertex_x_y_z(x, y, 0.0);
            }
            glEnd ();
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
                draw_vertex_x_y_z(x, y, 0.0);
            }
            for (i = 0; i < NSEG+1; i++) 
            {
                y = 1.5*MU - dy*(double)i;
                x = 0.25*y*y/MU - (MU + width) - LAMBDA;
                draw_vertex_x_y_z(x, y, 0.0);
            }
            glEnd ();
            
            glBegin(GL_LINE_LOOP);
            for (i = 0; i < NSEG+1; i++) 
            {
                y = -1.5*MU + dy*(double)i;
                x = LAMBDA + MU - 0.25*y*y/MU;
                draw_vertex_x_y_z(x, y, 0.0);
            }
            for (i = 0; i < NSEG+1; i++) 
            {
                y = 1.5*MU - dy*(double)i;
                x = LAMBDA + (MU + width) - 0.25*y*y/MU;
                draw_vertex_x_y_z(x, y, 0.0);
            }
            glEnd ();
            
//             if (FOCI)
//             {
//                 glColor3f(0.3, 0.3, 0.3);
//                 draw_circle(-LAMBDA, 0.0, r, NSEG);
//                 draw_circle(LAMBDA, 0.0, r, NSEG);
//             }

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
                    draw_vertex_x_y_z(x, y, 0.0);
                }
            }
            glEnd ();
            
//             if (FOCI)
//             {
//                 glColor3f(0.3, 0.3, 0.3);
//                 for (k=0; k<NPOLY; k++) 
//                 {
//                     alpha = APOLY*PID + (2.0*(double)k+1.0)*omega;
//                     draw_circle(LAMBDA*cos(alpha), LAMBDA*sin(alpha), r, NSEG);
//                 }
//             }
            
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
                draw_rectangle_3d(x, x, -x, -x);
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
                            draw_rectangle_3d(x, y, x+l, y+l);
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
                            draw_rectangle_3d(x, y, x+l, y+l);
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
                draw_vertex_x_y_z(x, y, 0.0);
            }
            
            /* straight parts */
            draw_vertex_x_y_z(-LAMBDA, width, 0.0);
            draw_vertex_x_y_z(-c, width, 0.0);
            draw_vertex_x_y_z(-c, MU, 0.0);
                        
            /* left half ellipse */
            for (i=0; i<=NSEG; i++)
            {
                phi = (double)i*dphi;
                x = -c + 0.5*MU*sin(phi);
                y = MU*cos(phi);
                draw_vertex_x_y_z(x, y, 0.0);
            }
            
            /* straight parts */
            draw_vertex_x_y_z(-c, -width, 0.0);
            draw_vertex_x_y_z(-LAMBDA, -width, 0.0);
            draw_vertex_x_y_z(-LAMBDA, -MU, 0.0);
            
            /* lower half ellipse */
            for (i=0; i<=NSEG; i++)
            {
                phi = (double)i*dphi;
                x = -LAMBDA*cos(phi);
                y = -MU - (1.0-MU)*sin(phi);
                draw_vertex_x_y_z(x, y, 0.0);
            }
            
            /* straight parts */
            draw_vertex_x_y_z(LAMBDA, -width, 0.0);
            draw_vertex_x_y_z(c, -width, 0.0);
            draw_vertex_x_y_z(c, -MU, 0.0);
            
            /* right half ellipse */
            for (i=0; i<=NSEG; i++)
            {
                phi = (double)i*dphi;
                x = c - 0.5*MU*sin(phi);
                y = -MU*cos(phi);
                draw_vertex_x_y_z(x, y, 0.0);
            }
            
            /* straight parts */
            draw_vertex_x_y_z(c, width, 0.0);
            draw_vertex_x_y_z(LAMBDA, width, 0.0);
            draw_vertex_x_y_z(LAMBDA, MU, 0.0);            
            glEnd ();
            break; 
        }
        case (D_TOKA_PRIME):
        {
            glBegin(GL_LINE_LOOP);
            tvertex_lineto_3d(polyline[0]);
            for (i=4; i<43; i++) tvertex_lineto_3d(polyline[i]);
            tvertex_lineto_3d(polyline[3]);
            tvertex_lineto_3d(polyline[2]);
            tvertex_lineto_3d(polyline[1]);

            tvertex_lineto_3d(polyline[44]);
            tvertex_lineto_3d(polyline[45]);
            for (i=84; i>45; i--) tvertex_lineto_3d(polyline[i]);
            glEnd();
            
            /* inner lines */ 
//             glLineWidth(BOUNDARY_WIDTH/2);
            glLineWidth(1);
            glColor3f(0.75, 0.75, 0.75);
            glBegin(GL_LINE_STRIP);
            tvertex_lineto_3d(polyline[0]);
            tvertex_lineto_3d(polyline[1]);
            tvertex_lineto_3d(polyline[2]);
            tvertex_lineto_3d(polyline[0]);
            tvertex_lineto_3d(polyline[3]);
            tvertex_lineto_3d(polyline[4]);
            glEnd();
            
            glBegin(GL_LINE_STRIP);
            tvertex_lineto_3d(polyline[0]);
            tvertex_lineto_3d(polyline[44]);
            tvertex_lineto_3d(polyline[45]);
            tvertex_lineto_3d(polyline[0]);
            tvertex_lineto_3d(polyline[46]);
            tvertex_lineto_3d(polyline[45]);
            glEnd();
            
            for (i=3; i<43; i++)
            {
                glBegin(GL_LINE_STRIP);
                tvertex_lineto_3d(polyline[i]);
                tvertex_lineto_3d(polyline[43]);
                glEnd();
                glBegin(GL_LINE_STRIP);
                tvertex_lineto_3d(polyline[i+42]);
                tvertex_lineto_3d(polyline[85]);
                glEnd();
            }
            
            break;
        }
        case (D_ISOSPECTRAL):
        {
            /* 1st triangle */
            glBegin(GL_LINE_LOOP);
            tvertex_lineto_3d(polyline[0]);
            tvertex_lineto_3d(polyline[4]);
            tvertex_lineto_3d(polyline[7]);
            tvertex_lineto_3d(polyline[1]);
            tvertex_lineto_3d(polyline[5]);
            tvertex_lineto_3d(polyline[8]);
            tvertex_lineto_3d(polyline[2]);
            tvertex_lineto_3d(polyline[3]);
            tvertex_lineto_3d(polyline[6]);
            glEnd();
            
            /* inner lines */
            glBegin(GL_LINE_LOOP);
            tvertex_lineto_3d(polyline[0]);
            tvertex_lineto_3d(polyline[1]);
            tvertex_lineto_3d(polyline[2]);
            tvertex_lineto_3d(polyline[0]);
            tvertex_lineto_3d(polyline[3]);
            tvertex_lineto_3d(polyline[2]);
            tvertex_lineto_3d(polyline[5]);
            tvertex_lineto_3d(polyline[1]);
            tvertex_lineto_3d(polyline[4]);
            glEnd();
            
            /* 2nd triangle */
            glBegin(GL_LINE_LOOP);
            tvertex_lineto_3d( polyline[9]);
            tvertex_lineto_3d(polyline[16]);
            tvertex_lineto_3d(polyline[13]);
            tvertex_lineto_3d(polyline[10]);
            tvertex_lineto_3d(polyline[17]);
            tvertex_lineto_3d(polyline[14]);
            tvertex_lineto_3d(polyline[11]);
            tvertex_lineto_3d(polyline[15]);
            tvertex_lineto_3d(polyline[12]);
            glEnd();
            
            /* inner lines */
            glBegin(GL_LINE_LOOP);
            tvertex_lineto_3d( polyline[9]);
            tvertex_lineto_3d(polyline[10]);
            tvertex_lineto_3d(polyline[11]);
            tvertex_lineto_3d( polyline[9]);
            tvertex_lineto_3d(polyline[13]);
            tvertex_lineto_3d(polyline[10]);
            tvertex_lineto_3d(polyline[14]);
            tvertex_lineto_3d(polyline[11]);
            tvertex_lineto_3d(polyline[12]);
            glEnd();
            break;
        }
        case (D_HOMOPHONIC):
        {
            /* 1st triangle */
            glBegin(GL_LINE_LOOP);
            tvertex_lineto_3d(polyline[1]);
            tvertex_lineto_3d(polyline[3]);
            tvertex_lineto_3d(polyline[4]);
            tvertex_lineto_3d(polyline[5]);
            tvertex_lineto_3d(polyline[6]);
            tvertex_lineto_3d(polyline[8]);
            tvertex_lineto_3d(polyline[9]);
            tvertex_lineto_3d(polyline[10]);
            tvertex_lineto_3d(polyline[12]);
            tvertex_lineto_3d(polyline[13]);
            tvertex_lineto_3d(polyline[15]);
            tvertex_lineto_3d(polyline[16]);
            tvertex_lineto_3d(polyline[17]);
            tvertex_lineto_3d(polyline[18]);
            tvertex_lineto_3d(polyline[20]);
            glEnd();
            
            /* inner lines */
            glLineWidth(BOUNDARY_WIDTH/2);
            glBegin(GL_LINE_STRIP);
            tvertex_lineto_3d(polyline[9]);
            tvertex_lineto_3d(polyline[1]);
            tvertex_lineto_3d(polyline[2]);
            tvertex_lineto_3d(polyline[5]);
            tvertex_lineto_3d(polyline[7]);
            tvertex_lineto_3d(polyline[2]);
            tvertex_lineto_3d(polyline[8]);
            tvertex_lineto_3d(polyline[21]);
            tvertex_lineto_3d(polyline[10]);
            tvertex_lineto_3d(polyline[2]);
            tvertex_lineto_3d(polyline[21]);
            tvertex_lineto_3d(polyline[11]);
            tvertex_lineto_3d(polyline[13]);
            tvertex_lineto_3d(polyline[21]);
            tvertex_lineto_3d(polyline[14]);
            tvertex_lineto_3d(polyline[20]);
            tvertex_lineto_3d(polyline[15]);
            tvertex_lineto_3d(polyline[19]);
            tvertex_lineto_3d(polyline[16]);
            tvertex_lineto_3d(polyline[18]);
            glEnd();
            
            /* 2nd triangle */
            glLineWidth(BOUNDARY_WIDTH);
            glBegin(GL_LINE_LOOP);
            tvertex_lineto_3d(polyline[22+10]);
            tvertex_lineto_3d(polyline[22+16]);
            tvertex_lineto_3d(polyline[22+17]);
            tvertex_lineto_3d(polyline[22+18]);
            tvertex_lineto_3d(polyline[22+12]);
            tvertex_lineto_3d(polyline[22+13]);
            tvertex_lineto_3d(polyline[22+15]);
            tvertex_lineto_3d(polyline[22+19]);
            tvertex_lineto_3d(polyline[22+20]);
            tvertex_lineto_3d(polyline[22+1]);
            tvertex_lineto_3d(polyline[22+4]);
            tvertex_lineto_3d(polyline[22+5]);
            tvertex_lineto_3d(polyline[22+7]);
            tvertex_lineto_3d(polyline[22+8]);
            tvertex_lineto_3d(polyline[22+9]);
            glEnd();
            
            /* inner lines */
            glLineWidth(BOUNDARY_WIDTH/2);
            glBegin(GL_LINE_STRIP);
            tvertex_lineto_3d(polyline[22+2]);
            tvertex_lineto_3d(polyline[22+6]);
            tvertex_lineto_3d(polyline[22+8]);
            tvertex_lineto_3d(polyline[22+2]);
            tvertex_lineto_3d(polyline[22+5]);
            tvertex_lineto_3d(polyline[22+3]);
            tvertex_lineto_3d(polyline[22+2]);
            tvertex_lineto_3d(polyline[22+1]);
            tvertex_lineto_3d(polyline[22+0]);
            tvertex_lineto_3d(polyline[22+21]);
            tvertex_lineto_3d(polyline[22+18]);
            tvertex_lineto_3d(polyline[22+16]);
            tvertex_lineto_3d(polyline[22+13]);
            tvertex_lineto_3d(polyline[22+21]);
            tvertex_lineto_3d(polyline[22+10]);
            tvertex_lineto_3d(polyline[22+12]);
            tvertex_lineto_3d(polyline[22+21]);
            tvertex_lineto_3d(polyline[22+14]);
            tvertex_lineto_3d(polyline[22+20]);
            tvertex_lineto_3d(polyline[22+15]);
            glEnd();
            break;
        }
        default:
        {
            break;
        }   
    }
}

void draw_polyline_visible(int j, int k, double margin)
/* hack to draw the billiard boundary in front of the wave */
/* only parts of the boundary having a small enough angle with the observer vector are drawn */
{
    double x, y, x1, y1, length, length1;
    static int first = 1;
    static double olength;
    
    if (first)
    {
        olength = module2(observer[0], observer[1]);
        first = 0;
    }
    
    x = polyline[j].x;
    y = polyline[j].y;
    x1 = polyline[k].x;
    y1 = polyline[k].y;
    length = module2(x,y);
    length1 = module2(x1,y1);
    if ((x*observer[0] + y*observer[1] > margin*length*olength)&&(x1*observer[0] + y1*observer[1] > margin*length1*olength))
    {
        glBegin(GL_LINE_STRIP);
        tvertex_lineto_3d(polyline[j]);
        tvertex_lineto_3d(polyline[k]);
        glEnd();
    }
}

void draw_vertex_visible(double x, double y, double margin)
/* hack to draw the billiard boundary in front of the wave */
/* only parts of the boundary having a small enough angle with the observer vector are drawn */
{
    static int first = 1;
    static double olength;
    
    if (first)
    {
        olength = module2(observer[0], observer[1]);
        first = 0;
    }
    
    if (x*observer[0] + y*observer[1] > margin*olength*module2(x,y))
        draw_vertex_x_y_z(x, y, 0.0);
    else 
    {
        glEnd();
        glBegin(GL_LINE_STRIP);
    }
}

void draw_billiard_3d_front(int fade, double fade_value)      
/* hack to draw the billiard boundary in front of the wave */
/* only parts of the boundary having a small enough angle with the observer vector are drawn */
{
    double x0, x, y, x1, y1, dx, dy, phi, r = 0.01, pos[2], pos1[2], alpha, dphi, omega, z, l, width, a, b, c, ymax, length, length1;
    int i, j, k, k1, k2, mr2;
    static int first = 1, nsides;
    static double olength;
    
    if (first)
    {
        olength = module2(observer[0], observer[1]);
        first = 0;
    }

    if (BLACK) 
    {
        if (fade) glColor3f(fade_value, fade_value, fade_value);
        else glColor3f(1.0, 1.0, 1.0);
    }
    else 
    {
        if (fade) glColor3f(1.0 - fade_value, 1.0 - fade_value, 1.0 - fade_value);
        else glColor3f(0.0, 0.0, 0.0);
    }
    glLineWidth(BOUNDARY_WIDTH);

    glEnable(GL_LINE_SMOOTH);

    switch (B_DOMAIN) {
        case (D_RECTANGLE):
        {
            glBegin(GL_LINE_STRIP);
            draw_vertex_x_y_z(LAMBDA, -1.0, 0.0);
            draw_vertex_x_y_z(LAMBDA, 1.0, 0.0);
            draw_vertex_x_y_z(-LAMBDA, 1.0, 0.0);
            glEnd();
            break;
        }
        case (D_STADIUM):
        {
            glBegin(GL_LINE_STRIP);
            for (i=0; i<=NSEG; i++)
            {
                phi = -PID + (double)i*PI/(double)NSEG;
                x = 0.5*LAMBDA + cos(phi);
                y = sin(phi);
                draw_vertex_visible(x, y, 0.0);
            }
            for (i=0; i<=NSEG; i++)
            {
                phi = PID + (double)i*PI/(double)NSEG;
                x = -0.5*LAMBDA + cos(phi);
                y = sin(phi);
                draw_vertex_visible(x, y, 0.0);
            }
            glEnd();
            break;
        }
        case (D_POLYGON):
        {
            omega = DPI/((double)NPOLY);
            glBegin(GL_LINE_STRIP);
            for (i=0; i<=NPOLY; i++)
            {
                x = cos(i*omega + APOLY*PID);
                y = sin(i*omega + APOLY*PID);
                draw_vertex_visible(x, y, 0.0);
            }
            glEnd ();
            break;
        }
        case (D_YOUNG):
        {
            glBegin(GL_LINE_STRIP);
            draw_vertex_x_y_z(-MU, YMIN, 0.0);
            draw_vertex_x_y_z(-MU, -LAMBDA-MU, 0.0);
            glEnd();
            
            glBegin(GL_LINE_STRIP);
            draw_vertex_x_y_z(-MU, YMAX, 0.0);
            draw_vertex_x_y_z(-MU, LAMBDA+MU, 0.0);
            glEnd();

            glBegin(GL_LINE_LOOP);
            draw_vertex_x_y_z(-MU, -LAMBDA+MU, 0.0);
            draw_vertex_x_y_z(-MU, LAMBDA-MU, 0.0);
            glEnd();
            break;
        }
        case (D_TOKA_PRIME):
        {
            draw_polyline_visible(0, 4, 0.2);
            for (i=4; i<41; i++) draw_polyline_visible(i, i+1, -0.1);
//             draw_polyline_visible(42, 3, 0.2);
            for (i=84; i>46; i--) draw_polyline_visible(i, i-1, -0.1);
            draw_polyline_visible(46, 0, 0.2);
            break;
        }
        case (D_STAR):
        {
            glLineWidth(BOUNDARY_WIDTH);
            for (i=0; i<npolyline; i++)
                draw_polyline_visible(i%npolyline, (i+1)%npolyline, 0.2);
            break;
        }
        default:
        {
            break;
        }   
    }
}

void compute_energy_field(double phi[NX*NY], double psi[NX*NY], short int xy_in[NX*NY], double energies[NX*NY])
/* computes cosine of angle between normal vector and vector light */
{
    int i, j;
    
    #pragma omp parallel for private(i,j)
    for (i=0; i<NX; i++)
        for (j=0; j<NY; j++)
        {
            if ((TWOSPEEDS)||(xy_in[i*NY+j])) energies[i*NY+j] = PLOT_SCALE_ENERGY*compute_energy_mod(phi, psi, xy_in, i, j);
            else energies[i*NY+j] = 0.0;
        }
}


void compute_phase_field(double phi[NX*NY], double psi[NX*NY], short int xy_in[NX*NY], double phase[NX*NY])
/* computes cosine of angle between normal vector and vector light */
{
    int i, j;
    double angle;
    
    #pragma omp parallel for private(i,j,angle)
    for (i=0; i<NX; i++)
        for (j=0; j<NY; j++)
        {
            if ((TWOSPEEDS)||(xy_in[i*NY+j])) 
            {
                angle = compute_phase(phi, psi, xy_in, i, j)/DPI + PHASE_SHIFT;
                if (angle >= 1.0) angle -= 1.0;
                else if (angle < 0.0) angle += 1.0;
                phase[i*NY+j] = angle;
            }
            else phase[i*NY+j] = 0.0;
        }
}


// void compute_log_energy_field(double *phi[NX], double *psi[NX], short int xy_in[NX*NY], double *energies[NX])
// /* computes cosine of angle between normal vector and vector light */
// {
//     int i, j;
//     
//     for (i=0; i<NX; i++)
//         for (j=0; j<NY; j++)
//         {
//             if ((TWOSPEEDS)||(xy_in[i*NY+j])) energies[i][j] = PLOT_SCALE_ENERGY*compute_energy(phi, psi, xy_in, i, j);
//             else energies[i][j] = 0.0;
//         }
// }


void compute_light_angle(double field[NX*NY], short int xy_in[NX*NY], double cos_angle[NX*NY])
/* computes cosine of angle between normal vector and vector light */
{
    int i, j;
    double gradx, grady, norm, pscal;
    static double dx, dy;
    static int first = 1;
    
    if (first)
    {
        dx = 2.0*(XMAX - XMIN)/(double)NX;
        dy = 2.0*(YMAX - YMIN)/(double)NY;
        first = 0;
    }
    
    #pragma omp parallel for private(i,j,gradx, grady, norm, pscal)
    for (i=1; i<NX-2; i++)
        for (j=1; j<NY-2; j++)
        {
            if ((TWOSPEEDS)||(xy_in[i*NY+j]))
            {
                gradx = (field[(i+1)*NY+j] - field[(i-1)*NY+j])/dx;
                grady = (field[i*NY+j+1] - field[i*NY+j-1])/dy;
                norm = sqrt(1.0 + gradx*gradx + grady*grady);
                pscal = -gradx*light[0] - grady*light[1] + 1.0;
                
                cos_angle[i*NY+j] = pscal/norm;
            }
        }
}


void energy_color_scheme(int palette, double energy, double rgb[3])
{
    if (COLOR_PALETTE >= COL_TURBO) color_scheme_asym_palette(COLOR_SCHEME, palette, energy, 1.0, 0, rgb);
    else color_scheme_palette(COLOR_SCHEME, palette, energy, 1.0, 0, rgb);
}

void log_energy_color_scheme(int palette, double energy, double rgb[3])
{
    color_scheme_palette(COLOR_SCHEME, palette, LOG_SHIFT + LOG_SCALE*log(energy), 1.0, 0, rgb);
//     if (energy > 1.0e-8) printf("energy = %.3lg, log energy = %.3lg\n", energy, LOG_SHIFT + LOG_SCALE*log(energy));
}

void phase_color_scheme(int palette, double phase, double rgb[3])
{
//     color_scheme_palette(C_ONEDIM_LINEAR, palette, phase/DPI, 1.0, 0, rgb);
    amp_to_rgb_palette(phase, rgb, palette);
}


void compute_interpolated_colors(int i, int j, double field[NX*NY], double palette, int plot, 
                                 double *z_sw, double *z_se, double *z_nw, double *z_ne, double *z_mid, 
                                 double rgb_e[3], double rgb_w[3], double rgb_n[3], double rgb_s[3])
{
    double zw, ze, zn, zs;
    
    *z_sw = field[i*NY+j];
    *z_se = field[(i+1)*NY+j];
    *z_nw = field[i*NY+j+1];
    *z_ne = field[(i+1)*NY+j+1];
                            
    *z_mid = 0.25*(*z_sw + *z_se + *z_nw + *z_ne);
                            
    zw = (*z_sw + *z_nw + *z_mid)/3.0;
    ze = (*z_se + *z_ne + *z_mid)/3.0;
    zs = (*z_sw + *z_se + *z_mid)/3.0;
    zn = (*z_nw + *z_ne + *z_mid)/3.0;
    
    if (plot == P_3D_ENERGY)
    {
        energy_color_scheme(palette, VSCALE_ENERGY*ze, rgb_e);
        energy_color_scheme(palette, VSCALE_ENERGY*zw, rgb_w);
        energy_color_scheme(palette, VSCALE_ENERGY*zn, rgb_n);
        energy_color_scheme(palette, VSCALE_ENERGY*zs, rgb_s);         
    }
    else if (plot == P_3D_LOG_ENERGY)
    {
        log_energy_color_scheme(palette, ze, rgb_e);
        log_energy_color_scheme(palette, zw, rgb_w);
        log_energy_color_scheme(palette, zn, rgb_n);
        log_energy_color_scheme(palette, zs, rgb_s);         
    }
    else
    {
        color_scheme_palette(COLOR_SCHEME, palette, VSCALE_AMPLITUDE*ze, 1.0, 0, rgb_e);
        color_scheme_palette(COLOR_SCHEME, palette, VSCALE_AMPLITUDE*zw, 1.0, 0, rgb_w);
        color_scheme_palette(COLOR_SCHEME, palette, VSCALE_AMPLITUDE*zn, 1.0, 0, rgb_n);
        color_scheme_palette(COLOR_SCHEME, palette, VSCALE_AMPLITUDE*zs, 1.0, 0, rgb_s); 
    }
}

void compute_interpolated_colors_new(int i, int j, double zfield[NX*NY], double colorfield[NX*NY], double palette, int cplot, 
                                 double *z_sw, double *z_se, double *z_nw, double *z_ne, double *z_mid, 
                                 double rgb_e[3], double rgb_w[3], double rgb_n[3], double rgb_s[3])
{
    double zw, ze, zn, zs, c_sw, c_se, c_nw, c_ne, c_mid;
    
    *z_sw = zfield[i*NY+j];
    *z_se = zfield[(i+1)*NY+j];
    *z_nw = zfield[i*NY+j+1];
    *z_ne = zfield[(i+1)*NY+j+1];
                            
    *z_mid = 0.25*(*z_sw + *z_se + *z_nw + *z_ne);
                            
    c_sw = colorfield[i*NY+j];
    c_se = colorfield[(i+1)*NY+j];
    c_nw = colorfield[i*NY+j+1];
    c_ne = colorfield[(i+1)*NY+j+1];
                            
    c_mid = 0.25*(c_sw + c_se + c_nw + c_ne);

    zw = (c_sw + c_nw + c_mid)/3.0;
    ze = (c_se + c_ne + c_mid)/3.0;
    zs = (c_sw + c_se + c_mid)/3.0;
    zn = (c_nw + c_ne + c_mid)/3.0;

    switch (cplot){
        case (P_3D_ENERGY):
        {
            energy_color_scheme(palette, VSCALE_ENERGY*ze, rgb_e);
            energy_color_scheme(palette, VSCALE_ENERGY*zw, rgb_w);
            energy_color_scheme(palette, VSCALE_ENERGY*zn, rgb_n);
            energy_color_scheme(palette, VSCALE_ENERGY*zs, rgb_s);
            break;
        }
        case (P_3D_LOG_ENERGY):
        {
            log_energy_color_scheme(palette, ze, rgb_e);
            log_energy_color_scheme(palette, zw, rgb_w);
            log_energy_color_scheme(palette, zn, rgb_n);
            log_energy_color_scheme(palette, zs, rgb_s);         
            break;
        }
        case (P_3D_PHASE):
        {
            phase_color_scheme(palette, ze, rgb_e);
            phase_color_scheme(palette, zw, rgb_w);
            phase_color_scheme(palette, zn, rgb_n);
            phase_color_scheme(palette, zs, rgb_s);         
            break;
        }
        default:
        {
//         hsl_to_rgb_palette(VSCALE_AMPLITUDE*ze, 0.9, 0.5, rgb_e, palette);
//         hsl_to_rgb_palette(VSCALE_AMPLITUDE*zw, 0.9, 0.5, rgb_w, palette);
//         hsl_to_rgb_palette(VSCALE_AMPLITUDE*zn, 0.9, 0.5, rgb_n, palette);
//         hsl_to_rgb_palette(VSCALE_AMPLITUDE*zs, 0.9, 0.5, rgb_s, palette);
            color_scheme_palette(COLOR_SCHEME, palette, VSCALE_AMPLITUDE*ze, 1.0, 0, rgb_e);
            color_scheme_palette(COLOR_SCHEME, palette, VSCALE_AMPLITUDE*zw, 1.0, 0, rgb_w);
            color_scheme_palette(COLOR_SCHEME, palette, VSCALE_AMPLITUDE*zn, 1.0, 0, rgb_n);
            color_scheme_palette(COLOR_SCHEME, palette, VSCALE_AMPLITUDE*zs, 1.0, 0, rgb_s); 
        }
    }
}

void compute_interpolated_colors_rde_temp(int i, int j, double zfield[NX*NY], double colorfield[NX*NY], double palette, int cplot, 
                                 double *z_sw, double *z_se, double *z_nw, double *z_ne, double *z_mid, 
                                 double rgb_e[3], double rgb_w[3], double rgb_n[3], double rgb_s[3])
{
    double zw, ze, zn, zs, c_sw, c_se, c_nw, c_ne, c_mid;
    
    *z_sw = zfield[i*NY+j];
    *z_se = zfield[(i+1)*NY+j];
    *z_nw = zfield[i*NY+j+1];
    *z_ne = zfield[(i+1)*NY+j+1];
                            
    *z_mid = 0.25*(*z_sw + *z_se + *z_nw + *z_ne);
                            
    c_sw = colorfield[i*NY+j];
    c_se = colorfield[(i+1)*NY+j];
    c_nw = colorfield[i*NY+j+1];
    c_ne = colorfield[(i+1)*NY+j+1];
                            
    c_mid = 0.25*(c_sw + c_se + c_nw + c_ne);

    zw = (c_sw + c_nw + c_mid)/3.0;
    ze = (c_se + c_ne + c_mid)/3.0;
    zs = (c_sw + c_se + c_mid)/3.0;
    zn = (c_nw + c_ne + c_mid)/3.0;

    switch (cplot){
        case (P_3D_ENERGY):
        {
            energy_color_scheme(palette, VSCALE_ENERGY*ze, rgb_e);
            energy_color_scheme(palette, VSCALE_ENERGY*zw, rgb_w);
            energy_color_scheme(palette, VSCALE_ENERGY*zn, rgb_n);
            energy_color_scheme(palette, VSCALE_ENERGY*zs, rgb_s);
            break;
        }
        case (P_3D_LOG_ENERGY):
        {
            log_energy_color_scheme(palette, ze, rgb_e);
            log_energy_color_scheme(palette, zw, rgb_w);
            log_energy_color_scheme(palette, zn, rgb_n);
            log_energy_color_scheme(palette, zs, rgb_s);         
            break;
        }
        case (P_3D_PHASE):
        {
            phase_color_scheme(palette, ze, rgb_e);
            phase_color_scheme(palette, zw, rgb_w);
            phase_color_scheme(palette, zn, rgb_n);
            phase_color_scheme(palette, zs, rgb_s);         
            break;
        }
        default:
        {
        hsl_to_rgb_palette(VSCALE_AMPLITUDE*ze, 0.9, 0.5, rgb_e, palette);
        hsl_to_rgb_palette(VSCALE_AMPLITUDE*zw, 0.9, 0.5, rgb_w, palette);
        hsl_to_rgb_palette(VSCALE_AMPLITUDE*zn, 0.9, 0.5, rgb_n, palette);
        hsl_to_rgb_palette(VSCALE_AMPLITUDE*zs, 0.9, 0.5, rgb_s, palette);
//             color_scheme_palette(COLOR_SCHEME, palette, VSCALE_AMPLITUDE*ze, 1.0, 0, rgb_e);
//             color_scheme_palette(COLOR_SCHEME, palette, VSCALE_AMPLITUDE*zw, 1.0, 0, rgb_w);
//             color_scheme_palette(COLOR_SCHEME, palette, VSCALE_AMPLITUDE*zn, 1.0, 0, rgb_n);
//             color_scheme_palette(COLOR_SCHEME, palette, VSCALE_AMPLITUDE*zs, 1.0, 0, rgb_s); 
        }
    }
}

void draw_wave_3d(double phi[NX*NY], double psi[NX*NY], short int xy_in[NX*NY], int plot, int cplot, int palette, int fade, double fade_value)
{
    int i, j, k, l, draw = 1;
    double xy[2], xy_screen[2], rgb[3], pos[2], ca, rgb_e[3], rgb_w[3], rgb_n[3], rgb_s[3]; 
    double z_sw, z_se, z_nw, z_ne, z_mid, zw, ze, zn, zs, min = 1000.0, max = 0.0;
    double xy_sw[2], xy_se[2], xy_nw[2], xy_ne[2], xy_mid[2];
    double energy;
    double *cos_angle, *energies, *phases;
    
    blank();
    draw_billiard_3d(fade, fade_value);
    cos_angle = (double *)malloc(NX*NY*sizeof(double));
    energies = (double *)malloc(NX*NY*sizeof(double));
    phases = (double *)malloc(NX*NY*sizeof(double));
    
    
    if ((plot == P_3D_ENERGY)||(plot == P_3D_LOG_ENERGY))
//     if (plot == P_3D_ENERGY)
    {
        compute_energy_field(phi, psi, xy_in, energies);
        compute_light_angle(energies, xy_in, cos_angle);
    }
//     else if (plot == P_3D_LOG_ENERGY)
//     {
//         compute_log_energy_field(phi, psi, xy_in, energies);
//         compute_light_angle(energies, xy_in, cos_angle);
//     }
    else if (plot != P_3D_AMPLITUDE) compute_light_angle(phi, xy_in, cos_angle);
    
    if (cplot == P_3D_PHASE) 
    {
        compute_phase_field(phi, psi, xy_in, phases);
//         compute_energy_field(phi, psi, xy_in, energies);
    }
    
//     if ((plot == P_3D_ANGLE)||(plot == P_3D_AMP_ANGLE)) compute_light_angle(phi, xy_in, cos_angle);
    
//     for (i=0; i<NX-2; i++)
//         for (j=0; j<NY-2; j++)
    for (i=1; i<NX-2; i++)
        for (j=1; j<NY-2; j++)
        {
            if ((TWOSPEEDS)||(xy_in[i*NY+j]))
            {
                switch (plot) {
                    case (P_3D_AMPLITUDE):
                    {
                        if (AMPLITUDE_HIGH_RES)
                        {
                            compute_interpolated_colors(i, j, phi, palette, plot, 
                                 &z_sw, &z_se, &z_nw, &z_ne, &z_mid, rgb_e, rgb_w, rgb_n, rgb_s);
                        }
                        else color_scheme_palette(COLOR_SCHEME, palette, VSCALE_AMPLITUDE*phi[i*NY+j], 1.0, 0, rgb);
                        break;
                    }
                    case (P_3D_ANGLE):
                    {
                        ca = cos_angle[i*NY+j];
                        if (ca < 0.0) ca = 0.0;
                        color_scheme_asym_palette(COLOR_SCHEME, palette, ca, 1.0, 0, rgb);
                        break;
                    }
                    case (P_3D_AMP_ANGLE):
                    {
                        ca = cos_angle[i*NY+j];
                        ca = (ca + 1.0)*0.4 + 0.2;
                        if ((FADE_IN_OBSTACLE)&&(!xy_in[i*NY+j])) ca *= 1.6;
                        if (AMPLITUDE_HIGH_RES)
                        {
                            compute_interpolated_colors_new(i, j, phi, phi, palette, cplot, 
                                 &z_sw, &z_se, &z_nw, &z_ne, &z_mid, rgb_e, rgb_w, rgb_n, rgb_s);
                            for (k=0; k<3; k++) rgb_e[k] *= ca;
                            for (k=0; k<3; k++) rgb_w[k] *= ca;
                            for (k=0; k<3; k++) rgb_n[k] *= ca;
                            for (k=0; k<3; k++) rgb_s[k] *= ca;
                        }
                        else 
                        {
                            color_scheme_palette(COLOR_SCHEME, palette, VSCALE_AMPLITUDE*phi[i*NY+j], 1.0, 0, rgb);
                            for (k=0; k<3; k++) rgb[k] *= ca;
                        }
                        break;
                    }
                    case (P_3D_ENERGY):
                    {
                        ca = cos_angle[i*NY+j];
                        ca = (ca + 1.0)*0.4 + 0.19;
                        if ((FADE_IN_OBSTACLE)&&(!xy_in[i*NY+j])) ca *= 1.6;
                        
                        if (AMPLITUDE_HIGH_RES)
                        {
                            compute_interpolated_colors(i, j, energies, palette, plot, 
                                 &z_sw, &z_se, &z_nw, &z_ne, &z_mid, rgb_e, rgb_w, rgb_n, rgb_s);
                            for (k=0; k<3; k++) 
                            {
                                rgb_e[k] *= ca;
                                rgb_w[k] *= ca;
                                rgb_n[k] *= ca;
                                rgb_s[k] *= ca;
                            }
                        }
                        else 
                        {
                            energy_color_scheme(palette, energies[i*NY+j], rgb);
                            for (k=0; k<3; k++) rgb[k] *= ca;
                        }
                        break;
                    }
                    case (P_3D_LOG_ENERGY):
                    {
                        ca = cos_angle[i*NY+j];
                        ca = (ca + 1.0)*0.4 + 0.19;
                        if ((FADE_IN_OBSTACLE)&&(!xy_in[i*NY+j])) ca *= 1.6;
                        
                        if (AMPLITUDE_HIGH_RES)
                        {
                            compute_interpolated_colors(i, j, energies, palette, plot, 
                                 &z_sw, &z_se, &z_nw, &z_ne, &z_mid, rgb_e, rgb_w, rgb_n, rgb_s);
                            for (k=0; k<3; k++) 
                            {
                                rgb_e[k] *= ca;
                                rgb_w[k] *= ca;
                                rgb_n[k] *= ca;
                                rgb_s[k] *= ca;
                            }
                        }
                        else 
                        {
                            log_energy_color_scheme(palette, energies[i*NY+j], rgb);
                            for (k=0; k<3; k++) rgb[k] *= ca;
                        }
                        break;
                    }
                }
                
                if (fade)
                {
                    for (k=0; k<3; k++)
                    {
                        rgb[k] *= fade_value;
                        rgb_n[k] *= fade_value;
                        rgb_s[k] *= fade_value;
                        rgb_e[k] *= fade_value;
                        rgb_w[k] *= fade_value;
                    }
                    
                }
                
                if ((plot == P_3D_ENERGY)||(plot == P_3D_LOG_ENERGY))
                {                    
                    draw = (xy_in[i*NY+j])&&(xy_in[(i+1)*NY+j])&&(xy_in[i*NY+j+1])&&(xy_in[(i+1)*NY+j+1]);
                }
                
                if ((plot != P_3D_AMPLITUDE)&&(AMPLITUDE_HIGH_RES)&&(draw))
                {
                    ij_to_xy(i, j, xy_sw);
                    ij_to_xy(i+1, j, xy_se);
                    ij_to_xy(i, j+1, xy_nw);
                    ij_to_xy(i+1, j+1, xy_ne);
                    
                    for (k=0; k<2; k++) xy_mid[k] = 0.25*(xy_sw[k] + xy_se[k] + xy_nw[k] + xy_ne[k]);
                    
                    glBegin(GL_TRIANGLE_FAN);
                    glColor3f(rgb_w[0], rgb_w[1], rgb_w[2]);
                    draw_vertex_xyz(xy_mid, z_mid);
                    draw_vertex_xyz(xy_nw, z_nw);
                    draw_vertex_xyz(xy_sw, z_sw);
                    
                    glColor3f(rgb_s[0], rgb_s[1], rgb_s[2]);
                    draw_vertex_xyz(xy_se, z_se);
                    
                    glColor3f(rgb_e[0], rgb_e[1], rgb_e[2]);
                    draw_vertex_xyz(xy_ne, z_ne);
                    
                    glColor3f(rgb_n[0], rgb_n[1], rgb_n[2]);
                    draw_vertex_xyz(xy_nw, z_nw);
                    glEnd ();
                }
                else
                {
                    glColor3f(rgb[0], rgb[1], rgb[2]);
    
                    glBegin(GL_TRIANGLE_FAN);
                    ij_to_xy(i, j, xy);
                    draw_vertex_xyz(xy, phi[i*NY+j]);
                    ij_to_xy(i+1, j, xy);
                    draw_vertex_xyz(xy, phi[(i+1)*NY+j]);
                    ij_to_xy(i+1, j+1, xy);
                    draw_vertex_xyz(xy, phi[(i+1)*NY+j+1]);
                    ij_to_xy(i, j+1, xy);
                    draw_vertex_xyz(xy, phi[i*NY+j+1]);
                    glEnd ();
                }
            }
            else if (DRAW_OUTSIDE_GRAY)     /* experimental */
            {
                printf("Drawing outside at (%i,%i)\n", i, j);
                glColor3f(0.5, 0.5, 0.5);
                glBegin(GL_TRIANGLE_FAN);
                ij_to_xy(i, j, xy);
                draw_vertex_xyz(xy, 0.0);
                ij_to_xy(i+1, j, xy);
                draw_vertex_xyz(xy, 0.0);
                ij_to_xy(i+1, j+1, xy);
                draw_vertex_xyz(xy, 0.0);
                ij_to_xy(i, j+1, xy);
                draw_vertex_xyz(xy, 0.0);
                glEnd ();
            }
        }

    free(cos_angle);
    free(energies);
    free(phases);
}

void draw_color_scheme_palette_3d(double x1, double y1, double x2, double y2, int plot, double min, double max, int palette, int fade, double fade_value)
{
    int j, k, ij_botleft[2], ij_topright[2], imin, imax, jmin, jmax;
    double y, dy, dy_e, dy_phase, rgb[3], value, lum, amp;
    
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
            case (P_3D_AMPLITUDE):
            {
                value = min + 1.0*dy*(double)(j - jmin);
                color_scheme_palette(COLOR_SCHEME, palette, 0.7*value, 1.0, 0, rgb);
                break;
            }
            case (P_3D_ANGLE):
            {
                value = 1.0*dy*(double)(j - jmin);
                color_scheme_asym_palette(COLOR_SCHEME, palette, value, 1.0, 0, rgb);
                break;
            }
            case (P_3D_AMP_ANGLE):
            {
                value = min + 1.0*dy*(double)(j - jmin);
                color_scheme_palette(COLOR_SCHEME, palette, 0.7*value, 1.0, 0, rgb);
                break;
            }
            case (P_3D_ENERGY):
            {
                value = dy_e*(double)(j - jmin)*100.0/E_SCALE;
                if (COLOR_PALETTE >= COL_TURBO) color_scheme_asym_palette(COLOR_SCHEME, palette, value, 1.0, 1, rgb);
                else color_scheme_palette(COLOR_SCHEME, palette, value, 1.0, 1, rgb);
                break;
            }
            case (P_3D_LOG_ENERGY):
            {
//                 value = LOG_SHIFT + LOG_SCALE*dy_e*(double)(j - jmin)*100.0/E_SCALE;
                value = LOG_SCALE*dy_e*(double)(j - jmin)*100.0/E_SCALE;
                color_scheme_asym_palette(COLOR_SCHEME, palette, value, 1.0, 1, rgb);
                break;
            }
            case (P_3D_PHASE):
            {
                value = dy_phase*(double)(j - jmin);
                color_scheme_palette(C_ONEDIM_LINEAR, palette, value, 1.0, 1, rgb);
                break;
            }
            case (Z_MODULE):
            {
                value = min + 1.0*dy*(double)(j - jmin);
                color_scheme_palette(COLOR_SCHEME, palette, 0.7*value, 1.0, 0, rgb);
                break;
            }
            case (Z_ARGUMENT):
            {
                value = dy_phase*(double)(j - jmin);
//                 hsl_to_rgb_palette(value, 0.9, 0.5, rgb, palette);
                color_scheme_palette(C_ONEDIM_LINEAR, palette, value, 1.0, 1, rgb);
                break;
            }
            case (Z_REALPART):
            {
                value = min + 1.0*dy*(double)(j - jmin);
                color_scheme_palette(COLOR_SCHEME, palette, value, 1.0, 1, rgb);
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
            draw_vertex_ij(j, imin);
            draw_vertex_ij(j, imax);
            draw_vertex_ij(j+1, imax);
            draw_vertex_ij(j+1, imin);            
        }
        else
        {
            draw_vertex_ij(imin, j);
            draw_vertex_ij(imax, j);
            draw_vertex_ij(imax, j+1);
            draw_vertex_ij(imin, j+1);
        }
    }
    glEnd ();
    
    if (fade) glColor3f(fade_value, fade_value, fade_value);
    else glColor3f(1.0, 1.0, 1.0);
    glLineWidth(BOUNDARY_WIDTH);
    draw_rectangle_noscale(x1, y1, x2, y2);
}


void draw_circular_color_scheme_palette_3d(double x1, double y1, double radius, int plot, double min, double max, int palette, int fade, double fade_value)
{
    int j, k, ij_center[2], ij_right[2], ic, jc, ir;
    double x, y, dy, dy_e, dy_phase, rgb[3], value, lum, amp, dphi, pos[2], phi, xy[2], zscale = 0.95;
    
//     printf("Drawing color bar\n");
    
    xy_to_ij(x1, y1, ij_center);
    xy_to_ij(x1 + radius, y1, ij_right);
    
//     rgb[0] = 0.0;   rgb[1] = 0.0;   rgb[2] = 0.0;
//     erase_area_rgb(0.5*(x1 + x2), x2 - x1, 0.5*(y1 + y2), y2 - y1, rgb);

    ic = ij_center[0];
    jc = ij_center[1];
    ir = ij_right[0] - ij_center[0];
//     imax = ij_topright[0];
//     jmax = ij_topright[1];    

        
//     glBegin(GL_TRIANGLE_FAN);
//     draw_vertex_ij(ic, jc);
    dy = (max - min)/360.0;
    dy_e = max/360.0;
    dy_phase = 1.0/360.0;
    dphi = DPI/360.0;
    
    for (j = 0; j < 361; j++)
    {
        switch (plot) {
            case (P_3D_AMPLITUDE):
            {
                value = min + 1.0*dy*(double)(j);
                color_scheme_palette(COLOR_SCHEME, palette, 0.7*value, 1.0, 0, rgb);
                break;
            }
            case (P_3D_ANGLE):
            {
                value = 1.0*dy*(double)(j);
                color_scheme_asym_palette(COLOR_SCHEME, palette, value, 1.0, 0, rgb);
                break;
            }
            case (P_3D_AMP_ANGLE):
            {
                value = min + 1.0*dy*(double)(j);
                color_scheme_palette(COLOR_SCHEME, palette, 0.7*value, 1.0, 0, rgb);
                break;
            }
            case (P_3D_ENERGY):
            {
                value = dy_e*(double)(j)*100.0/E_SCALE;
                if (COLOR_PALETTE >= COL_TURBO) color_scheme_asym_palette(COLOR_SCHEME, palette, value, 1.0, 1, rgb);
                else color_scheme_palette(COLOR_SCHEME, palette, value, 1.0, 1, rgb);
                break;
            }
            case (P_3D_LOG_ENERGY):
            {
                value = LOG_SCALE*dy_e*(double)(j)*100.0/E_SCALE;
                color_scheme_asym_palette(COLOR_SCHEME, palette, value, 1.0, 1, rgb);
                break;
            }
            case (P_3D_TOTAL_ENERGY):
            {
                value = dy_e*(double)(j)*100.0/E_SCALE;
                if (COLOR_PALETTE >= COL_TURBO) color_scheme_asym_palette(COLOR_SCHEME, palette, value, 1.0, 1, rgb);
                else color_scheme_palette(COLOR_SCHEME, palette, value, 1.0, 1, rgb);
                break;
            }
            case (P_3D_LOG_TOTAL_ENERGY):
            {
                value = LOG_SCALE*dy_e*(double)(j)*100.0/E_SCALE;
                color_scheme_asym_palette(COLOR_SCHEME, palette, value, 1.0, 1, rgb);
                break;
            }
            case (P_3D_MEAN_ENERGY):
            {
                value = dy_e*(double)(j)*100.0/E_SCALE;
                if (COLOR_PALETTE >= COL_TURBO) color_scheme_asym_palette(COLOR_SCHEME, palette, value, 1.0, 1, rgb);
                else color_scheme_palette(COLOR_SCHEME, palette, value, 1.0, 1, rgb);
                break;
            }
            case (P_3D_LOG_MEAN_ENERGY):
            {
                value = LOG_SCALE*dy_e*(double)(j)*100.0/E_SCALE;
                color_scheme_asym_palette(COLOR_SCHEME, palette, value, 1.0, 1, rgb);
                break;
            }
            case (P_3D_PHASE):
            {
                value = dy_phase*(double)(j);
                color_scheme_palette(C_ONEDIM_LINEAR, palette, value, 1.0, 1, rgb);
                break;
            }
            case (P_3D_FLUX_INTENSITY):
            {
                value = dy_e*(double)(j)*100.0/E_SCALE;
                if (COLOR_PALETTE >= COL_TURBO) color_scheme_asym_palette(COLOR_SCHEME, palette, value, 1.0, 1, rgb);
                else color_scheme_palette(COLOR_SCHEME, palette, value, 1.0, 1, rgb);
                break;
            }
            case (P_3D_FLUX_DIRECTION):
            {
                value = dy_phase*(double)(j);
                color_scheme_palette(C_ONEDIM_LINEAR, palette, value, 1.0, 1, rgb);
                break;
            }
            case (Z_POLAR):
            {
                value = dy_phase*(double)(j);
                color_scheme_palette(C_ONEDIM_LINEAR, palette, value, 1.0, 1, rgb);
                break;
            }
            case (Z_ARGUMENT):
            {
                value = dy_phase*(double)(j);
                color_scheme_palette(C_ONEDIM_LINEAR, palette, value, 1.0, 1, rgb);
                break;
            }
            case (Z_EULER_DIRECTION):
            {
                value = dy_phase*(double)(j);
                color_scheme_palette(C_ONEDIM_LINEAR, palette, value, 1.0, 1, rgb);
                break;
            }
            case (Z_EULER_DIRECTION_SPEED):
            {
                value = 0.5 - dy_phase*(double)(j);
                color_scheme_palette(C_ONEDIM_LINEAR, palette, value, 1.0, 1, rgb);
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
            case (Z_EULERC_VORTICITY):
             {
                value = min + 1.0*dy*(double)(j);
                printf("Palette value %.3lg\n", value);
                color_scheme_palette(COLOR_SCHEME, palette, 0.7*value, 1.0, 0, rgb);
                break;
            }
            case (Z_SWATER_DIRECTION_SPEED):
            {
                value = dy_phase*(double)(j) + 0.5 - 0.5*PHASE_SHIFT;
                if (value > 1.0) value -= 1.0;
                if (value < 0.0) value += 1.0;
                color_scheme_palette(C_ONEDIM_LINEAR, palette, value, 1.0, 1, rgb);
                break;
            }
               
        }
        if (fade) for (k=0; k<3; k++) rgb[k] *= fade_value;
        glColor3f(rgb[0], rgb[1], rgb[2]);
        
        glBegin(GL_TRIANGLE_FAN);
        draw_vertex_ij(ic, jc);
        x = cos(dphi*(double)j)*(double)ir;
        y = zscale*sin(dphi*(double)j)*(double)ir;
        xy[0] = XMIN + ((double)ic + x)*(XMAX-XMIN)/((double)NX);
        xy[1] = YMIN + ((double)jc + y)*(YMAX-YMIN)/((double)NY);
//         ij_to_xy(ic + x, jc + y, xy);
        glVertex2d(xy[0], xy[1]);
        x = cos(dphi*(double)(j+1))*(double)ir;
        y = zscale*sin(dphi*(double)(j+1))*(double)ir;
        xy[0] = XMIN + ((double)ic + x)*(XMAX-XMIN)/((double)NX);
        xy[1] = YMIN + ((double)jc + y)*(YMAX-YMIN)/((double)NY);
//         ij_to_xy(ic + x, jc + y, xy);
        glVertex2d(xy[0], xy[1]);
        glEnd ();
    }
//     draw_vertex_ij(ic + ir, jc);
//     draw_vertex_ij(ic, jc);
//     glEnd ();
    
    if (fade) glColor3f(fade_value, fade_value, fade_value);
    else glColor3f(1.0, 1.0, 1.0);
    glLineWidth(BOUNDARY_WIDTH*3/2);
    glEnable(GL_LINE_SMOOTH);
    
    dphi = DPI/(double)NSEG;
    glBegin(GL_LINE_LOOP);
    for (j = 0; j < NSEG; j++)
    {        
        x = cos(dphi*(double)j)*(double)ir;
        y = zscale*sin(dphi*(double)j)*(double)ir;
        xy[0] = XMIN + ((double)ic + x)*(XMAX-XMIN)/((double)NX);
        xy[1] = YMIN + ((double)jc + y)*(YMAX-YMIN)/((double)NY);
        glVertex2d(xy[0], xy[1]);
    }
    glEnd ();

}

