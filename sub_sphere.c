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
            
            ij_to_xy(NX-1-i,j,xy);
//             xy[0] = XMIN + ((double)(NX-i-1))*(XMAX-XMIN)/((double)NX);
//             xy[1] = YMIN + ((double)(j-DPOLE))*(YMAX-YMIN)/((double)(NY-2*DPOLE));
            
            xy[0] += phishift;
            if (xy[0] > XMAX) xy[0] += XMIN - XMAX;
            
            xy[1] *= (double)NY/(double)(NY-2*DPOLE);
            
            wsphere[i*NY+j].x2d = xy[0]; 
            wsphere[i*NY+j].y2d = xy[1]; 
        }
        
        /* cotangent, taking care of not dividing by zero */
        for (j=1; j<NY-1; j++) wsphere[i*NY+j].cottheta = wsphere[i*NY+j].ctheta/wsphere[i*NY+j].stheta;
        wsphere[i*NY].cottheta = wsphere[i*NY+1].cottheta;
        wsphere[i*NY + NY-1].cottheta = wsphere[i*NY+1].cottheta;
    }
}


void init_earth_map(t_wave_sphere wsphere[NX*NY])
/* init file from earth map */
{
    int i, j, ii, jj, k, nx, ny, maxrgb, nmaxpixels = 4915200, scan, rgbval, diff, sshift, nshift;
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
        printf("Image too large, increase nmaxpixels in choose_colors()\n");
        exit(0);
    }
    
    /* shift due to min/max latitudes of image */
    sshift = 0 + DPOLE;
    nshift = 0 + DPOLE; 
    
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
            ii = (int)(rx*(double)(NX-1 - i)) + nx/2;
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
            wsphere[i*NY+j].indomain = (diff < 10);
        }
    
    free(rgb_values);
    fclose(image_file);
}


int ij_to_sphere(int i, int j, double r, t_wave_sphere wsphere[NX*NY], double xyz[3])
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
    
    newr = wsphere[i*NY+j].radius;
    
    xyz[0] = wsphere[i*NY+j].x;
    xyz[1] = wsphere[i*NY+j].y;
    xyz[2] = wsphere[i*NY+j].z;
    
    pscal = xyz[0]*observer[0] + xyz[1]*observer[1] + xyz[2]*observer[2];
    
    xyz[0] *= newr;
    xyz[1] *= newr;
    xyz[2] *= newr;
    
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
        else if ((B_DOMAIN == D_SPHERE_JULIA)||(B_DOMAIN == D_SPHERE_JULIA_INV))
        {
            cos_rot = cos(JULIA_ROT*DPI/360.0);
            sin_rot = sin(JULIA_ROT*DPI/360.0);
        }
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
        case (D_SPHERE_EARTH):
        {
            return(wsphere[i*NY+j].indomain);
        }
        default:
        {
            printf("Function xy_in_billiard_sphere not defined for this billiard \n");
            return(0);
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
    int i, j;
    double xy[2], pscal, dist, dist2, stheta, ctheta;

    ctheta = cos(theta0 + PID);
    stheta = sin(theta0 + PID);
    
    printf("Initializing wave\n"); 
//     #pragma omp parallel for private(i,j,xy,dist2)
    for (i=0; i<NX; i++)
    {
        if (i%100 == 0) printf("Initializing column %i of %i\n", i, NX);
        for (j=0; j<NY; j++)
        {
            pscal = cos(wsphere[i*NY+j].phi - phi0)*wsphere[i*NY+j].stheta*stheta;
            pscal += wsphere[i*NY+j].ctheta*ctheta;
            
            dist = acos(pscal);
            dist2 = dist*dist;
            
	    xy_in[i*NY+j] = xy_in_billiard_sphere(i, j, wsphere);
            
	    if ((xy_in[i*NY+j])||(TWOSPEEDS)) 
                phi[i*NY+j] = INITIAL_AMP*exp(-dist2/INITIAL_VARIANCE)*cos(-sqrt(dist2)/INITIAL_WAVELENGTH);
            else phi[i*NY+j] = 0.0;
            psi[i*NY+j] = 0.0;
        }
    }
}

void add_circular_wave_sphere(double amp, double phi0, double theta0, double phi[NX*NY], double psi[NX*NY], short int xy_in[NX*NY], t_wave_sphere wsphere[NX*NY])
/* initialise field with drop - phi is wave height, psi is phi at time t-1 */
/* phi0 is longidude, theta0 is latitude */
{
    int i, j;
    double xy[2], pscal, dist, dist2, stheta, ctheta;

    ctheta = cos(theta0 + PID);
    stheta = sin(theta0 + PID);
    
    printf("Initializing wave\n"); 
//     #pragma omp parallel for private(i,j,xy,dist2)
    for (i=0; i<NX; i++)
    {
        if (i%100 == 0) printf("Initializing column %i of %i\n", i, NX);
        for (j=0; j<NY; j++)
        {
            pscal = cos(wsphere[i*NY+j].phi - phi0)*wsphere[i*NY+j].stheta*stheta;
            pscal += wsphere[i*NY+j].ctheta*ctheta;
            
            dist = acos(pscal);
            dist2 = dist*dist;
            
	    if ((xy_in[i*NY+j])||(TWOSPEEDS)) 
                phi[i*NY+j] += amp*INITIAL_AMP*exp(-dist2/INITIAL_VARIANCE)*cos(-sqrt(dist2)/INITIAL_WAVELENGTH);
            else phi[i*NY+j] = 0.0;
            psi[i*NY+j] = 0.0;
        }
    }
}


void compute_light_angle_sphere(short int xy_in[NX*NY], t_wave wave[NX*NY], t_wave_sphere wsphere[NX*NY], int movie)
/* computes cosine of angle between normal vector and vector light */
{
    int i, j;
    double x, y, z, norm, pscal, deltai[3], deltaj[3], deltar, n[3], r;
    static double dphi, dtheta;
    static int first = 1;
    
    if (first)
    {
        dphi = DPI/(double)NX;
        dtheta = PI/(double)NY;
        first = 0;
    }
    
    #pragma omp parallel for private(i,j)
    for (i=0; i<NX; i++)
        for (j=0; j<NY; j++)
        {
            r = 1.0 + RSCALE*(*wave[i*NY+j].p_zfield[movie]);
            if (r > RMAX) r = RMAX;
            wsphere[i*NY+j].radius = r;
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
                    
                    /* computation of normal vector */
                    n[0] = deltai[1]*deltaj[2] - deltai[2]*deltaj[1];
                    n[1] = deltai[2]*deltaj[0] - deltai[0]*deltaj[2];
                    n[2] = deltai[0]*deltaj[1] - deltai[1]*deltaj[0];
                                        
                    norm = sqrt(n[0]*n[0] + n[1]*n[1] + n[2]*n[2]);
                    
                    pscal = n[0]*light[0] + n[1]*light[1] + n[2]*light[2];
                
                    wave[i*NY+j].cos_angle = pscal/norm;
                }
                else 
                {
                    pscal = wsphere[i*NY+j].x*light[0] + wsphere[i*NY+j].y*light[1] + wsphere[i*NY+j].z*light[2];
                
                    wave[i*NY+j].cos_angle = pscal;
                }
            }
            
        for (i=0; i<NX-1; i++) wave[i*NY+NY-1].cos_angle = wave[i*NY+NY-2].cos_angle;
        for (j=0; j<NY-1; j++) wave[(NX-1)*NY+j].cos_angle = wave[(NX-2)*NY+j].cos_angle;
        wave[(NX-1)*NY+NY-1].cos_angle = wave[(NX-1)*NY+NY-2].cos_angle;
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

void compute_light_angle_sphere_2d(short int xy_in[NX*NY], t_wave wave[NX*NY], int movie)
/* computes cosine of angle between normal vector and vector light */
{
    int i, j;
    double gradx, grady, norm, pscal;
    static double dx, dy, vscale2;
    static int first = 1;
    
    if (first)
    {
        dx = 2.0*(XMAX - XMIN)/(double)NX;
        dy = 2.0*(YMAX - YMIN)/(double)NY;
        vscale2 = SHADE_SCALE_2D*SHADE_SCALE_2D;
        first = 0;
    }
    
    #pragma omp parallel for private(i,j,gradx, grady, norm, pscal)
    for (i=1; i<NX-1; i++)
        for (j=1; j<NY-1; j++)
        {
            if ((TWOSPEEDS)||(xy_in[i*NY+j]))
            {
                gradx = (*wave[(i+1)*NY+j].p_zfield[movie] - *wave[(i-1)*NY+j].p_zfield[movie])/dx;
                grady = (*wave[i*NY+j+1].p_zfield[movie] - *wave[i*NY+j-1].p_zfield[movie])/dy;
                
                norm = sqrt(vscale2 + gradx*gradx + grady*grady);
                pscal = -gradx*light[0] - grady*light[1] + SHADE_SCALE_2D;
                
                wave[i*NY+j].cos_angle = pscal/norm;
            }
        }
    
    /* i=0 */
    for (j=1; j<NY-1; j++)
    {
        if ((TWOSPEEDS)||(xy_in[j]))
        {
            gradx = (*wave[NY+j].p_zfield[movie] - *wave[(NY-1)*NY+j].p_zfield[movie])/dx;
            grady = (*wave[j+1].p_zfield[movie] - *wave[j-1].p_zfield[movie])/dy;
                
            norm = sqrt(vscale2 + gradx*gradx + grady*grady);
            pscal = -gradx*light[0] - grady*light[1] + SHADE_SCALE_2D;
                
            wave[j].cos_angle = pscal/norm;
        }
    }
    
    /* i=NX-1 */
    for (j=1; j<NY-1; j++)
    {
        if ((TWOSPEEDS)||(xy_in[(NX-1)*NY+j]))
        {
            gradx = (*wave[j].p_zfield[movie] - *wave[(NX-2)*NY+j].p_zfield[movie])/dx;
            grady = (*wave[(NX-1)*NY+j+1].p_zfield[movie] - *wave[(NX-1)*NY+j-1].p_zfield[movie])/dy;
                
            norm = sqrt(vscale2 + gradx*gradx + grady*grady);
            pscal = -gradx*light[0] - grady*light[1] + SHADE_SCALE_2D;
                
            wave[(NX-1)*NY+j].cos_angle = pscal/norm;
        }
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
    double x, y, ca;
    static int ishift;
    static double dx, dy;
    
    if (refresh)
    {
        compute_wave_fields(phi, psi, xy_in, zplot, cplot, wave);
        if (SHADE_3D) compute_light_angle_sphere(xy_in, wave, wsphere, movie);
        else if (SHADE_2D) compute_light_angle_sphere_2d(xy_in, wave, movie);
        compute_cfield_sphere(xy_in, cplot, palette, wave, fade, fade_value, movie);
        dx = (XMAX - XMIN)/(double)NX;
        dy = (YMAX - YMIN)/(double)(NY-2*DPOLE);
        ishift = (int)((double)NX*PHISHIFT/360.0);
    }
    
    glBegin(GL_QUADS);
    
    for (i=0; i<NX; i++)
        for (j=0; j<NY; j++)
        {
            if ((TWOSPEEDS)||(xy_in[i*NY+j]))
                glColor3f(wave[i*NY+j].rgb[0], wave[i*NY+j].rgb[1], wave[i*NY+j].rgb[2]);
            else
            {
                ca = wave[i*NY+j].cos_angle;
                ca = (ca + 1.0)*0.4 + 0.2;
                if (fade) ca *= fade_value;
                if (B_DOMAIN == D_SPHERE_EARTH)
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
}

void draw_wave_sphere_ij(int i, int iplus, int j, int jplus, int jcolor, int movie, double phi[NX*NY], double psi[NX*NY], short int xy_in[NX*NY], t_wave wave[NX*NY], t_wave_sphere wsphere[NX*NY], int zplot, int cplot, int palette, int fade, double fade_value)
/* draw wave at simulation grid point (i,j) */
{
    int k, l;
    double xyz[3], ca;
    
    if ((TWOSPEEDS)||(xy_in[i*NY+j]))
        glColor3f(wave[i*NY+jcolor].rgb[0], wave[i*NY+jcolor].rgb[1], wave[i*NY+jcolor].rgb[2]);
    else
    {
        ca = wave[i*NY+j].cos_angle;
        ca = (ca + 1.0)*0.4 + 0.2;
        if (fade) ca *= fade_value;
        if (B_DOMAIN == D_SPHERE_EARTH)
            glColor3f(wsphere[i*NY+j].r*ca, wsphere[i*NY+j].g*ca, wsphere[i*NY+j].b*ca);
        else glColor3f(COLOR_OUT_R*ca, COLOR_OUT_G*ca, COLOR_OUT_B*ca);
    }
    glBegin(GL_TRIANGLE_FAN);
    if (ij_to_sphere(i, j, *wave[i*NY+j].p_zfield[movie], wsphere, xyz))
        draw_vertex_sphere(xyz);
    if (ij_to_sphere(iplus, j, *wave[iplus*NY+j].p_zfield[movie], wsphere, xyz))
        draw_vertex_sphere(xyz);
    if (ij_to_sphere(iplus, jplus, *wave[iplus*NY+j+1].p_zfield[movie], wsphere, xyz))
        draw_vertex_sphere(xyz);
    if (ij_to_sphere(i, jplus, *wave[i*NY+j+1].p_zfield[movie], wsphere, xyz))
        draw_vertex_sphere(xyz);
    glEnd ();
}


void draw_wave_sphere_3D(int movie, double phi[NX*NY], double psi[NX*NY], short int xy_in[NX*NY], t_wave wave[NX*NY], t_wave_sphere wsphere[NX*NY], int zplot, int cplot, int palette, int fade, double fade_value, int refresh)
{
    int i, j, imax, imin;
    double observer_angle, angle2;
    
    blank();
            
    if (refresh)
    {
        compute_wave_fields(phi, psi, xy_in, zplot, cplot, wave);
        if (SHADE_3D) compute_light_angle_sphere(xy_in, wave, wsphere, movie);
        else if (SHADE_2D) compute_light_angle_sphere_2d(xy_in, wave, movie);
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
    
//     printf("Angle = %.5lg, angle2 = %.5lg, imin = %i, imax = %i\n", observer_angle, angle2, imin, imax);
    
    if (observer[2] > 0.0)
    {
        if (imin < imax)
        {
            for (i=imax; i>imin; i--)
                for (j=0; j<NY-2; j++)
                    draw_wave_sphere_ij(i, i+1, j, j+1, j+1, movie, phi, psi, xy_in, wave, wsphere, zplot, cplot, palette, fade, fade_value);
        
            for (i=imax+1; i<NX-1; i++)
                for (j=0; j<NY-2; j++)
                    draw_wave_sphere_ij(i, i+1, j, j+1, j+1, movie, phi, psi, xy_in, wave, wsphere, zplot, cplot, palette, fade, fade_value);
        
            for (j=0; j<NY-2; j++)
                draw_wave_sphere_ij(NX-1, 0, j, j+1, j+1, movie, phi, psi, xy_in, wave, wsphere, zplot, cplot, palette, fade, fade_value);
        
            for (i=0; i<=imin; i++)
                for (j=0; j<NY-2; j++)
                    draw_wave_sphere_ij(i, i+1, j, j+1, j+1, movie, phi, psi, xy_in, wave, wsphere, zplot, cplot, palette, fade, fade_value);
        }
        else
        {
            for (i=imax; i<imin; i++)
                for (j=0; j<NY-2; j++)
                    draw_wave_sphere_ij(i, i+1, j, j+1, j+1, movie, phi, psi, xy_in, wave, wsphere, zplot, cplot, palette, fade, fade_value);
        
            for (i=imax-1; i>=0; i--)
                for (j=0; j<NY-2; j++)
                    draw_wave_sphere_ij(i, i+1, j, j+1, j+1, movie, phi, psi, xy_in, wave, wsphere, zplot, cplot, palette, fade, fade_value);
        
            for (j=0; j<NY-2; j++)
                draw_wave_sphere_ij(NX-1, 0, j, j+1, j+1, movie, phi, psi, xy_in, wave, wsphere, zplot, cplot, palette, fade, fade_value);
        
            for (i=NX-2; i>=imin; i--)
                for (j=0; j<NY-2; j++)
                    draw_wave_sphere_ij(i, i+1, j, j+1, j+1, movie, phi, psi, xy_in, wave, wsphere, zplot, cplot, palette, fade, fade_value);
        }
    
        /* North pole */
        for (i=0; i<NX-1; i++) for (j=NY-2; j<NY-1; j++)
            draw_wave_sphere_ij(i, i+1, j-1, j, j, movie, phi, psi, xy_in, wave, wsphere, zplot, cplot, palette, fade, fade_value);
        
        for (j=NY-2; j<NY-1; j++)
            draw_wave_sphere_ij(NX-1, 0, j-1, j, NY-DPOLE, movie, phi, psi, xy_in, wave, wsphere, zplot, cplot, palette, fade, fade_value);
    }
    else
    {
        if (imin < imax)
        {
            for (i=imax; i>imin; i--)
                for (j=NY-3; j>=0; j--)
                    draw_wave_sphere_ij(i, i+1, j, j+1, j+1, movie, phi, psi, xy_in, wave, wsphere, zplot, cplot, palette, fade, fade_value);
        
            for (i=imax+1; i<NX-1; i++)
                for (j=NY-3; j>=0; j--)
                    draw_wave_sphere_ij(i, i+1, j, j+1, j+1, movie, phi, psi, xy_in, wave, wsphere, zplot, cplot, palette, fade, fade_value);
        
            for (j=NY-3; j>=0; j--)
                draw_wave_sphere_ij(NX-1, 0, j, j+1, j+1, movie, phi, psi, xy_in, wave, wsphere, zplot, cplot, palette, fade, fade_value);
        
            for (i=0; i<=imin; i++)
                for (j=NY-3; j>=0; j--)
                    draw_wave_sphere_ij(i, i+1, j, j+1, j+1, movie, phi, psi, xy_in, wave, wsphere, zplot, cplot, palette, fade, fade_value);
        }
        else
        {
            for (i=imax; i<imin; i++)
                for (j=NY-3; j>=0; j--)
                    draw_wave_sphere_ij(i, i+1, j, j+1, j+1, movie, phi, psi, xy_in, wave, wsphere, zplot, cplot, palette, fade, fade_value);
        
            for (i=imax-1; i>=0; i--)
                for (j=NY-3; j>=0; j--)
                    draw_wave_sphere_ij(i, i+1, j, j+1, j+1, movie, phi, psi, xy_in, wave, wsphere, zplot, cplot, palette, fade, fade_value);
        
            for (j=NY-3; j>=0; j--)
                draw_wave_sphere_ij(NX-1, 0, j, j+1, j+1, movie, phi, psi, xy_in, wave, wsphere, zplot, cplot, palette, fade, fade_value);
        
            for (i=NX-2; i>=imin; i--)
                for (j=NY-3; j>=0; j--)
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
