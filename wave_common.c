/* functions common to wave_billiard.c, wave_comparison.c, mangrove.c */

void init_xyin(short int * xy_in[NX])
/* initialise table xy_in, needed when obstacles are killed */
// short int * xy_in[NX];
//  
{
    int i, j;
    double xy[2];

    for (i=0; i<NX; i++)
        for (j=0; j<NY; j++)
        {
            ij_to_xy(i, j, xy);
	    xy_in[i][j] = xy_in_billiard(xy[0],xy[1]);
        }
}

void init_xyin_xrange(short int * xy_in[NX], int imin, int imax)
/* initialise table xy_in, needed when obstacles are killed */
// short int * xy_in[NX];
//  
{
    int i, j;
    double xy[2];

    for (i=imin; i<imax; i++)
        for (j=0; j<NY; j++)
        {
            ij_to_xy(i, j, xy);
	    xy_in[i][j] = xy_in_billiard(xy[0],xy[1]);
        }
}


void init_wave(double x, double y, double *phi[NX], double *psi[NX], short int * xy_in[NX])
/* initialise field with drop at (x,y) - phi is wave height, psi is phi at time t-1 */
{
    int i, j;
    double xy[2], dist2;

    for (i=0; i<NX; i++)
        for (j=0; j<NY; j++)
        {
            ij_to_xy(i, j, xy);
            dist2 = (xy[0]-x)*(xy[0]-x) + (xy[1]-y)*(xy[1]-y);
	    xy_in[i][j] = xy_in_billiard(xy[0],xy[1]);
            
	    if ((xy_in[i][j])||(TWOSPEEDS)) phi[i][j] = 0.2*exp(-dist2/0.005)*cos(-sqrt(dist2)/0.1);
// 	    if ((xy_in[i][j])||(TWOSPEEDS)) phi[i][j] = 0.2*exp(-dist2/0.001)*cos(-sqrt(dist2)/0.01);
// 	    if ((xy_in[i][j])||(TWOSPEEDS)) phi[i][j] = 0.2*exp(-dist2/0.00025)*cos(-sqrt(dist2)/0.005);
            else phi[i][j] = 0.0;
            psi[i][j] = 0.0;
        }
}

void init_circular_wave(double x, double y, double *phi[NX], double *psi[NX], short int * xy_in[NX])
/* initialise field with drop at (x,y) - phi is wave height, psi is phi at time t-1 */
{
    int i, j;
    double xy[2], dist2;

    printf("Initializing wave\n"); 
    for (i=0; i<NX; i++)
    {
        if (i%100 == 0) printf("Initializing column %i of %i\n", i, NX);
        for (j=0; j<NY; j++)
        {
            ij_to_xy(i, j, xy);
            dist2 = (xy[0]-x)*(xy[0]-x) + (xy[1]-y)*(xy[1]-y);
	    xy_in[i][j] = xy_in_billiard(xy[0],xy[1]);
            
	    if ((xy_in[i][j])||(TWOSPEEDS)) 
                phi[i][j] = INITIAL_AMP*exp(-dist2/INITIAL_VARIANCE)*cos(-sqrt(dist2)/INITIAL_WAVELENGTH);
            else phi[i][j] = 0.0;
            psi[i][j] = 0.0;
        }
    }
}

void init_wave_plus(double x, double y, double *phi[NX], double *psi[NX], short int * xy_in[NX])
/* initialise field with drop at (x,y) for y > 0 - phi is wave height, psi is phi at time t-1 */
{
    int i, j;
    double xy[2], dist2;

    for (i=0; i<NX; i++)
        for (j=0; j<NY; j++)
        {
            ij_to_xy(i, j, xy);
            dist2 = (xy[0]-x)*(xy[0]-x) + (xy[1]-y)*(xy[1]-y);
	    xy_in[i][j] = xy_in_billiard(xy[0],xy[1]);
            
	    if ((xy[1] > 0.0)&&((xy_in[i][j])||(TWOSPEEDS))) 
                phi[i][j] = INITIAL_AMP*exp(-dist2/INITIAL_VARIANCE)*cos(-sqrt(dist2)/INITIAL_WAVELENGTH);
            else phi[i][j] = 0.0;
            psi[i][j] = 0.0;
        }
}

void init_circular_wave_xplusminus(double xleft, double yleft, double xright, double yright, 
                                   double *phi[NX], double *psi[NX], short int * xy_in[NX])
/* initialise field with two drops, x > 0 and x < 0 do not communicate - phi is wave height, psi is phi at time t-1 */
{
    int i, j;
    double xy[2], dist2;

    for (i=0; i<NX; i++)
        for (j=0; j<NY; j++)
        {
            ij_to_xy(i, j, xy);
            if (xy[0] < 0.0) dist2 = (xy[0]-xleft)*(xy[0]-xleft) + (xy[1]-yleft)*(xy[1]-yleft);
            else dist2 = (xy[0]-xright)*(xy[0]-xright) + (xy[1]-yright)*(xy[1]-yright);
	    xy_in[i][j] = xy_in_billiard(xy[0],xy[1]);
            
	    if ((xy_in[i][j])||(TWOSPEEDS)) 
                phi[i][j] = INITIAL_AMP*exp(-dist2/INITIAL_VARIANCE)*cos(-sqrt(dist2)/INITIAL_WAVELENGTH);
            else phi[i][j] = 0.0;
            psi[i][j] = 0.0;
        }
}

void init_planar_wave(double x, double y, double *phi[NX], double *psi[NX], short int * xy_in[NX])
/* initialise field with drop at (x,y) - phi is wave height, psi is phi at time t-1 */
/* beta version, works for vertical planar wave only so far */
{
    int i, j;
    double xy[2], dist2;

    for (i=0; i<NX; i++)
    {
//         if (i%100 == 0) printf("Initializing column %i\n", i);
        for (j=0; j<NY; j++)
        {
            ij_to_xy(i, j, xy);
            dist2 = (xy[0]-x)*(xy[0]-x);
	    xy_in[i][j] = xy_in_billiard(xy[0],xy[1]);
            
            if ((xy_in[i][j])||(TWOSPEEDS)) 
                phi[i][j] = INITIAL_AMP*exp(-dist2/INITIAL_VARIANCE)*cos(-sqrt(dist2)/INITIAL_WAVELENGTH);
// 	    if ((xy_in[i][j])||(TWOSPEEDS)) phi[i][j] = 0.01*exp(-dist2/0.0005)*cos(-sqrt(dist2)/0.01);
// 	    if ((xy_in[i][j])||(TWOSPEEDS)) phi[i][j] = 0.01*exp(-dist2/0.0005)*cos(-sqrt(dist2)/0.02);
// 	    if ((xy_in[i][j])||(TWOSPEEDS)) phi[i][j] = 0.01*exp(-dist2/0.005)*cos(-sqrt(dist2)/0.02);
// 	    if ((xy_in[i][j])||(TWOSPEEDS)) phi[i][j] = 0.01*exp(-dist2/0.01)*cos(-sqrt(dist2)/0.02);
// 	    if ((xy_in[i][j])||(TWOSPEEDS)) phi[i][j] = 0.01*exp(-dist2/0.0005)*cos(-sqrt(dist2)/0.01);
//             if ((xy_in[i][j])||(TWOSPEEDS)) phi[i][j] = 0.01*exp(-dist2/0.05)*cos(-sqrt(dist2)/0.025);
            else phi[i][j] = 0.0;
            psi[i][j] = 0.0;
        }
    }
}


void init_wave_flat( double *phi[NX], double *psi[NX], short int * xy_in[NX])
/* initialise flat field - phi is wave height, psi is phi at time t-1 */
{
    int i, j;
    double xy[2], dist2;

    for (i=0; i<NX; i++) {
        if (i%100 == 0) printf("Wave and table xy_in - Initialising column %i of %i\n", i, NX);
        for (j=0; j<NY; j++)
        {
            ij_to_xy(i, j, xy);
	    xy_in[i][j] = xy_in_billiard(xy[0],xy[1]);
//             if ((i%10 == 0)&&(j == NY/2)) printf ("xy_in[%i][%i] = %i\n", i, j, xy_in[i][j]);
	    phi[i][j] = 0.0;
            psi[i][j] = 0.0;
        }
    }
}


void add_drop_to_wave(double factor, double x, double y, double *phi[NX], double *psi[NX])
/* OLD VERSION - add drop at (x,y) to the field with given prefactor */
{
    int i, j;
    double xy[2], dist2;

    for (i=0; i<NX; i++)
        for (j=0; j<NY; j++)
        {
            ij_to_xy(i, j, xy);
            dist2 = (xy[0]-x)*(xy[0]-x) + (xy[1]-y)*(xy[1]-y);
            phi[i][j] += 0.2*factor*exp(-dist2/0.001)*cos(-sqrt(dist2)/0.01);
        }
}

void add_circular_wave(double factor, double x, double y, double *phi[NX], double *psi[NX], short int * xy_in[NX])
/* add drop at (x,y) to the field with given prefactor */
{
    int i, j;
    double xy[2], dist2;

    for (i=0; i<NX; i++)
        for (j=0; j<NY; j++)
        {
            ij_to_xy(i, j, xy);
            dist2 = (xy[0]-x)*(xy[0]-x) + (xy[1]-y)*(xy[1]-y);
            if ((xy_in[i][j])||(TWOSPEEDS)) 
                phi[i][j] += INITIAL_AMP*factor*exp(-dist2/INITIAL_VARIANCE)*cos(-sqrt(dist2)/INITIAL_WAVELENGTH);
        }
}


void oscillate_linear_wave(double amplitude, double t, double x1, double y1, double x2, double y2, double *phi[NX],
                           double *psi[NX])
/* oscillating boundary condition at (x,y), beta version */
{
    int i, j, ij1[2], ij2[2], imin, imax, jmin, jmax, d = 5;
    double xy[2], dist2;
    
    xy_to_ij(x1, y1, ij1);
    xy_to_ij(x2, y2, ij2);
    imin = ij1[0] - d;     if (imin < 0) imin = 0;
    imax = ij2[0] + d;     if (imax >= NX) imax = NX-1;
    jmin = ij1[1] - d;     if (jmin < 0) jmin = 0;
    jmax = ij2[1] + d;     if (jmax >= NY) jmax = NY-1;

    for (i = imin; i < imax; i++)
        for (j = jmin; j < jmax; j++)
        {
            ij_to_xy(i, j, xy);
            dist2 = (xy[0]-x1)*(xy[0]-x1);      /* to be improved */
//             dist2 = (xy[0]-x)*(xy[0]-x) + (xy[1]-y)*(xy[1]-y);
//             if (dist2 < 0.01)
            if (dist2 < 0.001)
                phi[i][j] = amplitude*exp(-dist2/0.001)*cos(-sqrt(dist2)/0.01)*cos(t*OMEGA);
//                 phi[i][j] += 0.2*exp(-dist2/0.001)*cos(-sqrt(dist2)/0.01)*cos(t*OMEGA);
        }
}

void compute_gradient(double *phi[NX], double *psi[NX], double x, double y, double gradient[2])
{
    double velocity;
    int iplus, iminus, jplus, jminus, ij[2], i, j;
    
    xy_to_ij(x, y, ij);
    i = ij[0];
    j = ij[1];
                    
    iplus = (i+1);   if (iplus == NX) iplus = NX-1;
    iminus = (i-1);  if (iminus == -1) iminus = 0;
    jplus = (j+1);   if (jplus == NY) jplus = NY-1;
    jminus = (j-1);  if (jminus == -1) jminus = 0;
                        
    gradient[0] = (phi[iplus][j]-phi[i][j])*(phi[iplus][j]-phi[i][j]) 
        + (phi[i][j] - phi[iminus][j])*(phi[i][j] - phi[iminus][j]);
    gradient[1] = (phi[i][jplus]-phi[i][j])*(phi[i][jplus]-phi[i][j]) 
        + (phi[i][j] - phi[i][jminus])*(phi[i][j] - phi[i][jminus]);
}


double compute_energy(double *phi[NX], double *psi[NX], short int *xy_in[NX], int i, int j)
{
    double velocity, energy, gradientx2, gradienty2;
    int iplus, iminus, jplus, jminus;
    
    velocity = (phi[i][j] - psi[i][j]);
                    
    iplus = (i+1);   if (iplus == NX) iplus = NX-1;
    iminus = (i-1);  if (iminus == -1) iminus = 0;
    jplus = (j+1);   if (jplus == NY) jplus = NY-1;
    jminus = (j-1);  if (jminus == -1) jminus = 0;
                        
    gradientx2 = (phi[iplus][j]-phi[i][j])*(phi[iplus][j]-phi[i][j]) 
        + (phi[i][j] - phi[iminus][j])*(phi[i][j] - phi[iminus][j]);
    gradienty2 = (phi[i][jplus]-phi[i][j])*(phi[i][jplus]-phi[i][j]) 
        + (phi[i][j] - phi[i][jminus])*(phi[i][j] - phi[i][jminus]);
    if (xy_in[i][j]) return(E_SCALE*E_SCALE*(velocity*velocity + 0.5*COURANT*COURANT*(gradientx2+gradienty2)));
    else if (TWOSPEEDS) return(E_SCALE*E_SCALE*(velocity*velocity + 0.5*COURANTB*COURANTB*(gradientx2+gradienty2)));
    else return(0.0);
}


void compute_energy_flux(double *phi[NX], double *psi[NX], short int *xy_in[NX], int i, int j, double *gx, double *gy, double *arg, double *module)
/* computes energy flux given by c^2 norm(nabla u) du/dt*/
{
    double velocity, energy, gradientx, gradienty;
    int iplus, iminus, jplus, jminus;
    
    velocity = vabs(phi[i][j] - psi[i][j]);
                    
    iplus = (i+1);   if (iplus == NX) iplus = NX-1;
    iminus = (i-1);  if (iminus == -1) iminus = 0;
    jplus = (j+1);   if (jplus == NY) jplus = NY-1;
    jminus = (j-1);  if (jminus == -1) jminus = 0;
                        
    gradientx = (phi[iplus][j] - phi[iminus][j]);
    gradienty = (phi[i][jplus] - phi[i][jminus]);
    *arg = argument(gradientx,gradienty);
    if (*arg < 0.0) *arg += DPI;
    if (*arg > DPI) *arg -= DPI;
    
    if ((xy_in[i][j])||(TWOSPEEDS)) 
    {
        *module = velocity*module2(gradientx, gradienty);
        *gx = velocity*gradientx;
        *gy = velocity*gradienty;
    }
    else 
    {
        *module = 0.0;
        *gx = 0.0;
        *gy = 0.0;
    }
    
    //     if (xy_in[i][j]) return(E_SCALE*E_SCALE*(velocity*COURANT*module2(gradientx,gradienty)));
//     else if (TWOSPEEDS) return(E_SCALE*E_SCALE*(velocity*COURANTB*module2(gradientx,gradienty)));

}



double compute_variance(double *phi[NX], double *psi[NX], short int *xy_in[NX])
/* compute the variance of the field, to adjust color scheme */
{
    int i, j, n = 0;
    double variance = 0.0;

    for (i=1; i<NX; i++)
        for (j=1; j<NY; j++)
        {
            if (xy_in[i][j])
            {
                n++;
                variance += phi[i][j]*phi[i][j];
            }
        }
    if (n==0) n=1;
    return(variance/(double)n);
}


double compute_dissipation(double *phi[NX], double *psi[NX], short int *xy_in[NX], double x, double y)
{
    double velocity;
    int ij[2];
    
    xy_to_ij(x, y, ij);
    
    velocity = (phi[ij[0]][ij[1]] - psi[ij[0]][ij[1]]);
                    
    return(velocity*velocity);
}


void draw_wave(double *phi[NX], double *psi[NX], short int *xy_in[NX], double scale, int time, int plot)
/* draw the field */
{
    int i, j, iplus, iminus, jplus, jminus;
    double rgb[3], xy[2], x1, y1, x2, y2, velocity, energy, gradientx2, gradienty2;
    static double dtinverse = ((double)NX)/(COURANT*(XMAX-XMIN)), dx = (XMAX-XMIN)/((double)NX);

    glBegin(GL_QUADS);
    
//     printf("dtinverse = %.5lg\n", dtinverse);

    for (i=0; i<NX; i++)
        for (j=0; j<NY; j++)
        {
            if ((TWOSPEEDS)||(xy_in[i][j]))
            {
                if (plot == P_AMPLITUDE)
                {
                    /* make wave luminosity larger inside obstacles */
                    if (!(xy_in[i][j])) color_scheme_lum(COLOR_SCHEME, phi[i][j], scale, time, 0.7, rgb);
                    else color_scheme(COLOR_SCHEME, phi[i][j], scale, time, rgb);
                }
                else if (plot == P_ENERGY)
                    color_scheme(COLOR_SCHEME, compute_energy(phi, psi, xy_in, i, j), scale, time, rgb);
                else if (plot == P_MIXED)
                {
                    if (j > NY/2) color_scheme(COLOR_SCHEME, phi[i][j], scale, time, rgb);
                    else color_scheme(COLOR_SCHEME, compute_energy(phi, psi, xy_in, i, j), scale, time, rgb);
                }
                glColor3f(rgb[0], rgb[1], rgb[2]);

                glVertex2i(i, j);
                glVertex2i(i+1, j);
                glVertex2i(i+1, j+1);
                glVertex2i(i, j+1);
            }
        }

    glEnd ();
}

void draw_wave_e(double *phi[NX], double *psi[NX], double *total_energy[NX], double *color_scale[NX], short int *xy_in[NX], double scale, int time, int plot)
/* draw the field, new version with total energy option */
{
    int i, j, iplus, iminus, jplus, jminus;
    double rgb[3], xy[2], x1, y1, x2, y2, velocity, field_value, energy, gradientx2, gradienty2, r2;
    static double dtinverse = ((double)NX)/(COURANT*(XMAX-XMIN)), dx = (XMAX-XMIN)/((double)NX);

    glBegin(GL_QUADS);
    
//     printf("dtinverse = %.5lg\n", dtinverse);

    for (i=0; i<NX; i++)
        for (j=0; j<NY; j++)
        {
            if ((TWOSPEEDS)||(xy_in[i][j]))
            {
                switch (plot) {
                    case (P_AMPLITUDE):
                    {
                        /* make wave luminosity larger inside obstacles */
                        field_value = phi[i][j];
                        if (RESCALE_COLOR_IN_CENTER)
                        {
                            field_value *= color_scale[i][j];
                        }
//                         if (!(xy_in[i][j])) color_scheme_lum(COLOR_SCHEME, field_value, scale, time, 0.7, rgb);
//                         else 
                        color_scheme(COLOR_SCHEME, field_value, scale, time, rgb);
                        break;
                    }
                    case (P_ENERGY):
                    {
                        energy = compute_energy(phi, psi, xy_in, i, j);
                        if (RESCALE_COLOR_IN_CENTER)
                        {
                            energy *= color_scale[i][j];
                        }
                        /* adjust energy to color palette */
                        if (COLOR_PALETTE >= COL_TURBO) color_scheme_asym(COLOR_SCHEME, energy, scale, time, rgb);
                        else color_scheme(COLOR_SCHEME, energy, scale, time, rgb);
                        break;
                    }
                    case (P_MIXED):
                    {
                        if (j > NY/2) color_scheme(COLOR_SCHEME, phi[i][j], scale, time, rgb);
                        else color_scheme(COLOR_SCHEME, compute_energy(phi, psi, xy_in, i, j), scale, time, rgb);
                        break;
                    }
                    case (P_MEAN_ENERGY):
                    {
                        energy = compute_energy(phi, psi, xy_in, i, j);
                        total_energy[i][j] += energy;
                        if (COLOR_PALETTE >= COL_TURBO) 
                            color_scheme_asym(COLOR_SCHEME, total_energy[i][j]/(double)(time+1), scale, time, rgb);
                        else color_scheme(COLOR_SCHEME, total_energy[i][j]/(double)(time+1), scale, time, rgb);
                        break;
                    }
                    case (P_LOG_ENERGY):
                    {
                        energy = compute_energy(phi, psi, xy_in, i, j);
                        color_scheme(COLOR_SCHEME, LOG_SHIFT + LOG_SCALE*log(energy), scale, time, rgb);
                        break;
                    }
                    case (P_LOG_MEAN_ENERGY):
                    {
                        energy = compute_energy(phi, psi, xy_in, i, j);
                        if (energy == 0.0) energy = 1.0e-20;
                        total_energy[i][j] += energy;
                        energy = LOG_SHIFT + LOG_SCALE*log(total_energy[i][j]/(double)(time+1));
                        color_scheme(COLOR_SCHEME, energy, scale, time, rgb);
                        break;
                    }
                }
                glColor3f(rgb[0], rgb[1], rgb[2]);

                glVertex2i(i, j);
                glVertex2i(i+1, j);
                glVertex2i(i+1, j+1);
                glVertex2i(i, j+1);
            }
        }

    glEnd ();
}


void draw_wave_highres(int size, double *phi[NX], double *psi[NX], double *total_energy[NX], short int *xy_in[NX], double scale, int time, int plot)
/* draw the field, new version with total energy option */
{
    int i, j, iplus, iminus, jplus, jminus;
    double rgb[3], xy[2], x1, y1, x2, y2, velocity, energy, gradientx2, gradienty2;
    static double dtinverse = ((double)NX)/(COURANT*(XMAX-XMIN)), dx = (XMAX-XMIN)/((double)NX);

    glBegin(GL_QUADS);
    
//     printf("dtinverse = %.5lg\n", dtinverse);

    for (i=0; i<NX; i+=size)
        for (j=0; j<NY; j+=size)
        {
            if ((TWOSPEEDS)||(xy_in[i][j]))
            {
                switch (plot) {
                    case (P_AMPLITUDE):
                    {
//                         /* make wave luminosity larger inside obstacles */
//                         if (!(xy_in[i][j])) color_scheme_lum(COLOR_SCHEME, phi[i][j], scale, time, 0.7, rgb);
//                         else 
                        color_scheme(COLOR_SCHEME, phi[i][j], scale, time, rgb);
                        break;
                    }
                    case (P_ENERGY):
                    {
                        energy = compute_energy(phi, psi, xy_in, i, j);
                        /* adjust energy to color palette */
                        if (COLOR_PALETTE >= COL_TURBO) color_scheme_asym(COLOR_SCHEME, energy, scale, time, rgb);
                        else color_scheme(COLOR_SCHEME, energy, scale, time, rgb);
                        break;
                    }
                    case (P_MIXED):
                    {
                        if (j > NY/2) color_scheme(COLOR_SCHEME, phi[i][j], scale, time, rgb);
                        else color_scheme(COLOR_SCHEME, compute_energy(phi, psi, xy_in, i, j), scale, time, rgb);
                        break;
                    }
                    case (P_MEAN_ENERGY):
                    {
                        energy = compute_energy(phi, psi, xy_in, i, j);
                        total_energy[i][j] += energy;
                        if (COLOR_PALETTE >= COL_TURBO) 
                            color_scheme_asym(COLOR_SCHEME, total_energy[i][j]/(double)(time+1), scale, time, rgb);
                        else color_scheme(COLOR_SCHEME, total_energy[i][j]/(double)(time+1), scale, time, rgb);
                        break;
                    }
                    case (P_LOG_ENERGY):
                    {
                        energy = compute_energy(phi, psi, xy_in, i, j);
//                         energy = LOG_SHIFT + LOG_SCALE*log(energy);
//                         if (energy < 0.0) energy = 0.0;
                        color_scheme(COLOR_SCHEME, LOG_SHIFT + LOG_SCALE*log(energy), scale, time, rgb);
                        break;
                    }
                    case (P_LOG_MEAN_ENERGY):
                    {
                        energy = compute_energy(phi, psi, xy_in, i, j);
                        if (energy == 0.0) energy = 1.0e-20;
                        total_energy[i][j] += energy;
                        color_scheme(COLOR_SCHEME, LOG_SHIFT + LOG_SCALE*log(total_energy[i][j]/(double)(time+1)), scale, time, rgb);
                        break;
                    }
                    
                }
                glColor3f(rgb[0], rgb[1], rgb[2]);

                glVertex2i(i, j);
                glVertex2i(i+size, j);
                glVertex2i(i+size, j+size);
                glVertex2i(i, j+size);
            }
        }

    glEnd ();
}



void draw_wave_ediss(double *phi[NX], double *psi[NX], double *total_energy[NX], double *color_scale[NX], short int *xy_in[NX], 
                     double *gamma[NX], double gammamax, double scale, int time, int plot)
/* draw the field with luminosity depending on damping */
{
    int i, j, k, iplus, iminus, jplus, jminus;
    double rgb[3], xy[2], x1, y1, x2, y2, velocity, field_value, energy, gradientx2, gradienty2, r2;
    static double dtinverse = ((double)NX)/(COURANT*(XMAX-XMIN)), dx = (XMAX-XMIN)/((double)NX);

    glBegin(GL_QUADS);
    
//     printf("dtinverse = %.5lg\n", dtinverse);

    for (i=0; i<NX; i++)
        for (j=0; j<NY; j++)
        {
            if ((TWOSPEEDS)||(xy_in[i][j]))
            {
                switch (plot) {
                    case (P_AMPLITUDE):
                    {
                        /* make wave luminosity larger inside obstacles */
                        field_value = phi[i][j];
                        if (RESCALE_COLOR_IN_CENTER)
                        {
                            field_value *= color_scale[i][j];
                        }
                        if (!(xy_in[i][j])) color_scheme_lum(COLOR_SCHEME, field_value, scale, time, 0.7, rgb);
                        else color_scheme(COLOR_SCHEME, field_value, scale, time, rgb);
                        break;
                    }
                    case (P_ENERGY):
                    {
                        energy = compute_energy(phi, psi, xy_in, i, j);
                        if (RESCALE_COLOR_IN_CENTER)
                        {
                            energy *= color_scale[i][j];
                        }
                        /* adjust energy to color palette */
                        if (COLOR_PALETTE >= COL_TURBO) color_scheme_asym(COLOR_SCHEME, energy, scale, time, rgb);
                        else color_scheme(COLOR_SCHEME, energy, scale, time, rgb);
                        break;
                    }
                    case (P_MIXED):
                    {
                        if (j > NY/2) color_scheme(COLOR_SCHEME, phi[i][j], scale, time, rgb);
                        else color_scheme(COLOR_SCHEME, compute_energy(phi, psi, xy_in, i, j), scale, time, rgb);
                        break;
                    }
                    case (P_MEAN_ENERGY):
                    {
                        energy = compute_energy(phi, psi, xy_in, i, j);
                        total_energy[i][j] += energy;
                        if (COLOR_PALETTE >= COL_TURBO) 
                            color_scheme_asym(COLOR_SCHEME, total_energy[i][j]/(double)(time+1), scale, time, rgb);
                        else color_scheme(COLOR_SCHEME, total_energy[i][j]/(double)(time+1), scale, time, rgb);
                        break;
                    }
                    case (P_LOG_ENERGY):
                    {
                        energy = compute_energy(phi, psi, xy_in, i, j);
                        color_scheme(COLOR_SCHEME, LOG_SHIFT + LOG_SCALE*log(energy), scale, time, rgb);
                        break;
                    }
                    case (P_LOG_MEAN_ENERGY):
                    {
                        energy = compute_energy(phi, psi, xy_in, i, j);
                        if (energy == 0.0) energy = 1.0e-20;
                        total_energy[i][j] += energy;
                        color_scheme(COLOR_SCHEME, LOG_SHIFT + LOG_SCALE*log(total_energy[i][j]/(double)(time+1)), scale, time, rgb);
                        break;
                    }
                }
                for (k=0; k<3; k++) rgb[k]*= 1.0 - gamma[i][j]/gammamax;
                
                glColor3f(rgb[0], rgb[1], rgb[2]);

                glVertex2i(i, j);
                glVertex2i(i+1, j);
                glVertex2i(i+1, j+1);
                glVertex2i(i, j+1);
            }
        }

    glEnd ();
}


void draw_wave_highres_diss(int size, double *phi[NX], double *psi[NX], double *total_energy[NX], short int *xy_in[NX], 
                            double *gamma[NX], double gammamax, double scale, int time, int plot)
/* draw the field, new version with total energy option */
{
    int i, j, k, iplus, iminus, jplus, jminus;
    double rgb[3], xy[2], x1, y1, x2, y2, velocity, energy, gradientx2, gradienty2;
    static double dtinverse = ((double)NX)/(COURANT*(XMAX-XMIN)), dx = (XMAX-XMIN)/((double)NX);

    glBegin(GL_QUADS);
    
//     printf("dtinverse = %.5lg\n", dtinverse);

    for (i=0; i<NX; i+=size)
        for (j=0; j<NY; j+=size)
        {
            if ((TWOSPEEDS)||(xy_in[i][j]))
            {
                switch (plot) {
                    case (P_AMPLITUDE):
                    {
                        /* make wave luminosity larger inside obstacles */
                        if (!(xy_in[i][j])) color_scheme_lum(COLOR_SCHEME, phi[i][j], scale, time, 0.7, rgb);
                        else color_scheme(COLOR_SCHEME, phi[i][j], scale, time, rgb);
                        break;
                    }
                    case (P_ENERGY):
                    {
                        energy = compute_energy(phi, psi, xy_in, i, j);
                        /* adjust energy to color palette */
                        if (COLOR_PALETTE >= COL_TURBO) color_scheme_asym(COLOR_SCHEME, energy, scale, time, rgb);
                        else color_scheme(COLOR_SCHEME, energy, scale, time, rgb);
                        break;
                    }
                    case (P_MIXED):
                    {
                        if (j > NY/2) color_scheme(COLOR_SCHEME, phi[i][j], scale, time, rgb);
                        else color_scheme(COLOR_SCHEME, compute_energy(phi, psi, xy_in, i, j), scale, time, rgb);
                        break;
                    }
                    case (P_MEAN_ENERGY):
                    {
                        energy = compute_energy(phi, psi, xy_in, i, j);
                        total_energy[i][j] += energy;
                        if (COLOR_PALETTE >= COL_TURBO) 
                            color_scheme_asym(COLOR_SCHEME, total_energy[i][j]/(double)(time+1), scale, time, rgb);
                        else color_scheme(COLOR_SCHEME, total_energy[i][j]/(double)(time+1), scale, time, rgb);
                        break;
                    }
                    case (P_LOG_ENERGY):
                    {
                        energy = compute_energy(phi, psi, xy_in, i, j);
//                         energy = LOG_SHIFT + LOG_SCALE*log(energy);
//                         if (energy < 0.0) energy = 0.0;
                        color_scheme(COLOR_SCHEME, LOG_SHIFT + LOG_SCALE*log(energy), scale, time, rgb);
                        break;
                    }
                    case (P_LOG_MEAN_ENERGY):
                    {
                        energy = compute_energy(phi, psi, xy_in, i, j);
                        if (energy == 0.0) energy = 1.0e-20;
                        total_energy[i][j] += energy;
                        color_scheme(COLOR_SCHEME, LOG_SHIFT + LOG_SCALE*log(total_energy[i][j]/(double)(time+1)), scale, time, rgb);
                        break;
                    }
                    
                }
                for (k=0; k<3; k++) rgb[k]*= 1.0 - gamma[i][j]/gammamax;
                
                glColor3f(rgb[0], rgb[1], rgb[2]);

                glVertex2i(i, j);
                glVertex2i(i+size, j);
                glVertex2i(i+size, j+size);
                glVertex2i(i, j+size);
            }
        }

    glEnd ();
}


void draw_wave_epalette(double *phi[NX], double *psi[NX], double *total_energy[NX], double *total_flux,
                        double *color_scale[NX], 
                        short int *xy_in[NX], double scale, int time, int plot, int palette, int fade, double fade_value)
/* same as draw_wave_e, but with color scheme specification */
{
    int i, j, k, iplus, iminus, jplus, jminus;
    double rgb[3], xy[2], x1, y1, x2, y2, velocity, field_value, energy, gradientx2, gradienty2, r2, arg, mod, flux_factor, gx, gy, mgx, mgy, ffactor;
    static double dtinverse = ((double)NX)/(COURANT*(XMAX-XMIN)), dx = (XMAX-XMIN)/((double)NX);

    glBegin(GL_QUADS);
    
//     printf("dtinverse = %.5lg\n", dtinverse);

    for (i=0; i<NX; i++)
        for (j=0; j<NY; j++)
        {
            if ((TWOSPEEDS)||(xy_in[i][j]))
            {
                switch (plot) {
                    case (P_AMPLITUDE):
                    {
                        /* make wave luminosity larger inside obstacles */
                        field_value = phi[i][j];
                        if (RESCALE_COLOR_IN_CENTER)
                        {
                            field_value *= color_scale[i][j];
                        }
//                         if (!(xy_in[i][j])) color_scheme_lum(COLOR_SCHEME, field_value, scale, time, 0.7, rgb);
//                         else 
                        color_scheme_palette(COLOR_SCHEME, palette, field_value, scale, time, rgb);
                        break;
                    }
                    case (P_ENERGY):
                    {
                        energy = compute_energy(phi, psi, xy_in, i, j);
                        if (RESCALE_COLOR_IN_CENTER)
                        {
                            energy *= color_scale[i][j];
                        }
                        /* adjust energy to color palette */
                        if (COLOR_PALETTE >= COL_TURBO) color_scheme_asym_palette(COLOR_SCHEME, palette, energy, scale, time, rgb);
                        else color_scheme_palette(COLOR_SCHEME, palette, energy, scale, time, rgb);
                        break;
                    }
                    case (P_MIXED):
                    {
                        if (j > NY/2) color_scheme_palette(COLOR_SCHEME, palette, phi[i][j], scale, time, rgb);
                        else color_scheme_palette(COLOR_SCHEME, palette, compute_energy(phi, psi, xy_in, i, j), scale, time, rgb);
                        break;
                    }
                    case (P_MEAN_ENERGY):
                    {
                        energy = compute_energy(phi, psi, xy_in, i, j);
                        total_energy[i][j] += energy;
                        if (COLOR_PALETTE >= COL_TURBO) 
                            color_scheme_asym_palette(COLOR_SCHEME, palette, total_energy[i][j]/(double)(time+1), scale, time, rgb);
                        else color_scheme_palette(COLOR_SCHEME, palette, total_energy[i][j]/(double)(time+1), scale, time, rgb);
                        break;
                    }
                    case (P_LOG_ENERGY):
                    {
                        energy = compute_energy(phi, psi, xy_in, i, j);
                        color_scheme_palette(COLOR_SCHEME, palette, LOG_SHIFT + LOG_SCALE*log(energy), scale, time, rgb);
                        break;
                    }
                    case (P_LOG_MEAN_ENERGY):
                    {
                        energy = compute_energy(phi, psi, xy_in, i, j);
                        if (energy == 0.0) energy = 1.0e-20;
                        total_energy[i][j] += energy;
                        energy = LOG_SHIFT + LOG_SCALE*log(total_energy[i][j]/(double)(time+1));
                        color_scheme_palette(COLOR_SCHEME, palette, energy, scale, time, rgb);
                        break;
                    }
                    case (P_ENERGY_FLUX):
                    {
                        compute_energy_flux(phi, psi, xy_in, i, j, &gx, &gy, &arg, &mod);
//                         color_scheme_palette(C_ONEDIM_LINEAR, palette, arg/DPI, 1.0, 1, rgb);
//                         flux_factor = tanh(mod*E_SCALE);
//                         for (k=0; k<3; k++) rgb[k] *= flux_factor;
                        
                        color_scheme_asym_palette(COLOR_SCHEME, palette, mod*FLUX_SCALE, scale, time, rgb);
                        break;
                    }
                    case (P_TOTAL_ENERGY_FLUX):
                    {
//                         ffactor = 1.0;
                        compute_energy_flux(phi, psi, xy_in, i, j, &gx, &gy, &arg, &mod);
                        total_flux[2*NX*NY + 2*j*NX + 2*i] += gx;
                        total_flux[2*NX*NY + 2*j*NX + 2*i + 1] += gy;
                        total_flux[2*j*NX + 2*i] += total_flux[2*NX*NY + 2*j*NX + 2*i];
                        total_flux[2*j*NX + 2*i + 1] += total_flux[2*NX*NY + 2*j*NX + 2*i + 1];
//                         total_flux[2*j*NX + 2*i] *= 1.0 + 1.0/(double)(time+1);
//                         total_flux[2*j*NX + 2*i + 1] *= 1.0 + 1.0/(double)(time+1);
//                         total_flux[2*j*NX + 2*i] *= ffactor;
//                         total_flux[2*j*NX + 2*i + 1] *= ffactor;
//                         total_flux[2*j*NX + 2*i] += gx;
//                         total_flux[2*j*NX + 2*i + 1] += gy;
//                         total_flux[2*j*NX + 2*i] *= 1.0/ffactor;
//                         total_flux[2*j*NX + 2*i + 1] *= 1.0/ffactor;
//                         total_flux[2*j*NX + 2*i] += 0.1*gx;
//                         total_flux[2*j*NX + 2*i + 1] += 0.1*gy;
                        mgx = total_flux[2*j*NX + 2*i];
                        mgy = total_flux[2*j*NX + 2*i + 1];
//                         mgx = total_flux[2*j*NX + 2*i]/sqrt((double)(time+1));
//                         mgy = total_flux[2*j*NX + 2*i + 1]/sqrt((double)(time+1));
//                         mgx = total_flux[2*j*NX + 2*i]/(1.0 + 0.1*log((double)(time+2)));
//                         mgy = total_flux[2*j*NX + 2*i + 1]/(1.0 + 0.1*log((double)(time+2)));
                        mod = module2(mgx, mgy);
                        arg = argument(mgx, mgy);
                        if (arg < 0.0) arg += DPI;
                        color_scheme_palette(C_ONEDIM_LINEAR, palette, arg/DPI, 1.0, 1, rgb);
                        flux_factor = tanh(mod*FLUX_SCALE);
                        for (k=0; k<3; k++) rgb[k] *= flux_factor;
                        
//                         color_scheme_asym_palette(COLOR_SCHEME, palette, mod, scale, time, rgb);
                        break;
                    }
                }
                if (fade) for (k=0; k<3; k++) rgb[k] *= fade_value;
                glColor3f(rgb[0], rgb[1], rgb[2]);

                glVertex2i(i, j);
                glVertex2i(i+1, j);
                glVertex2i(i+1, j+1);
                glVertex2i(i, j+1);
            }
        }

    glEnd ();
}


void draw_wave_highres_palette(int size, double *phi[NX], double *psi[NX], double *total_energy[NX], double *total_flux, short int *xy_in[NX], double scale, int time, int plot, int palette, int fade, double fade_value)
/* same as draw_wave_highres, but with color scheme option */
{
    int i, j, k, iplus, iminus, jplus, jminus;
    double rgb[3], xy[2], x1, y1, x2, y2, velocity, energy, gradientx2, gradienty2, arg, mod, flux_factor, gx, gy, mgx, mgy;
    static double dtinverse = ((double)NX)/(COURANT*(XMAX-XMIN)), dx = (XMAX-XMIN)/((double)NX);

    glBegin(GL_QUADS);
    
//     printf("dtinverse = %.5lg\n", dtinverse);

    for (i=0; i<NX; i+=size)
        for (j=0; j<NY; j+=size)
        {
            if ((TWOSPEEDS)||(xy_in[i][j]))
            {
                switch (plot) {
                    case (P_AMPLITUDE):
                    {
//                         /* make wave luminosity larger inside obstacles */
//                         if (!(xy_in[i][j])) color_scheme_lum(COLOR_SCHEME, phi[i][j], scale, time, 0.7, rgb);
//                         else 
                        color_scheme_palette(COLOR_SCHEME, palette, phi[i][j], scale, time, rgb);
                        break;
                    }
                    case (P_ENERGY):
                    {
                        energy = compute_energy(phi, psi, xy_in, i, j);
                        /* adjust energy to color palette */
                        if (COLOR_PALETTE >= COL_TURBO) color_scheme_asym_palette(COLOR_SCHEME, palette, energy, scale, time, rgb);
                        else color_scheme_palette(COLOR_SCHEME, palette, energy, scale, time, rgb);
                        break;
                    }
                    case (P_MIXED):
                    {
                        if (j > NY/2) color_scheme_palette(COLOR_SCHEME, palette, phi[i][j], scale, time, rgb);
                        else color_scheme_palette(COLOR_SCHEME, palette, compute_energy(phi, psi, xy_in, i, j), scale, time, rgb);
                        break;
                    }
                    case (P_MEAN_ENERGY):
                    {
                        energy = compute_energy(phi, psi, xy_in, i, j);
                        total_energy[i][j] += energy;
                        if (COLOR_PALETTE >= COL_TURBO) 
                            color_scheme_asym_palette(COLOR_SCHEME, palette, total_energy[i][j]/(double)(time+1), scale, time, rgb);
                        else color_scheme_palette(COLOR_SCHEME, palette, total_energy[i][j]/(double)(time+1), scale, time, rgb);
                        break;
                    }
                    case (P_LOG_ENERGY):
                    {
                        energy = compute_energy(phi, psi, xy_in, i, j);
//                         energy = LOG_SHIFT + LOG_SCALE*log(energy);
//                         if (energy < 0.0) energy = 0.0;
                        color_scheme_palette(COLOR_SCHEME, palette, LOG_SHIFT + LOG_SCALE*log(energy), scale, time, rgb);
                        break;
                    }
                    case (P_LOG_MEAN_ENERGY):
                    {
                        energy = compute_energy(phi, psi, xy_in, i, j);
                        if (energy == 0.0) energy = 1.0e-20;
                        total_energy[i][j] += energy;
                        color_scheme_palette(COLOR_SCHEME, palette, LOG_SHIFT + LOG_SCALE*log(total_energy[i][j]/(double)(time+1)), scale, time, rgb);
                        break;
                    }
                    case (P_ENERGY_FLUX):
                    {
                        compute_energy_flux(phi, psi, xy_in, i, j, &gx, &gy, &arg, &mod);
//                         color_scheme_palette(C_ONEDIM_LINEAR, palette, arg/DPI, 1.0, 1, rgb);
//                         flux_factor = tanh(mod*E_SCALE);
//                         for (k=0; k<3; k++) rgb[k] *= flux_factor;
                        color_scheme_asym_palette(COLOR_SCHEME, palette, mod*FLUX_SCALE, scale, time, rgb);
                        break;
                    }
                    case (P_TOTAL_ENERGY_FLUX):
                    {
                        compute_energy_flux(phi, psi, xy_in, i, j, &gx, &gy, &arg, &mod);
                        total_flux[2*j*NX + 2*i] *= 0.99;
                        total_flux[2*j*NX + 2*i + 1] *= 0.99;
                        total_flux[2*j*NX + 2*i] += gx;
                        total_flux[2*j*NX + 2*i + 1] += gy;
//                         mgx = total_flux[2*j*NX + 2*i]/(double)(time+1);
//                         mgy = total_flux[2*j*NX + 2*i + 1]/(double)(time+1);
                        mgx = total_flux[2*j*NX + 2*i];
                        mgy = total_flux[2*j*NX + 2*i + 1];
//                         mgx = total_flux[2*j*NX + 2*i]/log((double)(time+2));
//                         mgy = total_flux[2*j*NX + 2*i + 1]/log((double)(time+2));
                        mod = module2(mgx, mgy);
                        arg = argument(mgx, mgy);
                        if (arg < 0.0) arg += DPI;
                        color_scheme_palette(C_ONEDIM_LINEAR, palette, arg/DPI, 1.0, 1, rgb);
                        flux_factor = tanh(mod*E_SCALE);
                        for (k=0; k<3; k++) rgb[k] *= flux_factor;
                        break;
                    }
                    
                }
                if (fade) for (k=0; k<3; k++) rgb[k] *= fade_value;
                glColor3f(rgb[0], rgb[1], rgb[2]);
                
                glVertex2i(i, j);
                glVertex2i(i+size, j);
                glVertex2i(i+size, j+size);
                glVertex2i(i, j+size);
            }
        }

    glEnd ();
}

/* modified function for "flattened" wave tables */

void init_circular_wave_mod(double x, double y, double phi[NX*NY], double psi[NX*NY], short int xy_in[NX*NY])
/* initialise field with drop at (x,y) - phi is wave height, psi is phi at time t-1 */
{
    int i, j;
    double xy[2], dist2;

    printf("Initializing wave\n"); 
//     #pragma omp parallel for private(i,j,xy,dist2)
    for (i=0; i<NX; i++)
    {
        if (i%100 == 0) printf("Initializing column %i of %i\n", i, NX);
        for (j=0; j<NY; j++)
        {
//             printf("i*NY+j = %i\n", i*NY+j);
            ij_to_xy(i, j, xy);
            dist2 = (xy[0]-x)*(xy[0]-x) + (xy[1]-y)*(xy[1]-y);
	    xy_in[i*NY+j] = xy_in_billiard(xy[0],xy[1]);
            
	    if ((xy_in[i*NY+j])||(TWOSPEEDS)) 
                phi[i*NY+j] = INITIAL_AMP*exp(-dist2/INITIAL_VARIANCE)*cos(-sqrt(dist2)/INITIAL_WAVELENGTH);
            else phi[i*NY+j] = 0.0;
            psi[i*NY+j] = 0.0;
        }
    }
}

void add_circular_wave_mod(double factor, double x, double y, double phi[NX*NY], double psi[NX*NY], short int xy_in[NX*NY])
/* add drop at (x,y) to the field with given prefactor */
{
    int i, j;
    double xy[2], dist2;

    #pragma omp parallel for private(i,j,xy,dist2)
    for (i=0; i<NX; i++)
        for (j=0; j<NY; j++)
        {
            ij_to_xy(i, j, xy);
            dist2 = (xy[0]-x)*(xy[0]-x) + (xy[1]-y)*(xy[1]-y);
            if ((xy_in[i*NY+j])||(TWOSPEEDS)) 
                phi[i*NY+j] += INITIAL_AMP*factor*exp(-dist2/INITIAL_VARIANCE)*cos(-sqrt(dist2)/INITIAL_WAVELENGTH);
        }
}

void init_wave_flat_mod(double phi[NX*NY], double psi[NX*NY], short int xy_in[NX*NY])
/* initialise flat field - phi is wave height, psi is phi at time t-1 */
{
    int i, j;
    double xy[2];

    #pragma omp parallel for private(i,j,xy)
    for (i=0; i<NX; i++) {
        if (i%100 == 0) printf("Wave and table xy_in - Initialising column %i of %i\n", i, NX);
        for (j=0; j<NY; j++)
        {
            ij_to_xy(i, j, xy);
	    xy_in[i*NY+j] = xy_in_billiard(xy[0],xy[1]);
	    phi[i*NY+j] = 0.0;
            psi[i*NY+j] = 0.0;
        }
    }
}


double compute_variance_mod(double phi[NX*NY], double psi[NX*NY], short int xy_in[NX*NY])
/* compute the variance of the field, to adjust color scheme */
{
    int i, j, n = 0;
    double variance = 0.0;

    #pragma omp parallel for private(i,j,variance)
    for (i=1; i<NX; i++)
        for (j=1; j<NY; j++)
        {
            if (xy_in[i*NY+j])
            {
                n++;
                variance += phi[i*NY+j]*phi[i*NY+j];
            }
        }
    if (n==0) n=1;
    return(variance/(double)n);
}

double compute_energy_mod(double phi[NX*NY], double psi[NX*NY], short int xy_in[NX*NY], int i, int j)
{
    double velocity, energy, gradientx2, gradienty2;
    int iplus, iminus, jplus, jminus;
    static int i1, j1, i2, j2, i3, j3, ij[2];
    static int first = 1;
    
    if (first)
    {
        xy_to_ij(-LAMBDA, -1.0, ij);
        i1 = ij[0] + 1;     j1 = ij[1] + 1;
        xy_to_ij(0.0, 0.0, ij);
        i2 = ij[0] - 1;     j2 = ij[1] - 1;
        xy_to_ij(LAMBDA, 1.0, ij);
        i3 = ij[0] - 1;     j3 = ij[1] - 1;
        first = 0;
    }
    
    velocity = (phi[i*NY+j] - psi[i*NY+j]);
    
/* avoid computing gradient across boundary of compared regions */
    if (COMPARISON)
    {
        iplus = (i+1);   if (iplus == NX) iplus = NX-1;
        iminus = (i-1);  if (iminus == -1) iminus = 0;
        jplus = (j+1);   if (jplus == NY) jplus = NY-1;     else if (jplus == NY/2) jplus = NY/2-1;
        jminus = (j-1);  if (jminus == -1) jminus = 0;      else if (jminus == NY/2-1) jminus = NY/2;
    }
    else
    {
        iplus = (i+1);   if (iplus == NX) iplus = NX-1;
        iminus = (i-1);  if (iminus == -1) iminus = 0;
        jplus = (j+1);   if (jplus == NY) jplus = NY-1;
        jminus = (j-1);  if (jminus == -1) jminus = 0;
    }
//     else
//     {
//         iplus = i+1;   
//         if ((j<j2)&&(iplus>i3-1)) iplus = i1;
//         else if ((j>=j2)&&(iplus>i2-1)) iplus = i1;
//         
//         iminus = i-1;   
//         if ((j<j2)&&(iminus<i1)) iminus = i3-1;
//         else if ((j>=j2)&&(iminus<i1)) iminus = i2-1;
//         
//         jplus = j+1;
//         if ((i<i2)&&(jplus>j3-1)) jplus = j1;
//         else if ((i>=i2)&&(jplus>j2-1)) jplus = j1;
//         
//         jminus = j-1;   
//         if ((i<i2)&&(jminus<j1)) jminus = j3-1;
//         else if ((i>=i2)&&(jminus<j1)) jminus = j2-1;
//         
//         jminus = j-1;
//     }
                
    if (B_COND == BC_LSHAPE)
    {
        if ((i<=i1)||(j<=j1)) return(0.0);
        if ((i>=i2-1)&&(j>=j2-1)) return(0.0);
        if ((i>=i3-1)||(j>=j3-1)) return(0.0);
    }
    
    
    gradientx2 = (phi[iplus*NY+j]-phi[i*NY+j])*(phi[iplus*NY+j]-phi[i*NY+j]) 
        + (phi[i*NY+j] - phi[iminus*NY+j])*(phi[i*NY+j] - phi[iminus*NY+j]);
    gradienty2 = (phi[i*NY+jplus]-phi[i*NY+j])*(phi[i*NY+jplus]-phi[i*NY+j]) 
        + (phi[i*NY+j] - phi[i*NY+jminus])*(phi[i*NY+j] - phi[i*NY+jminus]);
    if (xy_in[i*NY+j]) return(E_SCALE*E_SCALE*(velocity*velocity + 0.5*COURANT*COURANT*(gradientx2+gradienty2)));
    else if (TWOSPEEDS) return(E_SCALE*E_SCALE*(velocity*velocity + 0.5*COURANTB*COURANTB*(gradientx2+gradienty2)));
    else return(0.0);
}

void compute_energy_flux_mod(double phi[NX*NY], double psi[NX*NY], short int xy_in[NX*NY], int i, int j, double *gx, double *gy, double *arg, double *module)
/* computes energy flux given by c^2 norm(nabla u) du/dt*/
{
    double velocity, energy, gradientx, gradienty, max = 1.0e5, current_mod, current_arg;
    int iplus, iminus, jplus, jminus;
    
    if ((i == 0)||(i == NX-1)||(j == 0)||(j == NY-1))
    {
        current_mod = 0.0;
        current_arg = PI;
        *gx = 0.0;
        *gy = 0.0;
    }
    else if ((xy_in[i*NY+j])||(TWOSPEEDS)) 
    {
        velocity = vabs(phi[i*NY+j] - psi[i*NY+j]);
                    
//         iplus = (i+1);   /*if (iplus == NX) iplus = NX-1;*/
//         iminus = (i-1);  /*if (iminus == -1) iminus = 0;*/
//         jplus = (j+1);   /*if (jplus == NY) jplus = NY-1;*/
//         jminus = (j-1);  /*if (jminus == -1) jminus = 0;*/
                        
        gradientx = (phi[(i+1)*NY+j] - phi[(i-1)*NY+j]);
        gradienty = (phi[i*NY+j+1] - phi[i*NY+j-1]);
    
        if (gradientx > max) gradientx = max;
        else if (gradientx < -max) gradientx = -max;
        if (gradienty > max) gradienty = max;
        else if (gradienty < -max) gradienty = -max;
    
        current_mod = velocity*module2(gradientx, gradienty);
        if (current_mod > 1.0e-10)
        {
            current_arg = argument(gradientx,gradienty);
            if (current_arg < 0.0) current_arg += DPI;
            if (current_arg >= DPI) current_arg -= DPI;
        }
        else current_arg = PI;
        *gx = velocity*gradientx;
        *gy = velocity*gradienty;
    }
    else 
    {
        current_mod = 0.0;
        current_arg = PI;
        *gx = 0.0;
        *gy = 0.0;
    }
    
    *module = current_mod;
    *arg = current_arg;
}


double compute_phase(double phi[NX*NY], double psi[NX*NY], short int xy_in[NX*NY], int i, int j)
{
    double velocity, angle;
    
    velocity = (phi[i*NY+j] - psi[i*NY+j]);
            
    if (module2(phi[i*NY+j], velocity) < 1.0e-10) return(0.0); 
    else if (xy_in[i*NY+j]) 
    {
        angle = argument(phi[i*NY+j], PHASE_FACTOR*velocity/COURANT);
        if (angle < 0.0) angle += DPI;
        
        if ((i==NY/2)&&(j==NY/2)) printf("Phase = %.3lg Pi\n", angle/PI);
        return(angle);
    }
    else if (TWOSPEEDS) 
    {
        angle = argument(phi[i*NY+j], PHASE_FACTOR*velocity/COURANTB);
        if (angle < 0.0) angle += DPI;
        return(angle);
    }
    else return(0.0);
}

