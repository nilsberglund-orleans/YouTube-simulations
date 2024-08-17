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

void init_xy_tcc_in_xrange(short int * xy_in[NX], double *tcc_table[NX], double *tgamma_table[NX], int imin, int imax)
/* initialise table xy_in, needed when obstacles are killed */
// short int * xy_in[NX];
//  
{
    int i, j;
    double xy[2];

//     for (i=imin; i<imax; i++)
    for (i=NX/2; i<NX; i++)
        for (j=0; j<NY; j++)
        {
            ij_to_xy(i, j, xy);
	    xy_in[i][j] = xy_in_billiard(xy[0],xy[1]);
//             xy_in[i][j] = 1;
            if (xy_in[i][j] == 0) 
            {
                tcc_table[i][j] = COURANT*COURANT;
                tgamma_table[i][j] = GAMMA;
            }
                else if (xy_in[i][j] == 2)
            {
                tcc_table[i][j] = 0.0;
                tgamma_table[i][j] = 1.0;                            
            }
                else
            {
                tcc_table[i][j] = COURANTB*COURANTB;
                tgamma_table[i][j] = GAMMAB;
            }
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

void damp_circular_wave(double factor, int range, double x, double y, double *phi[NX], double *psi[NX])
/* damp field around (x,y) - phi is wave height, psi is phi at time t-1 */
{
    int i, j, ij[2], imin, imax, jmin, jmax;
    double r2, v;

//     printf("Damping wave\n");
    xy_to_ij(x, y, ij);
    imin = ij[0] - range;
    if (imin < 0) imin = 0;
    imax = ij[0] + range;
    if (imax > NX-1) imax = NX-1;
    jmin = ij[1] - range;
    if (jmin < 0) jmin = 0;
    jmax = ij[1] + range;
    if (jmax > NY-1) jmax = NY-1;
    v = (double)(range*range);
    
    for (i=imin; i<imax; i++)
    {
        for (j=jmin; j<jmax; j++)
        {
            r2 = (double)((i - ij[0])*(i - ij[0]) + (j - ij[1])*(j - ij[1]));
            factor *= 1.0 - exp(-0.5*r2/v);
            phi[i][j] *= factor;
            psi[i][j] *= factor;
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


void init_planar_travelling_wave(double x, double y, double vx, double *phi[NX], double *psi[NX], short int * xy_in[NX])
/* initialise field with drop at (x,y) - phi is wave height, psi is phi at time t-1 */
/* beta version, works for vertical planar wave only so far */
{
    int i, j, ij[2], jmin, jmax;
    double xy[2], dist2, dist2_shifted, xshift, y1;

    xy_to_ij(x, y, ij);
    jmax = ij[1];
    jmin = NY - jmax;

    xshift = vx*COURANT/(double)NX;
    for (i=0; i<NX; i++)
    {
        for (j=0; j<NY; j++)
//         for (j=jmin; j<jmax; j++)
        {
            ij_to_xy(i, j, xy);
            dist2 = (xy[0]-x)*(xy[0]-x);
	    dist2_shifted = (xy[0]-x-xshift)*(xy[0]-x-xshift);
	    xy_in[i][j] = xy_in_billiard(xy[0],xy[1]);
            
            if ((xy_in[i][j])||(TWOSPEEDS)) 
            {
                phi[i][j] = INITIAL_AMP*exp(-dist2/INITIAL_VARIANCE)*cos(-sqrt(dist2)/INITIAL_WAVELENGTH);
                psi[i][j] = INITIAL_AMP*exp(-dist2_shifted/INITIAL_VARIANCE)*cos(-sqrt(dist2_shifted)/INITIAL_WAVELENGTH);
            }
            else 
            {   
                phi[i][j] = 0.0;
                psi[i][j] = 0.0;
            }
        }
        for (j=0; j<jmin; j++)
        {   
            dist2 = (double)(jmin-1-j)*(YMAX-YMIN)/(double)NY;
            dist2 = 500.0*dist2*dist2;
            phi[i][j] *= exp(-dist2/INITIAL_VARIANCE);
            psi[i][j] *= exp(-dist2/INITIAL_VARIANCE);
        }
        for (j=jmax; j<NY; j++)
        {   
            phi[i][j] = phi[i][NY-1-j];
            psi[i][j] = psi[i][NY-1-j];
        }
    }
}

void add_planar_travelling_wave(double x, double y, double vx, double wavelength_factor, double *phi[NX], double *psi[NX], short int * xy_in[NX])
/* initialise field with drop at (x,y) - phi is wave height, psi is phi at time t-1 */
/* beta version, works for vertical planar wave only so far */
{
    int i, j;
    double xy[2], dist2, dist2_shifted, xshift, wavelength;

    xshift = vx*COURANT/(double)NX;
    wavelength = INITIAL_WAVELENGTH*wavelength_factor;
    for (i=0; i<NX; i++)
    {
        for (j=0; j<NY; j++)
        {
            ij_to_xy(i, j, xy);
            dist2 = (xy[0]-x)*(xy[0]-x);
	    dist2_shifted = (xy[0]-x-xshift)*(xy[0]-x-xshift);
	    
            if ((xy_in[i][j])||(TWOSPEEDS)) 
            {
                phi[i][j] += INITIAL_AMP*exp(-dist2/INITIAL_VARIANCE)*cos(-sqrt(dist2)/wavelength);
                psi[i][j] += INITIAL_AMP*exp(-dist2_shifted/INITIAL_VARIANCE)*cos(-sqrt(dist2_shifted)/wavelength);
            }
            else 
            {   
                phi[i][j] = 0.0;
                psi[i][j] = 0.0;
            }
        }
    }
}
void init_wave_flat( double *phi[NX], double *psi[NX], short int * xy_in[NX])
/* initialise flat field - phi is wave height, psi is phi at time t-1 */
{
    int i, j;
    double xy[2], dist2;

    if (!XYIN_INITIALISED) 
        for (i=0; i<NX; i++) 
        {
            if (i%100 == 0) printf("Table xy_in - Initialising column %i of %i\n", i, NX);
            for (j=0; j<NY; j++)
            {
                ij_to_xy(i, j, xy);
                xy_in[i][j] = xy_in_billiard(xy[0],xy[1]);
            }
        }

    for (i=0; i<NX; i++) {
        if (i%100 == 0) printf("Wave - Initialising column %i of %i\n", i, NX);
        for (j=0; j<NY; j++)
        {
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


void add_vertical_wave(double factor, double x, double y1, double y2, double *phi[NX], double *psi[NX], short int * xy_in[NX])
/* add wave around vertical segment (x,y1)-(x,y2) to the field with given prefactor */
{
    int i, j;
    double xy[2], dist2;

    for (i=0; i<NX; i++)
        for (j=0; j<NY; j++)
        {
            ij_to_xy(i, j, xy);
            if (xy[1] < y1) dist2 = (xy[0]-x)*(xy[0]-x) + (xy[1]-y1)*(xy[1]-y1);
            else if (xy[1] > y2) dist2 = (xy[0]-x)*(xy[0]-x) + (xy[1]-y2)*(xy[1]-y2);
            else dist2 = (xy[0]-x)*(xy[0]-x);
            if ((xy_in[i][j])||(TWOSPEEDS)) 
                phi[i][j] += INITIAL_AMP*factor*exp(-dist2/INITIAL_VARIANCE)*cos(-sqrt(dist2)/INITIAL_WAVELENGTH);
        }
}


void damp_vertical_wave(double factor, int range, double x, double y1, double y2, double *phi[NX], double *psi[NX])
/* damp field around (x,y) - phi is wave height, psi is phi at time t-1 */
{
    int i, j, ij1[2], ij2[2], imin, imax, jmin, jmax;
    double r2, v;

    xy_to_ij(x, y1, ij1);
    xy_to_ij(x, y2, ij2);
    imin = ij1[0] - range;
    if (imin < 0) imin = 0;
    imax = ij2[0] + range;
    if (imax > NX-1) imax = NX-1;
    jmin = ij1[1] - range;
    if (jmin < 0) jmin = 0;
    jmax = ij2[1] + range;
    if (jmax > NY-1) jmax = NY-1;
    v = (double)(range*range);
    
    for (i=imin; i<imax; i++)
    {
        for (j=jmin; j<jmax; j++)
        {
            if (j < ij1[1]) r2 = (double)((i - ij1[0])*(i - ij1[0]) + (j - ij1[1])*(j - ij1[1]));
            else if (j > ij2[1]) r2 = (double)((i - ij2[0])*(i - ij2[0]) + (j - ij2[1])*(j - ij2[1]));
            else r2 = (double)((i - ij1[0])*(i - ij1[0]));
            factor *= 1.0 - exp(-0.5*r2/v);
            phi[i][j] *= factor;
            psi[i][j] *= factor;
        }
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


void draw_wave_epalette(double *phi[NX], double *psi[NX], double *total_energy[NX], double *average_energy[NX], double *total_flux, double *color_scale[NX], short int *xy_in[NX], double scale, int time, int plot, int palette, int fade, double fade_value)
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


double wave_value(int i, int j, double *phi[NX], double *psi[NX], double *total_energy[NX], double *average_energy[NX], double *total_flux, double *color_scale[NX], short int *xy_in[NX], double scale, int time, int plot, int palette, double rgb[3])
/* compute value of wave and color */
{
    int k;
    double value, velocity, energy, gradientx2, gradienty2, arg, mod, flux_factor, gx, gy, mgx, mgy;
    
    switch (plot) {
        case (P_AMPLITUDE):
        {
            value = phi[i][j];
            color_scheme_palette(COLOR_SCHEME, palette, VSHIFT_AMPLITUDE + VSCALE_AMPLITUDE*value, scale, time, rgb);
            break;
        }
        case (P_ENERGY):
        {
            value = compute_energy(phi, psi, xy_in, i, j);
            /* adjust energy to color palette */
            if (COLOR_PALETTE >= COL_TURBO) color_scheme_asym_palette(COLOR_SCHEME, palette, value, scale, time, rgb);
            else color_scheme_palette(COLOR_SCHEME, palette, value, scale, time, rgb);
            break;
        }
        case (P_MIXED):
        {
            if (j > NY/2) value = phi[i][j];
            else value = compute_energy(phi, psi, xy_in, i, j);
            color_scheme_palette(COLOR_SCHEME, palette, value, scale, time, rgb);
            break;
        }
        case (P_MEAN_ENERGY):
        {
            energy = compute_energy(phi, psi, xy_in, i, j);
            total_energy[i][j] += energy;
            value = total_energy[i][j]/(double)(time+1);
            if (COLOR_PALETTE >= COL_TURBO) 
                color_scheme_asym_palette(COLOR_SCHEME, palette, value, scale, time, rgb);
            else color_scheme_palette(COLOR_SCHEME, palette, value, scale, time, rgb);
            break;
        }
        case (P_LOG_ENERGY):
        {
            energy = compute_energy(phi, psi, xy_in, i, j);
            if (energy < 1.0e-20) energy = 1.0e-20;
            value = LOG_SHIFT + LOG_SCALE*log(energy);
            color_scheme_palette(COLOR_SCHEME, palette, value, scale, time, rgb);
            break;
        }
        case (P_LOG_MEAN_ENERGY):
        {
            energy = compute_energy(phi, psi, xy_in, i, j);
            if (energy == 0.0) energy = 1.0e-20;
            total_energy[i][j] += energy;
            value = LOG_SHIFT + LOG_SCALE*log(total_energy[i][j]/(double)(time+1));
            color_scheme_palette(COLOR_SCHEME, palette, value, scale, time, rgb);
            break;
        }
        case (P_ENERGY_FLUX):
        {
            compute_energy_flux(phi, psi, xy_in, i, j, &gx, &gy, &arg, &mod);
//                         color_scheme_palette(C_ONEDIM_LINEAR, palette, arg/DPI, 1.0, 1, rgb);
//                         flux_factor = tanh(mod*E_SCALE);
//                         for (k=0; k<3; k++) rgb[k] *= flux_factor;
            value = mod*FLUX_SCALE;
            color_scheme_asym_palette(COLOR_SCHEME, palette, value, scale, time, rgb);
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
            value = arg/DPI;
            color_scheme_palette(C_ONEDIM_LINEAR, palette, value, 1.0, 1, rgb);
//             flux_factor = tanh(mod*E_SCALE);
            flux_factor = tanh(mod*FLUX_SCALE);
            for (k=0; k<3; k++) rgb[k] *= flux_factor;
            break;
        }
        case (P_AVERAGE_ENERGY):
        {
            energy = compute_energy(phi, psi, xy_in, i, j);
            average_energy[i][j] = AVRG_E_FACTOR*average_energy[i][j] + (1.0 - AVRG_E_FACTOR)*energy;
            value = average_energy[i][j];
            if (COLOR_PALETTE >= COL_TURBO) 
                color_scheme_asym_palette(COLOR_SCHEME, palette, value, scale, time, rgb);
            else color_scheme_palette(COLOR_SCHEME, palette, value, scale, time, rgb);
            break;
        }
        case (P_LOG_AVERAGE_ENERGY):
        {
            energy = compute_energy(phi, psi, xy_in, i, j);
            average_energy[i][j] = AVRG_E_FACTOR*average_energy[i][j] + (1.0 - AVRG_E_FACTOR)*energy;
            if (average_energy[i][j] == 0.0) value = LOG_SHIFT - 10.0*LOG_SCALE;
            else value = LOG_SHIFT + LOG_SCALE*log(average_energy[i][j]);
            color_scheme_palette(COLOR_SCHEME, palette, value, scale, time, rgb);
            break;
        }
    }
    
    if (RESCALE_COLOR_IN_CENTER)
    {
        value *= color_scale[i][j];
        for (k=0; k<3; k++) rgb[k] = 1.0 + color_scale[i][j]*(rgb[k]-1.0);
//         for (k=0; k<3; k++) rgb[k] = 0.5 + color_scale[i][j]*(rgb[k]-0.5);
    }
    
    return(value);
}


void draw_wave_profile_horizontal(double *values, int size, int average, int pnb, int fade, double fade_value)
/* draw a horizontal profile of the wave, if option DRAW_WAVE_PROFILE is active */
{
    int i, k, ij[2];
    double deltav, a, b, y;
    static int imin, imax, jmin, jmax, jmid, jval, d, first = 1;
    static double vmin[2], vmax[2], deltaj, ymin, average_vals[2*NX];
    
    if (first)
    {
        imin = 100;
        imax = NX - 250;
        if (PROFILE_AT_BOTTOM)
        {
            jmin = 50;
            jmax = 150;
        }
        else
        {
            jmin = NY - 150;
            jmax = NY - 50;
        }
        jmid = (jmin + jmax)/2;
        jmid -= (jmid%size);
        xy_to_ij(0.0, WAVE_PROFILE_Y, ij);
        jval = ij[1];
        jval -= (jval%size);
        d = (jmax - jmin)/10 + 1;
        deltaj = (double)(jmax - jmin - 2*d); 
        ymin = (double)(jmin + d);
        for (i=0; i<2*NX; i++) average_vals[i] = 0.0;
        for (k=0; k<2; k++)
        {
            vmin[k] = 1.0e10;
            vmax[k] = -vmin[k];
            if (average) vmin[k] = 0.0;
        }
        first = 0;
    }
    
    if ((average)&&(!fade))
    {
        for (i=imin; i<imax; i+=size)
            average_vals[pnb*NX+i] = 0.99*average_vals[pnb*NX+i] + 0.01*values[i*NY+jval]*values[i*NY+jval];
        for (i=imin; i<imax; i+=size)
            if (average_vals[pnb*NX+i] > vmax[pnb]) vmax[pnb] = average_vals[pnb*NX+i];
        
    }
    else if (!fade) 
    {
        for (i=imin; i<imax; i+=size)
        {
//         if (values[i*NY+jmin] > vmax) vmax = values[i*NY+jmin];
//         if (values[i*NY+jmin] < vmin) vmin = values[i*NY+jmin];        
            if (values[i*NY+jval] > vmax[pnb]) vmax[pnb] = values[i*NY+jval];
            if (values[i*NY+jval] < vmin[pnb]) vmin[pnb] = values[i*NY+jval];        
        }
    }
//     if (vmax <= vmin) vmax = vmin + 0.01;
//     if ((vmin[pnb] < 0.0)&&(vmax[pnb] > 0.0))
//     {
//         if (vmax[pnb] > -vmin[pnb]) vmin[pnb] = -vmax[pnb];
//         else if (vmax[pnb] < -vmin[pnb]) vmax[pnb] = -vmin[pnb];
//     }
//     printf("vmin = %.3lg, vmax = %.3lg\n", vmin, vmax);
    deltav = vmax[pnb]-vmin[pnb];
    if (deltav <= 0.0) deltav = 0.01;
    a = deltaj/deltav;
    b = ymin - a*vmin[pnb];
    
    erase_area_ij(imin, jmin, imax, jmax);
    
    if (fade) glColor3f(fade_value, fade_value, fade_value);
    else glColor3f(1.0, 1.0, 1.0);
    glBegin(GL_LINE_STRIP);
    for (i=imin; i<imax; i+=size)
    {
//         y = a*values[i*NY+jmin] + b;
        if (average) y = a*average_vals[pnb*NX+i] + b;
        else y = a*values[i*NY+jval] + b;
        if (y > ymin) glVertex2d((double)i, y);
    }
    glEnd();
    
    glBegin(GL_LINE_LOOP);
    glVertex2i(imin, jmin);
    glVertex2i(imax, jmin);
    glVertex2i(imax, jmax);
    glVertex2i(imin, jmax);
    glEnd();
    
    /* slightly decrease extremal values in case signal gets weaker */
    if (!fade)
    {
        vmax[pnb] *= 0.99;
        vmin[pnb] *= 0.99;
//         vmax[pnb] *= 0.95;
//         vmin[pnb] *= 0.95;
    }
}


void draw_wave_profile_vertical(double *values, int size, int average, int pnb, int fade, double fade_value)
/* draw a vertical profile of the wave, if option DRAW_WAVE_PROFILE is active */
{
    int j, k, ij[2];
    double deltav, a, b, x;
    static int imin, imax, ival, jmin, jmax, imid, d, first = 1;
    static double vmin[2], vmax[2], deltai, xmin, average_vals[2*NY];
    
    if (first)
    {        
        jmin = 50;
        jmax = NY - 50;
        imin = NX - 150;
        imax = NX - 50;
        imid = (imin + imax)/2;
        imid -= (imid%size);
        xy_to_ij(WAVE_PROFILE_X, 0.0, ij);
        ival = ij[0];
        ival -= (ival%size);
        d = (imax - imin)/10 + 1;
        deltai = (double)(imax - imin - 2*d); 
        xmin = (double)(imin + d);
        for (j=0; j<2*NY; j++) average_vals[j] = 0.0;
        for (k=0; k<2; k++)
        {
            vmin[k] = 1.0e10;
            vmax[k] = -vmin[k];
            if (average) vmin[k] = 0.0;
        }
        first = 0;
    }
    
    if (DRAW_WAVE_SOURCE == 2) 
    {
        xy_to_ij(focus_x, 0.0, ij);
        if (ij[0] < 0) ij[0] = 0;
        if (ij[0] > NX-1) ij[0] = NX-1;
        ival = ij[0];
        ival -= (ival%size);
    }
    
    if ((average)&&(!fade))
    {
        for (j=jmin; j<jmax; j+=size)
            average_vals[pnb*NY+j] = 0.99*average_vals[pnb*NY+j] + 0.01*values[ival*NY+j]*values[ival*NY+j];
        for (j=jmin; j<jmax; j+=size)
            if (average_vals[pnb*NY+j] > vmax[pnb]) vmax[pnb] = average_vals[pnb*NY+j];
        
    }
    
        
    else if (!fade)
    {
        for (j=jmin; j<jmax; j+=size)
        {
//             if (values[imin*NY+j] > vmax) vmax = values[imin*NY+j];
//             if (values[imin*NY+j] < vmin) vmin = values[imin*NY+j];        
            if (values[ival*NY+j] > vmax[pnb]) vmax[pnb] = values[ival*NY+j];
            if (values[ival*NY+j] < vmin[pnb]) vmin[pnb] = values[ival*NY+j];        
        }
        /* symmetrize output */
//         if ((vmin[pnb] < 0.0)&&(vmax[pnb] > 0.0))
//         {
//             if (vmax[pnb] > -vmin[pnb]) vmin[pnb] = -vmax[pnb];
//             else if (vmax[pnb] < -vmin[pnb]) vmax[pnb] = -vmin[pnb];
//         }

//         if (vmax > -vmin) vmin = -vmax;
//         else if (vmax < -vmin) vmax = -vmin;
    }
    
//     printf("vmin = %.3lg, vmax = %.3lg\n", vmin[pnb], vmax[pnb]);

//     if ((vmin[pnb] < 0.0)&&(vmax[pnb] > 0.0))
//     {
//         if (vmax[pnb] > -vmin[pnb]) vmin[pnb] = -vmax[pnb];
//         else if (vmax[pnb] < -vmin[pnb]) vmax[pnb] = -vmin[pnb];
//     }
    
    deltav = vmax[pnb]-vmin[pnb];
    if (deltav <= 0.0) deltav = 0.01;
    a = deltai/deltav;
    b = xmin - a*vmin[pnb];

    erase_area_ij(imin, jmin, imax, jmax);
    
    if (fade) glColor3f(fade_value, fade_value, fade_value);
    else glColor3f(1.0, 1.0, 1.0);
    glBegin(GL_LINE_STRIP);
    for (j=jmin; j<jmax; j+=size)
    {
        if (average) x = a*average_vals[pnb*NY+j] + b;
        else x = a*values[ival*NY+j] + b;
//         else x = a*values[imin*NY+j] + b;
        if (x >= xmin) glVertex2d(x, (double)j);
    }
    glEnd();
    
    glBegin(GL_LINE_LOOP);
    glVertex2i(imin, jmin);
    glVertex2i(imax, jmin);
    glVertex2i(imax, jmax);
    glVertex2i(imin, jmax);
    glEnd();
    
    /* slightly decrease extramal values in case signal gets weaker */
    if (!fade)
    {
        vmax[pnb] *= 0.99;
        vmin[pnb] *= 0.99;
    }
}


void draw_wave_profile(double *values, int size, int p_number, int fade, double fade_value)
/* draw a profile of the wave, if option DRAW_WAVE_PROFILE is active */
{
    if (VERTICAL_WAVE_PROFILE) 
        draw_wave_profile_vertical(values, size, AVERAGE_WAVE_PROFILE, p_number, fade, fade_value);
    if (HORIZONTAL_WAVE_PROFILE) 
        draw_wave_profile_horizontal(values, size, AVERAGE_WAVE_PROFILE, p_number, fade, fade_value);
}

void draw_exit_timeseries(double xleft, double yleft, double xtest, double ytest, double *values, int size, int top, int symmetrize, int fade, double fade_value)
/* draw a profile of the wave, if option DRAW_WAVE_TIMESERIES is active */
{
    int t, t1, s, padding = 50, nvals = TIMESERIES_NVALUES;
    double a, b, deltav, x, y, wmin, wmax;
    static int ij[2], ij_test[2], imin, imax, jmintop, jmaxtop, jminbot, jmaxbot, jmin, jmax, itest, jtest, jtesttop, jtestbot, counter[2], first  = 1, window = 5;
    static double timeseries[TIMESERIES_NVALUES][2], average[TIMESERIES_NVALUES][2], vmin[2], vmax[2], ymin, deltaj, deltai;
    
    if (first)
    {
        xy_to_ij(xleft, yleft, ij);
        xy_to_ij(xtest, ytest, ij_test);
        jmaxtop = NY - padding;
        jmintop = 2*ij[1] - jmaxtop;
        jtesttop = ij_test[1];
        
//         printf("jmintop = %i, jmaxtop = %i\n", jmintop, jmaxtop);
        
        xy_to_ij(xleft, -yleft, ij);
        xy_to_ij(xtest, -ytest, ij_test);
        jminbot = padding;
        jmaxbot = 2*ij[1] - jminbot;            
        jtestbot = ij_test[1];
        
//         printf("jminbot = %i, jmaxbot = %i\n", jminbot, jmaxbot);
        
        imin = ij[0];
        imax = NX - padding;
        deltai = (double)(imax - imin)/(double)nvals;
        
        vmin[0] = 1.0e10;
        vmax[0] = -vmin[0];
        vmin[1] = 1.0e10;
        vmax[1] = -vmin[1];
        for (t=0; t<nvals; t++) 
        {
            timeseries[t][0] = 0.0;
            timeseries[t][1] = 0.0;
        }
        
        counter[0] = 0;
        counter[1] = 0;
        first = 0;
    }
    
    if (top)
    {
        jmin = jmintop;
        jmax = jmaxtop;
        jtest = jtesttop;
    }
    else
    {
        jmin = jminbot;
        jmax = jmaxbot;
        jtest = jtestbot;
    }
    ymin = (double)(jmin);
    itest = ij_test[0];
    deltaj = (double)(jmax - jmin); 
    
    timeseries[counter[top]][top] = values[itest*NY + jtest] + values[(itest-1)*NY + jtest] + values[(itest)*NY + jtest + 1] + values[(itest-1)*NY + jtest + 1];
    counter[top]++;
    if (counter[top] == nvals) counter[top] = 0;
    
//     for (t=0; t<nvals; t++) printf("val[%i] = %.3lg\n", t, timeseries[t]);
        
    for (t=0; t<nvals; t++)
    {
        if (timeseries[t][top] > vmax[top]) vmax[top] = timeseries[t][top];
        if (timeseries[t][top] < vmin[top]) vmin[top] = timeseries[t][top];        
    }
    if (symmetrize)
    {
        if (vmax[top] > -vmin[top]) vmin[top] = -vmax[top];
        else if (vmax[top] < -vmin[top]) vmax[top] = -vmin[top];
    }
    
//     printf("vmin = %.3lg, vmax = %.3lg\n", vmin[top], vmax[top]);
    deltav = vmax[top]-vmin[top];
    if (deltav <= 0.0) deltav = 0.01;
    a = deltaj/deltav;
    b = ymin - a*vmin[top];
    
    /* compute average */
    for (t=window; t<nvals-window; t++)
    {
        average[t][top] = 0.0;
        for (s=-window; s<window; s++) 
        {
            t1 = counter[top] - t + s;
            if (t1 < 0) t1 += nvals;
            if (t1 > nvals-1) t1 -= nvals;
            average[t][top] += timeseries[t1][top]*timeseries[t1][top];
        }
        average[t][top] = sqrt(average[t][top]*0.5/(double)window);
    }
        
    erase_area_ij(imin, jmin, imax, jmax);
    
    if (fade) glColor3f(fade_value, 0.0, 0.0);
    else glColor3f(1.0, 0.0, 0.0);
    
    
    glBegin(GL_LINE_STRIP);
    for (t=0; t<nvals; t++)
    {
        t1 = counter[top] - t;
        if (t1 < 0) t1 += nvals;
        x = (double)imin + deltai*(double)t;
        y = a*timeseries[t1][top] + b;
        glVertex2d(x, y);
    }
    glEnd();
//     printf("(x,y) = (%.2f, %.2f)\n", x, y);
        
    
    /* draw average */
    if (fade) glColor3f(fade_value, fade_value, fade_value);
    else glColor3f(1.0, 1.0, 1.0);
    glBegin(GL_LINE_STRIP);
    for (t1=window; t1<nvals-window; t1++)
    {
        x = (double)imin + deltai*(double)t1;
        y = -a*average[t1][top] + b;
        if (x >= imin) glVertex2d(x, y);
    }
    glEnd();
    if (symmetrize)
    {
        glBegin(GL_LINE_STRIP);
        for (t1=window; t1<nvals-window; t1++)
        {
            x = (double)imin + deltai*(double)t1;
            y = a*average[t1][top] + b;
            glVertex2d(x, y);
        }
        glEnd();
    }
    
    if (fade) glColor3f(fade_value, fade_value, fade_value);
    else glColor3f(1.0, 1.0, 1.0);
    
    glBegin(GL_LINE_LOOP);
    glVertex2i(imin, jmin);
    glVertex2i(imax, jmin);
    glVertex2i(imax, jmax);
    glVertex2i(imin, jmax);
    glEnd();
}


void draw_exit_timeseries_old(double *values, int size, int fade, double fade_value)
/* draw a profile of the wave, if option DRAW_WAVE_TIMESERIES is active */
{
    int t, t1, s, padding = 50, nvals = TIMESERIES_NVALUES;
    double a, b, deltav, x, y, wmin, wmax;
    static int ij[2], ij_test[2], imin, imax, jmin, jmax, jmid, itest, jtest, counter = 0, first  = 1, window = 5;
    static double timeseries[TIMESERIES_NVALUES], average[TIMESERIES_NVALUES], vmin, vmax, ymin, deltaj, deltai;
    
    if (first)
    {
        xy_to_ij(LAMBDA, -1.0 + MU, ij);
//         xy_to_ij(LAMBDA - 0.1*MU, -1.0 + MU, ij_test);
        xy_to_ij(LAMBDA, -1.0 + 2.0*MU, ij_test);
        imin = ij[0];
        imax = NX - padding;
        jmin = padding;
        jmax = 2*ij[1] - jmin;
        jmid = (jmin + jmax)/2;
        ymin = (double)(jmin);
        itest = ij_test[0];
        jtest = ij_test[1];
        deltai = (double)(imax - imin)/(double)nvals;
        deltaj = (double)(jmax - jmin); 
        vmin = 1.0e10;
        vmax = -vmin;
        for (t=0; t<nvals; t++) timeseries[t] = 0.0;
        first = 0;
    }
    
    timeseries[counter] = values[itest*NY + jtest] + values[(itest-1)*NY + jtest] + values[(itest)*NY + jtest + 1] + values[(itest-1)*NY + jtest + 1];
    counter++;
    if (counter == nvals) counter = 0;
    
//     for (t=0; t<nvals; t++) printf("val[%i] = %.3lg\n", t, timeseries[t]);
        
    for (t=0; t<nvals; t++)
    {
        if (timeseries[t] > vmax) vmax = timeseries[t];
        if (timeseries[t] < vmin) vmin = timeseries[t];        
    }
    if (vmax > -vmin) vmin = -vmax;
    else if (vmax < -vmin) vmax = -vmin;
    
    printf("vmin = %.3lg, vmax = %.3lg\n", vmin, vmax);
    deltav = vmax-vmin;
    if (deltav <= 0.0) deltav = 0.01;
    a = deltaj/deltav;
    b = ymin - a*vmin;
    
    /* compute average */
    for (t=window; t<nvals-window; t++)
    {
        average[t] = 0.0;
        for (s=-window; s<window; s++) 
        {
            t1 = counter - t + s;
            if (t1 < 0) t1 += nvals;
            if (t1 > nvals-1) t1 -= nvals;
            average[t] += timeseries[t1]*timeseries[t1];
        }
        average[t] = sqrt(average[t]*0.5/(double)window);
    }
        
    erase_area_ij(imin, jmin, imax, jmax);
    
    if (fade) glColor3f(fade_value, fade_value, fade_value);
    else glColor3f(1.0, 1.0, 1.0);
    
    glBegin(GL_LINE_STRIP);
    for (t=0; t<nvals; t++)
    {
        t1 = counter - t;
        if (t1 < 0) t1 += nvals;
        x = (double)imin + deltai*(double)t;
        y = a*timeseries[t1] + b;
        glVertex2d(x, y);
    }
    glEnd();
    printf("(x,y) = (%.2f, %.2f)\n", x, y);
    
    /* draw average */
    if (fade) glColor3f(fade_value, 0.0, 0.0);
    else glColor3f(1.0, 0.0, 0.0);
    glBegin(GL_LINE_STRIP);
    for (t1=window; t1<nvals-window; t1++)
    {
        x = (double)imin + deltai*(double)t1;
        y = a*average[t1] + b;
        if (x >= imin) glVertex2d(x, y);
    }
    glEnd();
    glBegin(GL_LINE_STRIP);
    for (t1=window; t1<nvals-window; t1++)
    {
        x = (double)imin + deltai*(double)t1;
        y = -a*average[t1] + b;
        glVertex2d(x, y);
    }
    glEnd();
    
    if (fade) glColor3f(fade_value, fade_value, fade_value);
    else glColor3f(1.0, 1.0, 1.0);
    
    glBegin(GL_LINE_LOOP);
    glVertex2i(imin, jmin);
    glVertex2i(imax, jmin);
    glVertex2i(imax, jmax);
    glVertex2i(imin, jmax);
    glEnd();
}


void draw_entrance_timeseries(double xleft, double yleft, double *values, int size, int fade, double fade_value)
/* draw a profile of the wave, if option DRAW_WAVE_TIMESERIES is active */
{
    int t, t1, padding = 50, nvals = TIMESERIES_NVALUES, ishift = 4;
    double vmin, vmax, deltav, x, y;
    static int ij[2], imin, imax, jmin, jmax, jmid, itest, jtest, counter = 0, first  = 1, 
    time = 0;
    static double ymin, deltaj, deltai, a, b;
    
    if (first)
    {
        xy_to_ij(-xleft, yleft, ij);
        imin = padding;
        imax = ij[0];
        jmax = NY - padding;
        jmin = 2*ij[1] - jmax;
        jmid = (jmin + jmax)/2;
        ymin = (double)(jmin);
        deltai = (double)(imax - imin)/(double)nvals;
        deltaj = (double)(jmax - jmin); 
        b = (double)jmin + 0.25*deltaj;
        a = 0.5*deltaj;
        
        first = 0;
    }
    
    counter++;
        
    erase_area_ij(imin, jmin, imax, jmax);
    
    if (fade) glColor3f(fade_value, fade_value, fade_value);
    else glColor3f(1.0, 1.0, 1.0);
    
//     glEnable(GL_SCISSOR_TEST);
//     glScissor(imin, jmin, imax-imin, jmax-jmin);
    glBegin(GL_LINE_STRIP);
    for (t=0; t<nvals + ishift; t++)
    {
        t1 = counter - t + nvals;
        if (t1 >= NSTEPS) t1 = NSTEPS;
//         x = (double)imin + deltai*(double)t;
        x = (double)imin + deltai*(double)(t-ishift);
        y = a*(double)input_signal[t1] + b;
        if (t > ishift) glVertex2d(x, y);
    }
    glEnd();
//     glDisable(GL_SCISSOR_TEST);
    
    glBegin(GL_LINE_LOOP);
    glVertex2i(imin, jmin);
    glVertex2i(imax, jmin);
    glVertex2i(imax, jmax);
    glVertex2i(imin, jmax);
    glEnd();
}

void draw_wave_timeseries(double xleft, double yleft, double xtest, double ytest, double *values, int size, int fade, double fade_value)
/* draw a profile of the wave, if option DRAW_WAVE_TIMESERIES is active */
{
    int symmetrize = 1;
    
    draw_exit_timeseries(xleft, yleft, xtest, ytest, values, size, 0, symmetrize, fade, fade_value);
    if (DRAW_WAVE_TIMESERIES == 2) 
        draw_exit_timeseries(xleft, yleft, xtest, ytest, values, size, 1, symmetrize, fade, fade_value);
    draw_entrance_timeseries(xleft, yleft, values, size, fade, fade_value);
}

void draw_wave_timeseries_old(double *values, int size, int fade, double fade_value)
/* draw a profile of the wave, if option DRAW_WAVE_TIMESERIES is active */
{
    draw_exit_timeseries_old(values, size, fade, fade_value);
    draw_entrance_timeseries(1.0, 1.0, values, size, fade, fade_value);
}


void init_fade_table(double *tcc_table[NX], double *fade_table)
{
    int i, j;
    
    for (i=0; i<NX; i++)
        for (j=0; j<NY; j++)
            fade_table[i*NY + j] = pow(tcc_table[i][j]/courant2, 0.25);
//             fade_table[i*NY + j] = tcc_table[i][j]/courant2;
}

void draw_wave_highres_palette(int size, double *phi[NX], double *psi[NX], double *total_energy[NX], double *average_energy[NX], double *fade_table, double *total_flux, double *color_scale[NX], short int *xy_in[NX], double scale, int time, int plot, int palette, int pnumber, int fade, double fade_value)
/* same as draw_wave_highres, but with color scheme option */
{
    int i, j, k, iplus, iminus, jplus, jminus;
    double value, rgb[3], xy[2], x1, y1, x2, y2, velocity, energy, gradientx, gradienty, arg, mod, flux_factor, gx, gy, mgx, mgy, factor, norm, pscal, ca, vscale2;
//     double vmin, vmax, deltav;
    static double dtinverse = ((double)NX)/(COURANT*(XMAX-XMIN)), dx = (XMAX-XMIN)/((double)NX), col_ratio;
    static int first = 1;
    double *values, *shade, *rgbvals;
        
    if (first)
    {
        col_ratio = 1.0/1.01;
        first = 0;
    }
    
    values = (double *)malloc(NX*NY*sizeof(double));
    rgbvals = (double *)malloc(3*NX*NY*sizeof(double));
    #pragma omp parallel for private(i,j,rgb)
    for (i=0; i<NX; i+=size)
        for (j=0; j<NY; j+=size)
        {
            values[i*NY+j] = wave_value(i, j, phi, psi, total_energy, average_energy, total_flux, color_scale, xy_in, scale, time, plot, palette, rgb);
            rgbvals[i*NY+j] = rgb[0];
            rgbvals[NX*NY+i*NY+j] = rgb[1];
            rgbvals[2*NX*NY+i*NY+j] = rgb[2];
        }
        
    if (SHADE_2D)
    {
        vscale2 = SHADE_SCALE_2D*SHADE_SCALE_2D;
        shade = (double *)malloc(NX*NY*sizeof(double));
        #pragma omp parallel for private(i,j,gradientx,gradienty,norm,pscal,ca)
        for (i=0; i<NX-size; i+=size)
        {
            for (j=0; j<NY-size; j+=size)
            {
                gradientx = values[(i+1)*NY+j] - values[i*NY+j];
                gradienty = values[i*NY+j+1] - values[i*NY+j];
                norm = sqrt(vscale2 + gradientx*gradientx + gradienty*gradienty);
                pscal = -gradientx*light[0] - gradienty*light[1] + SHADE_SCALE_2D;
                ca = pscal/norm;
                ca = (ca + 1.0)*0.4 + 0.2;
                if (RESCALE_COLOR_IN_CENTER) for (k=0; k<3; k++) 
                    rgbvals[k*NX*NY+i*NY+j] = 1.0 + (ca*rgbvals[k*NX*NY+i*NY+j] - 1.0)*(color_scale[i][j]+0.01)*col_ratio;
//                     rgbvals[k*NX*NY+i*NY+j] = 0.5 + (ca*rgbvals[k*NX*NY+i*NY+j] - 0.5)*(color_scale[i][j]+0.01)*col_ratio;
                else for (k=0; k<3; k++) rgbvals[k*NX*NY+i*NY+j] *= ca;
            }
        }
    }
    
    glBegin(GL_QUADS);
    
//     printf("dtinverse = %.5lg\n", dtinverse);

    for (i=0; i<NX; i+=size)
        for (j=0; j<NY; j+=size)
        {
            if ((TWOSPEEDS)||(xy_in[i][j]))
            {
//                 value = wave_value(i, j, phi, psi, total_energy, average_energy, total_flux, xy_in, scale, time, plot, palette, rgb);
                
//                 if ((DRAW_WAVE_PROFILE)||(DRAW_WAVE_TIMESERIES)) values[i*NY+j] = value;
//                 rgb[0] = r[i*NY+j];
//                 rgb[1] = g[i*NY+j];
//                 rgb[2] = b[i*NY+j];
                
                if (FADE_IN_OBSTACLE)
                    for (k=0; k<3; k++) rgbvals[k*NX*NY+i*NY+j] *= fade_table[i*NY+j];
                        
                if (fade) for (k=0; k<3; k++) rgbvals[k*NX*NY+i*NY+j] *= fade_value;
                glColor3f(rgbvals[i*NY+j], rgbvals[NX*NY+i*NY+j], rgbvals[2*NX*NY+i*NY+j]);
                
                glVertex2i(i, j);
                glVertex2i(i+size, j);
                glVertex2i(i+size, j+size);
                glVertex2i(i, j+size);
            }
        }

    glEnd ();
    
    draw_billiard(fade, fade_value);
    
    if (DRAW_WAVE_TIMESERIES) draw_wave_timeseries(1.0, 1.0, 1.0, 0.5, values, size, fade, fade_value);
    if (DRAW_WAVE_PROFILE) draw_wave_profile(values, size, pnumber, fade, fade_value);
    if (SHADE_2D) free(shade);
    if (DRAW_WAVE_SOURCE >= 1)
    {
        glColor3f(0.0, 0.0, 0.0);
        for (k=0; k<N_SOURCES; k++) draw_circle(wave_source_x[k], wave_source_y[k], 0.005, 100);
        if (DRAW_WAVE_SOURCE >= 2)
            draw_circle(focus_x, 0.0, 0.005, 100);
    }
    
    free(values);
    free(rgbvals);
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

