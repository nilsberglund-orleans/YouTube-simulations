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
            
	    if ((xy_in[i][j])||(TWOSPEEDS)) phi[i][j] = 0.2*exp(-dist2/0.001)*cos(-sqrt(dist2)/0.01);
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
        for (j=0; j<NY; j++)
        {
            ij_to_xy(i, j, xy);
            dist2 = (xy[0]-x)*(xy[0]-x);
	    xy_in[i][j] = xy_in_billiard(xy[0],xy[1]);
            
// 	    if ((xy_in[i][j])||(TWOSPEEDS)) phi[i][j] = 0.01*exp(-dist2/0.0005)*cos(-sqrt(dist2)/0.01);
// 	    if ((xy_in[i][j])||(TWOSPEEDS)) phi[i][j] = 0.01*exp(-dist2/0.0005)*cos(-sqrt(dist2)/0.02);
// 	    if ((xy_in[i][j])||(TWOSPEEDS)) phi[i][j] = 0.01*exp(-dist2/0.005)*cos(-sqrt(dist2)/0.02);
// 	    if ((xy_in[i][j])||(TWOSPEEDS)) phi[i][j] = 0.01*exp(-dist2/0.01)*cos(-sqrt(dist2)/0.02);
	    if ((xy_in[i][j])||(TWOSPEEDS)) phi[i][j] = 0.01*exp(-dist2/0.0005)*cos(-sqrt(dist2)/0.01);
//             if ((xy_in[i][j])||(TWOSPEEDS)) phi[i][j] = 0.01*exp(-dist2/0.05)*cos(-sqrt(dist2)/0.025);
            else phi[i][j] = 0.0;
            psi[i][j] = 0.0;
        }
}


void init_wave_flat( double *phi[NX], double *psi[NX], short int * xy_in[NX])
/* initialise flat field - phi is wave height, psi is phi at time t-1 */
{
    int i, j;
    double xy[2], dist2;

    for (i=0; i<NX; i++)
        for (j=0; j<NY; j++)
        {
            ij_to_xy(i, j, xy);
	    xy_in[i][j] = xy_in_billiard(xy[0],xy[1]);
	    phi[i][j] = 0.0;
            psi[i][j] = 0.0;
        }
}


void add_drop_to_wave(double factor, double x, double y, double *phi[NX], double *psi[NX])
/* add drop at (x,y) to the field with given prefactor */
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
    else return(0);
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
                    color_scheme(COLOR_SCHEME, phi[i][j], scale, time, rgb);
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





