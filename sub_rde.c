double courant2;  /* Courant parameter squared */
double dx2;       /* spatial step size squared */
double intstep;   /* integration step */
double intstep1;  /* integration step used in absorbing boundary conditions */
double viscosity; /* viscosity (parameter in front of Laplacian) */
double rpslzb;    /* second parameter in Rock-Paper-Scissors-Lizard-Spock equation */

double gaussian()
/* returns standard normal random variable, using Box-Mueller algorithm */
{
    static double V1, V2, S;
    static int phase = 0;
    double X;

    if (phase == 0) 
    {
        do 
        {
        double U1 = (double)rand() / RAND_MAX;
        double U2 = (double)rand() / RAND_MAX;
        V1 = 2 * U1 - 1;
        V2 = 2 * U2 - 1;
        S = V1 * V1 + V2 * V2;
        } 
        while(S >= 1 || S == 0);
        X = V1 * sqrt(-2 * log(S) / S);
    } 
    else X = V2 * sqrt(-2 * log(S) / S);

    phase = 1 - phase;

    return X;
}

void init_random(double mean, double amplitude, double *phi[NFIELDS], short int xy_in[NX*NY])
/* initialise field with gaussian at position (x,y) */
{
    int i, j, k, in;
    double xy[2], dist2, module, phase, scale2;    

    printf("Initialising xy_in\n");
    #pragma omp parallel for private(i,j)
    for (i=0; i<NX; i++)
    {
        if (i%100 == 0) printf("Column %i of %i\n", i, NX);
        for (j=0; j<NY; j++)
        {
            ij_to_xy(i, j, xy);
            xy_in[i*NY+j] = xy_in_billiard(xy[0],xy[1]);
        }
    }
    
    printf("Initialising fields\n");
    for (k=0; k<NFIELDS; k++)
        for (i=0; i<NX; i++)
        {
            if (i%100 == 0) printf("Field %i of %i, column %i of %i\n", k, NFIELDS, i, NX);
            for (j=0; j<NY; j++)
            {
                if (xy_in[i*NY+j]) phi[k][i*NY+j] = mean + amplitude*(2.0*(double)rand()/RAND_MAX - 1.0);
                else phi[k][i*NY+j] = 0.0;
            }
        }
    printf("Fields initialised\n");
}



void init_gaussian(double x, double y, double mean, double amplitude, double scalex, 
                   double *phi[NX], short int xy_in[NX*NY])
/* initialise field with gaussian at position (x,y) */
{
    int i, j, in;
    double xy[2], dist2, module, phase, scale2;    

    scale2 = scalex*scalex;
    printf("Initialising field\n");
    for (i=0; i<NX; i++)
        for (j=0; j<NY; j++)
        {
            ij_to_xy(i, j, xy);
	    xy_in[i*NY+j] = xy_in_billiard(xy[0],xy[1]);

            in = xy_in[i*NY+j];
            if (in == 1)
            {
                dist2 = (xy[0]-x)*(xy[0]-x) + (xy[1]-y)*(xy[1]-y);
                module = amplitude*exp(-dist2/scale2);
                if (module < 1.0e-15) module = 1.0e-15;

                phi[i][j] = mean + module/scalex;
            }   /* boundary temperatures */
            else if (in >= 2) phi[i][j] = T_IN*pow(0.75, (double)(in-2));
//             else if (in >= 2) phi[i][j] = T_IN*pow(1.0 - 0.5*(double)(in-2), (double)(in-2));
//             else if (in >= 2) phi[i][j] = T_IN*(1.0 - (double)(in-2)/((double)MDEPTH))*(1.0 - (double)(in-2)/((double)MDEPTH));
            else phi[i][j] = T_OUT;
        }
}

void init_coherent_state(double x, double y, double px, double py, double scalex, double *phi[NFIELDS], short int xy_in[NX*NY])
/* initialise field with coherent state of position (x,y) and momentum (px, py) */
/* phi[0] is real part, phi[1] is imaginary part */
{
    int i, j;
    double xy[2], dist2, module, phase, scale2;    

    scale2 = scalex*scalex;
    for (i=0; i<NX; i++)
        for (j=0; j<NY; j++)
        {
            ij_to_xy(i, j, xy);
	    xy_in[i*NY+j] = xy_in_billiard(xy[0],xy[1]);

            if (xy_in[i*NY+j])
            {
                dist2 = (xy[0]-x)*(xy[0]-x) + (xy[1]-y)*(xy[1]-y);
                module = exp(-dist2/scale2);
                if (module < 1.0e-15) module = 1.0e-15;
                phase = (px*(xy[0]-x) + py*(xy[1]-y))/scalex;

                phi[0][i*NY+j] = module*cos(phase);
                phi[1][i*NY+j] = module*sin(phase);
            }
            else
            {
                phi[0][i*NY+j] = 0.0;
                phi[1][i*NY+j] = 0.0;
            }
        }
}

void init_fermion_state(double x, double y, double px, double py, double scalex, double *phi[NFIELDS], short int xy_in[NX*NY])
/* initialise field with antisymmetric coherent state of position (x,y) and momentum (px, py) */
/* phi[0] is real part, phi[1] is imaginary part */
{
    int i, j;
    double xy[2], dist2, module, phase, scale2;    

    scale2 = scalex*scalex;
    for (i=0; i<NX; i++)
        for (j=0; j<NY; j++)
        {
            ij_to_xy(i, j, xy);
	    xy_in[i*NY+j] = xy_in_billiard(xy[0],xy[1]);

            if (xy_in[i*NY+j])
            {
                dist2 = (xy[0]-x)*(xy[0]-x) + (xy[1]-y)*(xy[1]-y);
                module = exp(-dist2/scale2);
                if (module < 1.0e-15) module = 1.0e-15;
                phase = (px*(xy[0]-x) + py*(xy[1]-y))/scalex;

                phi[0][i*NY+j] = module*cos(phase);
                phi[1][i*NY+j] = module*sin(phase);
                
                /* antisymmetrize wave function */
                dist2 = (xy[1]-x)*(xy[1]-x) + (xy[0]-y)*(xy[0]-y);
                module = exp(-dist2/scale2);
                if (module < 1.0e-15) module = 1.0e-15;
                phase = (px*(xy[1]-x) + py*(xy[0]-y))/scalex;

                phi[0][i*NY+j] -= module*cos(phase);
                phi[1][i*NY+j] -= module*sin(phase);
            }
            else
            {
                phi[0][i*NY+j] = 0.0;
                phi[1][i*NY+j] = 0.0;
            }
        }
}

void init_boson_state(double x, double y, double px, double py, double scalex, double *phi[NFIELDS], short int xy_in[NX*NY])
/* initialise field with symmetric coherent state of position (x,y) and momentum (px, py) */
/* phi[0] is real part, phi[1] is imaginary part */
{
    int i, j;
    double xy[2], dist2, module, phase, scale2;    

    scale2 = scalex*scalex;
    for (i=0; i<NX; i++)
        for (j=0; j<NY; j++)
        {
            ij_to_xy(i, j, xy);
	    xy_in[i*NY+j] = xy_in_billiard(xy[0],xy[1]);

            if (xy_in[i*NY+j])
            {
                dist2 = (xy[0]-x)*(xy[0]-x) + (xy[1]-y)*(xy[1]-y);
                module = exp(-dist2/scale2);
                if (module < 1.0e-15) module = 1.0e-15;
                phase = (px*(xy[0]-x) + py*(xy[1]-y))/scalex;

                phi[0][i*NY+j] = module*cos(phase);
                phi[1][i*NY+j] = module*sin(phase);
                
                /* symmetrize wave function */
                dist2 = (xy[1]-x)*(xy[1]-x) + (xy[0]-y)*(xy[0]-y);
                module = exp(-dist2/scale2);
                if (module < 1.0e-15) module = 1.0e-15;
                phase = (px*(xy[1]-x) + py*(xy[0]-y))/scalex;

                phi[0][i*NY+j] += module*cos(phase);
                phi[1][i*NY+j] += module*sin(phase);
            }
            else
            {
                phi[0][i*NY+j] = 0.0;
                phi[1][i*NY+j] = 0.0;
            }
        }
}

void init_fermion_state2(double x, double y, double px, double py, double scalex, double *phi[NFIELDS], short int xy_in[NX*NY])
/* initialise field with antisymmetric coherent state of position (x,y) and momentum (px, py) */
/* phi[0] is real part, phi[1] is imaginary part */
{
    int i, j;
    double xy[2], dist2, module, phase, scale2, fx[2], fy[2], gx[2], gy[2];    

    scale2 = scalex*scalex;
    for (i=0; i<NX; i++)
        for (j=0; j<NY; j++)
        {
            ij_to_xy(i, j, xy);
	    xy_in[i*NY+j] = xy_in_billiard(xy[0],xy[1]);

            if (xy_in[i*NY+j])
            {
                dist2 = (xy[0]-x)*(xy[0]-x);
                module = exp(-dist2/scale2);
                if (module < 1.0e-15) module = 1.0e-15;
                phase = px*(xy[0]-x)/scalex;
                fx[0] = module*cos(phase);
                fx[1] = module*sin(phase);

                dist2 = (xy[0]-y)*(xy[0]-y);
                module = exp(-dist2/scale2);
                if (module < 1.0e-15) module = 1.0e-15;
                phase = px*(xy[0]-y)/scalex;
                fy[0] = module*cos(phase);
                fy[1] = module*sin(phase);

                dist2 = (xy[1]-x)*(xy[1]-x);
                module = exp(-dist2/scale2);
                if (module < 1.0e-15) module = 1.0e-15;
                phase = py*(xy[1]-x)/scalex;
                gx[0] = module*cos(phase);
                gx[1] = module*sin(phase);

                dist2 = (xy[1]-y)*(xy[1]-y);
                module = exp(-dist2/scale2);
                if (module < 1.0e-15) module = 1.0e-15;
                phase = py*(xy[1]-y)/scalex;
                gy[0] = module*cos(phase);
                gy[1] = module*sin(phase);

                phi[0][i*NY+j] = 0.5*(fx[0]*gy[0] - gx[0]*fy[0]);
                phi[1][i*NY+j] = 0.5*(fx[1]*gy[1] - gx[1]*fy[1]);
            }
            else
            {
                phi[0][i*NY+j] = 0.0;
                phi[1][i*NY+j] = 0.0;
            }
        }
}

void antisymmetrize_wave_function(double *phi[NFIELDS], short int xy_in[NX*NY])
/* force the wave function to be antisymmetric, for simulations of fermions */
{
    int i, j, k;
    
    #pragma omp parallel for private(i,j,k)
    for (i=0; i<NX; i++)
        for (j=0; j<=i; j++) if ((j<NY)&&(xy_in[i*NY+j]))
            for (k=0; k<2; k++)
            {
                phi[k][i*NY+j] = 0.5*(phi[k][i*NY+j] - phi[k][j*NY+i]);
                phi[k][j*NY+i] = -phi[k][i*NY+j];
            }
}

/*********************/
/* animation part    */
/*********************/

double wrap_angle(double angle)
{
    if (angle < 0.0) return(angle + DPI);
    else if (angle >= DPI) return (angle - DPI);
    else return(angle);
}

double unwrap_angle(double angle1, double angle2)
{
    if (angle2 - angle1 < -PI) return(angle2 + DPI);
    else if (angle2 - angle1 >= PI) return (angle2 - DPI);
    else return(angle2);
}


void compute_theta(double *phi[NFIELDS], short int xy_in[NX*NY], t_rde rde[NX*NY])
/* compute angle for rock-paper-scissors equation */
{
    int i, j;
    double u, v, w, rho, xc, yc, angle;
        
    #pragma omp parallel for private(i,j,u, v, w, rho, xc, yc, angle)
    for (i=0; i<NX; i++) for (j=0; j<NY; j++)
        if (xy_in[i*NY+j])
        {
            u = phi[0][i*NY+j];
            v = phi[1][i*NY+j];
            w = phi[2][i*NY+j];
            rho = u + v + w;
            u = u/rho;
            v = v/rho;
            yc = 0.5*(v - 1.0/3.0);
            xc = u - 1.0/3.0 + yc;
            yc *= sqrt(3.0);
            
            angle = argument(xc, yc);
            if (angle < 0.0) angle += DPI;
            if (angle >= DPI) angle -= DPI;
            rde[i*NY+j].theta = angle;
        }
        else rde[i*NY+j].theta = 0.0;
}

double colors_rpslz(int type)
{
    switch (type) {
        case (0): return(0.0);
        case (1): return(60.0);
        case (2): return(120.0);
        case (3): return(200.0);
        case (4): return(270.0);
    }
}

void compute_theta_rpslz(double *phi[NFIELDS], short int xy_in[NX*NY], t_rde rde[NX*NY], int cplot)
/* compute angle for rock-paper -scissors-lizard-Spock equation */
{
    int i, j, k, kmax;
    double rho, xc, yc, angle, max;
    static double sa[5], ca[5], shift;
    static int first = 1;
    
    if (first)
    {
//         shift = 0.0;
        for (i = 0; i < 5; i++)
        {
            ca[i] = cos(colors_rpslz(i)*DPI/360.0);
            sa[i] = sin(colors_rpslz(i)*DPI/360.0);
//             ca[i] = cos(shift + 0.2*DPI*(double)i);
//             sa[i] = sin(shift + 0.2*DPI*(double)i);
        }
        first = 0;
    }   
       
    switch (cplot) {
        case (Z_MAXTYPE_RPSLZ):
        {
            #pragma omp parallel for private(i,j,max,kmax)
            for (i=0; i<NX; i++) for (j=0; j<NY; j++)
            {
                if (xy_in[i*NY+j])
                {
                    max = 0.0;
                    kmax = 0;
                    for (k=0; k<5; k++)
                    {
                        if (phi[k][i*NY+j] > max)
                        {
                            max = phi[k][i*NY+j];
                            kmax = k;
                        }
                    }
                    angle = colors_rpslz(kmax);
                    rde[i*NY+j].theta = angle;
                }
                else rde[i*NY+j].theta = 0.0;
            }
            break;
        }
//         case (Z_THETA_RPSLZ):
        default: 
        {
            #pragma omp parallel for private(i,j,rho,xc,yc,angle)
            for (i=0; i<NX; i++) for (j=0; j<NY; j++)
            {
                if (xy_in[i*NY+j])
                {
                    xc = 0.0;
                    yc = 0.0;
                    for (k=0; k<5; k++)
                    {
                        xc += phi[k][i*NY+j]*ca[k];
                        yc += phi[k][i*NY+j]*sa[k];
                    }
            
                    angle = argument(xc, yc);
                    if (angle < 0.0) angle += DPI;
                    if (angle >= DPI) angle -= DPI;
                    rde[i*NY+j].theta = angle;
                }
                else rde[i*NY+j].theta = 0.0;
            }
            break;
        }
    }
}

void compute_gradient_xy(double phi[NX*NY], double gradient[2*NX*NY])
/* compute the gradient of the field */
{
    int i, j, iplus, iminus, jplus, jminus, padding = 0; 
    double deltaphi;
    double dx = (XMAX-XMIN)/((double)NX);
    
    dx = (XMAX-XMIN)/((double)NX);
    
    #pragma omp parallel for private(i,j,iplus,iminus,jplus,jminus,deltaphi)
    for (i=1; i<NX-1; i++)
        for (j=1; j<NY-1; j++)
        {
            deltaphi = phi[iplus*NY+j] - phi[iminus*NY+j];
            if (vabs(deltaphi) < 1.0e9) gradient[i*NY+j] = (deltaphi)/dx;
            else gradient[i*NY+j] = 0.0;
            
            deltaphi = phi[i*NY+jplus] - phi[i*NY+jminus];
            if (vabs(deltaphi) < 1.0e9) gradient[NX*NY+i*NY+j] = (deltaphi)/dx;
            else gradient[NX*NY+i*NY+j] = 0.0;
        }
        
    /* boundaries */
    for (i=0; i<NX; i++)
    {
        gradient[i*NY] = 0.0;
        gradient[NX*NY+i*NY] = 0.0;

        gradient[i*NY+NY-1] = 0.0;
        gradient[NX*NY+i*NY+NY-1] = 0.0;
    }
    
    for (j=1; j<NY-1; j++)
    {
        gradient[j] = 0.0;
        gradient[(NX-1)*NY+j] = 0.0;

        gradient[NX*NY+j] = 0.0;
        gradient[NX*NY+(NX-1)*NY+j] = 0.0;
    }
}


void compute_gradient_rde(double phi[NX*NY], t_rde rde[NX*NY])
/* compute the gradient of the field */
{
    int i, j, iplus, iminus, jplus, jminus; 
    double deltaphi;
    double dx = (XMAX-XMIN)/((double)NX);
    
    dx = (XMAX-XMIN)/((double)NX);
    
    #pragma omp parallel for private(i,j,iplus,iminus,jplus,jminus,deltaphi)
    for (i=0; i<NX; i++)
        for (j=0; j<NY; j++)
        {
            iplus = i+1;  if (iplus == NX) iplus = 0;
            iminus = i-1; if (iminus == -1) iminus = NX-1;
            jplus = j+1;  if (jplus == NX) jplus = 0;
            jminus = j-1; if (jminus == -1) jminus = NY-1;
            
            deltaphi = phi[iplus*NY+j] - phi[iminus*NY+j];
            if (deltaphi < -PI) deltaphi += DPI; 
            if (deltaphi > PI) deltaphi -= DPI; 
            if (vabs(deltaphi) < 1.0e9) rde[i*NY+j].nablax = (deltaphi)/dx;
            else rde[i*NY+j].nablax = 0.0;
            
            deltaphi = phi[i*NY+jplus] - phi[i*NY+jminus];
            if (deltaphi < -PI) deltaphi += DPI; 
            if (deltaphi > PI) deltaphi -= DPI; 
            if (vabs(deltaphi) < 1.0e9) rde[i*NY+j].nablay = (deltaphi)/dx;
            else rde[i*NY+j].nablay = 0.0;
            
            /* TO DO: improve this computation */
            rde[i*NY+j].field_norm = module2(rde[i*NY+j].nablax,rde[i*NY+j].nablay);
            rde[i*NY+j].field_arg = argument(rde[i*NY+j].nablax,rde[i*NY+j].nablay);
        }
}

void compute_gradient_theta(t_rde rde[NX*NY])
/* compute the gradient of the theta field */
{
    int i, j, iplus, iminus, jplus, jminus; 
    double deltaphi;
    double dx = (XMAX-XMIN)/((double)NX);
    
    dx = (XMAX-XMIN)/((double)NX);
    
    #pragma omp parallel for private(i,j,iplus,iminus,jplus,jminus,deltaphi)
    for (i=0; i<NX; i++)
        for (j=0; j<NY; j++)
        {
            iplus = i+1;  if (iplus == NX) iplus = 0;
            iminus = i-1; if (iminus == -1) iminus = NX-1;
            jplus = j+1;  if (jplus == NX) jplus = 0;
            jminus = j-1; if (jminus == -1) jminus = NY-1;
            
            deltaphi = rde[iplus*NY+j].theta - rde[iminus*NY+j].theta;
            if (deltaphi < -PI) deltaphi += DPI; 
            if (deltaphi > PI) deltaphi -= DPI; 
            if (vabs(deltaphi) < 1.0e9) rde[i*NY+j].nablax = (deltaphi)/dx;
            else rde[i*NY+j].nablax = 0.0;
            
            deltaphi = rde[i*NY+jplus].theta - rde[i*NY+jminus].theta;
            if (deltaphi < -PI) deltaphi += DPI; 
            if (deltaphi > PI) deltaphi -= DPI; 
            if (vabs(deltaphi) < 1.0e9) rde[i*NY+j].nablay = (deltaphi)/dx;
            else rde[i*NY+j].nablay = 0.0;
        }
}

void compute_curl(t_rde rde[NX*NY])
/* compute the curl of the field */
{
    int i, j, iplus, iminus, jplus, jminus; 
    double delta;
    double dx = (XMAX-XMIN)/((double)NX);
    
    dx = (XMAX-XMIN)/((double)NX);
    
    #pragma omp parallel for private(i,j,iplus,iminus,jplus,jminus,delta)
    for (i=0; i<NX; i++)
        for (j=0; j<NY; j++)
        {
            iplus = i+1;  if (iplus == NX) iplus = 0;
            iminus = i-1; if (iminus == -1) iminus = NX-1;
            jplus = j+1;  if (jplus == NX) jplus = 0;
            jminus = j-1; if (jminus == -1) jminus = NY-1;
            
            delta = (rde[i*NY+jplus].nablay - rde[i*NY+jminus].nablay - (rde[iplus*NY+j].nablax - rde[iminus*NY+j].nablax))/dx;
            
            if (vabs(delta)*CURL_SCALE < 1.0e8)  rde[i*NY+j].curl = CURL_SCALE*delta;
            else rde[i*NY+j].curl = 0.0;
        }
}

void compute_gradient_polar(t_rde rde[NX*NY], double factor)
/* compute the norm of the gradient field */
{
    int i, j; 
    double angle;
    
    #pragma omp parallel for private(i,j,angle)
    for (i=0; i<NX; i++)
        for (j=0; j<NY; j++)
        {
            rde[i*NY+j].field_norm = factor*module2(rde[i*NY+j].nablax, rde[i*NY+j].nablay);
            angle = argument(rde[i*NY+j].nablax, rde[i*NY+j].nablay) + PI*COLOR_PHASE_SHIFT;
            if (angle < 0.0) angle += DPI;
            if (angle >= DPI) angle -= DPI;
            rde[i*NY+j].field_arg = angle;
        }
}

void compute_field_module(double *phi[NFIELDS], t_rde rde[NX*NY], double factor)
/* compute the norm squared of first two fields */
{
    int i, j; 
    
    #pragma omp parallel for private(i,j)
    for (i=0; i<NX; i++)
        for (j=0; j<NY; j++)
        {
            rde[i*NY+j].field_norm = factor*(phi[0][i*NY+j]*phi[0][i*NY+j] + phi[1][i*NY+j]*phi[1][i*NY+j]);
        }
}

void compute_field_argument(double *phi[NFIELDS], t_rde rde[NX*NY])
/* compute the norm squared of first two fields */
{
    int i, j; 
    double arg;
    
    #pragma omp parallel for private(i,j,arg)
    for (i=0; i<NX; i++)
        for (j=0; j<NY; j++)
        {
            arg = argument(phi[0][i*NY+j], phi[1][i*NY+j]) + COLOR_PHASE_SHIFT*PI;
            while (arg < 0.0) arg += DPI;
            while (arg >= DPI) arg -= DPI;
            rde[i*NY+j].field_arg = arg;
        }
}

void compute_probabilities(t_rde rde[NX*NY], short int xy_in[NX*NY], double probas[2])
/* compute probabilities for Ehrenfest urns */
{
    int i, j;
    double pleft = 0.0, pright = 0.0, sum;
    
    #pragma omp parallel for private(j)
    for (j=0; j<NY; j++)
    {
        for (i=0; i<NX/2; i++) if (xy_in[i*NY+j]) pleft += rde[i*NY+j].field_norm;
        for (i=NX/2; i<NX; i++) if (xy_in[i*NY+j]) pright += rde[i*NY+j].field_norm;
    }
        
    sum = pleft + pright;
    probas[0] = pleft/sum;
    probas[1] = pright/sum;
}

// void compute_field_norm(double phi_x[NX*NY], double phi_y[NX*NY], double phi_norm[NX*NY], double factor)
// /* compute the norm of (phi_x, phi_y) */
// {
//     int i, j; 
//     
//     #pragma omp parallel for private(i,j)
//     for (i=0; i<NX; i++)
//         for (j=0; j<NY; j++)
//         {
//             phi_norm[i*NY+j] = factor*module2(phi_x[i*NY+j],phi_y[i*NY+j]);
//         }
// }


// void compute_field_argument(double phi_x[NX*NY], double phi_y[NX*NY], double phi_arg[NX*NY])
// /* compute the argument of (phi_x, phi_y) */
// {
//     int i, j; 
//     double angle;
//     
//     #pragma omp parallel for private(i,j, angle)
//     for (i=0; i<NX; i++)
//         for (j=0; j<NY; j++)
//         {
//             angle = argument(phi_x[i*NY+j],phi_y[i*NY+j]) + PI*COLOR_PHASE_SHIFT;
//             if (angle < 0.0) angle += DPI;
//             if (angle >= DPI) angle -= DPI;
//             phi_arg[i*NY+j] = 360.0*angle/DPI;
// //             phi_arg[i*NY+j] = 180.0*angle/DPI;
//         }
// }



void compute_laplacian_rde(double phi_in[NX*NY], double phi_out[NX*NY], short int xy_in[NX*NY])
/* computes the Laplacian of phi_in and stores it in phi_out */
{
    int i, j, iplus, iminus, jplus, jminus;
    
    #pragma omp parallel for private(i,j,iplus,iminus,jplus,jminus)
    for (i=1; i<NX-1; i++){
        for (j=1; j<NY-1; j++){
            if (xy_in[i*NY+j]){
                phi_out[i*NY+j] = phi_in[(i+1)*NY+j] + phi_in[(i-1)*NY+j] 
                + phi_in[i*NY+j+1] + phi_in[i*NY+j-1] - 4.0*phi_in[i*NY+j];
            }
        }
    }
    
    /* boundary conditions */
    switch (B_COND) {
        case (BC_PERIODIC):
        {
            /* left and right side */
            for (j = 0; j < NY; j++) 
            {
                jplus = j+1;  if (jplus == NX) jplus = 0;
                jminus = j-1; if (jminus == -1) jminus = NY-1;
                
                phi_out[j] = phi_in[jminus] + phi_in[jplus] + phi_in[(NX-1)*NY+j] + phi_in[NY+j] - 4.0*phi_in[j];
                phi_out[(NX-1)*NY+j] = phi_in[(NX-1)*NY+jminus] + phi_in[(NX-1)*NY+jplus] + phi_in[(NX-2)*NY+j] + phi_in[j] - 4.0*phi_in[(NX-1)*NY+j];
            }
            /* top and bottom side */
            for (i = 1; i < NX-1; i++) 
            {
                iplus = i+1;  /*if (iplus == NX) iplus = 0;*/
                iminus = i-1; /*if (iminus == -1) iminus = NX-1;*/
                
                phi_out[i*NY] = phi_in[iminus*NY] + phi_in[iplus*NY] + phi_in[i*NY+1] + phi_in[i*NY+NY-1] - 4.0*phi_in[i*NY];
                phi_out[i*NY+NY-1] = phi_in[iminus*NY+NY-1] + phi_in[iplus*NY+NY-1] + phi_in[i*NY] + phi_in[i*NY+NY-2] - 4.0*phi_in[i*NY+NY-1];
            }
        }
        case (BC_DIRICHLET):
        {
            /* TO DO */
            break;
        }
    }
}

void compute_light_angle_rde(short int xy_in[NX*NY], t_rde rde[NX*NY], double potential[NX*NY], int movie)
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
                gradx = (*rde[(i+1)*NY+j].p_zfield[movie] - *rde[(i-1)*NY+j].p_zfield[movie])/dx;
                grady = (*rde[i*NY+j+1].p_zfield[movie] - *rde[i*NY+j-1].p_zfield[movie])/dy;
                
                /* case where the potential is added to the z coordinate */
                if ((ADD_POTENTIAL)&&(ADD_POTENTIAL_TO_Z))
                {
                    gradx += ADD_POT_CONSTANT*(potential[(i+1)*NY+j] - potential[(i-1)*NY+j])/dx;
                    grady += ADD_POT_CONSTANT*(potential[i*NY+j+1] - potential[i*NY+j-1])/dx;
                }
                
                norm = sqrt(1.0 + gradx*gradx + grady*grady);
                pscal = -gradx*light[0] - grady*light[1] + 1.0;
                
                rde[i*NY+j].cos_angle = pscal/norm;
            }
        }
}


void draw_field_line(double x, double y, short int *xy_in[NX], double *nablax[NX], 
                     double *nablay[NX], double delta, int nsteps)
/* draw a field line of the gradient, starting in (x,y) - OLD VERSION */
{
    double x1, y1, x2, y2, pos[2], nabx, naby, norm2, norm;
    int i = 0, ij[2], cont = 1;
    
    glColor3f(1.0, 1.0, 1.0);
//     glColor3f(0.0, 0.0, 0.0);
    glLineWidth(FIELD_LINE_WIDTH);
    x1 = x;
    y1 = y;
    
//     printf("Drawing field line \n");

    glEnable(GL_LINE_SMOOTH);
    glBegin(GL_LINE_STRIP);
    xy_to_pos(x1, y1, pos);
    glVertex2d(pos[0], pos[1]);
    
    i = 0;
    while ((cont)&&(i < nsteps))
    {
        xy_to_ij(x1, y1, ij);
        
        if (ij[0] < 0) ij[0] = 0;
        if (ij[0] > NX-1) ij[0] = NX-1;
        if (ij[1] < 0) ij[1] = 0;
        if (ij[1] > NY-1) ij[1] = NY-1;
        
        nabx = nablax[ij[0]][ij[1]];
        naby = nablay[ij[0]][ij[1]];
        
        norm2 = nabx*nabx + naby*naby;
        
        if (norm2 > 1.0e-14)
        {
            /* avoid too large step size */
            if (norm2 < 1.0e-9) norm2 = 1.0e-9;
            norm = sqrt(norm2);
            x1 = x1 + delta*nabx/norm;
            y1 = y1 + delta*naby/norm;
        }
        else cont = 0;
        
        if (!xy_in[ij[0]*NY + ij[1]]) cont = 0;
        
        /* stop if the boundary is hit */
//         if (xy_in[ij[0]][ij[1]] != 1) cont = 0;
        
//         printf("x1 = %.3lg \t y1 = %.3lg \n", x1, y1);
                
        xy_to_pos(x1, y1, pos);
        glVertex2d(pos[0], pos[1]);
        
        i++;
    }
    glEnd();
}



void compute_field_color_rde(double value, int cplot, int palette, double rgb[3])
/* compute the color depending on the field value and color palette */
{
    switch (cplot) {
        case (Z_AMPLITUDE): 
        {
            color_scheme_palette(COLOR_SCHEME, palette, value, 1.0, 0, rgb);
            break;
        }
        case (Z_RGB): 
        {
            /* TO DO */
            break;
        }
        case (Z_POLAR): 
        {
//             hsl_to_rgb_palette(360.0*value/DPI, 0.9, 0.5, rgb, palette);
            color_scheme_palette(C_ONEDIM_LINEAR, palette, value/DPI, 1.0, 1, rgb);
            break;
        }
        case (Z_NORM_GRADIENT): 
        {
//             color_scheme_palette(COLOR_SCHEME, palette, value, 1.0, 0, rgb);
            color_scheme_asym_palette(COLOR_SCHEME, palette, value, 1.0, 0, rgb);
            break;
        }
        case (Z_ANGLE_GRADIENT): 
        {
            color_scheme_palette(C_ONEDIM_LINEAR, palette, value/DPI, 1.0, 1, rgb);
            break;
        }
        case (Z_NORM_GRADIENTX): 
        {
            hsl_to_rgb_palette(360.0*value/DPI, 0.9, 0.5, rgb, palette);
            break;
        }
        case (Z_ANGLE_GRADIENTX): 
        {
            color_scheme_palette(C_ONEDIM_LINEAR, palette, value/DPI, 1.0, 1, rgb);
            break;
        }
        case (Z_NORM_GRADIENT_INTENSITY): 
        {
            hsl_to_rgb_palette(360.0*value/DPI, 0.9, 0.5, rgb, palette);
            break;
        }
        case (Z_VORTICITY): 
        {
            hsl_to_rgb_palette(180.0*(1.0 - color_amplitude(value, 1.0, 0)), 0.9, 0.5, rgb, palette);
            break;
        }
        case (Z_VORTICITY_ABS): 
        {
            hsl_to_rgb_palette(360.0*(1.0 - vabs(color_amplitude(value, 1.0, 0))), 0.9, 0.5, rgb, palette);
            break;
        }
        case (Z_MAXTYPE_RPSLZ): 
        {
            hsl_to_rgb_palette(value, 0.9, 0.5, rgb, palette);
            break;
        }
        case (Z_THETA_RPSLZ): 
        {
            color_scheme_palette(C_ONEDIM_LINEAR, palette, value/DPI, 1.0, 1, rgb);
            break;
        }
        case (Z_NORM_GRADIENT_RPSLZ): 
        {
            color_scheme_palette(COLOR_SCHEME, palette, value, 1.0, 0, rgb);
            break;
        }
        case (Z_MODULE): 
        {
            color_scheme_asym_palette(COLOR_SCHEME, palette, VSCALE_AMPLITUDE*value, 1.0, 0, rgb);
            break;
        }
        case (Z_ARGUMENT): 
        {
//             color_scheme_palette(C_ONEDIM_LINEAR, palette, value/DPI, 1.0, 1, rgb);
            hsl_to_rgb_palette(360.0*value/DPI, 0.9, 0.5, rgb, palette);
            break;
        }
        case (Z_REALPART): 
        {
            color_scheme_palette(COLOR_SCHEME, palette, VSCALE_AMPLITUDE*value, 1.0, 0, rgb);
            break;
        }
    }
}

double adjust_field(double z, double pot)
/* add potential in case of option ADD_POTENTIAL_TO_Z */
{
    if ((ADD_POTENTIAL)&&(ADD_POTENTIAL_TO_Z)) return (z + ADD_POT_CONSTANT*pot);
    else return(z);
}  


double compute_interpolated_colors_rde(int i, int j, t_rde rde[NX*NY], double potential[NX*NY], double palette, int cplot, 
                                 double rgb_e[3], double rgb_w[3], double rgb_n[3], double rgb_s[3],
                                 double *z_sw, double *z_se, double *z_nw, double *z_ne, 
                                 int fade, double fade_value, int movie)
{
    int k;
    double cw, ce, cn, cs, c_sw, c_se, c_nw, c_ne, c_mid, ca, z_mid;
    double lum;
    
    *z_sw = *rde[i*NY+j].p_zfield[movie];
    *z_se = *rde[(i+1)*NY+j].p_zfield[movie];
    *z_nw = *rde[i*NY+j+1].p_zfield[movie];
    *z_ne = *rde[(i+1)*NY+j+1].p_zfield[movie];
    
    if ((ADD_POTENTIAL)&&(ADD_POTENTIAL_TO_Z))
    {
        *z_sw += ADD_POT_CONSTANT*potential[i*NY+j];
        *z_se += ADD_POT_CONSTANT*potential[(i+1)*NY+j];
        *z_nw += ADD_POT_CONSTANT*potential[i*NY+j+1];
        *z_ne += ADD_POT_CONSTANT*potential[(i+1)*NY+j+1];
    }
                            
    z_mid = 0.25*(*z_sw + *z_se + *z_nw + *z_ne);
    
    c_sw = *rde[i*NY+j].p_cfield[movie];
    c_se = *rde[(i+1)*NY+j].p_cfield[movie];
    c_nw = *rde[i*NY+j+1].p_cfield[movie];
    c_ne = *rde[(i+1)*NY+j+1].p_cfield[movie];
    
    if (COMPUTE_WRAP_ANGLE)
    {
        c_se = unwrap_angle(c_sw,c_se);
        c_nw = unwrap_angle(c_sw,c_nw);
        c_ne = unwrap_angle(c_sw,c_ne);
    }
    
    c_mid = 0.25*(c_sw + c_se + c_nw + c_ne);
                                    
    cw = (c_sw + c_nw + c_mid)/3.0;
    ce = (c_se + c_ne + c_mid)/3.0;
    cs = (c_sw + c_se + c_mid)/3.0;
    cn = (c_nw + c_ne + c_mid)/3.0;
    
    if (COMPUTE_WRAP_ANGLE)
    {
        cw = wrap_angle(cw);
        ce = wrap_angle(ce);
        cs = wrap_angle(cs);
        cn = wrap_angle(cn);
    }
    
    compute_field_color_rde(ce, cplot, palette, rgb_e);
    compute_field_color_rde(cw, cplot, palette, rgb_w);
    compute_field_color_rde(cn, cplot, palette, rgb_n);
    compute_field_color_rde(cs, cplot, palette, rgb_s);
    
//     if ((cplot == Z_ARGUMENT)||(cplot == Z_REALPART))
    if (cplot == Z_ARGUMENT)
    {
        lum = tanh(SLOPE_SCHROD_LUM*rde[i*NY+j].field_norm) + MIN_SCHROD_LUM;
        for (k=0; k<3; k++) 
        {
            rgb_e[k] *= lum;
            rgb_w[k] *= lum;
            rgb_n[k] *= lum;
            rgb_s[k] *= lum;
        }
    }
    if (SHADE_3D)
    {
        ca = rde[i*NY+j].cos_angle;
        ca = (ca + 1.0)*0.4 + 0.2;
        for (k=0; k<3; k++) 
        {
            rgb_e[k] *= ca;
            rgb_w[k] *= ca;
            rgb_n[k] *= ca;
            rgb_s[k] *= ca;
        }
    }
    if (fade)
        for (k=0; k<3; k++) 
        {
            rgb_e[k] *= fade_value;
            rgb_w[k] *= fade_value;
            rgb_n[k] *= fade_value;
            rgb_s[k] *= fade_value;
        }
    
    return(z_mid);
}


void compute_rde_fields(double *phi[NFIELDS], short int xy_in[NX*NY], int zplot, int cplot, t_rde rde[NX*NY])
/* compute the necessary auxiliary fields */
{
    int i, j;
    
    switch (RDE_EQUATION) {
        case (E_RPS):
        {
            if ((COMPUTE_THETA)||(COMPUTE_THETAZ))
                compute_theta(phi, xy_in, rde);
    
            if ((zplot == Z_NORM_GRADIENT)||(cplot == Z_NORM_GRADIENT))
            {
                compute_gradient_theta(rde);
//              compute_gradient_polar(rde, 1.0);
                compute_gradient_polar(rde, 0.003);
            }
    
            if ((zplot == Z_NORM_GRADIENTX)||(cplot == Z_NORM_GRADIENTX)||(zplot == Z_ANGLE_GRADIENTX)||(cplot == Z_ANGLE_GRADIENTX))
            {
                compute_gradient_rde(phi[0], rde);
                compute_gradient_polar(rde, 0.03);
            }
    
            if ((zplot == Z_VORTICITY)||(cplot == Z_VORTICITY)||(zplot == Z_VORTICITY_ABS)||(cplot == Z_VORTICITY_ABS))
            {
                compute_gradient_theta(rde);
                compute_curl(rde);
            }
            break;
        }
        case (E_RPSLZ):
        {
            compute_theta_rpslz(phi, xy_in, rde, cplot);
            if ((zplot == Z_NORM_GRADIENT_RPSLZ)||(cplot == Z_NORM_GRADIENT_RPSLZ))
            {
                compute_gradient_theta(rde);
                compute_gradient_polar(rde, 0.005);
            }
            break;
        }
        case (E_SCHRODINGER):
        {
            compute_field_module(phi, rde, 1.0);
            if ((zplot == Z_ARGUMENT)||(cplot == Z_ARGUMENT))
            {
                compute_field_argument(phi, rde);
            }
            break;
        }
        default : break;
    }
}

void init_zfield_rde(double *phi[NFIELDS], short int xy_in[NX*NY], int zplot, t_rde rde[NX*NY], int movie)
/* compute the necessary fields for the z coordinate */
{
    int i, j;
    
    switch(zplot) {
        case (Z_AMPLITUDE):
        {
            #pragma omp parallel for private(i,j)
            for (i=0; i<NX; i++) for (j=0; j<NY; j++) rde[i*NY+j].p_zfield[movie] = &phi[0][i*NY+j];
            break;
        }
        case (Z_POLAR):
        {
            #pragma omp parallel for private(i,j)
            for (i=0; i<NX; i++) for (j=0; j<NY; j++) rde[i*NY+j].p_zfield[movie] = &rde[i*NY+j].theta;
            break;
        }
        case (Z_NORM_GRADIENT):
        {
            #pragma omp parallel for private(i,j)
            for (i=0; i<NX; i++) for (j=0; j<NY; j++) rde[i*NY+j].p_zfield[movie] = &rde[i*NY+j].field_norm;
            break;
        }
        case (Z_ANGLE_GRADIENT):
        {
            #pragma omp parallel for private(i,j)
            for (i=0; i<NX; i++) for (j=0; j<NY; j++) rde[i*NY+j].p_zfield[movie] = &rde[i*NY+j].field_arg;
            break;
        }
        case (Z_NORM_GRADIENTX):
        {
            #pragma omp parallel for private(i,j)
            for (i=0; i<NX; i++) for (j=0; j<NY; j++) rde[i*NY+j].p_zfield[movie] = &rde[i*NY+j].field_norm;
            break;
        }
        case (Z_ANGLE_GRADIENTX):
        {
            #pragma omp parallel for private(i,j)
            for (i=0; i<NX; i++) for (j=0; j<NY; j++) rde[i*NY+j].p_zfield[movie] = &rde[i*NY+j].field_arg;
            break;
        }
        case (Z_NORM_GRADIENT_INTENSITY):
        {
            #pragma omp parallel for private(i,j)
            for (i=0; i<NX; i++) for (j=0; j<NY; j++) rde[i*NY+j].p_zfield[movie] = &rde[i*NY+j].field_norm;
            break;
        }
        case (Z_VORTICITY):
        {
            #pragma omp parallel for private(i,j)
            for (i=0; i<NX; i++) for (j=0; j<NY; j++) rde[i*NY+j].p_zfield[movie] = &rde[i*NY+j].curl;
            break;
        }
        case (Z_VORTICITY_ABS):
        {
            #pragma omp parallel for private(i,j)
            for (i=0; i<NX; i++) for (j=0; j<NY; j++) rde[i*NY+j].p_zfield[movie] = &rde[i*NY+j].curl;
            break;
        }
        case (Z_MAXTYPE_RPSLZ):
        {
            #pragma omp parallel for private(i,j)
            for (i=0; i<NX; i++) for (j=0; j<NY; j++) rde[i*NY+j].p_zfield[movie] = &rde[i*NY+j].theta;
            break;
        }
        case (Z_THETA_RPSLZ):
        {
            #pragma omp parallel for private(i,j)
            for (i=0; i<NX; i++) for (j=0; j<NY; j++) rde[i*NY+j].p_zfield[movie] = &rde[i*NY+j].theta;
            break;
        }
        case (Z_NORM_GRADIENT_RPSLZ):
        {
            #pragma omp parallel for private(i,j)
            for (i=0; i<NX; i++) for (j=0; j<NY; j++) rde[i*NY+j].p_zfield[movie] = &rde[i*NY+j].field_norm;
            break;
        }
        case (Z_MODULE):
        {
            #pragma omp parallel for private(i,j)
            for (i=0; i<NX; i++) for (j=0; j<NY; j++) rde[i*NY+j].p_zfield[movie] = &rde[i*NY+j].field_norm;
            break;
        }
        case (Z_ARGUMENT):
        {
            #pragma omp parallel for private(i,j)
            for (i=0; i<NX; i++) for (j=0; j<NY; j++) rde[i*NY+j].p_zfield[movie] = &rde[i*NY+j].field_arg;
            break;
        }
        case (Z_REALPART):
        {
            #pragma omp parallel for private(i,j)
            for (i=0; i<NX; i++) for (j=0; j<NY; j++) rde[i*NY+j].p_zfield[movie] = &phi[0][i*NY+j];
            break;
        }
    }
}

void init_cfield_rde(double *phi[NFIELDS], short int xy_in[NX*NY], int cplot, t_rde rde[NX*NY], int movie)
/* compute the colors */
{
    int i, j, k;
    
    switch(cplot) {
        case (Z_AMPLITUDE):
        {
            #pragma omp parallel for private(i,j)
            for (i=0; i<NX; i++) for (j=0; j<NY; j++) rde[i*NY+j].p_cfield[movie] = &phi[0][i*NY+j];
            break;
        }
        case (Z_POLAR):
        {
            #pragma omp parallel for private(i,j)
            for (i=0; i<NX; i++) for (j=0; j<NY; j++) rde[i*NY+j].p_cfield[movie] = &rde[i*NY+j].theta;
            break;
        }
        case (Z_NORM_GRADIENT):
        {
            #pragma omp parallel for private(i,j)
            for (i=0; i<NX; i++) for (j=0; j<NY; j++) rde[i*NY+j].p_cfield[movie] = &rde[i*NY+j].field_norm;
            break;
        }
        case (Z_ANGLE_GRADIENT):
        {
            #pragma omp parallel for private(i,j)
            for (i=0; i<NX; i++) for (j=0; j<NY; j++) rde[i*NY+j].p_cfield[movie] = &rde[i*NY+j].field_arg;
            break;
        }
        case (Z_NORM_GRADIENTX):
        {
            #pragma omp parallel for private(i,j)
            for (i=0; i<NX; i++) for (j=0; j<NY; j++) rde[i*NY+j].p_cfield[movie] = &rde[i*NY+j].field_norm;
            break;
        }
        case (Z_ANGLE_GRADIENTX):
        {
            #pragma omp parallel for private(i,j)
            for (i=0; i<NX; i++) for (j=0; j<NY; j++) rde[i*NY+j].p_cfield[movie] = &rde[i*NY+j].field_arg;
            break;
        }
        case (Z_NORM_GRADIENT_INTENSITY):
        {
            #pragma omp parallel for private(i,j)
            for (i=0; i<NX; i++) for (j=0; j<NY; j++) rde[i*NY+j].p_cfield[movie] = &rde[i*NY+j].field_norm;
            break;
        }
        case (Z_VORTICITY):
        {
            #pragma omp parallel for private(i,j)
            for (i=0; i<NX; i++) for (j=0; j<NY; j++) rde[i*NY+j].p_cfield[movie] = &rde[i*NY+j].curl;
            break;
        }
        case (Z_VORTICITY_ABS):
        {
            #pragma omp parallel for private(i,j)
            for (i=0; i<NX; i++) for (j=0; j<NY; j++) rde[i*NY+j].p_cfield[movie] = &rde[i*NY+j].curl;
            break;
        }
        case (Z_MAXTYPE_RPSLZ):
        {
            #pragma omp parallel for private(i,j)
            for (i=0; i<NX; i++) for (j=0; j<NY; j++) rde[i*NY+j].p_cfield[movie] = &rde[i*NY+j].theta;
            break;
        }
        case (Z_THETA_RPSLZ):
        {
            #pragma omp parallel for private(i,j)
            for (i=0; i<NX; i++) for (j=0; j<NY; j++) rde[i*NY+j].p_cfield[movie] = &rde[i*NY+j].theta;
            break;
        }
        case (Z_NORM_GRADIENT_RPSLZ):
        {
            #pragma omp parallel for private(i,j)
            for (i=0; i<NX; i++) for (j=0; j<NY; j++) rde[i*NY+j].p_cfield[movie] = &rde[i*NY+j].field_norm;
            break;
        }
        case (Z_MODULE):
        {
            #pragma omp parallel for private(i,j)
            for (i=0; i<NX; i++) for (j=0; j<NY; j++) rde[i*NY+j].p_cfield[movie] = &rde[i*NY+j].field_norm;
            break;
        }
        case (Z_ARGUMENT):
        {
            #pragma omp parallel for private(i,j)
            for (i=0; i<NX; i++) for (j=0; j<NY; j++) rde[i*NY+j].p_cfield[movie] = &rde[i*NY+j].field_arg;
            break;
        }
        case (Z_REALPART):
        {
            #pragma omp parallel for private(i,j)
            for (i=0; i<NX; i++) for (j=0; j<NY; j++) rde[i*NY+j].p_cfield[movie] = &phi[0][i*NY+j];
            break;
        }
    }
}

void compute_cfield_rde(short int xy_in[NX*NY], int cplot, int palette, t_rde rde[NX*NY], int fade, double fade_value, int movie)
/* compute the colors */
{
    int i, j, k;
    double ca, lum;
       
    #pragma omp parallel for private(i,j,k,ca,lum)
    for (i=0; i<NX; i++) for (j=0; j<NY; j++)
    {
        compute_field_color_rde(*rde[i*NY+j].p_cfield[movie], cplot, palette, rde[i*NY+j].rgb);
//         if ((cplot == Z_ARGUMENT)||(cplot == Z_REALPART))
        if (cplot == Z_ARGUMENT)
        {
            lum = tanh(SLOPE_SCHROD_LUM*rde[i*NY+j].field_norm);
            for (k=0; k<3; k++) rde[i*NY+j].rgb[k] *= lum;
        }
        if (SHADE_3D)
        {
            ca = rde[i*NY+j].cos_angle;
            ca = (ca + 1.0)*0.4 + 0.2;
            if ((FADE_IN_OBSTACLE)&&(!xy_in[i*NY+j])) ca *= 1.6;
            for (k=0; k<3; k++) rde[i*NY+j].rgb[k] *= ca;
        }
        if (fade) for (k=0; k<3; k++) rde[i*NY+j].rgb[k] *= fade_value;
    }
}



void draw_wave_2d_rde(short int xy_in[NX*NY], t_rde rde[NX*NY])
{
    int i, j;
    
    glBegin(GL_QUADS);
    
    for (i=0; i<NX; i++)
        for (j=0; j<NY; j++)
        {
            if (xy_in[i*NY+j])
            {
                glColor3f(rde[i*NY+j].rgb[0], rde[i*NY+j].rgb[1], rde[i*NY+j].rgb[2]);
                
                glVertex2i(i, j);
                glVertex2i(i+1, j);
                glVertex2i(i+1, j+1);
                glVertex2i(i, j+1);
            }
        }
    glEnd ();
    if (DRAW_BILLIARD) draw_billiard(0, 1.0);
}


void draw_wave_3d_ij_rde(int i, int j, int movie, double *phi[NFIELDS], short int xy_in[NX*NY], t_rde rde[NX*NY], 
                         double potential[NX*NY], int zplot, int cplot, int palette, int fade, double fade_value)
{
    int k, l, draw = 1;
    double xy[2], xy_screen[2], rgb[3], pos[2], ca, rgb_e[3], rgb_w[3], rgb_n[3], rgb_s[3]; 
    double z, z_sw, z_se, z_nw, z_ne, z_mid, zw, ze, zn, zs, min = 1000.0, max = 0.0;
    double xy_sw[2], xy_se[2], xy_nw[2], xy_ne[2], xy_mid[2];
    double energy;
    

    if (NON_DIRICHLET_BC) 
        draw = (xy_in[i*NY+j])&&(xy_in[(i+1)*NY+j])&&(xy_in[i*NY+j+1])&&(xy_in[(i+1)*NY+j+1]);
    else draw = (TWOSPEEDS)||(xy_in[i*NY+j]);
            
    if (draw)
    {
        if (AMPLITUDE_HIGH_RES > 0)
        {
            z_mid = compute_interpolated_colors_rde(i, j, rde, potential, palette, cplot, 
                                                        rgb_e, rgb_w, rgb_n, rgb_s, &z_sw, &z_se, &z_nw, &z_ne, 
                                                        fade, fade_value, movie);
            ij_to_xy(i, j, xy_sw);
            ij_to_xy(i+1, j, xy_se);
            ij_to_xy(i, j+1, xy_nw);
            ij_to_xy(i+1, j+1, xy_ne);
                    
            for (k=0; k<2; k++) xy_mid[k] = 0.25*(xy_sw[k] + xy_se[k] + xy_nw[k] + xy_ne[k]);
                       
            if (AMPLITUDE_HIGH_RES == 1)
            {                        
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
        }
        else
        {
            glColor3f(rde[i*NY+j].rgb[0], rde[i*NY+j].rgb[1], rde[i*NY+j].rgb[2]);
            
            glBegin(GL_TRIANGLE_FAN);
            ij_to_xy(i, j, xy);
            draw_vertex_xyz(xy, adjust_field(*rde[i*NY+j].p_zfield[movie], potential[i*NY+j]));
            ij_to_xy(i+1, j, xy);
            draw_vertex_xyz(xy, adjust_field(*rde[(i+1)*NY+j].p_zfield[movie], potential[(i+1)*NY+j]));
            ij_to_xy(i+1, j+1, xy);
            draw_vertex_xyz(xy, adjust_field(*rde[(i+1)*NY+j+1].p_zfield[movie], potential[(i+1)*NY+j+1]));
            ij_to_xy(i, j+1, xy);
            draw_vertex_xyz(xy, adjust_field(*rde[i*NY+j+1].p_zfield[movie], potential[i*NY+j+1]));
            glEnd ();
        }
    }
}



void draw_wave_3d_rde(int movie, double *phi[NFIELDS], short int xy_in[NX*NY], t_rde rde[NX*NY], double potential[NX*NY], 
                  int zplot, int cplot, int palette, int fade, double fade_value)
{
    int i, j;
    double observer_angle;
    
    blank();
    if (DRAW_BILLIARD) draw_billiard_3d(fade, fade_value);
    
    if (!ROTATE_VIEW)
    {
        for (i=0; i<NX-2; i++)
            for (j=0; j<NY-2; j++)
                draw_wave_3d_ij_rde(i, j, movie, phi, xy_in, rde, potential, zplot, cplot, palette, fade, fade_value);
    }
    else    /* draw facets in an order depending on the position of the observer */
    {
        observer_angle = argument(observer[0], observer[1]);
        observer_angle -= 0.5*PID;
        if (observer_angle < 0.0) observer_angle += DPI;
        printf("Observer_angle = %.3lg\n", observer_angle*360.0/DPI); 
        
        if ((observer_angle > 0.0)&&(observer_angle < PID))
        {
            for (j=0; j<NY-2; j++)
                for (i=0; i<NX-2; i++)
                    draw_wave_3d_ij_rde(i, j, movie, phi, xy_in, rde, potential, zplot, cplot, palette, fade, fade_value);
        }
        else if (observer_angle < PI)
        {
            for (i=NX-3; i>0; i--)
                for (j=0; j<NY-2; j++)
                    draw_wave_3d_ij_rde(i, j, movie, phi, xy_in, rde, potential, zplot, cplot, palette, fade, fade_value);
        }
        else if (observer_angle < 1.5*PI)
        {
             for (j=NY-3; j>0; j--)
                for (i=0; i<NX-2; i++)
                    draw_wave_3d_ij_rde(i, j, movie, phi, xy_in, rde, potential, zplot, cplot, palette, fade, fade_value);   
        }
        else
        {
            for (i=0; i<NX-2; i++)
                for (j=0; j<NY-2; j++)
                    draw_wave_3d_ij_rde(i, j, movie, phi, xy_in, rde, potential, zplot, cplot, palette, fade, fade_value);
        }
    }
            
    if (DRAW_BILLIARD_FRONT) draw_billiard_3d_front(fade, fade_value);
}


void draw_wave_rde(int movie, double *phi[NFIELDS], short int xy_in[NX*NY], t_rde rde[NX*NY], double potential[NX*NY], 
                   int zplot, int cplot, int palette, int fade, double fade_value, int refresh)
{
    int i, j, k, l, draw = 1;
    double xy[2], xy_screen[2], rgb[3], pos[2], ca, rgb_e[3], rgb_w[3], rgb_n[3], rgb_s[3]; 
    double z_sw, z_se, z_nw, z_ne, z_mid, zw, ze, zn, zs, min = 1000.0, max = 0.0;
    double xy_sw[2], xy_se[2], xy_nw[2], xy_ne[2], xy_mid[2];
    double energy;

    if (refresh)
    {
//         printf("Computing fields\n");
        compute_rde_fields(phi, xy_in, zplot, cplot, rde);      
//         printf("Computed fields\n");
        if ((PLOT_3D)&&(SHADE_3D)) compute_light_angle_rde(xy_in, rde, potential, movie);       
    }
    compute_cfield_rde(xy_in, cplot, palette, rde, fade, fade_value, movie);
    
    if (PLOT_3D) draw_wave_3d_rde(movie, phi, xy_in, rde, potential, zplot, cplot, palette, fade, fade_value);
    else draw_wave_2d_rde(xy_in, rde);
}


