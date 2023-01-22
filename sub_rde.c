double courant2;  /* Courant parameter squared */
double dx2;       /* spatial step size squared */
double intstep;   /* integration step */
double intstep1;  /* integration step used in absorbing boundary conditions */
double viscosity; /* viscosity (parameter in front of Laplacian) */
double rpslzb;    /* second parameter in Rock-Paper-Scissors-Lizard-Spock equation */
double flow_speed;  /* flow speed for laminar boundary conditions in Euler equation */

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

void add_coherent_state(double x, double y, double px, double py, double scalex, double *phi[NFIELDS], short int xy_in[NX*NY])
/* add to the field a coherent state of position (x,y) and momentum (px, py) */
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

void init_vortex_state(double amp, double x, double y, double scale, double density_mod, double *phi[NFIELDS], short int xy_in[NX*NY])
/* initialise field with vortex at position (x,y) with amplitude and variance scale */
/* for incompressible Euler, phi[0] is stream function, phi[1] is vorticity */
/* for compressible Euler, phi[0] is the density, phi[1] and phi[2] are velocity components */
{
    int i, j;
    double xy[2], dist2, module, phase, scale2, sign;    

//     scale2 = scale*scale;
    switch (RDE_EQUATION) {
        case (E_EULER_INCOMP): 
        {
            for (i=0; i<NX; i++)
                for (j=0; j<NY; j++)
                {
                    ij_to_xy(i, j, xy);
                    xy_in[i*NY+j] = xy_in_billiard(xy[0],xy[1]);

                    if (xy_in[i*NY+j])
                    {
                        dist2 = (xy[0]-x)*(xy[0]-x) + (xy[1]-y)*(xy[1]-y);
                        module = amp*exp(-dist2/vabs(scale));
                
                        if (scale > 0.0) sign = 1.0;
                        else sign = -1.0;
                
                        phi[1][i*NY+j] += sign*module;
                        phi[0][i*NY+j] -= sign*module;    /* approximate, stream function should solve Poisson equation */
                    }
                    else
                    {
                        phi[0][i*NY+j] = 0.0;
                        phi[1][i*NY+j] = 0.0;
                    }
                }
            break;
        }
        case (E_EULER_COMP):
        {
            for (i=0; i<NX; i++)
                for (j=0; j<NY; j++)
                {
                    ij_to_xy(i, j, xy);
                    xy_in[i*NY+j] = xy_in_billiard(xy[0],xy[1]);

                    if (xy_in[i*NY+j])
                    {
                        dist2 = (xy[0]-x)*(xy[0]-x) + (xy[1]-y)*(xy[1]-y);
                        module = amp*exp(-dist2/vabs(scale));
                
                        if (scale > 0.0) sign = 1.0;
                        else sign = -1.0;
                
//                         phi[0][i*NY+j] = 1.0;
                        /* nonconstant density to make things more interesting */
                        phi[0][i*NY+j] = 0.5 + density_mod*module/amp;
                        phi[1][i*NY+j] = -sign*module*(xy[1]-y)/vabs(scale);
                        phi[2][i*NY+j] = sign*module*(xy[0]-x)/vabs(scale);   
                    }
                    else
                    {
                        phi[0][i*NY+j] = 1.0;
                        phi[1][i*NY+j] = 0.0;
                        phi[2][i*NY+j] = 0.0;
                    }
                }
            break;
        }
    }
}

void add_vortex_state(double amp, double x, double y, double scale, double density_mod, double *phi[NFIELDS], short int xy_in[NX*NY])
/* add vortex at position (x,y) with variance scale to field */
/* for incompressible Euler, phi[0] is stream function, phi[1] is vorticity */
/* for compressible Euler, phi[0] is the density, phi[1] and phi[2] are velocity components */
{
    int i, j;
    double xy[2], dist2, module, phase, scale2, sign;    

//     scale2 = scale*scale;
    switch (RDE_EQUATION) {
        case (E_EULER_INCOMP): 
        {
            for (i=0; i<NX; i++)
                for (j=0; j<NY; j++)
                {
                    ij_to_xy(i, j, xy);
                    xy_in[i*NY+j] = xy_in_billiard(xy[0],xy[1]);

                    if (xy_in[i*NY+j])
                    {
                        dist2 = (xy[0]-x)*(xy[0]-x) + (xy[1]-y)*(xy[1]-y);
                        module = amp*exp(-dist2/vabs(scale));
                
                        if (scale > 0.0) sign = 1.0;
                        else sign = -1.0;
                
                        phi[1][i*NY+j] -= sign*module;
                        phi[0][i*NY+j] = sign*module;    /* approximate, stream function should solve Poisson equation */
                    }
                    else
                    {
                        phi[0][i*NY+j] = 0.0;
                        phi[1][i*NY+j] = 0.0;
                    }
                }
            break;
        }
        case (E_EULER_COMP):
        {
            for (i=0; i<NX; i++)
                for (j=0; j<NY; j++)
                {
                    ij_to_xy(i, j, xy);
                    xy_in[i*NY+j] = xy_in_billiard(xy[0],xy[1]);

                    if (xy_in[i*NY+j])
                    {
                        dist2 = (xy[0]-x)*(xy[0]-x) + (xy[1]-y)*(xy[1]-y);
                        module = amp*exp(-dist2/vabs(scale));
                
                        if (scale > 0.0) sign = 1.0;
                        else sign = -1.0;
                
//                         phi[0][i*NY+j] = 1.0;
                        /* nonconstant density to make things more interesting */
                        phi[0][i*NY+j] += 0.5 + density_mod*sign*module/amp;
                        phi[1][i*NY+j] += sign*module*(xy[1]-y)/vabs(scale);
                        phi[2][i*NY+j] -= sign*module*(xy[0]-x)/vabs(scale);   
                    }
                    else
                    {
                        phi[0][i*NY+j] = 1.0;
                        phi[1][i*NY+j] = 0.0;
                        phi[2][i*NY+j] = 0.0;
                    }
                }
            break;
        }
    }
}


void init_shear_flow(double amp, double delta, double rho, int nx, int ny, double yshift, double *phi[NFIELDS], short int xy_in[NX*NY])
/* initialise field with a shear flow */
/* phi[0] is stream function, phi[1] is vorticity */
/* amp is global amplitude */
/* delta is the amplitude of the periodic perturbation in the x direction */
/* rho controls the size of the transition zone in the y direction */
/* nx is number of oscillations in x direction */
/* ny is number of layers in y direction */
{
    int i, j;
    double xy[2], y1, a, b, f, cplus, cminus;    

    a = (double)nx*DPI/(XMAX - XMIN);
    b = 0.5;
    
    for (i=0; i<NX; i++)
        for (j=0; j<NY; j++)
        {
            ij_to_xy(i, j, xy);
	    xy_in[i*NY+j] = xy_in_billiard(xy[0],xy[1]);
            
            y1 = xy[1]*(double)ny/YMAX - yshift; 
            while (y1 > 1.0) y1 -= 2.0;
            while (y1 < -1.0) y1 += 2.0;

            if (xy_in[i*NY+j])
            {
                f = delta*cos(a*xy[0]);
                cplus = cosh((y1 + b)/rho);
                cminus = cosh((y1 - b)/rho);
                
                if (y1 > 0.0)
                {
                    phi[1][i*NY+j] = amp*(f + 1.0/(rho*cminus*cminus));
                    phi[0][i*NY+j] = amp*(f/(a*a) + rho/(cminus*cminus));    
                }
                else
                {
                    phi[1][i*NY+j] = amp*(f - 1.0/(rho*cplus*cplus));
                    phi[0][i*NY+j] = amp*(f/(a*a) - rho/(cplus*cplus));    
                }
                
                phi[2][i*NY+j] = 0.0;
            }
            else
            {
                phi[0][i*NY+j] = 0.0;
                phi[1][i*NY+j] = 0.0;
                phi[2][i*NY+j] = 0.0;
            }
        }
}

void set_boundary_laminar_flow(double amp, double modulation, double period, double yshift, double *phi[NFIELDS], short int xy_in[NX*NY], int imin, int imax, int jmin, int jmax)
/* enfoce laminar flow in x direction on top and bottom boundaries */
/* phi[0] is stream function, phi[1] is vorticity */
/* amp is global amplitude */
{
    int i, j;
    double xy[2], y1, a, b, f, cplus, cminus;    
    
    a = period*PI/YMAX;
    
    for (i=imin; i<imax; i++)
        for (j=jmin; j<jmax; j++)
        {
            ij_to_xy(i, j, xy);
	    xy_in[i*NY+j] = xy_in_billiard(xy[0],xy[1]);
            y1 = xy[1] + yshift;
            if (xy_in[i*NY+j])
            {
                phi[0][i*NY+j] = amp*(y1 + modulation*sin(a*y1)/a);
                phi[1][i*NY+j] = amp*modulation*a*sin(a*y1);
            }
        }
}


void init_laminar_flow(double amp, double xmodulation, double ymodulation, double xperiod, double yperiod, double yshift, double density_mod, double *phi[NFIELDS], short int xy_in[NX*NY])
/* initialise field with a laminar flow in x direction */
/* phi[0] is stream function, phi[1] is vorticity */
/* amp is global amplitude */
{
    int i, j;
    double xy[2], y1, a, b, f, cplus, cminus;    
    
    a = xperiod*PI/YMAX;
    b = yperiod*PI/XMAX;
    
    switch (RDE_EQUATION) {
        case (E_EULER_INCOMP): 
        {
            for (i=0; i<NX; i++)
                for (j=0; j<NY; j++)
                {
                    ij_to_xy(i, j, xy);
                    xy_in[i*NY+j] = xy_in_billiard(xy[0],xy[1]);
                    y1 = xy[1] + yshift;

                    if (xy_in[i*NY+j])
                    {
                        phi[0][i*NY+j] = amp*(y1 + xmodulation*sin(a*y1)/a);
                        phi[1][i*NY+j] = amp*xmodulation*a*sin(a*y1);
                    }
                    else
                    {
                        phi[0][i*NY+j] = 0.0;
                        phi[1][i*NY+j] = 0.0;
                    }
                }
            break;
        }
        case (E_EULER_COMP):
        {
            for (i=0; i<NX; i++)
                for (j=0; j<NY; j++)
                {
                    ij_to_xy(i, j, xy);
                    xy_in[i*NY+j] = xy_in_billiard(xy[0],xy[1]);
                    y1 = xy[1] + yshift;

                    if (xy_in[i*NY+j])
                    {
                        phi[0][i*NY+j] = 1.0 + density_mod*cos(a*y1);
                        phi[1][i*NY+j] = amp*(1.0 + xmodulation*cos(a*y1));
                        phi[2][i*NY+j] = ymodulation*sin(b*xy[0]);
                    }
                    else
                    {
                        phi[0][i*NY+j] = 1.0;
                        phi[1][i*NY+j] = 0.0;
                        phi[2][i*NY+j] = 0.0;
                    }
                }
            break;
        }
    }
}

void initialize_bcfield(double bc_field[NX*NY])
/* apply smooth modulation to adapt initial state to obstacles */ 
{
    int i, j;
    double xy[2], r, f;
    
    switch (OBSTACLE_GEOMETRY) {
        case (D_SINAI): 
        {
            for (i=0; i<NX; i++)
                for (j=0; j<NY; j++)
                {
                    ij_to_xy(i, j, xy);
                    r = module2(xy[0], xy[1]);
                    f = 0.5*(1.0 + tanh(BC_STIFFNESS*(r - 1.0*LAMBDA))); 
                    bc_field[i*NY+j] = f;
                }
            break;
        }
    }
}


void adapt_state_to_bc(double *phi[NFIELDS], double bc_field[NX*NY], short int xy_in[NX*NY])
/* apply smooth modulation to adapt initial state to obstacles */ 
{
    int i, j, field;
    double xy[2], r, f;
    
    #pragma omp parallel for private(i,j)
    for (i=0; i<NX; i++)
        for (j=0; j<NY; j++) if (xy_in[i*NY+j])
            {
                for (field = 1; field < NFIELDS; field++) 
                    phi[field][i*NY+j] *= bc_field[i*NY+j];
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
    double deltaphi, maxgradient = 1.0e10;
    double dx = (XMAX-XMIN)/((double)NX);
    
    dx = (XMAX-XMIN)/((double)NX);
    
    #pragma omp parallel for private(i,j,iplus,iminus,jplus,jminus,deltaphi)
    for (i=1; i<NX-1; i++)
        for (j=1; j<NY-1; j++)
        {
            iplus = i+1;
            iminus = i-1;
            jplus = j+1;
            jminus = j-1;
            
            deltaphi = phi[iplus*NY+j] - phi[iminus*NY+j];
            if (vabs(deltaphi) < maxgradient) gradient[i*NY+j] = (deltaphi)/dx;
            else gradient[i*NY+j] = 0.0;
            
            deltaphi = phi[i*NY+jplus] - phi[i*NY+jminus];
            if (vabs(deltaphi) < maxgradient) gradient[NX*NY+i*NY+j] = (deltaphi)/dx;
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

void compute_gradient_euler(double phi[NX*NY], double gradient[2*NX*NY], double yshift)
/* compute the gradient of the field */
{
    int i, j, iplus, iminus, jplus, jminus, padding = 0; 
    double deltaphi, maxgradient = 1.0e10;
    double dx = (XMAX-XMIN)/((double)NX);
    
    dx = (XMAX-XMIN)/((double)NX);
    
//     #pragma omp parallel for private(i)
//     for (i=0; i<2*NX*NY; i++) gradient[i] = 0.0;
    
    #pragma omp parallel for private(i,j,iplus,iminus,jplus,jminus)
    for (i=1; i<NX-1; i++)
        for (j=1; j<NY-1; j++)
        {
            iplus = i+1;
            iminus = i-1;
            jplus = j+1;
            jminus = j-1;
            
            gradient[i*NY+j] = 0.5*(phi[iplus*NY+j] - phi[iminus*NY+j]);
            gradient[NX*NY+i*NY+j] = 0.5*(phi[i*NY+jplus] - phi[i*NY+jminus]);
        }
        
    /* boundaries */
    for (i=1; i<NX-1; i++)
    {
        iplus = i+1;    
        iminus = i-1;   
        
        j = 0;
        jplus = 1;
        jminus = NY-1;
            
        gradient[i*NY+j] = 0.5*(phi[iplus*NY+j] - phi[iminus*NY+j]);
        gradient[NX*NY+i*NY+j] = 0.5*(phi[i*NY+jplus] - phi[i*NY+jminus] - yshift);
        
//         if (i == 1) printf("psi+ = %.5lg, psi- = %.5lg, gradient = %.5lg\n", phi[i*NY+jplus], phi[i*NY+jminus], gradient[NX*NY+i*NY+j]);
        
        j = NY-1;
        jplus = 0;
        jminus = NY-2;
        
        gradient[i*NY+j] = 0.5*(phi[iplus*NY+j] - phi[iminus*NY+j]);
        gradient[NX*NY+i*NY+j] = 0.5*(phi[i*NY+jplus] - phi[i*NY+jminus] + yshift);
    }
    
    for (j=1; j<NY-1; j++)
    {
        jplus = j+1;
        jminus = j-1;

        i = 0;
        iplus = 1; 
        iminus = NX-1; 

        gradient[i*NY+j] = 0.5*(phi[iplus*NY+j] - phi[iminus*NY+j]);
        gradient[NX*NY+i*NY+j] = 0.5*(phi[i*NY+jplus] - phi[i*NY+jminus]);
        
//         printf("j = %i, psi+ = %.5lg, psi- = %.5lg, gradient = %.5lg\n", j, phi[i*NY+jplus], phi[i*NY+jminus], gradient[NX*NY+i*NY+j]);

            
        i = NX-1;
        iplus = 0;
        iminus = NX-2;
        
        gradient[i*NY+j] = 0.5*(phi[iplus*NY+j] - phi[iminus*NY+j]);
        gradient[NX*NY+i*NY+j] = 0.5*(phi[i*NY+jplus] - phi[i*NY+jminus]);
    }
    
    /* corners */
    i = 0;  iplus = 1;  iminus = NX-1;
    
    j = 0;  jplus = 1;  jminus = NY-1;
    gradient[i*NY+j] = 0.5*(phi[iplus*NY+j] - phi[iminus*NY+j]);
    gradient[NX*NY+i*NY+j] = 0.5*(phi[i*NY+jplus] - phi[i*NY+jminus] + yshift);
    
    j = NY-1;  jplus = 0;  jminus = NY-2;
    gradient[i*NY+j] = 0.5*(phi[iplus*NY+j] - phi[iminus*NY+j]);
    gradient[NX*NY+i*NY+j] = 0.5*(phi[i*NY+jplus] - phi[i*NY+jminus] - yshift);
    
    i = NX-1;  iplus = 0;  iminus = NX-2;
    
    j = 0;  jplus = 1;  jminus = NY-1;
    gradient[i*NY+j] = 0.5*(phi[iplus*NY+j] - phi[iminus*NY+j]);
    gradient[NX*NY+i*NY+j] = 0.5*(phi[i*NY+jplus] - phi[i*NY+jminus] + yshift);
    
    j = NY-1;  jplus = 0;  jminus = NY-2;
    gradient[i*NY+j] = 0.5*(phi[iplus*NY+j] - phi[iminus*NY+j]);
    gradient[NX*NY+i*NY+j] = 0.5*(phi[i*NY+jplus] - phi[i*NY+jminus] - yshift);
}


void compute_gradient_euler_test(double phi[NX*NY], double gradient[2*NX*NY], short int xy_in[NX*NY])
/* compute the gradient of the field */
{
    int i, j, iplus, iminus, jplus, jminus, padding = 0; 
    double deltaphi, maxgradient = 1.0e10;
    double dx = (XMAX-XMIN)/((double)NX);
    
    dx = (XMAX-XMIN)/((double)NX);
    
    #pragma omp parallel for private(i)
    for (i=0; i<2*NX*NY; i++) gradient[i] = 0.0;
    
    #pragma omp parallel for private(i,j,iplus,iminus,jplus,jminus)
    for (i=1; i<NX-1; i++)
    {
        for (j=1; j<NY-1; j++) 
        {
            iplus = i+1;
            iminus = i-1;
            jplus = j+1;
            jminus = j-1;
            
            if ((xy_in[iplus*NY+j])&&(xy_in[iminus*NY+j]))
                gradient[i*NY+j] = 0.5*(phi[iplus*NY+j] - phi[iminus*NY+j]);
            if ((xy_in[i*NY+jplus])&&(xy_in[i*NY+jminus]))
                gradient[NX*NY+i*NY+j] = 0.5*(phi[i*NY+jplus] - phi[i*NY+jminus]);
        }
    }
        
    /* boundaries */
    for (i=1; i<NX-1; i++)
    {
        iplus = i+1;    
        iminus = i-1;   
        
        j = 0;
        jplus = 1;
        jminus = NY-1;
            
        gradient[i*NY+j] = 0.5*(phi[iplus*NY+j] - phi[iminus*NY+j]);
        gradient[NX*NY+i*NY+j] = 0.5*(phi[i*NY+jplus] - phi[i*NY+jminus]);
        
//         if (i == 1) printf("psi+ = %.5lg, psi- = %.5lg, gradient = %.5lg\n", phi[i*NY+jplus], phi[i*NY+jminus], gradient[NX*NY+i*NY+j]);
        
        j = NY-1;
        jplus = 0;
        jminus = NY-2;
        
        gradient[i*NY+j] = 0.5*(phi[iplus*NY+j] - phi[iminus*NY+j]);
        gradient[NX*NY+i*NY+j] = 0.5*(phi[i*NY+jplus] - phi[i*NY+jminus]);
    }
    
    for (j=1; j<NY-1; j++)
    {
        jplus = j+1;
        jminus = j-1;

        i = 0;
        iplus = 1; 
        iminus = NX-1; 

        gradient[i*NY+j] = 0.5*(phi[iplus*NY+j] - phi[iminus*NY+j]);
        gradient[NX*NY+i*NY+j] = 0.5*(phi[i*NY+jplus] - phi[i*NY+jminus]);
        
//         printf("j = %i, psi+ = %.5lg, psi- = %.5lg, gradient = %.5lg\n", j, phi[i*NY+jplus], phi[i*NY+jminus], gradient[NX*NY+i*NY+j]);

            
        i = NX-1;
        iplus = 0;
        iminus = NX-2;
        
        gradient[i*NY+j] = 0.5*(phi[iplus*NY+j] - phi[iminus*NY+j]);
        gradient[NX*NY+i*NY+j] = 0.5*(phi[i*NY+jplus] - phi[i*NY+jminus]);
    }
    
    /* corners */
    i = 0;  iplus = 1;  iminus = NX-1;
    
    j = 0;  jplus = 1;  jminus = NY-1;
    gradient[i*NY+j] = 0.5*(phi[iplus*NY+j] - phi[iminus*NY+j]);
    gradient[NX*NY+i*NY+j] = 0.5*(phi[i*NY+jplus] - phi[i*NY+jminus]);
    
    j = NY-1;  jplus = 0;  jminus = NY-2;
    gradient[i*NY+j] = 0.5*(phi[iplus*NY+j] - phi[iminus*NY+j]);
    gradient[NX*NY+i*NY+j] = 0.5*(phi[i*NY+jplus] - phi[i*NY+jminus]);
    
    i = NX-1;  iplus = 0;  iminus = NX-2;
    
    j = 0;  jplus = 1;  jminus = NY-1;
    gradient[i*NY+j] = 0.5*(phi[iplus*NY+j] - phi[iminus*NY+j]);
    gradient[NX*NY+i*NY+j] = 0.5*(phi[i*NY+jplus] - phi[i*NY+jminus]);
    
    j = NY-1;  jplus = 0;  jminus = NY-2;
    gradient[i*NY+j] = 0.5*(phi[iplus*NY+j] - phi[iminus*NY+j]);
    gradient[NX*NY+i*NY+j] = 0.5*(phi[i*NY+jplus] - phi[i*NY+jminus]);
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
            jplus = j+1;  if (jplus == NY) jplus = 0;
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
            jplus = j+1;  if (jplus == NY) jplus = 0;
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
            jplus = j+1;  if (jplus == NY) jplus = 0;
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

void compute_field_log(double *phi[NFIELDS], t_rde rde[NX*NY])
/* compute the log of a field */
{
    int i, j; 
    double value;
    
    #pragma omp parallel for private(i,j,value)
    for (i=0; i<NX; i++)
        for (j=0; j<NY; j++)
        {
            value = vabs(phi[1][i*NY+j]);
            if (value < LOG_MIN) value = LOG_MIN;
            rde[i*NY+j].log_vorticity = LOG_SHIFT + LOG_SCALE*log(value);
        }
}

void compute_pressure(double *phi[NFIELDS], t_rde rde[NX*NY])
/* compute the Laplacian of the pressure */
{
    int i, j, iplus, iminus, jplus, jminus; 
    double value, dx, psixx, psiyy, psixy;
    static double dx4;
    static int first = 1;
    
    if (first)
    {
        dx = (XMAX-XMIN)/((double)NX);
        dx4 = dx*dx;
        dx4 = dx4*dx4;
        first = 0;
    }
    
    #pragma omp parallel for private(i,j,iplus,iminus,jplus,jminus,psixx,psiyy,psixy)
    for (i=0; i<NX; i++)
        for (j=0; j<NY; j++)
        {
            iplus = i+1;  if (iplus == NX) iplus = 0;
            iminus = i-1; if (iminus == -1) iminus = NX-1;
            jplus = j+1;  if (jplus == NY) jplus = 0;
            jminus = j-1; if (jminus == -1) jminus = NY-1;
            
            psixx = phi[0][iplus*NY+j] - 2.0*phi[0][i*NY+j] + phi[0][iminus*NY+j];
            psiyy = phi[0][i*NY+jplus] - 2.0*phi[0][i*NY+j] + phi[0][i*NY+jminus];
            psixy = phi[0][iplus*NY+jplus] - phi[0][iplus*NY+jminus] - phi[0][iminus*NY+jplus] + phi[0][iminus*NY+jminus];
            
            rde[i*NY+j].Lpressure = (psixx*psiyy - psixy*psixy)/dx4;
        }
}


void compute_pressure_laplacian(double *phi[NFIELDS], double *l_pressure)
/* compute the Laplacian of the pressure */
{
    int i, j, iplus, iminus, jplus, jminus; 
    double value, dx, psixx, psiyy, psixy, dx4;
    static double scaling;
    static int first = 1;
    
    if (first)
    {
        dx = (XMAX-XMIN)/((double)NX);
        dx4 = dx*dx;
        scaling = 1.0/(2.0*dx4*dx4);
        first = 0;
    }
    
    #pragma omp parallel for private(i,j,iplus,iminus,jplus,jminus,psixx,psiyy,psixy)
    for (i=0; i<NX; i++)
        for (j=0; j<NY; j++)
        {
            iplus = i+1;  if (iplus == NX) iplus = 0;
            iminus = i-1; if (iminus == -1) iminus = NX-1;
            jplus = j+1;  if (jplus == NY) jplus = 0;
            jminus = j-1; if (jminus == -1) jminus = NY-1;
            
            psixx = phi[0][iplus*NY+j] - 2.0*phi[0][i*NY+j] + phi[0][iminus*NY+j];
            psiyy = phi[0][i*NY+jplus] - 2.0*phi[0][i*NY+j] + phi[0][i*NY+jminus];
            psixy = phi[0][iplus*NY+jplus] - phi[0][iplus*NY+jminus] - phi[0][iminus*NY+jplus] + phi[0][iminus*NY+jminus];
            
            l_pressure[i*NY+j] = (psixx*psiyy - psixy*psixy)*scaling;
        }
}


void compute_speed(double *phi[NFIELDS], t_rde rde[NX*NY])
/* compute the log of a field */
{
    int i, j; 
    double value;
    
    #pragma omp parallel for private(i,j,value)
    for (i=0; i<NX; i++)
        for (j=0; j<NY; j++)
        {
            value = module2(phi[1][i*NY+j], phi[2][i*NY+j]);
            rde[i*NY+j].field_norm = value;
        }
}

void compute_vorticity(t_rde rde[NX*NY])
/* compute the log of a field */
{
    int i, j; 
    double value;
    
    #pragma omp parallel for private(i,j,value)
    for (i=0; i<NX; i++)
        for (j=0; j<NY; j++)
        {
            rde[i*NY+j].curl = rde[i*NY+j].dxv - rde[i*NY+j].dyu + VORTICITY_SHIFT;
        }
}

void compute_velocity_gradients(double *phi[NFIELDS], t_rde rde[NX*NY])
/* compute the gradients of the velocity field (for Euler equation) */
{
    int i, j, k, iplus, iminus, jplus, jminus, padding = 0; 
    double deltaphi, maxgradient = 1.0e10;
//     double dx = (XMAX-XMIN)/((double)NX);
    
//     dx = (XMAX-XMIN)/((double)NX);
        
    #pragma omp parallel for private(i,j,iplus,iminus,jplus,jminus)
    for (i=1; i<NX-1; i++)
        for (j=1; j<NY-1; j++) 
        {
            iplus = i+1;
            iminus = i-1;
            jplus = j+1;
            jminus = j-1;
            
            rde[i*NY+j].dxu = 0.5*(phi[1][iplus*NY+j] - phi[1][iminus*NY+j]);
            rde[i*NY+j].dyu = 0.5*(phi[1][i*NY+jplus] - phi[1][i*NY+jminus]);
           
            rde[i*NY+j].dxv = 0.5*(phi[2][iplus*NY+j] - phi[2][iminus*NY+j]);
            rde[i*NY+j].dyv = 0.5*(phi[2][i*NY+jplus] - phi[2][i*NY+jminus]);
        }
        
    /* boundaries */
    for (i=1; i<NX-1; i++)
    {
        iplus = i+1;    
        iminus = i-1;   
        
        j = 0;
        jplus = 1;
        jminus = NY-1;
            
        rde[i*NY+j].dxu = 0.5*(phi[1][iplus*NY+j] - phi[1][iminus*NY+j]);
        rde[i*NY+j].dyu = 0.5*(phi[1][i*NY+jplus] - phi[1][i*NY+jminus]);
        
        rde[i*NY+j].dxv = 0.5*(phi[2][iplus*NY+j] - phi[2][iminus*NY+j]);
        rde[i*NY+j].dyv = 0.5*(phi[2][i*NY+jplus] - phi[2][i*NY+jminus]);
        
        j = NY-1;
        jplus = 0;
        jminus = NY-2;
        
        rde[i*NY+j].dxu = 0.5*(phi[1][iplus*NY+j] - phi[1][iminus*NY+j]);
        rde[i*NY+j].dyu = 0.5*(phi[1][i*NY+jplus] - phi[1][i*NY+jminus]);
        
        rde[i*NY+j].dxv = 0.5*(phi[2][iplus*NY+j] - phi[2][iminus*NY+j]);
        rde[i*NY+j].dyv = 0.5*(phi[2][i*NY+jplus] - phi[2][i*NY+jminus]);
    }
    
    for (j=1; j<NY-1; j++)
    {
        jplus = j+1;
        jminus = j-1;

        i = 0;
        iplus = 1; 
        iminus = NX-1; 

        rde[i*NY+j].dxu = 0.5*(phi[1][iplus*NY+j] - phi[1][iminus*NY+j]);
        rde[i*NY+j].dyu = 0.5*(phi[1][i*NY+jplus] - phi[1][i*NY+jminus]);
        
        rde[i*NY+j].dxv = 0.5*(phi[2][iplus*NY+j] - phi[2][iminus*NY+j]);
        rde[i*NY+j].dyv = 0.5*(phi[2][i*NY+jplus] - phi[2][i*NY+jminus]);
                    
        i = NX-1;
        iplus = 0;
        iminus = NX-2;
        
        rde[i*NY+j].dxu = 0.5*(phi[1][iplus*NY+j] - phi[1][iminus*NY+j]);
        rde[i*NY+j].dyu = 0.5*(phi[1][i*NY+jplus] - phi[1][i*NY+jminus]);
        
        rde[i*NY+j].dxv = 0.5*(phi[2][iplus*NY+j] - phi[2][iminus*NY+j]);
        rde[i*NY+j].dyv = 0.5*(phi[2][i*NY+jplus] - phi[2][i*NY+jminus]);
    }
    
    /* corners */
    i = 0;  iplus = 1;  iminus = NX-1;
    
    j = 0;  jplus = 1;  jminus = NY-1;
    rde[i*NY+j].dxu = 0.5*(phi[1][iplus*NY+j] - phi[1][iminus*NY+j]);
    rde[i*NY+j].dyu = 0.5*(phi[1][i*NY+jplus] - phi[1][i*NY+jminus]);
        
    rde[i*NY+j].dxv = 0.5*(phi[2][iplus*NY+j] - phi[2][iminus*NY+j]);
    rde[i*NY+j].dyv = 0.5*(phi[2][i*NY+jplus] - phi[2][i*NY+jminus]);
    
    j = NY-1;  jplus = 0;  jminus = NY-2;
    rde[i*NY+j].dxu = 0.5*(phi[1][iplus*NY+j] - phi[1][iminus*NY+j]);
    rde[i*NY+j].dyu = 0.5*(phi[1][i*NY+jplus] - phi[1][i*NY+jminus]);
        
    rde[i*NY+j].dxv = 0.5*(phi[2][iplus*NY+j] - phi[2][iminus*NY+j]);
    rde[i*NY+j].dyv = 0.5*(phi[2][i*NY+jplus] - phi[2][i*NY+jminus]);
    
    i = NX-1;  iplus = 0;  iminus = NX-2;
    
    j = 0;  jplus = 1;  jminus = NY-1;
    rde[i*NY+j].dxu = 0.5*(phi[1][iplus*NY+j] - phi[1][iminus*NY+j]);
    rde[i*NY+j].dyu = 0.5*(phi[1][i*NY+jplus] - phi[1][i*NY+jminus]);
        
    rde[i*NY+j].dxv = 0.5*(phi[2][iplus*NY+j] - phi[2][iminus*NY+j]);
    rde[i*NY+j].dyv = 0.5*(phi[2][i*NY+jplus] - phi[2][i*NY+jminus]);
    
    j = NY-1;  jplus = 0;  jminus = NY-2;
    rde[i*NY+j].dxu = 0.5*(phi[1][iplus*NY+j] - phi[1][iminus*NY+j]);
    rde[i*NY+j].dyu = 0.5*(phi[1][i*NY+jplus] - phi[1][i*NY+jminus]);
        
    rde[i*NY+j].dxv = 0.5*(phi[2][iplus*NY+j] - phi[2][iminus*NY+j]);
    rde[i*NY+j].dyv = 0.5*(phi[2][i*NY+jplus] - phi[2][i*NY+jminus]);
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
    
    #pragma omp parallel for private(i,j)
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
                jplus = j+1;  if (jplus == NY) jplus = 0;
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
            break;
        }
        case (BC_DIRICHLET):
        {
            /* left and right side */
            for (j = 1; j < NY-1; j++) 
            {
                phi_out[j] = phi_in[j-1] + phi_in[j+1] + phi_in[NY+j] - 3.0*phi_in[j];
                phi_out[(NX-1)*NY+j] = phi_in[(NX-1)*NY+j-1] + phi_in[(NX-1)*NY+j+1] + phi_in[(NX-2)*NY+j] - 3.0*phi_in[(NX-1)*NY+j];
            }
            /* top and bottom side */
            for (i = 1; i < NX-1; i++) 
            {
                phi_out[i*NY] = phi_in[(i-1)*NY] + phi_in[(i+1)*NY] + phi_in[i*NY+1] - 3.0*phi_in[i*NY];
                phi_out[i*NY+NY-1] = phi_in[(i-1)*NY+NY-1] + phi_in[(i+1)*NY+NY-1] + phi_in[i*NY+NY-2] - 3.0*phi_in[i*NY+NY-1];
            }
            /* corners */
            phi_out[0] = phi_in[1] + phi_in[NY] - 2.0*phi_in[0];
            phi_out[NY-1] = phi_in[NY-2] + phi_in[NY+NY-1] - 2.0*phi_in[NY-1];
            phi_out[(NX-1)*NY] = phi_in[(NX-2)*NY] + phi_in[(NX-1)*NY+1] - 2.0*phi_in[(NX-1)*NY];
            phi_out[(NX-1)*NY+NY-1] = phi_in[(NX-2)*NY+NY-1] + phi_in[(NX-1)*NY-2] - 2.0*phi_in[(NX-1)*NY+NY-1];
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
                if (((ADD_POTENTIAL)||(ADD_MAGNETIC_FIELD))&&(ADD_POTENTIAL_TO_Z))
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
        case (Z_EULER_VORTICITY): 
        {
            if (value < 0.0) value = -value;
            color_scheme_asym_palette(COLOR_SCHEME, palette, VSCALE_AMPLITUDE*value, 1.0, 0, rgb);
            break;
        }
        case (Z_EULER_LOG_VORTICITY): 
        {
//             if (value < 0.0) value = -value;
//             if (value < 1.0e-10) value = 1.0e-10;
//             color_scheme_palette(COLOR_SCHEME, palette, LOG_SCALE*value + LOG_SHIFT, 1.0, 0, rgb);
            color_scheme_palette(COLOR_SCHEME, palette, LOG_SCALE*value, 1.0, 0, rgb);
            break;
        }
        case (Z_EULER_VORTICITY_ASYM): 
        {
            color_scheme_palette(COLOR_SCHEME, palette, VSCALE_AMPLITUDE*value, 1.0, 0, rgb);
            break;
        }
        case (Z_EULER_LPRESSURE): 
        {
            color_scheme_palette(COLOR_SCHEME, palette, VSCALE_AMPLITUDE*value, 1.0, 0, rgb);
            break;
        }
        case (Z_EULER_PRESSURE): 
        {
            if (value + PRESSURE_SHIFT < 1.0e-10) value = 1.0e-10 - PRESSURE_SHIFT;
            color_scheme_palette(COLOR_SCHEME, palette, VSCALE_PRESSURE*log(value + PRESSURE_SHIFT) + PRESSURE_LOG_SHIFT, 1.0, 0, rgb);
//             color_scheme_palette(COLOR_SCHEME, palette, VSCALE_PRESSURE*(value - AVERAGE_PRESSURE), 1.0, 0, rgb);
            break;
        }
        case (Z_EULER_DENSITY): 
        {
            color_scheme_palette(COLOR_SCHEME, palette, VSCALE_DENSITY*(value-1.0), 1.0, 0, rgb);
            break;
        }
        case (Z_EULER_SPEED): 
        {
            if (ASYM_SPEED_COLOR) 
                color_scheme_asym_palette(COLOR_SCHEME, palette, VSCALE_SPEED*value, 1.0, 0, rgb);
            else 
                color_scheme_palette(COLOR_SCHEME, palette, VSCALE_SPEED*(value-VMEAN_SPEED), 1.0, 0, rgb);
            break;
        }
        case (Z_EULERC_VORTICITY): 
        {
            color_scheme_palette(COLOR_SCHEME, palette, VSCALE_VORTICITY*(value-VORTICITY_SHIFT), 1.0, 0, rgb);
            break;
        }
    }
}

double adjust_field(double z, double pot)
/* add potential in case of option ADD_POTENTIAL_TO_Z */
{
    if (((ADD_POTENTIAL)||(ADD_MAGNETIC_FIELD))&&(ADD_POTENTIAL_TO_Z)) return (z + ADD_POT_CONSTANT*pot);
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
    
    if (((ADD_POTENTIAL)||(ADD_MAGNETIC_FIELD))&&(ADD_POTENTIAL_TO_Z))
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
        case (E_EULER_INCOMP):
        {
            if ((zplot == Z_EULER_LOG_VORTICITY)||(cplot == Z_EULER_LOG_VORTICITY))
                compute_field_log(phi, rde);
            if ((zplot == Z_EULER_LPRESSURE)||(cplot == Z_EULER_LPRESSURE))
                compute_pressure(phi, rde);
            break;
        }
        case (E_EULER_COMP):
        {
            if ((zplot == Z_EULER_SPEED)||(cplot == Z_EULER_SPEED))
                compute_speed(phi, rde);
            if ((zplot == Z_EULERC_VORTICITY)||(cplot == Z_EULERC_VORTICITY))
                compute_vorticity(rde);
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
        case (Z_EULER_VORTICITY):
        {
            #pragma omp parallel for private(i,j)
            for (i=0; i<NX; i++) for (j=0; j<NY; j++) rde[i*NY+j].p_zfield[movie] = &phi[1][i*NY+j];
            break;
        }
        case (Z_EULER_LOG_VORTICITY):
        {
            #pragma omp parallel for private(i,j)
            for (i=0; i<NX; i++) for (j=0; j<NY; j++) rde[i*NY+j].p_zfield[movie] = &rde[i*NY+j].log_vorticity;
            break;
        }
        case (Z_EULER_VORTICITY_ASYM):
        {
            #pragma omp parallel for private(i,j)
            for (i=0; i<NX; i++) for (j=0; j<NY; j++) rde[i*NY+j].p_zfield[movie] = &phi[1][i*NY+j];
            break;
        }
        case (Z_EULER_LPRESSURE):
        {
            #pragma omp parallel for private(i,j)
            for (i=0; i<NX; i++) for (j=0; j<NY; j++) rde[i*NY+j].p_zfield[movie] = &rde[i*NY+j].Lpressure;
            break;
        }
        case (Z_EULER_PRESSURE):
        {
            #pragma omp parallel for private(i,j)
            for (i=0; i<NX; i++) for (j=0; j<NY; j++) rde[i*NY+j].p_zfield[movie] = &phi[2][i*NY+j];
            break;
        }
        case (Z_EULER_DENSITY):
        {
            #pragma omp parallel for private(i,j)
            for (i=0; i<NX; i++) for (j=0; j<NY; j++) rde[i*NY+j].p_zfield[movie] = &phi[0][i*NY+j];
            break;
        }
        case (Z_EULER_SPEED):
        {
            #pragma omp parallel for private(i,j)
            for (i=0; i<NX; i++) for (j=0; j<NY; j++) rde[i*NY+j].p_zfield[movie] = &rde[i*NY+j].field_norm;
            break;
        }
        case (Z_EULERC_VORTICITY):
        {
            #pragma omp parallel for private(i,j)
            for (i=0; i<NX; i++) for (j=0; j<NY; j++) rde[i*NY+j].p_zfield[movie] = &rde[i*NY+j].curl;
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
        case (Z_EULER_VORTICITY):
        {
            #pragma omp parallel for private(i,j)
            for (i=0; i<NX; i++) for (j=0; j<NY; j++) rde[i*NY+j].p_cfield[movie] = &phi[1][i*NY+j];
            break;
        }
        case (Z_EULER_LOG_VORTICITY):
        {
            #pragma omp parallel for private(i,j)
            for (i=0; i<NX; i++) for (j=0; j<NY; j++) rde[i*NY+j].p_cfield[movie] = &rde[i*NY+j].log_vorticity;
            break;
        }
        case (Z_EULER_VORTICITY_ASYM):
        {
            #pragma omp parallel for private(i,j)
            for (i=0; i<NX; i++) for (j=0; j<NY; j++) rde[i*NY+j].p_cfield[movie] = &phi[1][i*NY+j];
            break;
        }
        case (Z_EULER_LPRESSURE):
        {
            #pragma omp parallel for private(i,j)
            for (i=0; i<NX; i++) for (j=0; j<NY; j++) rde[i*NY+j].p_cfield[movie] = &rde[i*NY+j].Lpressure;
            break;
        }
        case (Z_EULER_PRESSURE):
        {
            #pragma omp parallel for private(i,j)
            for (i=0; i<NX; i++) for (j=0; j<NY; j++) rde[i*NY+j].p_cfield[movie] = &phi[2][i*NY+j];
            break;
        }
        case (Z_EULER_DENSITY):
        {
            #pragma omp parallel for private(i,j)
            for (i=0; i<NX; i++) for (j=0; j<NY; j++) rde[i*NY+j].p_cfield[movie] = &phi[0][i*NY+j];
            break;
        }
        case (Z_EULER_SPEED):
        {
            #pragma omp parallel for private(i,j)
            for (i=0; i<NX; i++) for (j=0; j<NY; j++) rde[i*NY+j].p_cfield[movie] = &rde[i*NY+j].field_norm;
            break;
        }
        case (Z_EULERC_VORTICITY):
        {
            #pragma omp parallel for private(i,j)
            for (i=0; i<NX; i++) for (j=0; j<NY; j++) rde[i*NY+j].p_cfield[movie] = &rde[i*NY+j].curl;
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


void draw_wave_3d_ij_rde(int i, int j, int movie, double *phi[NFIELDS], short int xy_in[NX*NY], t_rde rde[NX*NY], double potential[NX*NY], int zplot, int cplot, int palette, int fade, double fade_value)
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

void draw_wave_3d_ij_rde_periodic(int shiftx, int shifty, int i, int j, int movie, double *phi[NFIELDS], 
                         short int xy_in[NX*NY], t_rde rde[NX*NY], double potential[NX*NY], int zplot, int cplot, int palette, int fade, double fade_value)
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
                draw_vertex_xyz_shift(xy_mid, z_mid, shiftx, shifty);
                draw_vertex_xyz_shift(xy_nw, z_nw, shiftx, shifty);
                draw_vertex_xyz_shift(xy_sw, z_sw, shiftx, shifty);
                    
                glColor3f(rgb_s[0], rgb_s[1], rgb_s[2]);
                draw_vertex_xyz_shift(xy_se, z_se, shiftx, shifty);
                    
                glColor3f(rgb_e[0], rgb_e[1], rgb_e[2]);
                draw_vertex_xyz_shift(xy_ne, z_ne, shiftx, shifty);
                    
                glColor3f(rgb_n[0], rgb_n[1], rgb_n[2]);
                draw_vertex_xyz_shift(xy_nw, z_nw, shiftx, shifty);
                glEnd ();
            }
        }
        else
        {
            glColor3f(rde[i*NY+j].rgb[0], rde[i*NY+j].rgb[1], rde[i*NY+j].rgb[2]);
            
            glBegin(GL_TRIANGLE_FAN);
            ij_to_xy(i, j, xy);
            draw_vertex_xyz_shift(xy, adjust_field(*rde[i*NY+j].p_zfield[movie], potential[i*NY+j]), shiftx, shifty);
            ij_to_xy(i+1, j, xy);
            draw_vertex_xyz_shift(xy, adjust_field(*rde[(i+1)*NY+j].p_zfield[movie], potential[(i+1)*NY+j]), shiftx, shifty);
            ij_to_xy(i+1, j+1, xy);
            draw_vertex_xyz_shift(xy, adjust_field(*rde[(i+1)*NY+j+1].p_zfield[movie], potential[(i+1)*NY+j+1]), shiftx, shifty);
            ij_to_xy(i, j+1, xy);
            draw_vertex_xyz_shift(xy, adjust_field(*rde[i*NY+j+1].p_zfield[movie], potential[i*NY+j+1]), shiftx, shifty);
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
        for (i=BORDER_PADDING; i<NX-2-BORDER_PADDING; i++)
            for (j=BORDER_PADDING; j<NY-2-BORDER_PADDING; j++)
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


void draw_periodicised_wave_3d(int movie, double *phi[NFIELDS], short int xy_in[NX*NY], t_rde rde[NX*NY], double potential[NX*NY], 
                  int zplot, int cplot, int palette, int fade, double fade_value)
/* a variant where the wave is repeated periodically in x and y directions (beta) */
{
    int i, j, shiftx, shifty;
    double observer_angle;
    
    blank();
    if (DRAW_BILLIARD) draw_billiard_3d(fade, fade_value);
    
    if (!ROTATE_VIEW)
    {
        for (shiftx = -1; shiftx < 2; shiftx++)
            for (shifty = -2; shifty < 2; shifty++)
                for (i=BORDER_PADDING; i<NX-1-BORDER_PADDING; i++)
                    for (j=BORDER_PADDING; j<NY-1-BORDER_PADDING; j++)
                        draw_wave_3d_ij_rde_periodic(shiftx, shifty, i, j, movie, phi, xy_in, rde, potential, zplot, cplot, palette, fade, fade_value);
    }
}


void draw_wave_rde(int movie, double *phi[NFIELDS], short int xy_in[NX*NY], t_rde rde[NX*NY],
                   double potential[NX*NY], int zplot, int cplot, int palette, int fade, double fade_value, int refresh)
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
    
    if (PLOT_3D) 
    {
        if (DRAW_PERIODICISED) 
            draw_periodicised_wave_3d(movie, phi, xy_in, rde, potential, zplot, cplot, palette, fade, fade_value);
        else draw_wave_3d_rde(movie, phi, xy_in, rde, potential, zplot, cplot, palette, fade, fade_value);
    }
    else draw_wave_2d_rde(xy_in, rde);
}

void draw_tracers(double *phi[NFIELDS], double tracers[2*N_TRACERS*NSTEPS], int time, int fade, double fade_value)
/* draw trajectories of tracers */
{
    int tracer, t, t1, length = 50;
    double x1, y1, x2, y2, lum;
    
    glColor3f(1.0, 1.0, 1.0);
    glLineWidth(1);
    
    t1 = time - length;
    if (t1 < 1) t1 = 1;
    
    for (t = t1 + 1; t < time; t++) 
        for (tracer = 0; tracer < N_TRACERS; tracer++)
        {
            x1 = tracers[2*(t-1)*N_TRACERS + 2*tracer];
            y1 = tracers[2*(t-1)*N_TRACERS + 2*tracer + 1];
        
            x2 = tracers[2*t*N_TRACERS + 2*tracer];
            y2 = tracers[2*t*N_TRACERS + 2*tracer + 1];
            
            lum = 1.0 - 0.75*(double)(time - t)/(double)length;
            if (fade) lum *= fade_value;
            
            glColor3f(lum, lum, lum);
            if (module2(x2 - x1, y2 - y1) < 0.2) draw_line(x1, y1, x2, y2);
            
//             printf("time = %i, tracer = %i, coord = %i, x1 = %.2lg, y1 = %.2lg, x2 = %.2lg, y2 = %.2lg\n", t, tracer,2*t*N_TRACERS + 2*tracer, x1, y1, x2, y2);
        }
    
}
