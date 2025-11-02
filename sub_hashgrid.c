/* routines dealing with hashgrid for lennardjones.c */
/* Note: a possible improvement would be to use pointers instead of integers */
/* when referencing other hashgrid cells or particles */ 

int hashx_sphere[HASHY];    /* number of hash cells in x direction for sphere bc */
double dx_sphere[HASHY];    /* width of hash cells in x direction for sphere bc */

 double module2(double x, double y)   /* Euclidean norm */
 {
	double m;

	m = sqrt(x*x + y*y);
	return(m);
 }

int bc_grouped(int bc)
/* regroup boundary conditions by type: rectangular = 0, periodic = 1, Klein = 2, Boy = 3, L shape = 4, ... */
{
    switch (bc) {
        case (BC_SCREEN): return(0);
        case (BC_RECTANGLE): return(0);
        case (BC_CIRCLE): return(0);
        case (BC_PERIODIC): return(1);
        case (BC_PERIODIC_CIRCLE): return(1);
        case (BC_EHRENFEST): return(0);
        case (BC_PERIODIC_FUNNEL): return(1);
        case (BC_RECTANGLE_LID): return(0);
        case (BC_RECTANGLE_WALL): return(0);
        case (BC_PERIODIC_TRIANGLE): return(1);
        case (BC_KLEIN): return(2);
        case (BC_SCREEN_BINS): return(0);
        case (BC_BOY): return(3);
        case (BC_GENUS_TWO): return(4);
        case (BC_ABSORBING): return(0);
        case (BC_REFLECT_ABS): return(0);
        case (BC_REFLECT_ABS_BOTTOM): return(0);
        case (BC_REFLECT_ABS_RIGHT): return(0);
        case (BC_SPHERE): return(5);
        default: 
        {
            printf("Warning: Hashgrid will not be properly initialised, update bc_grouped()\n\n");
            return(-1);
        }
    }
}


int mhash(int i, int j)
{
    return(i*HASHY + j);
}


int hash_cell(double x, double y)
/* compute hash grid position of particle at (x,y) */
/* returns number of hash cell */
{
    static int first = 1;
    static double lx, ly, padding, dy_inverse;
    int i, j;
    
    if (first)
    {
        if (!NO_WRAP_BC) padding = 0.0;
        else padding = HASHGRID_PADDING;
        lx = BCXMAX - BCXMIN + 2.0*padding;
        ly = BCYMAX - (BCYMIN) + 2.0*padding;
        dy_inverse = (double)HASHY/(PI - 2.0*POLAR_PADDING);
        first = 0;
    }
    
    if (bc_grouped(BOUNDARY_COND) == 5)     /* spherical bc */
    {
        j = (int)((y-POLAR_PADDING)*dy_inverse);
        if (j<0) j = 0;
        else if (j>=HASHY) j = HASHY-1;
        
        i = (int)(x/dx_sphere[j]);
        if (i<0) i = 0;
        else if (i>=hashx_sphere[j]) i = hashx_sphere[j]-1;
        
    }
    else 
    {
        if (CENTER_VIEW_ON_OBSTACLE) x -= xshift;
    
        i = (int)((double)HASHX*(x - BCXMIN + padding)/lx);
        j = (int)((double)HASHY*(y - BCYMIN + padding)/ly);
    
        if (i<0) i = 0;
        else if (i>=HASHX) i = HASHX-1;
        if (j<0) j = 0;
        else if (j>=HASHY) j = HASHY-1;
    }
    
//     printf("Hash_cell (%.3lg, %.3lg): %i\n", x, y, mhash(i,j));
    return(mhash(i,j));
//     printf("Mapped (%.3f,%.3f) to (%i, %i)\n", x, y, ij[0], ij[1]);
}


void init_hashcell_coordinates(t_hashgrid hashgrid[HASHX*HASHY])
/* initialise coordinates of corners of hashcells, needed for option COLOR_BACKGROUND */
{
    static int first = 1;
    static double lx, ly, padding, dx, dy;
    int i, j, n;
    double x, y, dummy = -100.0;
    
    if (first)
    {
        if (!NO_WRAP_BC) padding = 0.0;
        else padding = HASHGRID_PADDING;
        lx = BCXMAX - BCXMIN + 2.0*padding;
        ly = BCYMAX - (BCYMIN) + 2.0*padding;
        dx = 1.0*lx/(double)HASHX;
        if (bc_grouped(BOUNDARY_COND) == 5)
        {
            dy = (PI - 2.0*POLAR_PADDING)/(double)(HASHY);
        }
        else dy = 1.0*ly/(double)HASHY;
        first = 0;
    }
    
    if (bc_grouped(BOUNDARY_COND) == 5)     /* spherical bc */
    {
        for (j=0; j<HASHY; j++)
        {
            for (i=0; i<hashx_sphere[j]; i++)
            {
                n = mhash(i, j);
                hashgrid[n].x1 = (double)i*dx_sphere[j];
                hashgrid[n].x2 = (double)(i+1)*dx_sphere[j];
                hashgrid[n].y1 = POLAR_PADDING + (double)j*dy;
                hashgrid[n].y2 = POLAR_PADDING + (double)(j+1)*dy;
                
                /* TODO: correct area */
                hashgrid[n].area = dx_sphere[j]*dy*sin(0.5*(hashgrid[n].y1 + hashgrid[n].y2));
                
                printf("Cell %i corners (%.3lg, %.3lg), (%.3lg, %.3lg)\n", n, hashgrid[n].x1, hashgrid[n].y1, hashgrid[n].x2, hashgrid[n].y2); 
            }
            for (i=hashx_sphere[j]; i<HASHX; i++)
            {
                n = mhash(i, j);
                hashgrid[n].x1 = dummy;
                hashgrid[n].x2 = dummy;
                hashgrid[n].y1 = dummy;
                hashgrid[n].y2 = dummy;
                hashgrid[n].area = 1.0;
            }
        }
        /* correction at poles */
        n = mhash(0,0);
        hashgrid[n].y1 = 0.0;
        hashgrid[n].area = DPI*(POLAR_PADDING + dy)*sin(0.5*hashgrid[n].y2);
        n = mhash(0, HASHY-1);
        hashgrid[n].y2 = PI;
        hashgrid[n].area = DPI*(POLAR_PADDING + dy)*sin(0.5*(PI + hashgrid[n].y1));
    }
    else for (i=0; i<HASHX; i++)
        for (j=0; j<HASHY; j++)
        {
            x = BCXMIN - padding + (double)i*dx;
            y = BCYMIN - padding + (double)j*dy;
            n = mhash(i, j);
            hashgrid[n].x1 = x;
            hashgrid[n].y1 = y;
            hashgrid[n].x2 = x + dx;
            hashgrid[n].y2 = y + dy;
            hashgrid[n].charge = 0.0;
        }
}


void init_hashgrid(t_hashgrid hashgrid[HASHX*HASHY])
/* initialise table of neighbouring cells for each hashgrid cell, depending on boundary condition */
{
    int i, j, k, p, q, m, i1, j1, m1, pmin, pmax;
    double dy, x1, x2, xx1, xx2, padding;
    short int sym;
    
    printf("Initializing hash grid\n");
    
    /* initialize hashx_sphere[] and dx_sphere[] in case of spherical bc */
    if (bc_grouped(BOUNDARY_COND) == 5) 
    {
        dy = PI/((double)HASHY);
        for (j=0; j<HASHY/2 + 1; j++)
        {
            hashx_sphere[j] = (int)((double)HASHX*(sin(dy*(double)j)));
            if (hashx_sphere[j] == 0) hashx_sphere[j] = 1;
            if (hashx_sphere[j] > HASHX-1) hashx_sphere[j] = HASHX-1;
            dx_sphere[j] = DPI/(double)hashx_sphere[j];
            
            printf("hashx_sphere[%i] = %i, dx_sphere[%i] = %.3lg\n", j, hashx_sphere[j], j, dx_sphere[j]);
            fprintf(lj_log, "hashx_sphere[%i] = %i, dx_sphere[%i] = %.3lg\n", j, hashx_sphere[j], j, dx_sphere[j]);
        }
        for (j=HASHY/2 + 1; j<HASHY; j++)
        {
            hashx_sphere[j] = hashx_sphere[HASHY-1-j];
            dx_sphere[j] = dx_sphere[HASHY-1-j];
            
            printf("hashx_sphere[%i] = %i, dx_sphere[%i] = %.3lg\n", j, hashx_sphere[j], j, dx_sphere[j]);
        }
    }
    
    if ((COLOR_BACKGROUND)||(bc_grouped(BOUNDARY_COND) == 5)) 
        init_hashcell_coordinates(hashgrid);
    
    /* bulk of the table */
    if (bc_grouped(BOUNDARY_COND) == 5)     /* spherical bc */
    {
        /* dummy values for safety */
        for (i=0; i<HASHX*HASHY; i++) hashgrid[i].nneighb = 0;
        
        for (j=1; j<HASHY-1; j++)
        {
            padding = 1.0*dx_sphere[j];
            for (i=1; i<hashx_sphere[j]-1; i++)
            {
                m = mhash(i, j);
//                 printf("m = %i\n", m); 
                hashgrid[m].nneighb = 3;
                hashgrid[m].neighbour[0] = mhash(i-1,j);
                hashgrid[m].neighbour[1] = mhash(i,j);
                hashgrid[m].neighbour[2] = mhash(i+1,j);
                
                /* neighbours in layer above and below */
                if (j==1) pmin = 1;
                else pmin = -1;
                if (j==HASHY-2) pmax = 0;
                else pmax = 2;
                for (p=pmin; p<pmax; p+=2)
                    for (i1=0; i1<hashx_sphere[j+p]; i1++)
                    {
                        m1 = mhash(i1, j+p);
                        x1 = hashgrid[m].x1 - padding;
                        x2 = hashgrid[m].x2 + padding;
                        xx1 = hashgrid[m1].x1;
                        xx2 = hashgrid[m1].x2;
                        if (((xx2 >= x1)&&(xx2 <= x2))||((xx1 >= x1)&&(xx1 <= x2)))
                        {
                            hashgrid[m].neighbour[hashgrid[m].nneighb] = m1;
                            hashgrid[m].nneighb++;
                        }
                    }
                    
                /* extra treatment for pole neighbours */
                if (j==1)
                {
                    hashgrid[m].neighbour[hashgrid[m].nneighb] = 0;
                    hashgrid[m].nneighb++;
                }
                if (j==HASHY-2)
                {
                    hashgrid[m].neighbour[hashgrid[m].nneighb] = HASHY-1;
                    hashgrid[m].nneighb++;
                }
                    
//                 printf("hashgrid[%i].nneighb = %i\n", m, hashgrid[m].nneighb);
                
                /* check there are not too many neighbours */
                if (hashgrid[m].nneighb >= HASHMAXNEIGH)
                {
                    printf("(i,j) = (%i,%i) Error: HASHMAXNEIGH should be at least %i\n", i, j, hashgrid[m].nneighb + 1);
                    exit(1);
                }
            }
        }
    }
    else for (i=0; i<HASHX-1; i++)
        for (j=0; j<HASHY-1; j++)
        {
            m = mhash(i, j);
            hashgrid[m].nneighb = 9;
            
            for (p=-1; p<2; p++)
                for (q=-1; q<2; q++)
                    hashgrid[m].neighbour[3*p+q+4] = mhash(i+p, j+q);
        }
        
    /* different boundary conditions */
    switch (bc_grouped(BOUNDARY_COND)) {
        case (0): /* rectangular b.c. */
        {
            /* left/right boundaries */
            for (j=1; j<HASHY-1; j++)
            {
                for (i=0; i<HASHX; i+=HASHX-1) hashgrid[mhash(i, j)].nneighb = 6;
                
                for (q=-1; q<2; q++)
                {
                    for (p=0; p<2; p++) hashgrid[mhash(0, j)].neighbour[2*(q+1)+p] = mhash(p, j+q);
                    for (p=-1; p<1; p++) hashgrid[mhash(HASHX-1, j)].neighbour[2*(q+1)+p+1] = mhash(HASHX-1+p, j+q);
                }
            }
            
            /* top/down boundaries */
            for (i=1; i<HASHX-1; i++)
            {
                for (j=0; j<HASHY; j+=HASHY-1) hashgrid[mhash(i, j)].nneighb = 6;
                
                for (p=-1; p<2; p++)
                {
                    for (q=0; q<2; q++) hashgrid[mhash(i, 0)].neighbour[2*(p+1)+q] = mhash(i+p, q);
                    for (q=-1; q<1; q++) hashgrid[mhash(i, HASHY-1)].neighbour[2*(p+1)+q+1] = mhash(i+p, HASHY-1+q);
                }
            }
            
            /* corners */
            for (i=0; i<HASHX; i+=HASHX-1) for (j=0; j<HASHY; j+=HASHY-1) hashgrid[mhash(i, j)].nneighb = 4;
            for (p=0; p<2; p++) for (q=0; q<2; q++) 
            {
                hashgrid[mhash(0,0)].neighbour[2*p+q] = mhash(p, q);
                hashgrid[mhash(HASHX-1,0)].neighbour[2*p+q] = mhash(HASHX-1-p, q);
                hashgrid[mhash(0,HASHY-1)].neighbour[2*p+q] = mhash(p, HASHY-1-q);
                hashgrid[mhash(HASHX-1,HASHY-1)].neighbour[2*p+q] = mhash(HASHX-1-p, HASHY-1-q);                
            }
            break;
        }
        
        case(1): /* periodic b.c. */
        {
            /* left/right boundaries */
            for (j=0; j<HASHY; j++) for (i=0; i<HASHX; i+=HASHX-1) 
            {
                hashgrid[mhash(i, j)].nneighb = 9;
                
                for (p=-1; p<2; p++) for (q=-1; q<2; q++)
                {
                    i1 = (i+p+HASHX)%HASHX;
                    j1 = (j+q+HASHY)%HASHY;
                    hashgrid[mhash(i,j)].neighbour[3*(p+1)+q+1] = mhash(i1, j1);
                }
            }
            
            /* top/bottom boundaries */
            for (i=1; i<HASHX-1; i++) for (j=0; j<HASHY; j+=HASHY-1) 
            {
                hashgrid[mhash(i, j)].nneighb = 9;
                
                for (p=-1; p<2; p++) for (q=-1; q<2; q++)
                {
                    i1 = (i+p+HASHX)%HASHX;
                    j1 = (j+q+HASHY)%HASHY;
                    hashgrid[mhash(i,j)].neighbour[3*(p+1)+q+1] = mhash(i1, j1);
                }
            }
            break;
        }
        
        case(2): /* Klein bottle b.c. */
        {
            /* left/right boundaries */
            for (j=0; j<HASHY; j++) for (i=0; i<HASHX; i+=HASHX-1) 
            {
                hashgrid[mhash(i, j)].nneighb = 9;
                
                for (p=-1; p<2; p++) for (q=-1; q<2; q++)
                {
                    i1 = (i+p+HASHX)%HASHX;
                    if (((i == 0)&&(p >= 0))||((i == HASHX-1)&&(p <= 0))) j1 = (j+q+HASHY)%HASHY;
                    else j1 = (2*HASHY-1-j-q)%HASHY;
                    hashgrid[mhash(i,j)].neighbour[3*(p+1)+q+1] = mhash(i1, j1);
                }
            }
            
            /* top/bottom boundaries */
            for (i=1; i<HASHX-1; i++) for (j=0; j<HASHY; j+=HASHY-1) 
            {
                hashgrid[mhash(i, j)].nneighb = 9;
                
                for (p=-1; p<2; p++) for (q=-1; q<2; q++)
                {
                    i1 = (i+p+HASHX)%HASHX;
                    j1 = (j+q+HASHY)%HASHY;
                    hashgrid[mhash(i,j)].neighbour[3*(p+1)+q+1] = mhash(i1, j1);
                }
            }
            break;
        }
        
        case(3): /* Boy surface/projective plane b.c. */
        {
            /* left/right boundaries */
            for (j=1; j<HASHY-1; j++) for (i=0; i<HASHX; i+=HASHX-1) 
            {
                hashgrid[mhash(i, j)].nneighb = 9;
                
                for (p=-1; p<2; p++) for (q=-1; q<2; q++)
                {
                    i1 = (i+p+HASHX)%HASHX;
                    if (((i == 0)&&(p >= 0))||((i == HASHX-1)&&(p <= 0))) j1 = (j+q+HASHY)%HASHY;
                    else j1 = (2*HASHY-1-j-q)%HASHY;
                    hashgrid[mhash(i,j)].neighbour[3*(p+1)+q+1] = mhash(i1, j1);
                }
            }
            
            /* top/bottom boundaries */
            for (i=1; i<HASHX-1; i++) for (j=0; j<HASHY; j+=HASHY-1) 
            {
                hashgrid[mhash(i, j)].nneighb = 9;
                
                for (p=-1; p<2; p++) for (q=-1; q<2; q++)
                {
                    if (((j == 0)&&(q >= 0))||((j == HASHY-1)&&(q <= 0))) i1 = (i+p+HASHX)%HASHX;
                    else i1 = (2*HASHX-1-i-p)%HASHX;
                    j1 = (j+q+HASHY)%HASHY;
                    hashgrid[mhash(i,j)].neighbour[3*(p+1)+q+1] = mhash(i1, j1);
                }
            }
            
            /* corners */
            for (i=0; i<HASHX; i+=HASHX-1) for (j=0; j<HASHY; j+=HASHY-1) hashgrid[mhash(i, j)].nneighb = 7;
            
            for (p=0; p<2; p++) for (q=0; q<2; q++) hashgrid[mhash(0,0)].neighbour[2*p+q] = mhash(p,q);
            hashgrid[mhash(0,0)].neighbour[4] = mhash(HASHX-1,HASHY-1);
            hashgrid[mhash(0,0)].neighbour[5] = mhash(HASHX-2,HASHY-1);
            hashgrid[mhash(0,0)].neighbour[6] = mhash(HASHX-1,HASHY-2);

            for (p=0; p<2; p++) for (q=0; q<2; q++) hashgrid[mhash(HASHX-1,0)].neighbour[2*p+q] = mhash(HASHX-1-p,q);
            hashgrid[mhash(HASHX-1,0)].neighbour[4] = mhash(0,HASHY-1);
            hashgrid[mhash(HASHX-1,0)].neighbour[5] = mhash(1,HASHY-1);
            hashgrid[mhash(HASHX-1,0)].neighbour[6] = mhash(0,HASHY-2);
            
            for (p=0; p<2; p++) for (q=0; q<2; q++) hashgrid[mhash(0,HASHY-1)].neighbour[2*p+q] = mhash(p,HASHY-1-q);
            hashgrid[mhash(0,HASHY-1)].neighbour[4] = mhash(HASHX-1,0);
            hashgrid[mhash(0,HASHY-1)].neighbour[5] = mhash(HASHX-2,1);
            hashgrid[mhash(0,HASHY-1)].neighbour[6] = mhash(HASHX-1,2);

            for (p=0; p<2; p++) for (q=0; q<2; q++) hashgrid[mhash(HASHX-1,HASHY-1)].neighbour[2*p+q] = mhash(HASHX-1-p,HASHY-1-q);
            hashgrid[mhash(HASHX-1,HASHY-1)].neighbour[4] = mhash(0,1);
            hashgrid[mhash(HASHX-1,HASHY-1)].neighbour[5] = mhash(1,1);
            hashgrid[mhash(HASHX-1,HASHY-1)].neighbour[6] = mhash(0,2);
            break;
        }
        
        case(4): /* genus 2 (L-shape) b.c. */
        {
            if ((HASHX%2 == 1)||(HASHY%2 == 1))
            {
                printf("Error: HASHX and HASHY must be even for this boundary condition\n");
                exit(0);
            }
            
            /* dummy value for unused cells */
            for (i=HASHX/2; i<HASHX; i++) for (j=HASHY/2; j<HASHY; j++) hashgrid[mhash(i, j)].nneighb = 0;
            
            /* left boundary */
            for (j=0; j<HASHY; j++) 
            {
                hashgrid[mhash(0, j)].nneighb = 9;
                
                for (p=-1; p<2; p++) for (q=-1; q<2; q++)
                {
                    j1 = (j+q+HASHY)%HASHY;
                    if (j1 < HASHY/2) i1 = (p+HASHX)%HASHX;
                    else i1 = (p+HASHX/2)%(HASHX/2);
                    hashgrid[mhash(0, j)].neighbour[3*(p+1)+q+1] = mhash(i1, j1);
                }
            }
            
            /* right boundary */
            for (j=0; j<HASHY; j++) 
            {
                if (j < HASHY/2)
                {
                    hashgrid[mhash(HASHX-1, j)].nneighb = 9;
                
                    for (p=-1; p<2; p++) for (q=-1; q<2; q++)
                    {
                        j1 = (j+q+HASHY)%HASHY;
                        if (j1 < HASHY/2) i1 = (HASHX-1+p)%HASHX;
                        else i1 = (HASHX/2-1+p)%(HASHX/2);
                        hashgrid[mhash(HASHX-1,j)].neighbour[3*(p+1)+q+1] = mhash(i1, j1);
                    }
                }
                else 
                {
                    hashgrid[mhash(HASHX/2-1, j)].nneighb = 9;
//                     hashgrid[mhash(HASHX-1, j)].nneighb = 0;        /* dummy value for unused boundary */
                
                    for (p=-1; p<2; p++) for (q=-1; q<2; q++)
                    {
                        j1 = (j+q+HASHY)%HASHY;
                        if (j1 < HASHY/2) i1 = (HASHX-1+p)%HASHX;
                        else i1 = (HASHX/2-1+p)%(HASHX/2);
                        hashgrid[mhash(HASHX/2-1, j)].neighbour[3*(p+1)+q+1] = mhash(i1, j1);
                    }
                }
                
            }
            
            /* bottom boundary */
            for (i=1; i<HASHX-1; i++) 
            {
                hashgrid[mhash(i, 0)].nneighb = 9;
                
                for (p=-1; p<2; p++) for (q=-1; q<2; q++)
                {
                    i1 = (i+p+HASHX)%HASHX;
                    if (i1 < HASHX/2) j1 = (q+HASHY)%HASHY;
                    else j1 = (q+HASHY/2)%(HASHY/2);
                    hashgrid[mhash(i,0)].neighbour[3*(p+1)+q+1] = mhash(i1, j1);
                }
            }
            
            /* top boundary */
            for (i=1; i<HASHX-1; i++) 
            {
                if (i < HASHX/2)
                {
                    hashgrid[mhash(i, HASHY-1)].nneighb = 9;
                
                    for (p=-1; p<2; p++) for (q=-1; q<2; q++)
                    {
                        i1 = (i+p+HASHX/2)%(HASHX/2);
                        if (i1 < HASHX/2) j1 = (HASHY-1+q)%HASHY;
                        else j1 = (HASHY/2-1+q)%(HASHY/2);
                        hashgrid[mhash(i, HASHY-1)].neighbour[3*(p+1)+q+1] = mhash(i1, j1);
                    }
                }
                else
                {
                    hashgrid[mhash(i, HASHY/2-1)].nneighb = 9;
//                     hashgrid[mhash(i, HASHY-1)].nneighb = 0;    /* dummy value for unused boundary */
                
                    for (p=-1; p<2; p++) for (q=-1; q<2; q++)
                    {
                        i1 = (i+p+HASHX)%HASHX;
                        if (i1 < HASHX/2) j1 = (HASHY-1+q)%HASHY;
                        else j1 = (HASHY/2-1+q)%(HASHY/2);
                        hashgrid[mhash(i, HASHY/2-1)].neighbour[3*(p+1)+q+1] = mhash(i1, j1);
                    }
                }                    
            }
            
            /* TO DO : add more cells for "corners" ? */
            break;
        }
        case (5):   /* spherical boundary conditions */
        {
            printf("Left border\n");
            /* left border */
            for (j=1; j<HASHY-1; j++)
            {
                padding = dx_sphere[j];
                m = mhash(0, j);
//                 printf("j = %i, m = %i\n", j, m);
                if (j==1)
                {
                    hashgrid[m].nneighb = 4;
                    hashgrid[m].neighbour[0] = mhash(hashx_sphere[j]-1,j);
                    hashgrid[m].neighbour[1] = mhash(0,j);
                    hashgrid[m].neighbour[2] = mhash(1,j);
                    hashgrid[m].neighbour[3] = mhash(hashx_sphere[j+1]-1,j+1);
                }
                else if (j==HASHY-2)
                {
                    hashgrid[m].nneighb = 4;
                    hashgrid[m].neighbour[0] = mhash(hashx_sphere[j]-1,j);
                    hashgrid[m].neighbour[1] = mhash(0,j);
                    hashgrid[m].neighbour[2] = mhash(1,j);
                    hashgrid[m].neighbour[3] = mhash(hashx_sphere[j-1]-1,j-1);
                }
                else
                {
                    hashgrid[m].nneighb = 5;
                    hashgrid[m].neighbour[0] = mhash(hashx_sphere[j]-1,j);
                    hashgrid[m].neighbour[1] = mhash(0,j);
                    hashgrid[m].neighbour[2] = mhash(1,j);
                    hashgrid[m].neighbour[3] = mhash(hashx_sphere[j+1]-1,j+1);
                    hashgrid[m].neighbour[4] = mhash(hashx_sphere[j-1]-1,j-1);
                    
//                     hashgrid[m].nneighb = 7;
//                     hashgrid[m].neighbour[5] = mhash(2,j);
//                     hashgrid[m].neighbour[6] = mhash(hashx_sphere[j]-2,j);                    
                }
                
                
                /* neighbours in layer above and below */
                for (p=-1; p<2; p+=2)
                    for (i1=0; i1<hashx_sphere[j+p]; i1++)
                    {
                        m1 = mhash(i1, j+p);
                        x1 = hashgrid[m].x1 - padding;
                        x2 = hashgrid[m].x2 + padding;
                        xx1 = hashgrid[m1].x1;
                        xx2 = hashgrid[m1].x2;
                        if (((xx2 >= x1)&&(xx2 <= x2))||((xx1 >= x1)&&(xx1 <= x2)))
                        {
                            hashgrid[m].neighbour[hashgrid[m].nneighb] = m1;
                            hashgrid[m].nneighb++;
                        }
                    }
                    
                /* check there are not too many neighbours */
                if (hashgrid[m].nneighb >= HASHMAXNEIGH)
                {
                    printf("(i,j) = (0,%i) Error: HASHMAXNEIGH should be at least %i\n", j, hashgrid[m].nneighb + 1);
                    exit(1);
                }
            }
            
            printf("Right border\n");
            /* right border */
            for (j=1; j<HASHY-1; j++)
            {
                padding = dx_sphere[j];
                m = mhash(hashx_sphere[j]-1, j);
                if (j==1)
                {
                    hashgrid[m].nneighb = 4;
                    hashgrid[m].neighbour[0] = mhash(hashx_sphere[j]-2,j);
                    hashgrid[m].neighbour[1] = mhash(hashx_sphere[j]-1,j);
                    hashgrid[m].neighbour[2] = mhash(0,j);
                    hashgrid[m].neighbour[3] = mhash(0,j+1);
                }
                else if (j==HASHY-2)
                {
                    hashgrid[m].nneighb = 4;
                    hashgrid[m].neighbour[0] = mhash(hashx_sphere[j]-2,j);
                    hashgrid[m].neighbour[1] = mhash(hashx_sphere[j]-1,j);
                    hashgrid[m].neighbour[2] = mhash(0,j);
                    hashgrid[m].neighbour[3] = mhash(0,j-1);
                }
                else
                {
                    hashgrid[m].nneighb = 5;
                    hashgrid[m].neighbour[0] = mhash(hashx_sphere[j]-2,j);
                    hashgrid[m].neighbour[1] = mhash(hashx_sphere[j]-1,j);
                    hashgrid[m].neighbour[2] = mhash(0,j);
                    hashgrid[m].neighbour[3] = mhash(0,j+1);
                    hashgrid[m].neighbour[4] = mhash(0,j-1);
                    
//                     hashgrid[m].nneighb = 7;
//                     hashgrid[m].neighbour[5] = mhash(1,j);
//                     hashgrid[m].neighbour[6] = mhash(hashx_sphere[j]-3,j);
                }
                
                /* neighbours in layer above and below */
                for (p=-1; p<2; p+=2)
                    for (i1=0; i1<hashx_sphere[j+p]; i1++)
                    {
                        m1 = mhash(i1, j+p);
                        x1 = hashgrid[m].x1 - padding;
                        x2 = hashgrid[m].x2 + padding;
                        xx1 = hashgrid[m1].x1;
                        xx2 = hashgrid[m1].x2;
                        if (((xx2 >= x1)&&(xx2 <= x2))||((xx1 >= x1)&&(xx1 <= x2)))
                        {
                            hashgrid[m].neighbour[hashgrid[m].nneighb] = m1;
                            hashgrid[m].nneighb++;
                        }
                    }
                    
                /* check there are not too many neighbours */
                if (hashgrid[m].nneighb >= HASHMAXNEIGH)
                {
                    printf("(i,j) = (%i,%i) Error: HASHMAXNEIGH should be at least %i\n", hashx_sphere[j]-1, j, hashgrid[m].nneighb + 1);
                    exit(1);
                }
            }
            
            printf("Top border\n");
            /* top border/North pole */
            m = mhash(0, HASHY-1);
            hashgrid[m].nneighb = hashx_sphere[HASHY-2] + 1;
            hashgrid[m].neighbour[0] = m;
            fprintf(lj_log, "Top border: m = %i, %i neighbours\n", m, hashgrid[m].nneighb);
            for (i1=0; i1<hashx_sphere[HASHY-2]; i1++)
                hashgrid[m].neighbour[i1+1] = mhash(i1,HASHY-2);
            
            /* check there are not too many neighbours */
            if (hashgrid[m].nneighb >= HASHMAXNEIGH)
            {
                printf("Error at North Pole: HASHMAXNEIGH should be at least %i\n", hashgrid[m].nneighb + 1);
                exit(1);
            }
            
            printf("Bottom border\n");
            /* bottom border/South pole */
            m = mhash(0, 0);
            hashgrid[m].nneighb = hashx_sphere[1] + 1;
            hashgrid[m].neighbour[0] = m;
            for (i1=0; i1<hashx_sphere[1]; i1++)
                hashgrid[m].neighbour[i1+1] = mhash(i1,1);
            
            /* check there are not too many neighbours */
            if (hashgrid[m].nneighb >= HASHMAXNEIGH)
            {
                printf("Error at South Pole: HASHMAXNEIGH should be at least %i\n", hashgrid[m].nneighb + 1);
                exit(1);
            }
            
            /* symmetrize hashgrid */
            for (i=0; i<HASHX*HASHY; i++) 
            {
                for (j=0; j<hashgrid[i].nneighb; j++)
                {
                    m = hashgrid[i].neighbour[j];
                    sym = 0;
                    for (k=0; k<hashgrid[m].nneighb; k++)
                        if (hashgrid[m].neighbour[k] == i) sym = 1;
                    if (!sym) 
                    {
                        fprintf(lj_log, "Symmetrising hashcells %i and %i\n", i, m);
                        hashgrid[m].neighbour[hashgrid[m].nneighb] = i;
                        hashgrid[m].nneighb++;
                    }
                }
            }
            
            for (i=0; i<HASHX*HASHY; i++)
            {
                p = hashgrid[i].nneighb;
                printf("Hash cell %i has %i neighbours: ", i, p);
                fprintf(lj_log, "Hash cell %i has %i neighbours: ", i, p);
                for (j=0; j<p; j++)
                {
                    printf("%i ", hashgrid[i].neighbour[j]);
                    fprintf(lj_log, "%i ", hashgrid[i].neighbour[j]);
                }
                printf("\n");
                fprintf(lj_log, "\n");
            }
            
            
            
//             sleep(5);
            
            printf("Hashgrid initialised\n");
            
            break;
        }
        
        default: /* do nothing */;
    }
    
//     for (i=0; i<HASHX; i++) 
//     {
//         for (j=0; j<HASHY; j++)
//         {
//             for (k=0; k<hashgrid[mhash(i,j)].nneighb; k++)
//             {
//                 m = hashgrid[mhash(i,j)].neighbour[k];
//                 p = m/HASHY;
//                 q = m%HASHY;
//                 if (vabs((double)(p-i)) + vabs((double)(q-j)) > 2.0)
//                 printf("Grid cell (%i, %i) - neighbour %i = %i = (%i, %i)\n", i, j, k, m, p, q);
//             }
// //         sleep(1);
//         }
// //         sleep(1);
//     }

//     if (COLOR_BACKGROUND) init_hashcell_coordinates(hashgrid);
    
//     sleep(1);
}

void update_hashgrid(t_particle* particle, t_obstacle *obstacle, t_hashgrid* hashgrid, int verbose)
{
    int i, j, k, n, m, max = 0, hashcell, maxcell = -1;
    
//     printf("Updating hashgrid_number\n");
    for (i=0; i<HASHX*HASHY; i++) hashgrid[i].number = 0;
//     printf("Updated hashgrid_number\n");
        
    /* place each particle in hash grid */
    for (k=0; k<ncircles; k++)
        if (particle[k].active)
        {
//             printf("placing circle %i\t", k);
            hashcell = hash_cell(particle[k].xc, particle[k].yc);
            n = hashgrid[hashcell].number;
            if (n < HASHMAX) hashgrid[hashcell].particles[n] = k;
            else printf("Too many particles in hash cell (%i, %i), try increasing HASHMAX\n", i, j);
            hashgrid[hashcell].number++;
            particle[k].hashcell = hashcell;
            
            if (n > max) 
            {
                max = n;
                maxcell = hashcell;
            }
        }
        
    if (OSCILLATE_OBSTACLES) 
    {
        for (i=0; i<HASHX*HASHY; i++) hashgrid[i].nobs = 0;
        for (k=0; k<nobstacles; k++)
            if (obstacle[k].active)
            {
                hashcell = hash_cell(obstacle[k].xc, obstacle[k].yc);
                /* note: only one obstacle per cell is counted */
                hashgrid[hashcell].nobs = 1;
                hashgrid[hashcell].obstacle = k;
            }
    }
    
    if(verbose) printf("Maximal number of particles per hash cell: %i, in cell %i\n", max, maxcell);
}


int wrap_particle(t_particle* particle, double *px, double *py)
/* relocate particles in case of periodic and similar boundary conditions */
{
    double x, y, x1, y1, x2, y2;
    int move = 0;
    
    x = particle->xc;
    y = particle->yc;
    x1 = x;
    y1 = y;
    
    switch (bc_grouped(BOUNDARY_COND)) {
        case (0): /* rectangular b.c. */
        {
            /* do nothing */
            return(0);
            break;
        }
        case (1): /* periodic b.c. */
        {
            if (x < BCXMIN) 
            {
                x1 += BCXMAX - BCXMIN;
                move++;
            }
            else if (x > BCXMAX) 
            {
                x1 += BCXMIN - BCXMAX;
                move++;
            }
            if (y > BCYMAX) 
            {
                y1 += BCYMIN - BCYMAX;
                move++;
            }
            else if (y < BCYMIN) 
            {
                y1 += BCYMAX - BCYMIN;
                move++;
            }
            particle->xc = x1;
            particle->yc = y1;
            return(move);
            break;
        }
        case (2): /* Klein bottle b.c. */
        {
            if (y > BCYMAX) 
            {
                y1 += BCYMIN - BCYMAX;
                move++;
            }
            else if (y < BCYMIN) 
            {
                y1 += BCYMAX - BCYMIN;
                move++;
            }
            if (x < BCXMIN) 
            {
                x1 += BCXMAX - BCXMIN;
                y1 = -y1;
                *py *= -1.0;
                particle->angle = DPI - particle->angle;
                move++;
            }
            else if (x > BCXMAX) 
            {
                x1 += BCXMIN - BCXMAX;
                y1 = -y1;
                particle->angle = DPI - particle->angle;
                *py *= -1.0;
                move++;
            }
            particle->xc = x1;
            particle->yc = y1;
            return(move);
            break;
        }
        case (3): /* Boy surface b.c. */
        {
            if ((y < BCYMAX)&&(y > BCYMIN))
            {
                if (x > BCXMAX) 
                {
                    x1 += BCXMIN - BCXMAX;
                    y1 = -y1;
                    particle->angle *= -1.0;
                    particle->vy *= -1.0;
                    *py *= -1.0;
                    move++;
                }
                else if (x < BCXMIN) 
                {
                    x1 += BCXMAX - BCXMIN;
                    y1 = -y1;
                    particle->angle *= -1.0;
                    particle->vy *= -1.0;
                    *py *= -1.0;
                    move++;
                }
            }
//             x = x1;
//             y = y1;
            if ((x < BCXMAX)&&(x > BCXMIN))
            {
                if (y > BCYMAX) 
                {
                    y1 += BCYMIN - BCYMAX;
                    x1 = -x1;
                    particle->angle = PI - particle->angle;
                    particle->vx *= -1.0;
                    *px *= -1.0;
                    move++;
                }
                else if (y < BCYMIN) 
                {
                    y1 += BCYMAX - BCYMIN;
                    x1 = -x1;
                    particle->angle = PI - particle->angle;
                    particle->vx *= -1.0;
                    *px *= -1.0;
                    move++;
                }
            }
            if (((x >= BCXMAX)||(x <= BCXMIN))&&((y >= BCYMAX)||(y <= BCYMIN)))
            {
                /* This case can lead to numerical instabilities, and can be avoided by putting obstacles in the corners */
                printf("Double wrap!\n");
                if (x >= BCXMAX) x1 = BCXMAX - BCXMIN - x1;
                else x1 = -BCXMAX + BCXMIN - x1;
                if (y >= BCXMAX) y1 = BCYMAX - BCYMIN - y1;
                else y1 = -BCYMAX + BCYMIN - y1;
                
                particle->vx *= -1.0;
                *px *= -1.0;
                particle->vy *= -1.0;
                *py *= -1.0;
                
                
//                 if (x1 >= BCXMAX) x1 = BCXMAX - 1.0e-5;
//                 else if (x1 <= BCXMIN) x1 = BCXMIN + 1.0e-5;
//                 if (y1 >= BCXMAX) y1 = BCYMAX - 1.0e-5;
//                 else if (y1 <= BCYMIN) y1 = BCYMIN + 1.0e-5;
                move++;
            }
                            
            particle->xc = x1;
            particle->yc = y1;
            return(move);
            break;
        }
        case (4): /* genus two (L-shaped domain) b.c. */
        {
            if ((x > 0.0)&&(y > 0.0)) 
            {
                if (x > y) y1 -= 0.5*(BCYMAX - BCYMIN);
                else x1 -= 0.5*(BCXMAX - BCXMIN);
                move++;
            }
            
            else
            {
            
                if (x < BCXMIN) 
                {
                    if (y < 0.0) x1 += BCXMAX - BCXMIN;
                    else x1 += 0.5*(BCXMAX - BCXMIN);
                    move++;
                }
                else 
                {
                    if ((y < 0.0)&&(x > BCXMAX)) 
                    {
                        x1 += BCXMIN - BCXMAX;
                        move++;
                    }
                    else if ((y >= 0.0)&&(x > 0.0)&&(x < OBSTACLE_RADIUS)) 
                    {
                        x1 += 0.5*(BCXMIN - BCXMAX);
                        move++;
                    }
                }
            
                if (y < BCYMIN) 
                {
                    if (x1 < 0.0) y1 += BCYMAX - BCYMIN;
                    else y1 += 0.5*(BCYMAX - BCYMIN);
                    move++;
                }
                else 
                {
                    if ((x1 < 0.0)&&(y > BCYMAX)) 
                    {
                        y1 += BCYMIN - BCYMAX;
                        move++;
                    }
                    else if ((x1 >= 0.0)&&(y > 0.0)&&(y < OBSTACLE_RADIUS)) 
                    {
                        y1 += 0.5*(BCYMIN - BCYMAX);
                        move++;
                    }
                }
            }
            
//             if (move > 0) printf("Moved particle from (%.3lg, %.3lg) to (%.3lg, %.3lg)\n", x, y, x1, y1);
            
            particle->xc = x1;
            particle->yc = y1;
            return(move);
            break;
        }
        case (5):   /* spherical bc */
        {
            if (x < 0.0) 
            {
                x1 += DPI;
//                 printf("Moving particle by 2Pi\n");
                move++;
            }
            else if (x > DPI) 
            {
                x1 -= DPI;
//                 printf("Moving particle by -2Pi\n");
                move++;
            }
            if (y > PI) 
            {
                x1 += PI;
                if (x1 > DPI) x1 -= DPI;
                y1 = PI;
                particle->vy *= -1.0;
                *py *= -1.0;
                move++;
            }
            else if (y < 0.0) 
            {
                x1 += PI;
                if (x1 > DPI) x1 -= DPI;
                y1 = 0.0;
                particle->vy *= -1.0;
                *py *= -1.0;
                move++;
            }
            particle->xc = x1;
            particle->yc = y1;
            return(move);
            break;
        }
        default:
        {
            /* do nothing */
            break;
        }
    }
    
}


void wrap_relative_positions(double x1, double y1, double *x2, double *y2)
/* computes relative positions of particles, taking boundary conditions into account */
{
double dxhalf = 0.5*(BCXMAX - BCXMIN), dyhalf = 0.5*(BCYMAX - BCYMIN);
double dx, dy, x3, y3, x4, y4, dh, dv, wrap_x, wrap_y;
int verbose = 0;

    dx = *x2 - x1;
    dy = *y2 - y1;
        
    switch (bc_grouped(BOUNDARY_COND)) {
        case (0): /* rectangular b.c. */
        {
        /* do nothing */
            break;
        }
        case (1): /* periodic b.c. */
        {
            if (dx > dxhalf) *x2 -= BCXMAX - BCXMIN;
            else if (-dx > dxhalf) *x2 += BCXMAX - BCXMIN;
            if (dy > dyhalf) *y2 -= BCYMAX - BCYMIN;
            else if (-dy > dyhalf) *y2 += BCYMAX - BCYMIN;
//             printf("(x2,y2) = (%.3lg, %.3lg)\n", *x2, *y2);
            break;
        }
        case (2): /* Klein bottle b.c. */
        {
            if (dx > dxhalf) 
            {
                *x2 -= BCXMAX - BCXMIN;
                *y2 *= -1.0;
            }
            else if (-dx > dxhalf) 
            {
                *x2 += BCXMAX - BCXMIN;
                *y2 *= -1.0;
            }     
            dy = *y2 - y1;
            if (dy > dyhalf) *y2 -= BCYMAX - BCYMIN;
            else if (-dy > dyhalf) *y2 += BCYMAX - BCYMIN;
            break;
        }
        case (3): /* Boy surface - TO FIX */
        {
//             wrap_x = 2.0*(BCXMAX - BCXMIN)/(double)HASHX;
//             wrap_y = 2.0*(BCYMAX - BCYMIN)/(double)HASHY;
            wrap_x = dxhalf;
            wrap_y = dyhalf;
            
            /* find which wrapped point is closest to (x1, y1) */
            dh = 100.0;
            if (dx > wrap_x)
            {
                x3 = *x2 - BCXMAX + BCXMIN;
                y3 = -*y2;
                dh = (x3-x1)*(x3-x1) + (y3-y1)*(y3-y1);
            }
            else if (-dx > wrap_x)
            {
                x3 = *x2 + BCXMAX - BCXMIN;
                y3 = -*y2;
                dh = (x3-x1)*(x3-x1) + (y3-y1)*(y3-y1);
            }
            if ((verbose)&&(dh < 100.0)&&(dh > 0.0))
            {
                printf("Case 1: (x1,y1) = (%.3lg, %.3lg)\t", x1, y1);
                printf("(x2,y2) = (%.3lg, %.3lg) -> (%.3lg, %.3lg)\n", *x2, *y2, x3, y3);
            }

            dv = 100.0;
            if (dy > wrap_y)
            {
                x4 = -*x2;
                y4 = *y2 - BCYMAX + BCYMIN;
                dv = (x4-x1)*(x4-x1) + (y4-y1)*(y4-y1);
            }
            else if (-dy > wrap_y)
            {
                x4 = -*x2;
                y4 = *y2 + BCYMAX - BCYMIN;
                dv = (x4-x1)*(x4-x1) + (y4-y1)*(y4-y1);
            }
            if ((verbose)&&(dv < 100.0)&&(dv > 0.0))
            {
                printf("Case 2: (x1,y1) = (%.3lg, %.3lg)\t", x1, y1);
                printf("(x2,y2) = (%.3lg, %.3lg) -> (%.3lg, %.3lg)\n", *x2, *y2, x4, y4);
            }
            
            if (dh < 100.0)
            {
                *x2 = x3;
                *y2 = y3;
            }
            if (dv < dh)
            {
                *x2 = x4;
                *y2 = y4;
            }
            break;
        }
        case (4): /* genus 2 (L-shaped) b.c. */
        {
            x3 = *x2;
            y3 = *y2;
            if ((x1 < 0.0)&&(y1 < 0.0))
            {
                if (dx > dxhalf) *x2 -= (BCXMAX - BCXMIN);
                if (dy > dyhalf) *y2 -= (BCYMAX - BCYMIN);
            }
            else if ((x1 >= 0.0)&&(y1 < 0.0))
            {
                if (dx < -dxhalf) *x2 += (BCXMAX - BCXMIN);
                if (dy > 0.5*dyhalf) *y2 -= 0.5*(BCYMAX - BCYMIN);
                else if (dy < -0.5*dyhalf) *y2 += 0.5*(BCYMAX - BCYMIN);
            }
            else if ((x1 < 0.0)&&(y1 >= 0.0))
            {
                if (dy < -dyhalf) *y2 += BCYMAX - BCYMIN;
                if (dx > 0.5*dxhalf) *x2 -= 0.5*(BCXMAX - BCXMIN);
                else if (dx < -0.5*dxhalf) *x2 += 0.5*(BCXMAX - BCXMIN);
            }
            if ((verbose)&&(module2(*x2 - x1, *y2 - y1) > dyhalf))
            {
                printf("(x1,y1) = (%.3lg, %.3lg)\n", x1, y1);
                printf("(x2,y2) = (%.3lg, %.3lg) ", x3, y3);
                printf("-> (%.3lg, %.3lg)\n", *x2, *y2);
            }
//             printf("(x2,y2) = (%.3lg, %.3lg)\n", *x2, *y2);
            break;
        }
       case (5): /* spherical b.c. */
        {
            if (dx > dxhalf) *x2 -= BCXMAX - BCXMIN;
            else if (-dx > dxhalf) *x2 += BCXMAX - BCXMIN;
//             if (dy > dyhalf) *y2 -= BCYMAX - BCYMIN;
//             else if (-dy > dyhalf) *y2 += BCYMAX - BCYMIN;
            /* TODO: poles? */
//             printf("(x2,y2) = (%.3lg, %.3lg)\n", *x2, *y2);
            break;
        }
     }
}


void compute_relative_positions(t_particle particle[NMAXCIRCLES], t_hashgrid hashgrid[HASHX*HASHY])
{
    int j, i0, j0, m0, k, m, p, q, n = 0;
    double x1, x2, y1, y2, xtemp, ytemp;
    
    for (j=0; j < ncircles; j++) if (particle[j].active)
    {
//         i0 = particle[j].hashx;
//         j0 = particle[j].hashy;
//         m0 = mhash(i0, j0);
        m0 = particle[j].hashcell;
        x1 = particle[j].xc;
        y1 = particle[j].yc;
        n = 0;
    
        for (q=0; q<hashgrid[m0].nneighb; q++)
        {
            m = hashgrid[m0].neighbour[q];
            for (k=0; k<hashgrid[m].number; k++)
            if ((hashgrid[m].particles[k]!=j)&&(particle[hashgrid[m].particles[k]].active))
            {
                if (n < 9*HASHMAX)
                {
                    p = hashgrid[m].particles[k];
                
                    x2 = particle[p].xc;
                    y2 = particle[p].yc;
                    xtemp = x2;
                    ytemp = y2;
                
                    if (bc_grouped(BOUNDARY_COND) != 0) wrap_relative_positions(x1, y1, &x2, &y2);
                
                    particle[j].hashneighbour[n] = p;
                    particle[j].deltax[n] = x2 - x1;
                    particle[j].deltay[n] = y2 - y1;
//                     if ((j%50 == 0)&&((vabs(x2-x1)>1.0)||(vabs(y2-y1)>1.0))) 
//                     if (((vabs(x2-x1)>0.7)||(vabs(y2-y1)>0.7))) 
//                     {
//                         printf("(x1, y1) = (%.3lg, %.3lg), (x2, y2) = (%.3lg, %.3lg) -> (%.3lg, %.3lg)\n", x1, y1, xtemp, ytemp, x2, y2);
//                         printf("particle[%i].delta[%i] = (%.4lg, %.4lg)\n", j, n, particle[j].deltax[n], particle[j].deltay[n]);
//                     }
                    n++;
                }
                else printf("Not enough memory in particle.deltax, particle.deltay\n");
            }
        }
        particle[j].hash_nneighb = n;
    }
}

