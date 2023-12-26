/* routines dealing with hashgrid for lennardjones.c */
/* Note: a possible improvement would be to use pointers instead of integers */
/* when referencing other hashgrid cells or particles */ 

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
    static double lx, ly, padding;
    int i, j;
    
    if (first)
    {
        if (!NO_WRAP_BC) padding = 0.0;
        else padding = HASHGRID_PADDING;
        lx = BCXMAX - BCXMIN + 2.0*padding;
        ly = BCYMAX - (BCYMIN) + 2.0*padding;
        first = 0;
    }
    
    if (CENTER_VIEW_ON_OBSTACLE) x -= xshift;
    
    i = (int)((double)HASHX*(x - BCXMIN + padding)/lx);
    j = (int)((double)HASHY*(y - BCYMIN + padding)/ly);
    
    if (i<0) i = 0;
    else if (i>=HASHX) i = HASHX-1;
    if (j<0) j = 0;
    else if (j>=HASHY) j = HASHY-1;
    
    return(mhash(i,j));
//     printf("Mapped (%.3f,%.3f) to (%i, %i)\n", x, y, ij[0], ij[1]);
}


void init_hashcell_coordinates(t_hashgrid hashgrid[HASHX*HASHY])
/* initialise coordinates of corners of hashcells, needed for option COLOR_BACKGROUND */
{
    static int first = 1;
    static double lx, ly, padding, dx, dy;
    int i, j, n;
    double x, y;
    
    if (first)
    {
        if (!NO_WRAP_BC) padding = 0.0;
        else padding = HASHGRID_PADDING;
        lx = BCXMAX - BCXMIN + 2.0*padding;
        ly = BCYMAX - (BCYMIN) + 2.0*padding;
        dx = 1.0*lx/(double)HASHX;
        dy = 1.0*ly/(double)HASHY;
        first = 0;
    }
    
    for (i=0; i<HASHX; i++)
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
    int i, j, k, p, q, m, i1, j1;
    
    printf("Initializing hash grid\n");
    
    /* bulk of the table */
    for (i=0; i<HASHX-1; i++)
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

    if (COLOR_BACKGROUND) init_hashcell_coordinates(hashgrid);
    
//     sleep(1);
}

void update_hashgrid(t_particle* particle, t_hashgrid* hashgrid, int verbose)
{
    int i, j, k, n, m, max = 0, hashcell;
    
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
            
            if (n > max) max = n;
        }
    
    if(verbose) printf("Maximal number of particles per hash cell: %i\n", max);
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

