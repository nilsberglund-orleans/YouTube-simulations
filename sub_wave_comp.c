/* some changes in sub_wave.c required by wave_comparison.c */ 

short int circletop[NMAXCIRCLES];  /* set to 1 if circle is in top half */


void init_circle_config_half(int pattern, int top, t_circle circles[NMAXCIRCLES])
/* initialise the arrays circlex, circley, circlerad and circleactive */
/* for billiard shape D_CIRCLES */
{
    int i, j, k, n, ncirc0, n_p_active, ncandidates=5000, naccepted, nnew; 
    double dx, dy, p, phi, r, r0, ra[5], sa[5], height, x, y = 0.0, gamma, ymean, ytop, ybottom, dpoisson = 3.05*MU, xmax;
    short int active_poisson[NMAXCIRCLES], far;
   
    ymean = 0.5*(YMIN + YMAX);
    switch (pattern) {
        case (C_SQUARE):
        {
            dy = (YMAX - YMIN)/((double)NGRIDY);
            for (i = 0; i < NGRIDX; i++)
                for (j = 0; j < NGRIDY/2; j++)
                {
                    printf("i = %i, j = %i, n = %i\n", i, j, n);
                    n = ncircles + NGRIDY*i/2 + j;
                    circles[n].xc = ((double)(i-NGRIDX/2) + 0.5)*dy;
                    y = ((double)j + 0.5)*dy;
                    if (top) circles[n].yc = ymean + y;
                    else circles[n].yc = ymean - y;
                    if (top) circles[n].radius = MU;
                    else circles[n].radius = MUB;
                    circles[n].active = 1;
                    circletop[n] = top;
                }
            ncircles += NGRIDX*NGRIDY/2;
            break;
        }
        case (C_HEX):
        {
            dy = (YMAX - YMIN)/((double)NGRIDY);
            dx = dy*0.5*sqrt(3.0);
            for (i = 0; i < NGRIDX; i++)
                for (j = 0; j < NGRIDY/2+2; j++)
                {
                    n = ncircles + (NGRIDY+2)*i/2 + j;
                    circles[n].xc = ((double)(i-NGRIDX/2) + 0.5)*dy;
                    y = ((double)j - 0.5)*dy;
                    if ((i+NGRIDX)%2 == 1) y += 0.5*dy;
                    if (top) circles[n].yc = ymean + 0.5*dy + y;
                    else circles[n].yc = ymean - 0.5*dy - y;
                    if (top) circles[n].radius = MU;
                    else circles[n].radius = MUB;
                    circles[n].active = 1;
                    circletop[n] = top;
                }
            ncircles += NGRIDX*(NGRIDY+2)/2;
            break;
        }
        case (C_HEX_NONUNIF):
        {
            dy = (YMAX - YMIN)/((double)NGRIDY);
            dx = dy*0.5*sqrt(3.0);
            for (i = 0; i < NGRIDX; i++)
                for (j = 0; j < NGRIDY/2+2; j++)
                {
                    n = ncircles + (NGRIDY+2)*i/2 + j;
                    x = ((double)(i-NGRIDX/2) + 0.5)*dy;
                    xmax = ((double)(NGRIDX/2) - 0.5)*dy;
                    if (top) circles[n].xc = x - HEX_NONUNIF_COMPRESSSION*(x*x - xmax*xmax);
                    else circles[n].xc = x - HEX_NONUNIF_COMPRESSSION_B*(x*x - xmax*xmax);
                    y = ((double)j - 0.5)*dy;
                    if ((i+NGRIDX)%2 == 1) y += 0.5*dy;
                    if (top) circles[n].yc = ymean + 0.5*dy + y;
                    else circles[n].yc = ymean - 0.5*dy - y;
                    if (top) circles[n].radius = MU;
                    else circles[n].radius = MUB;
                    circles[n].active = 1;
                    circletop[n] = top;
                }
            ncircles += NGRIDX*(NGRIDY+2)/2;
            break;
        }
        case (C_RAND_DISPLACED):
        {
            dy = (YMAX - YMIN)/((double)NGRIDY);
            for (i = 0; i < NGRIDX; i++)
                for (j = 0; j < NGRIDY/2; j++)
                {
                    n = ncircles + NGRIDY*i/2 + j;
                    circles[n].xc = ((double)(i-NGRIDX/2) + 0.5*((double)rand()/RAND_MAX - 0.5))*dy;
                    if (NGRIDX%2 == 0) circles[n].xc += 0.5*dy;
                    y = ((double)j + 0.5 + 0.5*((double)rand()/RAND_MAX - 0.5))*dy;
                    if (top) circles[n].yc = ymean + y;
                    else circles[n].yc = ymean - y;
                    if (top) circles[n].radius = MU*sqrt(1.0 + 0.8*((double)rand()/RAND_MAX - 0.5));
                    else circles[n].radius = MUB*sqrt(1.0 + 0.8*((double)rand()/RAND_MAX - 0.5));
                    circles[n].active = 1;
                    circletop[n] = top;
                    printf("n = %i, x = %.3lg\n", n, circles[n].xc);
                }
            ncircles += NGRIDX*NGRIDY/2;
            break;
        }
        case (C_RAND_PERCOL):
        {
            dy = (YMAX - YMIN)/((double)NGRIDY);
            for (i = 0; i < NGRIDX; i++)
                for (j = 0; j < NGRIDY/2; j++)
                {
                    n = ncircles + NGRIDY*i/2 + j;
                    circles[n].xc = ((double)(i-NGRIDX/2) + 0.5)*dy;
                    y = YMIN + ((double)j + 0.5)*dy;
                    if ( ((top)&&(y > 0.0))||((!top)&&(y <= 0.0)) )
                        circles[n].yc = y;
                    if (top) circles[n].radius = MU;
                    else circles[n].radius = MUB;
                    p = (double)rand()/RAND_MAX;
                    if (p < P_PERCOL) circles[n].active = 1;
                    else circles[n].active = 0;
                    circletop[n] = top;
                }
            ncircles += NGRIDX*NGRIDY/2;
            break;
        }
        case (C_RAND_POISSON):
        {
            for (i = 0; i < NPOISSON/2; i++)
            {
                n = ncircles + i;
                circles[n].xc = LAMBDA*(2.0*(double)rand()/RAND_MAX - 1.0);
                if (top) y = ymean + (YMAX-ymean)*(double)rand()/RAND_MAX;
                else y = ymean + (YMIN-ymean)*(double)rand()/RAND_MAX;
                if ( ((top)&&(y > 0.0))||((!top)&&(y <= 0.0)) )
                    circles[n].yc = y;
                if (top) circles[n].radius = MU;
                else circles[n].radius = MUB;
                circles[n].active = 1;
                circletop[n] = top;
                printf("n = %i, x = %.3lg\n", n, circles[n].xc);
            }
            ncircles += NPOISSON/2;
            break;
        }
        case (C_CLOAK):
        {
            for (i = 0; i < 40; i++)
                for (j = 0; j < 5; j++)
                {
                    n = ncircles + 5*i + j;
                    phi = (double)i*DPI/40.0;
                    r = LAMBDA*0.5*(1.0 + (double)j/5.0);
                    circles[n].xc = r*cos(phi);
                    y = r*sin(phi);
                    if ( ((top)&&(y > 0.0))||((!top)&&(y <= 0.0)) )
                        circles[n].yc = y;
                    if (top) circles[n].radius = MU;
                    else circles[n].radius = MUB;
                    circles[n].active = 1;
                    circletop[n] = top;
                }
            ncircles += 200;
            break;
        }
        case (C_CLOAK_A):   /* optimized model A1 by C. Jo et al */
        {
            ra[0] = 0.0731;     sa[0] = 1.115;
            ra[1] = 0.0768;     sa[1] = 1.292;
            ra[2] = 0.0652;     sa[2] = 1.464;
            ra[3] = 0.056;      sa[3] = 1.633;
            ra[4] = 0.0375;     sa[4] = 1.794;
            for (i = 0; i < 21; i++)
                for (j = 0; j < 5; j++)
                {
                    n = ncircles + 5*i + j;
                    phi = (double)i*DPI/40.0;
                    r = LAMBDA*sa[j];
                    circles[n].xc = r*cos(phi);
                    circles[n].yc = r*sin(phi);
                    
                    if (top) y = r*sin(phi);
                    else y = -r*sin(phi);
                        
                    circles[n].yc = y;
                    
                    circles[n].radius = LAMBDA*ra[j];
                    circles[n].active = 1;
                    circletop[n] = top;
                }
            ncircles += 105;
            
            /* add circle in the center */
            circles[ncircles].xc = 0.0;
            circles[ncircles].yc = 0.0;
            if (top) circles[ncircles].radius = MU;
            else circles[ncircles].radius = MUB;
            circles[ncircles].active = 2;
            circletop[ncircles] = top; 
            ncircles += 1;
            break;
        }
        case (C_POISSON_DISC):
        {
            printf("Generating Poisson disc sample\n");
            /* generate first circle */
            n = ncircles;
            circles[n].xc = LAMBDA*(2.0*(double)rand()/RAND_MAX - 1.0);
            if (top) y = ymean + (YMAX-ymean)*(double)rand()/RAND_MAX;
            else y = ymean + (YMIN-ymean)*(double)rand()/RAND_MAX;
            circles[n].yc = y;
            if (top) circles[n].radius = MU;
            else circles[n].radius = MUB;
            circles[n].active = 1;
            circletop[n] = top;
            active_poisson[n] = 1;
            n_p_active = 1;
            ncirc0 = 1;
            
            while ((n_p_active > 0)&&(ncircles + ncirc0 < NMAXCIRCLES))
            {
                /* randomly select an active circle */
                i = rand()%(ncirc0);
                n = ncircles + i;
                while (!active_poisson[ncircles + i]) i = rand()%(ncirc0);                 
//                 printf("Starting from circle %i at (%.3f,%.3f)\n", i, circles[i].xc, circles[i].yc);
                /* generate new candidates */
                naccepted = 0;
                for (j=0; j<ncandidates; j++)
                {
                    r = dpoisson*(2.0*(double)rand()/RAND_MAX + 1.0);
                    phi = DPI*(double)rand()/RAND_MAX;
                    x = circles[n].xc + r*cos(phi);
                    y = circles[n].yc + r*sin(phi);
//                        printf("Testing new circle at (%.3f,%.3f)\t", x, y);
                    far = 1;
                    for (k=ncircles; k<ncircles + ncirc0; k++) if ((k!=n))
                    {
                        /* new circle is far away from circle k */
                        far = far*((x - circles[k].xc)*(x - circles[k].xc) + (y - circles[k].yc)*(y - circles[k].yc) >= dpoisson*dpoisson);
                        /* new circle is in domain */
                        if (top) far = far*(vabs(x) < LAMBDA)*(y < YMAX)*(y > 0.0);
                        else far = far*(vabs(x) < LAMBDA)*(y > YMIN)*(y < 0.0);
                    }
                    if (far)    /* accept new circle */
                    {
                        printf("New circle at (%.3f,%.3f) accepted\n", x, y);
                        nnew = ncircles + ncirc0;
                        circles[nnew].xc = x;
                        circles[nnew].yc = y;
                        if (top) circles[nnew].radius = MU;
                        else circles[nnew].radius = MUB;
                        circles[nnew].active = 1;
                        active_poisson[nnew] = 1;
//                         circleactive[nnew] = 1;
                        circletop[nnew] = top;
                        ncirc0++;
                        n_p_active++;
                        naccepted++;
                    }
//                     else printf("Rejected\n");
                }
                if (naccepted == 0)    /* inactivate circle i */ 
                {
//                     printf("No candidates work, inactivate circle %i\n", ncircles + i);
                    active_poisson[ncircles + i] = 0;
                    n_p_active--;
                }
                printf("%i active circles\n", n_p_active);
//                 sleep(1);
            }
            
            printf("Already existing: %i circles\n", ncircles);
            ncircles += ncirc0;
            printf("Generated %i circles\n", ncirc0);
            printf("Total: %i circles\n", ncircles);
            break;
        }
        case (C_GOLDEN_MEAN):
        {
            gamma = (sqrt(5.0) - 1.0)*0.5;    /* golden mean */
            height = YMAX - ymean;
            dx = 2.0*LAMBDA/150.0;
            if (top) y = ymean + 0.5*height;
            else y = ymean - 0.5*height; 
            for (n = 0; n < 150; n++)
            {
                circles[ncircles + n].xc = -LAMBDA + n*dx;
                circles[ncircles + n].yc = y;
                
                if (top)
                {
                    y += height*gamma; 
                    if (y > YMAX) y -= height;
                }
                else 
                {
                    y -= height*gamma; 
                    if (y < YMIN) y += height;
                   
                }
                
                if (top) circles[ncircles + n].radius = MU;
                else circles[ncircles + n].radius = MUB;
                circles[ncircles + n].active = 1;
                circletop[ncircles] = top; 
            }
            
            /* test for circles that overlap top or bottom boundary */
            ncirc0 = ncircles;
            ncircles += 150;
            if (top) ytop = YMAX;
            else ytop = ymean;
            
            if (top) ybottom = ymean;
            else ybottom = YMIN;
            
            for (n=0; n < 150; n++)
            {
                if (circles[ncirc0+n].yc + circles[ncirc0 + n].radius > ytop)
                {
                    circles[ncircles].xc = circles[ncirc0 + n].xc;
                    circles[ncircles].yc = circles[ncirc0+n].yc - height;
                    circles[ncircles].radius = MU;
                    circles[ncircles].active = 1;
                    ncircles ++;
                }
                else if (circles[ncirc0+n].yc - circles[ncirc0 + n].radius < ybottom)
                {
                    circles[ncircles].xc = circles[ncirc0 + n].xc;
                    circles[ncircles].yc = circles[ncirc0+n].yc + height;
                    circles[ncircles].radius = MU;
                    circles[ncircles].active = 1;
                    ncircles ++;
                }
            }
            
            break;
        }
        case (C_GOLDEN_SPIRAL):
        {
            circles[ncircles].xc = 0.0;
            circles[ncircles].yc = 0.0;
            if (top) circles[ncircles].radius = MU;
            else circles[ncircles].radius = MUB;
            circles[ncircles].active = 1;
            circletop[ncircles] = top; 
            
            gamma = (sqrt(5.0) - 1.0)*PI;    /* golden mean times 2Pi */
            phi = 0.0;
            r0 = 2.0*MU;
            r = r0 + MU;
            
            for (i=0; i<1000; i++) 
            {
                x = r*cos(phi);
                y = r*sin(phi);
                
                phi += gamma;
                r += MU*r0/r;
                
                if ((vabs(x) < LAMBDA)&&(vabs(y) < YMAX + MU))
                {
                    circles[ncircles].xc = x;
                    circles[ncircles].yc = y;
                    circles[ncircles].radius = MU;
                    if (((top)&&(circles[ncircles].yc < YMAX + MU)&&(circles[ncircles].yc > ymean - MU))
                        ||((!top)&&(circles[ncircles].yc < ymean + MU)&&(circles[ncircles].yc > YMIN - MU)))
                    {
                        circles[ncircles].active = 1;
                        circletop[ncircles] = top; 
                        ncircles++;
                    }
                }
            }
        break;
        }
        case (C_ONE):
        {
            circles[ncircles].xc = 0.0;
            circles[ncircles].yc = 0.0;
            if (top) circles[ncircles].radius = MU;
            else circles[ncircles].radius = MUB;
            circles[ncircles].active = 1;
            circletop[ncircles] = top; 
            ncircles += 1;
            break;
        }
        case (C_TWO):
        {
            circles[ncircles].xc = 0.0;
            circles[ncircles].yc = 0.0;
            if (top) circles[ncircles].radius = MU;
            else circles[ncircles].radius = MUB;
            circles[ncircles].active = 2;
            circletop[ncircles] = top; 
            ncircles += 1;

            circles[ncircles].xc = 0.0;
            circles[ncircles].yc = 0.0;
            if (top) circles[ncircles].radius = 2.0*MU;
            else circles[ncircles].radius = 2.0*MUB;
            circles[ncircles].active = 1;
            circletop[ncircles] = top; 
            ncircles += 1; 
            break;
        }
        case (C_NOTHING):
        {
            ncircles += 0;
            break;
        }
        default: 
        {
            printf("Function init_circle_config not defined for this pattern \n");
        }
    }
}


void init_circle_config_comp(t_circle circles[NMAXCIRCLES])
/* initialise the arrays circlex, circley, circlerad and circleactive */
/* for billiard shape D_CIRCLES */
{
    ncircles = 0;
    init_circle_config_half(CIRCLE_PATTERN, 1, circles);
    init_circle_config_half(CIRCLE_PATTERN_B, 0, circles);
}

void init_circle_config_energy(t_circle circles[NMAXCIRCLES])
/* initialise the arrays circlex, circley, circlerad and circleactive */
/* for billiard shape D_CIRCLES */
{
    ncircles = 0;
    init_circle_config_half(CIRCLE_PATTERN, 0, circles);
}


void init_polygon_config_half(int pattern, int top, int random_angle, int xdep_angle, t_polygon polygons[NMAXCIRCLES])
/* initialise the polygon configuration, for billiard shape D_CIRCLES */
/* uses init_circle_config, this is where C++ would be more elegant */
{
    int i;
    t_circle circle[NMAXCIRCLES];
    
//     ncircles = 0;
    init_circle_config_half(pattern, top, circle);
    for (i=0; i<NMAXCIRCLES; i++)
    {
        polygons[i].xc = circle[i].xc;
        polygons[i].yc = circle[i].yc;
        polygons[i].radius = circle[i].radius;
        polygons[i].active = circle[i].active;
        polygons[i].nsides = NPOLY;
        
        if ((top)&&(circletop[i])) 
        {
            if (random_angle) polygons[i].angle = DPI*(double)rand()/RAND_MAX;
            else if (xdep_angle) polygons[i].angle = APOLY + PID*POLY_ROTATION_ANGLE*polygons[i].xc;
            else polygons[i].angle = APOLY;
        }
        else if ((!top)&&(!circletop[i]))
        {
            if (random_angle) polygons[i].angle = DPI*(double)rand()/RAND_MAX;
            else if (xdep_angle) polygons[i].angle = APOLY_B + PID*POLY_ROTATION_ANGLE*polygons[i].xc;
            else polygons[i].angle = APOLY_B;            
        }
            
        if (i < ncircles) printf("(x,y) = (%.2f, %.2f), r = %.2f, angle = %.2f, sides = %i\n", polygons[i].xc, polygons[i].yc, polygons[i].radius, polygons[i].angle, polygons[i].nsides);
    }
}

void init_polygon_config_comp(t_polygon polygons[NMAXCIRCLES])
/* initialise polygon configuration for billiard shape D_POLYGONS */
{
    ncircles = 0;
    init_polygon_config_half(CIRCLE_PATTERN, 1, RANDOM_POLY_ANGLE, XDEP_POLY_ANGLE, polygons);
    init_polygon_config_half(CIRCLE_PATTERN_B, 0, RANDOM_POLY_ANGLE_B, XDEP_POLY_ANGLE_B, polygons);
}



int xy_in_billiard_half(double x, double y, int domain, int pattern, int top)
/* returns 1 if (x,y) represents a point in the billiard */
{
    double l2, r2, r2mu, omega, c, angle, z, x1, y1, x2, y2, u, v, u1, v1, dx, dy, width, a, b;
    int i, j, k, k1, k2, condition, type;

    switch (domain) {
        case (D_NOTHING):
        {
            return(1);
            break;
        }
        case (D_MENGER):       
        {
            x1 = 0.5*(x+1.0);
            y1 = 0.5*(y+1.0);
            for (k=0; k<MDEPTH; k++)
            {
                x1 = x1*(double)MRATIO;
                y1 = y1*(double)MRATIO;
                if ((top)&&(y > 0.0)&&(vabs(x)<1.0)&&(vabs(y)<1.0)&&(((int)x1 % MRATIO)==MRATIO/2)&&(((int)y1 % MRATIO)==MRATIO/2))                 return(0);
                else if ((!top)&&(y < 0.0)&&(vabs(x)<1.0)&&(vabs(y)<1.0)&&(((int)x1 % MRATIO)==MRATIO/2)&&(((int)y1 % MRATIO)==MRATIO/2)) return(0); 
            }
            return(1);
        }
        case (D_CIRCLES):
        {
            for (i = 0; i < ncircles; i++)
                if (circles[i].active != 0) 
                {
                    /* choose specific type according to value of circles[i].active */
                    if (circles[i].active == 1) type = 0;
                    else type = circles[i].active;
                    
                    x1 = circles[i].xc;
                    y1 = circles[i].yc;
                    r2 = circles[i].radius*circles[i].radius;
                    if ((top)&&(circletop[i])&&(y > 0.0)&&((x-x1)*(x-x1) + (y-y1)*(y-y1) < r2)) return(type); 
                    else if ((!top)&&(!circletop[i])&&(y < 0.0)&&((x-x1)*(x-x1) + (y-y1)*(y-y1) < r2)) return(type); 
                }
            return(1);
        }
        case (D_POLYGONS):
        {
            for (i = 0; i < ncircles; i++) 
                if (polygons[i].active)
                {
                    /* choose specific type according to value of circles[i].active */
                    if (polygons[i].active == 1) type = 0;
                    else type = polygons[i].active;
                    
                    if ((top)&&(circletop[i])&&(y > 0.0)&&(in_tpolygon(x, y, polygons[i]))) return(type); 
                    else if ((!top)&&(!circletop[i])&&(y < 0.0)&&(in_tpolygon(x, y, polygons[i]))) return(type); 
                }
            return(1);
        }
        case (D_FRESNEL):
        {
            if (vabs(y) > 0.9*vabs(LAMBDA)) return(1);
            if (vabs(x) > MU) return(1);
            
            x1 = sqrt(LAMBDA*LAMBDA - y*y) - vabs(LAMBDA);
            while (x1 <= 0.0) x1 += MU;
            if (LAMBDA > 0.0)
            {
                if (x < x1) return(0);
                else return(1);
            }
            else 
            {
                if (x > -x1) return(0);
                else return(1);
            }
        }
        case (D_CIRCLE_SEGMENT):
        {
            if (vabs(y) > 0.9*vabs(LAMBDA)) return(1);
            
            y1 = 0.9*LAMBDA;
            x1 = sqrt(LAMBDA*LAMBDA - y1*y1) - vabs(LAMBDA) + MU;
            if ((LAMBDA > 0.0)&&(x < x1)) return(1);
            else if ((LAMBDA < 0.0)&&(x > -x1)) return(1);
            
            x1 = sqrt(LAMBDA*LAMBDA - y*y) - vabs(LAMBDA) + MU;
            if (LAMBDA > 0.0)
            {
                if (x < x1) return(0);
                else return(1);
            }
            else
            {
                if (x > -x1) return(0);
                else return(1);                
            }
        }
        case (D_TWO_LENSES_WALL):
        {
            width = 0.2;
            a = 1.9;
            b = 1.9 - 2.0*LAMBDA + 2.0*width;
            x1 = vabs(x);
            if (module2(x1 - a, y) > LAMBDA) return(1);
            if (module2(x1 - b, y) > LAMBDA) return(1);
            return(0);
        }
        case (D_TWO_LENSES_OBSTACLE):
        {
            width = 0.2;
            a = 1.9;
            b = 1.9 - 2.0*LAMBDA + 2.0*width;
            x1 = vabs(x);
            if (module2(x1 - a, y) > LAMBDA) return(1);
            if (module2(x1 - b, y) > LAMBDA) return(1);
            return(0);
        }
        default:
        {
            printf("Function ij_in_billiard not defined for this billiard \n");
            return(0);
        }
    }
}


int xy_in_billiard_comp(double x, double y)
/* returns 1 if (x,y) represents a point in the billiard */
{
    if (y > 0) return(xy_in_billiard_half(x, y, B_DOMAIN, CIRCLE_PATTERN, 1));
    else return(xy_in_billiard_half(x, y, B_DOMAIN_B, CIRCLE_PATTERN_B, 0));
}


int ij_in_billiard_comp(int i, int j)
/* returns 1 if (i,j) represents a point in the billiard */
{
    double xy[2];

    ij_to_xy(i, j, xy);

    return(xy_in_billiard_comp(xy[0], xy[1]));
}


void draw_billiard_half(int domain, int pattern, int top, int fade, double fade_value)      
/* draws the billiard boundary */
/* two domain version implemented for D_CIRCLES */
{
    double x0, x, y, x1, y1, dx, dy, phi, r = 0.01, pos[2], pos1[2], alpha, dphi, omega, z, l, signtop, width, a, b, arcangle;
    int i, j, k, k1, k2, mr2;
    static int first = 1;
    
    glEnable(GL_SCISSOR_TEST);
    if (top) glScissor(0.0, YMID, NX, YMID);
//     if (top) glScissor(0.0, YMID, NX, YMAX);
    else glScissor(0.0, 0.0, NX, YMID);

    if (fade)
    {
        if (BLACK) glColor3f(fade_value, fade_value, fade_value);
        else glColor3f(1.0 - fade_value, 1.0 - fade_value, 1.0 - fade_value);        
    }
    else
    {
        if (BLACK) glColor3f(1.0, 1.0, 1.0);
        else glColor3f(0.0, 0.0, 0.0);
    }
    glLineWidth(5);
    
    if (top) signtop = 1.0;
    else signtop = -1.0;

    glEnable(GL_LINE_SMOOTH);

    switch (domain) {
        case (D_MENGER):
        {
            glLineWidth(3);
//             draw_rectangle(XMIN, -1.0, XMAX, 1.0);
            
            /* level 1 */
            if (MDEPTH > 0)
            {
                glLineWidth(2);
                x = 1.0/((double)MRATIO);
                draw_rectangle(-x, 0.0, x, signtop*x);
            }
            
            /* level 2 */
            if (MDEPTH > 1)
            {
                glLineWidth(1);
                mr2 = MRATIO*MRATIO;
                l = 2.0/((double)mr2);
                
                for (i=0; i<MRATIO; i++)
                    for (j=MRATIO/2; j<MRATIO; j++)
                        if ((i!=MRATIO/2)||(j!=MRATIO/2))
                        {
                            x = -1.0 - 0.5*l + 2.0*((double)i + 0.5)/((double)MRATIO);
                            y = -1.0 - 0.5*l + 2.0*((double)j + 0.5)/((double)MRATIO);
                            y1 = y+l;
                            if (y < 0.0) y = 0.0;
                            if (y1 < 0.0) y1 = 0.0;
                            draw_rectangle(x, signtop*y, x+l, signtop*y1);
                        }
            }
            
            /* level 3 */
            if (MDEPTH > 2)
            {
                glLineWidth(1);
                l = 2.0/((double)(mr2*MRATIO));
                
                for (i=0; i<mr2; i++)
                    for (j=mr2/2; j<mr2; j++)
                        if ( (((i%MRATIO!=MRATIO/2))||(j%MRATIO!=MRATIO/2)) && (((i/MRATIO!=MRATIO/2))||(j/MRATIO!=MRATIO/2)) )
                        {
                            x = -1.0 - 0.5*l + 2.0*((double)i + 0.5)/((double)mr2);
                            y = -1.0 - 0.5*l + 2.0*((double)j + 0.5)/((double)mr2);
                            y1 = y+l;
                            if (y < 0.0) y = 0.0;
                            if (y1 < 0.0) y1 = 0.0;
                            draw_rectangle(x, signtop*y, x+l, signtop*y1);
                        }
            }
            
            break;
        }        
        case (D_CIRCLES):
        {
            glLineWidth(2);
            for (i = 0; i < ncircles; i++) 
                if ((circles[i].active)&&(circletop[i] == top)) 
                {
                    glBegin(GL_LINE_STRIP);
                    for (k=0; k<=NSEG; k++)
                    {
                        phi = (double)k*DPI/(double)NSEG;
                        x = circles[i].xc + circles[i].radius*cos(phi);
                        y = circles[i].yc + circles[i].radius*sin(phi);
                        if ((top)&&(circletop[i])&&(y < 0.0)) y = 0.0;
                        if ((!top)&&(!circletop[i])&&(y > 0.0)) y = 0.0;
                        xy_to_pos(x, y, pos);
                        glVertex2d(pos[0], pos[1]);
                    }
                    glEnd ();
                }
            break;
        }
        case (D_POLYGONS):
        {
            glLineWidth(BOUNDARY_WIDTH);
            for (i = 0; i < ncircles; i++) 
                if ((polygons[i].active)&&(circletop[i] == top)) 
                {
                    if ((top)&&(circletop[i])) draw_tpolygon(polygons[i]);
                    else if ((!top)&&(!circletop[i])) draw_tpolygon(polygons[i]);
                }
            break;
        }
        case (D_FRESNEL):
        {
            glLineWidth(BOUNDARY_WIDTH);
            glBegin(GL_LINE_LOOP);
            for (i=0; i<npolyline; i++) tvertex_lineto(polyline[i]);
            glEnd();
            break;
        }
        case (D_CIRCLE_SEGMENT):
        {
            glLineWidth(BOUNDARY_WIDTH);
            glBegin(GL_LINE_LOOP);
            for (i=0; i<NSEG; i++) 
            {
                y = -0.9*LAMBDA + (double)i*1.8*LAMBDA/(double)NSEG;
                if (LAMBDA > 0.0) x = sqrt(LAMBDA*LAMBDA - y*y) - LAMBDA + MU;
                else x = -sqrt(LAMBDA*LAMBDA - y*y) - LAMBDA - MU;
                xy_to_pos(x, y, pos);
                glVertex2d(pos[0], pos[1]);
            }
            y = 0.9*LAMBDA;
            if (LAMBDA > 0.0) x = sqrt(LAMBDA*LAMBDA - y*y) - LAMBDA + MU;
            else x = -sqrt(LAMBDA*LAMBDA - y*y) - LAMBDA - MU;
            xy_to_pos(x, y, pos);
            glVertex2d(pos[0], pos[1]);
            glEnd();
            break;
        }
        case (D_TWO_LENSES_WALL):
        {
            width = 0.2;
            a = 1.9;
            b = 1.9 - 2.0*LAMBDA + 2.0*width;
            if (first)
            {
                arcangle = acos(1.0 - width/LAMBDA);
                first = 0;
            }
            draw_circle_arc(a, 0.0, LAMBDA, PI-arcangle, 2.0*arcangle, NSEG);
            draw_circle_arc(b, 0.0, LAMBDA, -arcangle, 2.0*arcangle, NSEG);
            draw_circle_arc(-a, 0.0, LAMBDA, -arcangle, 2.0*arcangle, NSEG);
            draw_circle_arc(-b, 0.0, LAMBDA, PI-arcangle, 2.0*arcangle, NSEG);
            
            width = 0.05;
            draw_rectangle(-width, MU, width, YMAX + 1.0);
            draw_rectangle(-width, YMIN-1.0, width, -MUB);
            break;
        }
        case (D_TWO_LENSES_OBSTACLE):
        {
            width = 0.2;
            a = 1.9;
            b = 1.9 - 2.0*LAMBDA + 2.0*width;
            if (first)
            {
                arcangle = acos(1.0 - width/LAMBDA);
                first = 0;
            }
            draw_circle_arc(a, 0.0, LAMBDA, PI-arcangle, 2.0*arcangle, NSEG);
            draw_circle_arc(b, 0.0, LAMBDA, -arcangle, 2.0*arcangle, NSEG);
            draw_circle_arc(-a, 0.0, LAMBDA, -arcangle, 2.0*arcangle, NSEG);
            draw_circle_arc(-b, 0.0, LAMBDA, PI-arcangle, 2.0*arcangle, NSEG);
            
            width = 0.05;
            draw_rectangle(-width, -MUB, width, MU);
            break;
        }
        case (D_MANDELBROT):
        {
            /* Do nothing */
            break;
        }
        case (D_MANDELBROT_CIRCLE):
        {
            /* Do nothing */
            break;
        }
        case (D_JULIA):
        {
            /* Do nothing */
            break;
        }
        case (D_NOTHING):
        {
            /* Do nothing */
            break;
        }   
        default:
        {
            printf("Function draw_billiard not defined for this billiard \n");
        }
    }
    
    glDisable(GL_SCISSOR_TEST);
}


void draw_billiard_comp(int fade, double fade_value)      /* draws the billiard boundary */
{
    draw_billiard_half(B_DOMAIN, CIRCLE_PATTERN, 1, fade, fade_value);
    draw_billiard_half(B_DOMAIN_B, CIRCLE_PATTERN_B, 0, fade, fade_value);
}


void int_planar_wave_comp(double x, double y, double *phi[NX], double *psi[NX], short int * xy_in[NX])
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
	    xy_in[i][j] = xy_in_billiard_comp(xy[0],xy[1]);
            
	    if ((xy_in[i][j])||(TWOSPEEDS)) phi[i][j] = 0.01*exp(-dist2/0.0005)*cos(-sqrt(dist2)/0.01);
            else phi[i][j] = 0.0;
            psi[i][j] = 0.0;
        }
}


void init_wave_flat_comp( double *phi[NX], double *psi[NX], short int * xy_in[NX])
/* initialise flat field - phi is wave height, psi is phi at time t-1 */
{
    int i, j;
    double xy[2], dist2;

    for (i=0; i<NX; i++)
        for (j=0; j<NY; j++)
        {
            ij_to_xy(i, j, xy);
	    xy_in[i][j] = xy_in_billiard_comp(xy[0],xy[1]);
	    phi[i][j] = 0.0;
            psi[i][j] = 0.0;
        }
}


void add_circular_wave_comp(double factor, double x, double y, double *phi[NX], double *psi[NX], short int * xy_in[NX], int half)
/* add drop at (x,y) to the field with given prefactor */
{
    int i, j;
    double xy[2], dist2;

    for (i=0; i<NX; i++)
    {
        if (half == 1) for (j=NY/2; j<NY; j++)
        {
            ij_to_xy(i, j, xy);
            dist2 = (xy[0]-x)*(xy[0]-x) + (xy[1]-y)*(xy[1]-y);
            if ((xy_in[i][j])||(TWOSPEEDS)) 
                phi[i][j] += INITIAL_AMP*factor*exp(-dist2/INITIAL_VARIANCE)*cos(-sqrt(dist2)/INITIAL_WAVELENGTH);
        }
        else for (j=0; j<NY/2; j++)
        {
            ij_to_xy(i, j, xy);
            dist2 = (xy[0]-x)*(xy[0]-x) + (xy[1]-y)*(xy[1]-y);
            if ((xy_in[i][j])||(TWOSPEEDS)) 
                phi[i][j] += INITIAL_AMP*factor*exp(-dist2/INITIAL_VARIANCE)*cos(-sqrt(dist2)/INITIAL_WAVELENGTH);
        }
    }
}

void draw_wave_comp(double *phi[NX], double *psi[NX], short int *xy_in[NX], double scale, int time, int plot)
/* draw the field */
{
    int i, j, iplus, iminus, jplus, jminus;
    double rgb[3], xy[2], x1, y1, x2, y2, velocity, energy, gradientx2, gradienty2, pos[2];
    static double dtinverse = ((double)NX)/(COURANT*(XMAX-XMIN)), dx = (XMAX-XMIN)/((double)NX);

    glBegin(GL_QUADS);
    
//     printf("dtinverse = %.5lg\n", dtinverse);

    for (i=0; i<NX; i++)
        for (j=0; j<NY; j++)
        {
            if (((TWOSPEEDS)&&(xy_in[i][j] != 2))||(xy_in[i][j] == 1))
            {
                switch (plot) {
                    case (P_AMPLITUDE):
                    {
                        color_scheme(COLOR_SCHEME, phi[i][j], scale, time, rgb);
                        break;
                    }
                    case (P_ENERGY):
                    {
                        energy = compute_energy(phi, psi, xy_in, i, j);
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
//                     case (P_MEAN_ENERGY):
//                     {
//                         energy = compute_energy(phi, psi, xy_in, i, j);
//                         total_energy[i][j] += energy;
//                         if (COLOR_PALETTE >= COL_TURBO) 
//                             color_scheme_asym(COLOR_SCHEME, total_energy[i][j]/(double)(time+1), scale, time, rgb);
//                         else color_scheme(COLOR_SCHEME, total_energy[i][j]/(double)(time+1), scale, time, rgb);
//                         break;
//                     }
                    case (P_LOG_ENERGY):
                    {
                        energy = compute_energy(phi, psi, xy_in, i, j);
                        color_scheme(COLOR_SCHEME, LOG_SHIFT + LOG_SCALE*log(energy), scale, time, rgb);
                        break;
                    }
//                     case (P_LOG_MEAN_ENERGY):
//                     {
//                         energy = compute_energy(phi, psi, xy_in, i, j);
//                         if (energy == 0.0) energy = 1.0e-20;
//                         total_energy[i][j] += energy;
//                         color_scheme(COLOR_SCHEME, LOG_SCALE*log(total_energy[i][j]/(double)(time+1)), scale, time, rgb);
//                         break;
//                     }
                }
                    
//                 if (PLOT == P_AMPLITUDE)
//                     color_scheme(COLOR_SCHEME, phi[i][j], scale, time, rgb);
//                 else if (PLOT == P_ENERGY)
//                 {
//                     energy = compute_energy(phi, psi, xy_in, i, j);
//                     if (COLOR_PALETTE >= COL_TURBO) color_scheme_asym(COLOR_SCHEME, energy, scale, time, rgb);
//                     else color_scheme(COLOR_SCHEME, energy, scale, time, rgb);
//                 }
//                 else if (PLOT == P_MIXED)
//                 {
//                     if (j > NY/2) color_scheme(COLOR_SCHEME, phi[i][j], scale, time, rgb);
//                     else color_scheme(COLOR_SCHEME, compute_energy(phi, psi, xy_in, i, j), scale, time, rgb);
//                 }
                glColor3f(rgb[0], rgb[1], rgb[2]);

                glVertex2i(i, j);
                glVertex2i(i+1, j);
                glVertex2i(i+1, j+1);
                glVertex2i(i, j+1);
            }
        }

    glEnd ();
    
    /* draw horizontal mid line */
    glColor3f(1.0, 1.0, 1.0);
    glBegin(GL_LINE_STRIP);
    xy_to_pos(XMIN, 0.5*(YMIN+YMAX), pos);
    glVertex2d(pos[0], pos[1]);
    xy_to_pos(XMAX, 0.5*(YMIN+YMAX), pos);
    glVertex2d(pos[0], pos[1]);
    glEnd();
}

void draw_wave_comp_highres_palette(int size, double *phi[NX], double *psi[NX], double *total_energy[NX], short int *xy_in[NX], double scale, int time, int plot, int palette, int fade, double fade_value)
/* draw the field */
{
    int i, j, iplus, iminus, jplus, jminus, k;
    double rgb[3], xy[2], x1, y1, x2, y2, velocity, energy, gradientx2, gradienty2, pos[2];
    static double dtinverse = ((double)NX)/(COURANT*(XMAX-XMIN)), dx = (XMAX-XMIN)/((double)NX);

    glBegin(GL_QUADS);
    
//     printf("dtinverse = %.5lg\n", dtinverse);

    for (i=0; i<NX; i++)
        for (j=0; j<NY; j++)
        {
            if (((TWOSPEEDS)&&(xy_in[i][j] != 2))||(xy_in[i][j] == 1))
            {
                switch (plot) {
                    case (P_AMPLITUDE):
                    {
                        color_scheme_palette(COLOR_SCHEME, palette, phi[i][j], scale, time, rgb);
                        break;
                    }
                    case (P_ENERGY):
                    {
                        energy = compute_energy(phi, psi, xy_in, i, j);
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
                        color_scheme(COLOR_SCHEME, LOG_SHIFT + LOG_SCALE*log(energy), scale, time, rgb);
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
                }
                    
//                 if (PLOT == P_AMPLITUDE)
//                     color_scheme(COLOR_SCHEME, phi[i][j], scale, time, rgb);
//                 else if (PLOT == P_ENERGY)
//                 {
//                     energy = compute_energy(phi, psi, xy_in, i, j);
//                     if (COLOR_PALETTE >= COL_TURBO) color_scheme_asym(COLOR_SCHEME, energy, scale, time, rgb);
//                     else color_scheme(COLOR_SCHEME, energy, scale, time, rgb);
//                 }
//                 else if (PLOT == P_MIXED)
//                 {
//                     if (j > NY/2) color_scheme(COLOR_SCHEME, phi[i][j], scale, time, rgb);
//                     else color_scheme(COLOR_SCHEME, compute_energy(phi, psi, xy_in, i, j), scale, time, rgb);
//                 }
                if (fade) for (k=0; k<3; k++) rgb[k] *= fade_value;
                glColor3f(rgb[0], rgb[1], rgb[2]);
                
                glVertex2i(i, j);
                glVertex2i(i+size, j);
                glVertex2i(i+size, j+size);
                glVertex2i(i, j+size);
            }
        }

    glEnd ();
    
    /* draw horizontal mid line */
    if (fade) glColor3f(fade_value, fade_value, fade_value);
    else glColor3f(1.0, 1.0, 1.0);
    glBegin(GL_LINE_STRIP);
    xy_to_pos(XMIN, 0.5*(YMIN+YMAX), pos);
    glVertex2d(pos[0], pos[1]);
    xy_to_pos(XMAX, 0.5*(YMIN+YMAX), pos);
    glVertex2d(pos[0], pos[1]);
    glEnd();
}

void compute_energy_tblr(double *phi[NX], double *psi[NX], short int *xy_in[NX], double energies[6])
/* compute total energy in top/bottom left/right boxes */ 
{
    int i, j, ij[2];
    double energy = 0.0, rescale, pos, xleft = XMAX, xright = XMIN, emax = 1.2, xc, r;
    double *energy_ij[NX];
    static short int first = 1;
    static int ileft, iright, jmid = NY/2; 
    static double sqremax;
    
    for (i=0; i<NX; i++) energy_ij[i] = (double *)malloc(NY*sizeof(double));
        
    
    if (first) /* compute box limits */
    {
        /* find leftmost and rightmost circle */
        for (i=0; i<ncircles; i++) 
        {
            if ((circles[i].active)&&(circles[i].xc - circles[i].radius < xleft)) xleft = circles[i].xc - circles[i].radius;
            if ((polygons[i].active)&&(polygons[i].xc - polygons[i].radius < xleft)) xleft = polygons[i].xc - polygons[i].radius;
            if ((circles[i].active)&&(circles[i].xc + circles[i].radius > xright)) xright = circles[i].xc + circles[i].radius; 
            if ((polygons[i].active)&&(polygons[i].xc + polygons[i].radius > xright)) xright = polygons[i].xc +         polygons[i].radius; 
        }
//         for (i=0; i<ncircles; i++) 
//             if ((circles[i].active)&&(circles[i].xc - circles[i].radius < xleft)) xleft = circles[i].xc - circles[i].radius; 
//         for (i=0; i<ncircles; i++) 
//             if ((circles[i].active)&&(circles[i].xc + circles[i].radius > xright)) xright = circles[i].xc + circles[i].radius; 
        
        xy_to_ij(xleft, 0.0, ij);
        ileft = ij[0];
        xy_to_ij(xright, 0.0, ij);
        iright = ij[0];
        first = 0;
        
        printf("xleft = %.3lg, xright = %.3lg", xleft, xright);
    }
    
    for (i=0; i<NX; i++)
        for (j=0; j<NY; j++)
            energy_ij[i][j] = compute_energy(phi, psi, xy_in, i, j);
    
    /* prevent local energy from growing too large */
    if (FLOOR)
    {
        for (i=10; i<NX; i++)
            for (j=0; j<NY; j++)
                if ((xy_in[i][j])&&(energy_ij[i][j] > emax))
                {
                    rescale = sqrt(emax/energy_ij[i][j]);
                    if (j%100 == 0) printf("Rescaling at (%i,%i) by %.5lg\n", i, j, rescale);
                    phi[i][j] = phi[i][j]*rescale;
                    psi[i][j] = psi[i][j]*rescale;
                }
                else if (energy_ij[i][j] > 0.1*emax)
                {
                    rescale = sqrt(0.1*emax/energy_ij[i][j]);
                    if (j%10 == 0) printf("Rescaling at (%i,%i) by %.5lg\n", i, j, rescale);
                    phi[i][j] = phi[i][j]*rescale;
                    psi[i][j] = psi[i][j]*rescale;
                }
    }
    
    /* top left box */
    for (i=0; i<ileft; i++)
        for (j=jmid; j<NY; j++)
            energy += energy_ij[i][j];
    energies[0] = energy;
    
    /* top middle box */
    energy = 0.0;
    for (i=ileft; i<iright; i++)
        for (j=jmid; j<NY; j++)
            energy += energy_ij[i][j];
    energies[1] = energy;
    
    /* top right box */
    energy = 0.0;
    for (i=iright; i<NX; i++)
        for (j=jmid; j<NY; j++)
            energy += energy_ij[i][j];
    energies[2] = energy;
    
    /* bottom left box */
    energy = 0.0;
    for (i=0; i<ileft; i++)
        for (j=0; j<jmid; j++)
            energy += energy_ij[i][j];
    energies[3] = energy;
    
    /* bottom middle box */
    energy = 0.0;
    for (i=ileft; i<iright; i++)
        for (j=0; j<jmid; j++)
            energy += energy_ij[i][j];
    energies[4] = energy;
    
    /* bottom right box */
    energy = 0.0;
    for (i=iright; i<NX; i++)
        for (j=0; j<jmid; j++)
            energy += energy_ij[i][j];
    energies[5] = energy;
    
//     printf("Energies: %.5lg, %.5lg, %.5lg\n %.5lg, %.5lg, %.5lg\n", energies[0], energies[1], energies[2], energies[3], energies[4], energies[5]);
    
    for (i=0; i<NX; i++) free(energy_ij[i]);
}

void print_energies(double energies[6], double top_energy, double bottom_energy)
{
    char message[50];
    double ytop, ybot, pos[2], centerx = -0.075, boxxright = XMAX - 0.17, textxright = XMAX - 0.28, boxwidth = 0.1, 
            boxheight = 0.05, leftboxshift = 0.185, centerboxshift = 0.085, text_color;

    if (BLACK_TEXT) text_color = 0.0;
    else text_color = 1.0;
    
    /* adapt sizes of text areas to high resolution */
    if (WINWIDTH > 1280)
    {
        boxheight = 0.035;
//         centerx -= 0.2;
        boxwidth = 0.08;
        leftboxshift = 0.16;
        centerboxshift = 0.06;
    }
    
    ytop = YMAX - 0.1;
    ybot = YMIN + 0.05;
    if (DRAW_COLOR_SCHEME) 
    {   
        boxxright -= 0.35;
        textxright -= 0.35;
    }
    
    erase_area(XMIN + leftboxshift, ytop + 0.025, boxwidth, boxheight);

    glColor3f(text_color, text_color, text_color);
    sprintf(message, "%.3f", energies[0]/top_energy);
    xy_to_pos(XMIN + 0.1, ytop, pos);
    write_text(pos[0], pos[1], message);
    
    erase_area(centerx + centerboxshift, ytop + 0.025, boxwidth, boxheight);

    glColor3f(text_color, text_color, text_color);
    sprintf(message, "%.3f", energies[1]/top_energy);
    xy_to_pos(centerx, ytop, pos);
    write_text(pos[0], pos[1], message);
    
    erase_area(boxxright, ytop + 0.025, boxwidth + 0.05, boxheight);

    glColor3f(text_color, text_color, text_color);
    sprintf(message, "%.5f", energies[2]/top_energy);
    xy_to_pos(textxright, ytop, pos);
    write_text(pos[0], pos[1], message);
    
    erase_area(XMIN + leftboxshift, ybot + 0.025, boxwidth, boxheight);

    glColor3f(text_color, text_color, text_color);
    sprintf(message, "%.3f", energies[3]/bottom_energy);
    xy_to_pos(XMIN + 0.1, ybot, pos);
    write_text(pos[0], pos[1], message);

    erase_area(centerx + centerboxshift, ybot + 0.025, boxwidth, boxheight);

    glColor3f(text_color, text_color, text_color);
    sprintf(message, "%.3f", energies[4]/bottom_energy);
    xy_to_pos(centerx, ybot, pos);
    write_text(pos[0], pos[1], message);

    erase_area(boxxright, ybot + 0.025, boxwidth + 0.05, boxheight);

    glColor3f(text_color, text_color, text_color);
    sprintf(message, "%.5f", energies[5]/bottom_energy);
    xy_to_pos(textxright, ybot, pos);
    write_text(pos[0], pos[1], message);
    
}

void init_ior_2d_comp(short int *xy_in[NX], double *tcc_table[NX], double ior_angle)
/* compute variable index of refraction */
/* should be at some point merged with 3D version in suv_wave_3d.c */
{
    int i, j, k, n, inlens, ncircles;
    double courant2 = COURANT*COURANT, courantb2 = COURANTB*COURANTB, lambda1, mu1;
    double u, v, u1, x, y, xy[2], norm2, speed, r2, c, salpha, h, ll, ca, sa, x1, y1, dx, dy, sum, sigma, x0, y0, rgb[3];
    double xc[NGRIDX*NGRIDY], yc[NGRIDX*NGRIDY], height[NGRIDX*NGRIDY];
    static double xc_stat[NGRIDX*NGRIDY], yc_stat[NGRIDX*NGRIDY], sigma_stat;
    static int first = 1;
    t_circle circles[NMAXCIRCLES];
    
    rgb[0] = 1.0;
    rgb[1] = 1.0;
    rgb[2] = 1.0;
    
    if (VARIABLE_IOR)
    {
        switch (IOR) {
            case (IOR_LENS_WALL):
            {
                printf("Initializing IOR_LENS_WALL\n");
                for (i=0; i<NX; i++){
                    for (j=0; j<NY/2; j++){
                        ij_to_xy(i, j, xy);
                        x = xy[0];
                        y = xy[1];
                        if ((vabs(x) < 0.05)&&(vabs(y) > MUB)) tcc_table[i][j] = 0.0;
                        else 
                        {
                            if (xy_in[i][j] != 0) tcc_table[i][j] = courant2;
                            else tcc_table[i][j] = courantb2;
                        }
                    }
                    for (j=NY/2; j<NY; j++){
                        ij_to_xy(i, j, xy);
                        x = xy[0];
                        y = xy[1];
                        if ((vabs(x) < 0.05)&&(vabs(y) > MU)) tcc_table[i][j] = 0.0;
                        else 
                        {
                            if (xy_in[i][j] != 0) tcc_table[i][j] = courant2;
                            else tcc_table[i][j] = courantb2;
                        }
                    }
                }
                
                break;
            }
            case (IOR_LENS_OBSTACLE):
            {
                printf("Initializing IOR_LENS_OBSTACLE\n");
                for (i=0; i<NX; i++){
                    for (j=0; j<NY/2; j++){
                        ij_to_xy(i, j, xy);
                        x = xy[0];
                        y = xy[1];
                        if ((vabs(x) < 0.05)&&(vabs(y) < MUB)) tcc_table[i][j] = 0.0;
                        else 
                        {
                            if (xy_in[i][j] != 0) tcc_table[i][j] = courant2;
                            else tcc_table[i][j] = courantb2;
                        }
                    }
                    for (j=NY/2; j<NY; j++){
                        ij_to_xy(i, j, xy);
                        x = xy[0];
                        y = xy[1];
                        if ((vabs(x) < 0.05)&&(vabs(y) < MU)) tcc_table[i][j] = 0.0;
                        else 
                        {
                            if (xy_in[i][j] != 0) tcc_table[i][j] = courant2;
                            else tcc_table[i][j] = courantb2;
                        }
                    }
                }
                
                break;
            }
            default:
            {
                for (i=0; i<NX; i++){
                    for (j=0; j<NY; j++){
                        tcc_table[i][j] = COURANT;
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
                    tcc_table[i][j] = courant2;
                }
                else if (TWOSPEEDS)
                {
                    tcc_table[i][j] = courantb2;
                }
            }
        }
    }
}
