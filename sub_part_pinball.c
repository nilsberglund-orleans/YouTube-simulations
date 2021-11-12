/* global variables needed for circle configuration with periodic boundary conditions */

// short int double_circle[NMAXCIRCLES];   /* set to 1 if a circle is a translate of another one on the boundary */
// int partner_circle[NMAXCIRCLES];        /* number of circle of which current circle is a copy */    



void init_circles_pinball(int circle_pattern, t_circle circles[NMAXCIRCLES])
{
    int i, j, k, n, ncirc0, n_p_active, ncandidates=5000, naccepted; 
    double dx, dy, xx[4], yy[4], x, y, xk, yk, hy, gamma, height, phi, r0, r, dpoisson = 3.25*MU;
    short int active_poisson[NMAXCIRCLES], far;
    
    switch (circle_pattern) {
        case (C_FOUR_CIRCLES):
        {
            ncircles = 4;
            
            circles[0].xc = 1.0;
            circles[0].yc = 0.0;
            circles[0].radius = 0.8;
                        
            circles[1].xc = -1.0;
            circles[1].yc = 0.0;
            circles[1].radius = 0.8;
            
            circles[2].xc = 0.0;
            circles[2].yc = 0.8;
            circles[2].radius = 0.4;
            
            circles[3].xc = 0.0;
            circles[3].yc = -0.8;
            circles[3].radius = 0.4;
                        
            for (i=0; i<4; i++) 
            {   
                circles[i].active = 1;
                circles[i].color = 0;
                circles[i].new = 0;
            }

            break;
        }
        case (C_SQUARE):
        {
            ncircles = NCX*NCY;
            dy = (BOXYMAX - BOXYMIN)/((double)NCY);
//             dy = (YMAX - YMIN)/((double)NCY);
            for (i = 0; i < NCX; i++)
                for (j = 0; j < NCY; j++)
                {
                    n = NCY*i + j;
                    circles[n].xc = ((double)(i-NCX/2) + 0.5)*dy;
                    circles[n].yc = BOXYMIN + ((double)j + 0.5)*dy;
                    circles[n].radius = MU;
                    circles[n].active = 1;
                    circles[n].color = 0;
                    circles[n].new = 0;
                }
            break;
        }
        case (C_HEX):
        {
            ncircles = NCX*(NCY+1);
            dy = (YMAX - YMIN)/((double)NCY);
            dx = dy*0.5*sqrt(3.0);
            for (i = 0; i < NCX; i++)
                for (j = 0; j < NCY+1; j++)
                {
                    n = (NCY+1)*i + j;
//                     circles[n].xc = ((double)(i-NCX/2) + 0.5)*dy;
                    circles[n].xc = ((double)(i-NCX/2))*dy;
                    if (NCX % 2 == 0) circles[n].xc += 0.5*dy;
                    circles[n].yc = YMIN + ((double)j - 0.5)*dy;
                    if ((i+NCX)%2 == 1) circles[n].yc += 0.5*dy;
                    circles[n].radius = MU;
                    circles[n].active = 1;
                    circles[n].color = 0;
                    circles[n].new = 0;
                }
            break;
        }
        case (C_TRI):
        {
            ncircles = NCX*(NCY+1);
            dy = (BOXYMAX - BOXYMIN)/((double)NCY);
            dx = dy*0.5*sqrt(3.0);
            for (i = 0; i < NCX; i++)
                for (j = 0; j < NCY+1; j++)
                {
                    n = (NCY+1)*i + j;
                    circles[n].xc = ((double)(i-NCX/2) + 0.5)*dx;
                    circles[n].yc = BOXYMIN + ((double)j)*dy;
                    if ((i+NCX)%2 == 1) circles[n].yc += 0.5*dy;
                    circles[n].radius = MU;
                    circles[n].active = 1;
                    circles[n].color = 0;
                    circles[n].new = 0;
                    
                    /* take care of periodic boundary conditions */
                    if (B_DOMAIN == D_CIRCLES_IN_TORUS)
                    {
                        if ((j == NCY)&&((i+NCX)%2 == 0)) 
                        {
                            circles[n].double_circle = 1;
                            circles[n].partner = (NCY+1)*i;
                        }
                        else 
                        {
                            circles[n].double_circle = 0;
                            if ((j == 0)&&((i+NCX)%2 == 0)) circles[n].partner = (NCY+1)*i + NCY;
                            else circles[n].partner = n;
                        }
                    }
                    else circles[n].double_circle = 0;
                }
            break;
        }
        case (C_GOLDEN_MEAN):
        {
            ncircles = NCX*NCY;
            gamma = (sqrt(5.0) - 1.0)*0.5;    /* golden mean */
            height = YMAX - YMIN;
            dx = 2.0*LAMBDA/((double)ncircles);
            for (n = 0; n < ncircles; n++)
            {
                circles[n].xc = -LAMBDA + n*dx;
                circles[n].yc = y;
                y += height*gamma; 
                if (y > YMAX) y -= height;
                circles[n].radius = MU;
                circles[n].active = 1;
                circles[n].color = 0;
                circles[n].new = 0;
            }
            break;
        }
        case (C_GOLDEN_SPIRAL):
        {
            ncircles = 1;
            circles[0].xc = 0.0;
            circles[0].yc = 0.0;
            
            gamma = (sqrt(5.0) - 1.0)*PI;    /* golden mean times 2Pi */
            phi = 0.0;
            r0 = 2.0*MU;
            r = r0 + MU;
            
            for (i=0; i<NGOLDENSPIRAL; i++) 
            {
                x = r*cos(phi);
                y = r*sin(phi);
                
                phi += gamma;
                r += MU*r0/r;
                
                if ((vabs(x) < LAMBDA)&&(vabs(y) < YMAX + MU))
                {
                    circles[ncircles].xc = x;
                    circles[ncircles].yc = y;
                    ncircles++;
                }
            }
            
            for (i=0; i<ncircles; i++)
            {
                circles[i].radius = MU;
                circles[i].color = 0;
                circles[i].new = 0;
                /* inactivate circles outside the domain */
                if ((circles[i].yc < YMAX + MU)&&(circles[i].yc > YMIN - MU)) circles[i].active = 1;
            }
            break;
        }
        case (C_RAND_DISPLACED):
        {
            ncircles = NCX*NCY;
            dy = (YMAX - YMIN)/((double)NCY);
            for (i = 0; i < NCX; i++)
                for (j = 0; j < NCY; j++)
                {
                    n = NCY*i + j;
                    circles[n].xc = ((double)(i-NCX/2) + 0.5*((double)rand()/RAND_MAX - 0.5))*dy;
                    circles[n].yc = YMIN + ((double)j + 0.5 + 0.5*((double)rand()/RAND_MAX - 0.5))*dy;
                    circles[n].radius = MU*sqrt(1.0 + 0.8*((double)rand()/RAND_MAX - 0.5));
                    circles[n].active = 1;
                    circles[n].color = 0;
                    circles[n].new = 0;
                }
            break;
        }
        case (C_RAND_POISSON):
        {
            ncircles = NPOISSON;
            for (n = 0; n < NPOISSON; n++)
            {
                circles[n].xc = LAMBDA*(2.0*(double)rand()/RAND_MAX - 1.0);
                circles[n].yc = (BOXYMAX - BOXYMIN)*(double)rand()/RAND_MAX + BOXYMIN;
                circles[n].radius = MU;
                circles[n].active = 1;
                circles[n].color = 0;
                circles[n].new = 0;
                circles[n].double_circle = 0;
                circles[n].partner = n;
                
                /* take care of periodic boundary conditions */
                if (B_DOMAIN == D_CIRCLES_IN_TORUS)
                {
                    /* inactivate circles in corners */
                    if ((vabs(circles[n].xc) > LAMBDA - MU)&&((circles[n].yc < BOXYMIN + MU)||(circles[n].yc > BOXYMAX - MU)))
                        circles[n].active = 0;
                    
                    if (circles[n].xc < - LAMBDA + MU) 
                    {
                        circles[ncircles].xc = circles[n].xc + 2.0*LAMBDA;
                        circles[ncircles].yc = circles[n].yc;
                        circles[ncircles].partner = n;
                        circles[n].partner = ncircles;
                        ncircles++;
                    }
                    else if (circles[n].xc > LAMBDA - MU)
                    {
                        circles[ncircles].xc = circles[n].xc - 2.0*LAMBDA;
                        circles[ncircles].yc = circles[n].yc;
                        circles[ncircles].partner = n;
                        circles[n].partner = ncircles;
                        ncircles++;                        
                    }
                    if (circles[n].yc < BOXYMIN + MU) 
                    {
                        circles[ncircles].xc = circles[n].xc;
                        circles[ncircles].yc = circles[n].yc + BOXYMAX - BOXYMIN;
                        circles[ncircles].partner = n;
                        circles[n].partner = ncircles;
                        ncircles++;
                    }
                    else if (circles[n].yc > BOXYMAX - MU) 
                    {
                        circles[ncircles].xc = circles[n].xc;
                        circles[ncircles].yc = circles[n].yc - BOXYMAX + BOXYMIN;
                        circles[ncircles].partner = n;
                        circles[n].partner = ncircles;
                        ncircles++;
                    }
                }
            }
            printf("%i circles\n", ncircles);
            if (B_DOMAIN == D_CIRCLES_IN_TORUS) for (n = NPOISSON; n < ncircles; n++)
            {
//                 printf("circle %i at (%.3f, %.3f)\n", n, circles[n].xc, circles[n].yc); 
                circles[n].radius = MU;
                if (circles[circles[n].partner].active) circles[n].active = 1;
                circles[n].color = 0;
                circles[n].new = 0;
                circles[n].double_circle = 1;
            }

            break;
        }
        case (C_POISSON_DISC):
        {
            printf("Generating Poisson disc sample\n");
            /* generate first circle */
            circles[0].xc = LAMBDA*(2.0*(double)rand()/RAND_MAX - 1.0);
            circles[0].yc = (BOXYMAX - BOXYMIN)*(double)rand()/RAND_MAX + BOXYMIN;
            active_poisson[0] = 1;
            circles[0].double_circle = 0;
            circles[0].partner = 0;
            n_p_active = 1;
            ncircles = 1;
            
            while ((n_p_active > 0)&&(ncircles < NMAXCIRCLES))
            {
                /* randomly select an active circle */
                i = rand()%(ncircles);
                while (!active_poisson[i]) i = rand()%(ncircles);                 
//                 printf("Starting from circle %i at (%.3f,%.3f)\n", i, circlex[i], circley[i]);
                /* generate new candidates */
                naccepted = 0;
                for (j=0; j<ncandidates; j++)
                {
                    r = dpoisson*(2.0*(double)rand()/RAND_MAX + 1.0);
                    phi = DPI*(double)rand()/RAND_MAX;
                    x = circles[i].xc + r*cos(phi);
                    y = circles[i].yc + r*sin(phi);
//                        printf("Testing new circle at (%.3f,%.3f)\t", x, y);
                    far = 1;
                    for (k=0; k<ncircles; k++) if ((k!=i))
                    {
                        xk = circles[k].xc; 
                        yk = circles[k].yc;
                        hy = BOXYMAX - BOXYMIN;
                        
                        /* new circle is far away from circle k */
                        far = far*((x - xk)*(x - xk) + (y - yk)*(y - yk) >= dpoisson*dpoisson);
                        far = far*((x - xk)*(x - xk) + (y - yk + hy)*(y - yk + hy) >= dpoisson*dpoisson);
                        far = far*((x - xk)*(x - xk) + (y - yk - hy)*(y - yk - hy) >= dpoisson*dpoisson);
                        far = far*((x - xk + 2.0*LAMBDA)*(x - xk + 2.0*LAMBDA) + (y - yk)*(y - yk) >= dpoisson*dpoisson);
                        far = far*((x - xk - 2.0*LAMBDA)*(x - xk - 2.0*LAMBDA) + (y - yk)*(y - yk) >= dpoisson*dpoisson);                        
                    }                           
                    /* new circle is in domain */
                    far = far*(vabs(x) < LAMBDA + 0.0)*(y < BOXYMAX + 0.0)*(y > BOXYMIN - 0.0);
                    
                    /* exclude circles in corners */
                    if ((x > LAMBDA - MU)&&((y < BOXYMIN + MU)||(y > BOXYMAX - MU)))
                        far = 0;
                    
                    if (far)    /* accept new circle */
                    {
                        printf("New circle at (%.3f,%.3f) accepted\n", x, y);
                        circles[ncircles].xc = x;
                        circles[ncircles].yc = y;
                        circles[ncircles].radius = MU;
                        circles[ncircles].active = 1;
                        circles[ncircles].color = 0;
                        circles[ncircles].new = 0;
                        active_poisson[ncircles] = 1;
                        circles[ncircles].double_circle = 0;
                        circles[ncircles].partner = ncircles;
                        ncircles++;
                        n_p_active++;
                        naccepted++;
                        
                        /* take care of periodic boundary conditions */
                        if (B_DOMAIN == D_CIRCLES_IN_TORUS) 
                        {
                            n = ncircles - 1;
                            if (x < - LAMBDA + MU)
                            {
                                circles[ncircles].xc = x + 2.0*LAMBDA;
                                circles[ncircles].yc = y;
                                circles[ncircles].radius = MU;
                                circles[ncircles].active = 1;
                                circles[ncircles].color = 0;
                                circles[ncircles].new = 0;
                                active_poisson[ncircles] = 1;
                                circles[n].double_circle = 1;
                                circles[ncircles].double_circle = 0;
                                circles[n].partner = ncircles;
                                circles[ncircles].partner = n;
                                ncircles++;
                                n_p_active++;
                                naccepted++;
                            }
                            else if (x > LAMBDA - MU)
                            {
                                circles[ncircles].xc = x - 2.0*LAMBDA;
                                circles[ncircles].yc = y;
                                circles[ncircles].radius = MU;
                                circles[ncircles].active = 1;
                                circles[ncircles].color = 0;
                                circles[ncircles].new = 0;
                                active_poisson[ncircles] = 1;
                                circles[n].double_circle = 1;
                                circles[ncircles].double_circle = 0;
                                circles[n].partner = ncircles;
                                circles[ncircles].partner = n;
                                ncircles++;
                                n_p_active++;
                                naccepted++;
                            }
                            if (y < BOXYMIN + MU)
                            {
                                circles[ncircles].xc = x;
                                circles[ncircles].yc = y + BOXYMAX - BOXYMIN;
                                circles[ncircles].radius = MU;
                                circles[ncircles].active = 1;
                                circles[ncircles].color = 0;
                                circles[ncircles].new = 0;
                                active_poisson[ncircles] = 1;
                                circles[n].double_circle = 1;
                                circles[ncircles].double_circle = 0;
                                circles[n].partner = ncircles;
                                circles[ncircles].partner = n;
                                ncircles++;
                                n_p_active++;
                                naccepted++;
                            }
                            else if (y > BOXYMAX - MU)
                            {
                                circles[ncircles].xc = x;
                                circles[ncircles].yc = y - BOXYMAX + BOXYMIN;
                                circles[ncircles].radius = MU;
                                circles[ncircles].active = 1;
                                circles[ncircles].color = 0;
                                circles[ncircles].new = 0;
                                active_poisson[ncircles] = 1;
                                circles[n].double_circle = 1;
                                circles[ncircles].double_circle = 0;
                                circles[n].partner = ncircles;
                                circles[ncircles].partner = n;
                                ncircles++;
                                n_p_active++;
                                naccepted++;
                            }
                        }
                    }
//                        else printf("Rejected\n");
                }
                if (naccepted == 0)    /* inactivate circle i */ 
                {
//                     printf("No candidates work, inactivate circle %i\n", i);
                    active_poisson[i] = 0;
                    n_p_active--;
                }
                printf("%i active circles\n", n_p_active);
            }
            
            printf("Generated %i circles\n", ncircles);
            
//             for (i=0; i<ncircles; i++)
//                 printf("Circle %i at (%.3f, %.3f), partner %i\n", i, circlex[i], circley[i], partner_circle[i]);
            break;
        }
//         case (C_LASER):
//         {
//             ncircles = 17;
//             
//             xx[0] = 0.5*(X_SHOOTER + X_TARGET);
//             xx[1] = LAMBDA - 0.5*(X_TARGET - X_SHOOTER);
//             xx[2] = -xx[0];
//             xx[3] = -xx[1];
//             
//             yy[0] = 0.5*(Y_SHOOTER + Y_TARGET);
//             yy[1] = 1.0 - 0.5*(Y_TARGET - Y_SHOOTER);
//             yy[2] = -yy[0];
//             yy[3] = -yy[1];
// 
//             for (i = 0; i < 4; i++)
//                 for (j = 0; j < 4; j++)
//                 {
//                     circlex[4*i + j] = xx[i];
//                     circley[4*i + j] = yy[j];
//                     
//                 }
//                 
//             circlex[ncircles - 1] = X_TARGET;
//             circley[ncircles - 1] = Y_TARGET;
//             
//             for (i=0; i<ncircles - 1; i++)
//             {
//                 circlerad[i] = MU;
//                 circleactive[i] = 1;
//             }
//             
//             circlerad[ncircles - 1] = 0.5*MU;
//             circleactive[ncircles - 1] = 2;
//             
//             break;
//         }
        default: 
        {
            printf("Function init_circle_config not defined for this pattern \n");
        }
    }
}
