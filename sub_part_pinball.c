/* global variables needed for circle configuration with periodic boundary conditions */

short int double_circle[NMAXCIRCLES];   /* set to 1 if a circle is a translate of another one on the boundary */
int partner_circle[NMAXCIRCLES];        /* number of circle of which current circle is a copy */    

void init_circle_config_pinball(int circle_pattern)
{
    int i, j, k, n, ncirc0, n_p_active, ncandidates=5000, naccepted; 
    double dx, dy, xx[4], yy[4], x, y, gamma, height, phi, r0, r, dpoisson = 3.25*MU;
    short int active_poisson[NMAXCIRCLES], far;
    
    switch (circle_pattern) {
        case (C_FOUR_CIRCLES):
        {
            ncircles = 4;
            
            circlex[0] = 1.0;
            circley[0] = 0.0;
            circlerad[0] = 0.8;
            
            circlex[1] = -1.0;
            circley[1] = 0.0;
            circlerad[1] = 0.8;
            
            circlex[2] = 0.0;
            circley[2] = 0.8;
            circlerad[2] = 0.4;
            
            circlex[3] = 0.0;
            circley[3] = -0.8;
            circlerad[3] = 0.4;
            
            for (i=0; i<4; i++) 
            {   
                circleactive[i] = 1;
                circlecolor[i] = 0;
                newcircle[i] = 0;
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
                    circlex[n] = ((double)(i-NCX/2) + 0.5)*dy;
                    circley[n] = BOXYMIN + ((double)j + 0.5)*dy;
                    circlerad[n] = MU;
                    circleactive[n] = 1;
                    circlecolor[n] = 0;
                    newcircle[n] = 0;
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
//                     circlex[n] = ((double)(i-NCX/2) + 0.5)*dy;
                    circlex[n] = ((double)(i-NCX/2))*dy;
                    if (NCX % 2 == 0) circlex[n] += 0.5*dy;
                    circley[n] = YMIN + ((double)j - 0.5)*dy;
                    if ((i+NCX)%2 == 1) circley[n] += 0.5*dy;
                    circlerad[n] = MU;
                    circleactive[n] = 1;
                    circlecolor[n] = 0;
                    newcircle[n] = 0;
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
                    circlex[n] = ((double)(i-NCX/2) + 0.5)*dx;
                    circley[n] = BOXYMIN + ((double)j)*dy;
                    if ((i+NCX)%2 == 1) circley[n] += 0.5*dy;
                    circlerad[n] = MU;
                    circleactive[n] = 1;
                    circlecolor[n] = 0;
                    newcircle[n] = 0;
                    
                    /* take care of periodic boundary conditions */
                    if (B_DOMAIN == D_CIRCLES_IN_TORUS)
                    {
                        if ((j == NCY)&&((i+NCX)%2 == 0)) 
                        {
                            double_circle[n] = 1;
                            partner_circle[n] = (NCY+1)*i;
                        }
                        else 
                        {
                            double_circle[n] = 0;
                            if ((j == 0)&&((i+NCX)%2 == 0)) partner_circle[n] = (NCY+1)*i + NCY;
                            else partner_circle[n] = n;
                        }
                    }
                    else double_circle[n] = 0;
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
                circlex[n] = -LAMBDA + n*dx;
                circley[n] = y;
                y += height*gamma; 
                if (y > YMAX) y -= height;
                circlerad[n] = MU;
                circleactive[n] = 1;
                circlecolor[n] = 0;
                newcircle[n] = 0;
            }
            break;
        }
        case (C_GOLDEN_SPIRAL):
        {
            ncircles = 1;
            circlex[0] = 0.0;
            circley[0] = 0.0;
            
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
                    circlex[ncircles] = x;
                    circley[ncircles] = y;
                    ncircles++;
                }
            }
            
            for (i=0; i<ncircles; i++)
            {
                circlerad[i] = MU;
                circlecolor[i] = 0;
                newcircle[i] = 0;
                /* inactivate circles outside the domain */
                if ((circley[i] < YMAX + MU)&&(circley[i] > YMIN - MU)) circleactive[i] = 1;
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
                    circlex[n] = ((double)(i-NCX/2) + 0.5*((double)rand()/RAND_MAX - 0.5))*dy;
                    circley[n] = YMIN + ((double)j + 0.5 + 0.5*((double)rand()/RAND_MAX - 0.5))*dy;
                    circlerad[n] = MU*sqrt(1.0 + 0.8*((double)rand()/RAND_MAX - 0.5));
                    circleactive[n] = 1;
                    circlecolor[n] = 0;
                    newcircle[n] = 0;
                }
            break;
        }
        case (C_RAND_POISSON):
        {
            ncircles = NPOISSON;
            for (n = 0; n < NPOISSON; n++)
            {
                circlex[n] = LAMBDA*(2.0*(double)rand()/RAND_MAX - 1.0);
                circley[n] = (BOXYMAX - BOXYMIN)*(double)rand()/RAND_MAX + BOXYMIN;
                circlerad[n] = MU;
                circleactive[n] = 1;
                circlecolor[n] = 0;
                newcircle[n] = 0;
                double_circle[n] = 0;
                partner_circle[n] = n;
                
                /* take care of periodic boundary conditions */
                if (B_DOMAIN == D_CIRCLES_IN_TORUS)
                {
                    /* inactivate circles in corners */
                    if ((vabs(circlex[n]) > LAMBDA - MU)&&((circley[n] < BOXYMIN + MU)||(circley[n] > BOXYMAX - MU)))
                        circleactive[n] = 0;
                    
                    if (circlex[n] < - LAMBDA + MU) 
                    {
                        circlex[ncircles] = circlex[n] + 2.0*LAMBDA;
                        circley[ncircles] = circley[n];
                        partner_circle[ncircles] = n;
                        partner_circle[n] = ncircles;
                        ncircles++;
                    }
                    else if (circlex[n] > LAMBDA - MU)
                    {
                        circlex[ncircles] = circlex[n] - 2.0*LAMBDA;
                        circley[ncircles] = circley[n];
                        partner_circle[ncircles] = n;
                        partner_circle[n] = ncircles;
                        ncircles++;                        
                    }
                    if (circley[n] < BOXYMIN + MU) 
                    {
                        circlex[ncircles] = circlex[n];
                        circley[ncircles] = circley[n] + BOXYMAX - BOXYMIN;
                        partner_circle[ncircles] = n;
                        partner_circle[n] = ncircles;
                        ncircles++;
                    }
                    else if (circley[n] > BOXYMAX - MU) 
                    {
                        circlex[ncircles] = circlex[n];
                        circley[ncircles] = circley[n] - BOXYMAX + BOXYMIN;
                        partner_circle[ncircles] = n;
                        partner_circle[n] = ncircles;
                        ncircles++;
                    }
                }
            }
            printf("%i circles\n", ncircles);
            if (B_DOMAIN == D_CIRCLES_IN_TORUS) for (n = NPOISSON; n < ncircles; n++)
            {
//                 printf("circle %i at (%.3f, %.3f)\n", n, circlex[n], circley[n]); 
                circlerad[n] = MU;
                if (circleactive[partner_circle[n]]) circleactive[n] = 1;
                circlecolor[n] = 0;
                newcircle[n] = 0;
                double_circle[n] = 1;
            }

            break;
        }
        case (C_POISSON_DISC):
        {
            printf("Generating Poisson disc sample\n");
            /* generate first circle */
            circlex[0] = LAMBDA*(2.0*(double)rand()/RAND_MAX - 1.0);
            circley[0] = (YMAX - YMIN)*(double)rand()/RAND_MAX + YMIN;
            active_poisson[0] = 1;
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
                    x = circlex[i] + r*cos(phi);
                    y = circley[i] + r*sin(phi);
//                        printf("Testing new circle at (%.3f,%.3f)\t", x, y);
                    far = 1;
                    for (k=0; k<ncircles; k++) if ((k!=i))
                    {
                        /* new circle is far away from circle k */
                        far = far*((x - circlex[k])*(x - circlex[k]) + (y - circley[k])*(y - circley[k]) >=     dpoisson*dpoisson);
                        /* new circle is in domain */
                        far = far*(vabs(x) < LAMBDA)*(y < YMAX)*(y > YMIN);
                    }
                    if (far)    /* accept new circle */
                    {
                        printf("New circle at (%.3f,%.3f) accepted\n", x, y);
                        circlex[ncircles] = x;
                        circley[ncircles] = y;
                        circlerad[ncircles] = MU;
                        circleactive[ncircles] = 1;
                        circlecolor[ncircles] = 0;
                        newcircle[ncircles] = 0;
                        active_poisson[ncircles] = 1;
                        ncircles++;
                        n_p_active++;
                        naccepted++;
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
