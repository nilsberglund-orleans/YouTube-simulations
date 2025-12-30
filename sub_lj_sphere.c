/* some graphic routines moved here from su_lj.c */

#include "colors_waves.c"

void blank()
{
    if (BLACK) glClearColor(0.0, 0.0, 0.0, 1.0);
    else glClearColor(1.0, 1.0, 1.0, 1.0);
    glClear(GL_COLOR_BUFFER_BIT);
}

/*********************/
/* some basic math   */
/*********************/

 double vabs(double x)     /* absolute value */
 {
	double res;

	if (x<0.0) res = -x;
	else res = x;
	return(res);
 }

  double ipow(double x, int n)
 {
    double y;
    int i;
    
    y = x;
    for (i=1; i<n; i++) y *= x;
    
    return(y);
 }
 

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

double argument(double x, double y)
 {
	double alph;

	if (x!=0.0)
	{
		alph = atan(y/x);
		if (x<0.0)
			alph += PI;
	}
	else
	{
		alph = PID;
		if (y<0.0)
			alph = PI*1.5;
	}
	return(alph);
 }
 
double dist_sphere(double phi1, double theta1, double phi2, double theta2)
/* returns angle between points given in spherical coordinates */
{
    double dist;
    
    dist = sin(theta1)*sin(theta2)*cos(phi1 - phi2);
    dist += cos(theta1)*cos(theta2);
                
    return(acos(dist));
}


void compute_midpoint_sphere(double phi1, double theta1, double phi2, double theta2, double *phi, double *theta)
/* compute midpoint between two points on the sphere */
{
    double x1, y1, z1, x2, y2, z2, x3, y3, z3, r;
    
    x1 = cos(phi1)*sin(theta1);
    y1 = sin(phi1)*sin(theta1);
    z1 = -cos(theta1);
    
    x2 = cos(phi2)*sin(theta2);
    y2 = sin(phi2)*sin(theta2);
    z2 = -cos(theta2);
    
    x3 = x1 + x2;
    y3 = y1 + y2;
    z3 = z1 + z2;
    
    r = sqrt(x3*x3 + y3*y3 + z3*z3);
    *theta = acos(-z3/r);
    *phi = argument(x3, y3);
    if (*phi < 0.0) *phi += DPI;
}


void compute_convex_combination_sphere(double phi1, double theta1, double phi2, double theta2, double *phia, double *thetaa, double *phib, double *thetab, double lambda)
/* compute convex combination of two points on the sphere */
{
    double x1, y1, z1, x2, y2, z2, x3, y3, z3, r, lam1;
    
    lam1 = 1.0 - lambda;
    
    x1 = cos(phi1)*sin(theta1);
    y1 = sin(phi1)*sin(theta1);
    z1 = -cos(theta1);
    
    x2 = cos(phi2)*sin(theta2);
    y2 = sin(phi2)*sin(theta2);
    z2 = -cos(theta2);
    
    x3 = lambda*x1 + lam1*x2;
    y3 = lambda*y1 + lam1*y2;
    z3 = lambda*z1 + lam1*z2;
    
    r = sqrt(x3*x3 + y3*y3 + z3*z3);
    *thetaa = acos(-z3/r);
    *phia = argument(x3, y3);
    if (*phia < 0.0) *phia += DPI;
    
    x3 = lam1*x1 + lambda*x2;
    y3 = lam1*y1 + lambda*y2;
    z3 = lam1*z1 + lambda*z2;
    
    r = sqrt(x3*x3 + y3*y3 + z3*z3);
    *thetab = acos(-z3/r);
    *phib = argument(x3, y3);
    if (*phib < 0.0) *phib += DPI;
}

int in_polygon(double x, double y, double r, int npoly, double apoly)
/* test whether (x,y) is in regular polygon of npoly sides inscribed in circle of radius r, turned by apoly Pi/2 */
{
    int condition = 1, k;
    double omega, cw, angle; 
    
    omega = DPI/((double)npoly);
    cw = cos(omega*0.5);
    for (k=0; k<npoly; k++)  
    {
        angle = -apoly*PID + ((double)k+0.5)*omega;
        condition = condition*(x*cos(angle) + y*sin(angle) < r*cw);
    }
    return(condition);
}

 
double neighbour_color(int nnbg)
{
    if (nnbg > 7) nnbg = 7;
    switch(nnbg){
        case (7): return(340.0);
        case (6): return(310.0);
        case (5): return(260.0);
        case (4): return(200.0);
        case (3): return(140.0);
        case (2): return(100.0);
        case (1): return(70.0);
        default:  return(30.0);
    }   
}


double partner_color(int np)
{
    switch(np){
        case (0): return(340.0);
        case (1): return(260.0);
        case (2): return(210.0);
        case (3): return(140.0);
        case (4): return(70.0);
        default:  return(20.0);
//         case (0): return(70.0);
//         case (1): return(200.0);
//         case (2): return(280.0);
//         case (3): return(140.0);
//         case (4): return(320.0);
//         default:  return(20.0);
    }   
}


double type_hue(int type)
{
    int hue;
    double t2;
    static double b, hmax;
    static int first = 1;
    
    if (first)
    {
        hmax = 360.0;
        b = 16.0*(hmax - HUE_TYPE3);
        first = 0;
    }
    
//     if ((RD_REACTION == CHEM_CATALYTIC_A2D)&&(type == 4)) return(HUE_TYPE3); 
    
    if ((RD_REACTION == CHEM_ABDACBE)&&(type == 4)) return(HUE_TYPE3); 
    if ((RD_REACTION == CHEM_ABDACBE)&&(type == 5)) return(280.0); 
    
    switch (type) {
        case (0): return(HUE_TYPE0);
        case (1): return(HUE_TYPE1);
        case (2): return(HUE_TYPE2);
        case (3): return(HUE_TYPE3);
        case (4): return(HUE_TYPE4);
        case (5): return(HUE_TYPE5);
        case (6): return(HUE_TYPE6);
        case (7):
        {
            if (RD_REACTION == CHEM_BZ) return(HUE_TYPE2);
            else return(HUE_TYPE7);
        }
        case (8): return(HUE_TYPE8);
        default:
        {
            if (RD_REACTION == CHEM_BZ)
            {
                if (type == 7) return(HUE_TYPE2);
                if (type == 8) type = 5;
            }
            else if ((RD_REACTION == CHEM_BRUSSELATOR)&&(type >= 5)) return(70.0);
            t2 = (double)(type*type);
            hue = (hmax*t2 - b)/t2;
            return(hue);
        }
    }
}

void compute_particle_colors(t_particle particle, t_cluster cluster[NMAXCIRCLES], int plot, double rgb[3], double rgbx[3], double rgby[3], t_particle other_particle[NMAXCIRCLES])
{
    double ej, angle, hue, huex, huey, lum, x, y, ccluster, xg, yg, vx, vy, amoment, dist, density, dx, dy;
    int i, k, p, q, cl, mol[2], mass[2];
    
    switch (plot) {
        case (P_KINETIC): 
        {
            ej = particle.energy;
            if (ej > 0.0) 
            {
                hue = ENERGY_HUE_MIN + (ENERGY_HUE_MAX - ENERGY_HUE_MIN)*ej/PARTICLE_EMAX;
                if (hue > ENERGY_HUE_MIN) hue = ENERGY_HUE_MIN;
                if (hue < ENERGY_HUE_MAX) hue = ENERGY_HUE_MAX;
            }
            break;
        }
        case (P_NEIGHBOURS): 
        {
            hue = neighbour_color(particle.neighb);
            break;
        }
        case (P_BONDS):
        {
            hue = neighbour_color(particle.neighb);
            break;
        }
        case (P_ANGLE):
        {
            angle = particle.angle;
            hue = angle*particle.spin_freq/DPI;
            hue -= (double)((int)hue);
            huex = (DPI - angle)*particle.spin_freq/DPI;
            huex -= (double)((int)huex);
            angle = PI - angle;
            if (angle < 0.0) angle += DPI;
            huey = angle*particle.spin_freq/DPI;
            huey -= (double)((int)huey);
            hue = PARTICLE_HUE_MIN + (PARTICLE_HUE_MAX - PARTICLE_HUE_MIN)*(hue);
            huex = PARTICLE_HUE_MIN + (PARTICLE_HUE_MAX - PARTICLE_HUE_MIN)*(huex);
            huey = PARTICLE_HUE_MIN + (PARTICLE_HUE_MAX - PARTICLE_HUE_MIN)*(huey);
            break;
        }
        case (P_TYPE):
        {
            hue = type_hue(particle.type);
            break;
        }
        case (P_DIRECTION): 
        {
            hue = argument(particle.vx, particle.vy);
            if (hue > DPI) hue -= DPI;
            if (hue < 0.0) hue += DPI;
            hue = PARTICLE_HUE_MIN + (PARTICLE_HUE_MAX - PARTICLE_HUE_MIN)*(hue)/DPI;
            break;
        }
        case (P_DIRECT_ENERGY): 
        {
            hue = argument(particle.vx, particle.vy);
            if (hue > DPI) hue -= DPI;
            if (hue < 0.0) hue += DPI;
            hue = PARTICLE_HUE_MIN + (PARTICLE_HUE_MAX - PARTICLE_HUE_MIN)*(hue)/DPI;
            if (particle.energy < 0.1*PARTICLE_EMAX) lum = 10.0*particle.energy/PARTICLE_EMAX;
            else lum = 1.0;
            break;
        }
        case (P_ANGULAR_SPEED): 
        {
            hue = 160.0*(1.0 + tanh(SLOPE*particle.omega));
//                printf("omega = %.3lg, hue = %.3lg\n", particle[j].omega, hue);
            break;
        }
        case (P_DIFF_NEIGHB): 
        {
            hue = (double)(particle.diff_neighb+1)/(double)(particle.neighb+1);
//                hue = PARTICLE_HUE_MIN + (PARTICLE_HUE_MAX - PARTICLE_HUE_MIN)*hue;
//             hue = 180.0*(1.0 + hue);
            hue = 20.0 + 320.0*hue;
            break;
        }
        case (P_THERMOSTAT):
        {
            if (particle.thermostat) hue = 30.0;
            else hue = 270.0;
            break;
        }
        case (P_INITIAL_POS):
        {
            hue = particle.color_hue;
            break;
        }
        case (P_NUMBER):
        {
            hue = particle.color_hue;
            break;
        }
        case (P_EMEAN): 
        {
            ej = particle.emean;
            if (ej > 0.0) 
            {
                hue = ENERGY_HUE_MIN + (ENERGY_HUE_MAX - ENERGY_HUE_MIN)*ej/PARTICLE_EMAX;
                if (hue > ENERGY_HUE_MIN) hue = ENERGY_HUE_MIN;
                if (hue < ENERGY_HUE_MAX) hue = ENERGY_HUE_MAX;
            }
            break;
        }
        case (P_EMEAN_DENSITY): 
        {
            ej = particle.emean;
            cl = particle.cluster;
            ej *= PARTICLE_MASS/cluster[cl].mass;
            if (ej > 0.0) 
            {
                hue = ENERGY_HUE_MIN + (ENERGY_HUE_MAX - ENERGY_HUE_MIN)*ej/PARTICLE_EMAX;
                if (hue > ENERGY_HUE_MIN) hue = ENERGY_HUE_MIN;
                if (hue < ENERGY_HUE_MAX) hue = ENERGY_HUE_MAX;
            }
            break;
        }
        case (P_LOG_EMEAN): 
        {
            ej = particle.emean;
            if (ej > 0.0) 
            {
                ej = log(ej/PARTICLE_EMIN)/log(PARTICLE_EMAX/PARTICLE_EMIN);
                if (ej < 0.0) ej = 0.0;
                else if (ej > 1.0) ej = 1.0;
                hue = ENERGY_HUE_MIN + (ENERGY_HUE_MAX - ENERGY_HUE_MIN)*ej;
//                 if (hue > ENERGY_HUE_MIN) hue = ENERGY_HUE_MIN;
//                 if (hue < ENERGY_HUE_MAX) hue = ENERGY_HUE_MAX;
            }
            break;
        }
        case (P_DIRECT_EMEAN): 
        {
            hue = particle.dirmean + COLOR_HUESHIFT*PI;
//             printf("dirmean = %.3lg\n", particle.dirmean);
            if (hue > DPI) hue -= DPI;
            if (hue < 0.0) hue += DPI;
            hue = PARTICLE_HUE_MIN + (PARTICLE_HUE_MAX - PARTICLE_HUE_MIN)*(hue)/DPI;
            ej = particle.emean;
            if (ej < 0.5*PARTICLE_EMAX) lum = 2.0*ej/PARTICLE_EMAX;
            else lum = 1.0;
            break;
        }
        case (P_NOPARTICLE): 
        {
            hue = 0.0;
            lum = 1.0;
            break;
        }
        case (P_NPARTNERS):
        {
            hue = partner_color(particle.npartners);
            break;
        }
        case (P_CHARGE):
        {
            hue = (-CHARGE_HUE_RANGE*tanh(BG_CHARGE_SLOPE*particle.charge)+1.0)*180.0;
            break;
        }
        case (P_MOL_ANGLE):
        {
            p = particle.p0;
            q = particle.p1;
            x = other_particle[q].xc - other_particle[p].xc;
            y = other_particle[q].yc - other_particle[p].yc;
            /* deal with periodic boundary conditions */
            if (x > 0.5*(XMAX - XMIN)) x -= (XMAX - XMIN);
            else if (x < 0.5*(XMIN - XMAX)) x += (XMAX - XMIN);
            if (y > 0.5*(YMAX - YMIN)) y -= (YMAX - YMIN);
            else if (y < 0.5*(YMIN - YMAX)) y += (YMAX - YMIN);
            angle = (argument(x, y) + COLOR_HUESHIFT*PI)*MOL_ANGLE_FACTOR;
//             printf("Particle p = %i, mol_angle = %i\n", p, particle.mol_angle);
//             angle = argument(x, y)*(double)particle.mol_angle;
//             angle = argument(x, y)*(double)(other_particle[particle.p0].npartners);
            while (angle > DPI) angle -= DPI;
            while (angle < 0.0) angle += DPI;
            hue = PARTICLE_HUE_MIN + (PARTICLE_HUE_MAX - PARTICLE_HUE_MIN)*(angle)/DPI;
            break;
        }
        case (P_CLUSTER):
        {
//             cluster = (double)(particle.cluster)/(double)(ncircles);
            ccluster = (double)(particle.cluster_color)/(double)(ncircles);
            ccluster -= (double)((int)ccluster);
            hue = PARTICLE_HUE_MIN + (PARTICLE_HUE_MAX - PARTICLE_HUE_MIN)*ccluster;
            break;
        }
        case (P_CLUSTER_SIZE):
        {
//             cluster = (double)(particle.cluster)/(double)(ncircles);
            ccluster = 1.0 - 5.0/((double)particle.cluster_size + 4.0);
            hue = PARTICLE_HUE_MIN + (PARTICLE_HUE_MAX - PARTICLE_HUE_MIN)*ccluster;
            break;
        }
        case (P_CLUSTER_SELECTED):
        {
            cl = particle.cluster;
            if (cluster[cl].selected) hue = COLOR_HUE_CLUSTER_SELECTED;
            else hue = COLOR_HUE_CLUSTER_NOT_SELECTED;
            break;
        }
        case (P_COLLISION):
        {
            hue = (double)particle.collision;
            if (hue > 0.0) hue = atan(0.25*(0.03*hue + 1.0))/PID;
//             {
//                 hue += 10.0;
//                 hue *= 1.0/(40.0 + hue);
//             }
            hue = PARTICLE_HUE_MIN + (PARTICLE_HUE_MAX - PARTICLE_HUE_MIN)*hue;
            break;
        }
        case (P_RADIUS):
        {
//             hue = atan(5.0*(particle.radius/MU - 0.75))/PID;
            hue = (particle.radius/MU - RANDOM_RADIUS_MIN)/RANDOM_RADIUS_RANGE;
//             hue = 0.5*(hue + 1.0);
            hue = PARTICLE_HUE_MIN + (PARTICLE_HUE_MAX - PARTICLE_HUE_MIN)*hue;
            break;
        }
        case (P_MOL_ANG_MOMENTUM):
        {
            /* move this computation to evolve_clusters? */
            mol[0] = other_particle[particle.p0].cluster;
            mol[1] = other_particle[particle.p1].cluster;
            mass[0] = 1.0/cluster[mol[0]].mass_inv;
            mass[1] = 1.0/cluster[mol[1]].mass_inv;
            xg = mass[0]*cluster[mol[0]].xg + mass[1]*cluster[mol[1]].xg;
            yg = mass[0]*cluster[mol[0]].yg + mass[1]*cluster[mol[1]].yg;
            amoment = 0.999*cluster[mol[0]].lmean;
            for (i=0; i<2; i++)
            {
                q = mol[i];
                x = cluster[q].xg - xg;
                y = cluster[q].yg - yg;
                if (x > 0.5*(XMAX - XMIN)) x -= (XMAX - XMIN);
                else if (x < 0.5*(XMIN - XMAX)) x += (XMAX - XMIN);
                if (y > 0.5*(YMAX - YMIN)) y -= (YMAX - YMIN);
                else if (y < 0.5*(YMIN - YMAX)) y += (YMAX - YMIN);
                vx = cluster[q].vx;
                vy = cluster[q].vy;
                amoment += 0.001*(x*vy - y*vx)*mass[i];
//                 printf("(x,y) = (%.3lg, %.3lg), (vx, vy) = (%.3lg, %.3lg)\n", x, y, vx, vy);
            }
            
            for (i=0; i<2; i++) cluster[mol[i]].lmean = amoment;
//             printf("cluster %i: p0 = %i, p1 = %i, L = %.3lg\n", cl, particle.p0, particle.p1, amoment);
            
            hue = 180.0*(1.0 + tanh(cluster[mol[0]].lmean/PARTICLE_LMAX));
//             printf("L/LMAX = %.3lg, hue = %.3lg\n", amoment/PARTICLE_LMAX, hue); 
            break;
        }
        case (P_DENSITY): 
        {
            density = 0.0;
            for (i=0; i<particle.hash_nneighb; i++)
            {
                dx = particle.deltax[i];
                dy = particle.deltay[i];
                if (BOUNDARY_COND == BC_SPHERE)
                {
                    dist = dist_sphere(particle.xc, particle.yc, particle.xc + dx, particle.yc + dy);
                    dist *= dist;
                }
                else dist = dx*dx + dy*dy; 
                density += 1.0/dist;
            }
            density -= PARTICLE_DMIN;
//             printf("Density = %.3lg\n", density); 
            hue = ENERGY_HUE_MAX + (ENERGY_HUE_MIN - ENERGY_HUE_MAX)*density/PARTICLE_DMAX;
            if (hue > ENERGY_HUE_MIN) hue = ENERGY_HUE_MIN;
            if (hue < ENERGY_HUE_MAX) hue = ENERGY_HUE_MAX;
            break;
        }
        
    }
    
    switch (plot) {
        case (P_KINETIC):  
        {
//             hsl_to_rgb_turbo(hue, 0.9, 0.5, rgb);
//             hsl_to_rgb_turbo(hue, 0.9, 0.5, rgbx);
//             hsl_to_rgb_turbo(hue, 0.9, 0.5, rgby);
            hsl_to_rgb_palette(hue, 0.9, 0.5, rgb, COLOR_PALETTE_EKIN);
            hsl_to_rgb_palette(hue, 0.9, 0.5, rgbx, COLOR_PALETTE_EKIN);
            hsl_to_rgb_palette(hue, 0.9, 0.5, rgby, COLOR_PALETTE_EKIN);
            break;
        }
        case (P_BONDS):  
        {
            hsl_to_rgb_turbo(hue, 0.9, 0.5, rgb);
            hsl_to_rgb_turbo(hue, 0.9, 0.5, rgbx);
            hsl_to_rgb_turbo(hue, 0.9, 0.5, rgby);
            break;
        }
        case (P_DIRECTION): 
        {
            hsl_to_rgb_palette(hue, 0.9, 0.5, rgb, COLOR_PALETTE_DIRECTION);
            hsl_to_rgb_palette(hue, 0.9, 0.5, rgbx, COLOR_PALETTE_DIRECTION);
            hsl_to_rgb_palette(hue, 0.9, 0.5, rgby, COLOR_PALETTE_DIRECTION);
            break;
        }
        case (P_ANGLE): 
        {
            hsl_to_rgb_palette(hue, 0.9, 0.5, rgb, COLOR_PALETTE_ANGLE);
            hsl_to_rgb_palette(hue, 0.9, 0.5, rgbx, COLOR_PALETTE_ANGLE);
            hsl_to_rgb_palette(hue, 0.9, 0.5, rgby, COLOR_PALETTE_ANGLE);
            break;
        }
        case (P_DIRECT_ENERGY): 
        {
            hsl_to_rgb_twilight(hue, 0.9, 0.5, rgb);
            hsl_to_rgb_twilight(hue, 0.9, 0.5, rgbx);
            hsl_to_rgb_twilight(hue, 0.9, 0.5, rgby);
            for (i=0; i<3; i++)
            {
                rgb[i] *= lum;
                rgbx[i] *= lum;
                rgby[i] *= lum;
            }
            break;
        }
        case (P_DIFF_NEIGHB):  
        {
            hsl_to_rgb_palette(hue, 0.9, 0.5, rgb, COLOR_PALETTE_DIFFNEIGH);
            hsl_to_rgb_palette(hue, 0.9, 0.5, rgbx, COLOR_PALETTE_DIFFNEIGH);
            hsl_to_rgb_palette(hue, 0.9, 0.5, rgby, COLOR_PALETTE_DIFFNEIGH);
//             hsl_to_rgb_twilight(hue, 0.9, 0.5, rgb);
//             hsl_to_rgb_twilight(hue, 0.9, 0.5, rgbx);
//             hsl_to_rgb_twilight(hue, 0.9, 0.5, rgby);
            break;
        }
        case (P_INITIAL_POS): 
        {
            hsl_to_rgb_palette(hue, 0.9, 0.5, rgb, COLOR_PALETTE_INITIAL_POS);
            hsl_to_rgb_palette(hue, 0.9, 0.5, rgbx, COLOR_PALETTE_INITIAL_POS);
            hsl_to_rgb_palette(hue, 0.9, 0.5, rgby, COLOR_PALETTE_INITIAL_POS);
            break;
        }
        case (P_EMEAN):  
        {
            hsl_to_rgb_palette(hue, 0.9, 0.5, rgb, COLOR_PALETTE_EKIN);
            hsl_to_rgb_palette(hue, 0.9, 0.5, rgbx, COLOR_PALETTE_EKIN);
            hsl_to_rgb_palette(hue, 0.9, 0.5, rgby, COLOR_PALETTE_EKIN);
            break;
        }
        case (P_LOG_EMEAN):  
        {
            hsl_to_rgb_palette(hue, 0.9, 0.5, rgb, COLOR_PALETTE_EKIN);
            hsl_to_rgb_palette(hue, 0.9, 0.5, rgbx, COLOR_PALETTE_EKIN);
            hsl_to_rgb_palette(hue, 0.9, 0.5, rgby, COLOR_PALETTE_EKIN);
            break;
        }
        case (P_DIRECT_EMEAN): 
        {
            hsl_to_rgb_palette(hue, 0.9, 0.5, rgb, COLOR_PALETTE_DIRECTION);
            hsl_to_rgb_palette(hue, 0.9, 0.5, rgbx, COLOR_PALETTE_DIRECTION);
            hsl_to_rgb_palette(hue, 0.9, 0.5, rgby, COLOR_PALETTE_DIRECTION);
            for (i=0; i<3; i++)
            {
                rgb[i] *= lum;
                rgbx[i] *= lum;
                rgby[i] *= lum;
            }
            break;
        }
        case (P_CHARGE):  
        {
            hsl_to_rgb_palette(hue, 0.9, 0.5, rgb, COLOR_PALETTE_CHARGE);
            hsl_to_rgb_palette(hue, 0.9, 0.5, rgbx, COLOR_PALETTE_CHARGE);
            hsl_to_rgb_palette(hue, 0.9, 0.5, rgby, COLOR_PALETTE_CHARGE);
            break;
        }
        case (P_MOL_ANGLE):  
        {
            hsl_to_rgb_palette(hue, 0.9, 0.5, rgb, COLOR_PALETTE_ANGLE);
            hsl_to_rgb_palette(hue, 0.9, 0.5, rgbx, COLOR_PALETTE_ANGLE);
            hsl_to_rgb_palette(hue, 0.9, 0.5, rgby, COLOR_PALETTE_ANGLE);
            break;
        }
        case (P_CLUSTER):  
        {
            hsl_to_rgb_palette(hue, 0.9, 0.5, rgb, COLOR_PALETTE_CLUSTER);
            hsl_to_rgb_palette(hue, 0.9, 0.5, rgbx, COLOR_PALETTE_CLUSTER);
            hsl_to_rgb_palette(hue, 0.9, 0.5, rgby, COLOR_PALETTE_CLUSTER);
            break;
        }
        case (P_CLUSTER_SIZE):  
        {
            hsl_to_rgb_palette(hue, 0.9, 0.5, rgb, COLOR_PALETTE_CLUSTER_SIZE);
            hsl_to_rgb_palette(hue, 0.9, 0.5, rgbx, COLOR_PALETTE_CLUSTER_SIZE);
            hsl_to_rgb_palette(hue, 0.9, 0.5, rgby, COLOR_PALETTE_CLUSTER_SIZE);
            break;
        }
        case (P_CLUSTER_SELECTED):  
        {
            hsl_to_rgb_palette(hue, 0.9, 0.5, rgb, COLOR_PALETTE_CLUSTER_SELECTED);
            hsl_to_rgb_palette(hue, 0.9, 0.5, rgbx, COLOR_PALETTE_CLUSTER_SELECTED);
            hsl_to_rgb_palette(hue, 0.9, 0.5, rgby, COLOR_PALETTE_CLUSTER_SELECTED);
            break;
        }
        case (P_MOL_ANG_MOMENTUM):  
        {
            hsl_to_rgb_palette(hue, 0.9, 0.5, rgb, COLOR_PALETTE_ANGULAR_MOMENTUM);
            hsl_to_rgb_palette(hue, 0.9, 0.5, rgbx, COLOR_PALETTE_ANGULAR_MOMENTUM);
            hsl_to_rgb_palette(hue, 0.9, 0.5, rgby, COLOR_PALETTE_ANGULAR_MOMENTUM);
            break;
        }
        default: 
        {
            hsl_to_rgb(hue, 0.9, 0.5, rgb);
            hsl_to_rgb(hue, 0.9, 0.5, rgbx);
            hsl_to_rgb(hue, 0.9, 0.5, rgby);
        }
    }
}


void compute_all_particle_colors(t_particle particle[NMAXCIRCLES], t_cluster cluster[NMAXCIRCLES], int plot)
/* compute the colors of all particles */
{
    int i, k;
    double rgb[3], rgbx[3], rgby[3];
    
    for (i=0; i<ncircles; i++) if (particle[i].active)
    {
        compute_particle_colors(particle[i], cluster, plot, rgb, rgbx, rgby, particle);
        for (k=0; k<3; k++)
        {
            particle[i].rgb[k] = rgb[k];
            particle[i].rgbx[k] = rgbx[k];
            particle[i].rgby[k] = rgby[k];            
        }
    }
    
}

void compute_background_color(t_particle particle[NMAXCIRCLES], t_segment segment[NMAXSEGMENTS], t_obstacle obstacle[NMAXOBSTACLES], int bg_color, t_hashgrid hashgrid[HASHX*HASHY])
/* color background according to particle properties */
{
    int i, j, k, n, p, q, m, nnb, number, avrg_fact, obs;
    double rgb[3], hue, value, p1, p2, pp1, pp2, oldhue, valx, valy, lum;
    static int first = 1, counter = 0;
    static double area_factor;
    
    if (first)
    {
        area_factor = hashgrid[mhash(0,HASHY/2)].area;
        first = 0;
    }
    
    p1 = 0.75;
    p2 = 1.0 - p1;
//     pp1 = 0.95;
    pp1 = 0.99;
    pp2 = 1.0 - pp1;
   
    for (i=0; i<HASHX; i++)
        for (j=0; j<HASHY; j++)
        {
            n = mhash(i, j);
            if (first) 
            {
                hashgrid[n].hue1 = 180.0;
                hashgrid[n].hue2 = 180.0;
            }
            /* set two old values for option DOUBLE_MOVIE */
            if (DOUBLE_MOVIE)
            {
                if (counter) oldhue = hashgrid[n].hue1;
                else oldhue = hashgrid[n].hue2;
            }
            else oldhue = hashgrid[n].hue1;
            
            switch (bg_color) {
                case (BG_DENSITY):
                {
                    nnb = hashgrid[n].nneighb;
                    number = 0;
                    for (q=0; q<nnb; q++)
                    {
                        m = hashgrid[n].neighbour[q];
                        number += hashgrid[m].number;
                    }
                    number += hashgrid[n].number;
                    
//                     hue = 50.0*(double)hashgrid[n].number;
                    hue = 75.0*(double)number/(double)(nnb + 1);
                    if (BOUNDARY_COND == BC_SPHERE) 
                    {
                        hue *= 5.0*area_factor/hashgrid[n].area;
                    }
                    hue = p1*oldhue + p2*hue;
                    rgb[0] = hue/360.0;
                    rgb[1] = hue/360.0;
                    rgb[2] = hue/360.0;
                    break;
                }
                case (BG_NEIGHBOURS):
                {
                    value = 0.0;
                    nnb = hashgrid[n].nneighb;
                    for (q=0; q<nnb; q++)
                    {
                        m = hashgrid[n].neighbour[q];
                        for (k=0; k<hashgrid[m].number; k++)
                        {
                            p = hashgrid[m].particles[k];
                            value += (double)particle[p].neighb;
                        }
                    }
                    /* hashcell n counts double */
                    for (k=0; k<hashgrid[n].number; k++)
                    {
                        p = hashgrid[n].particles[k];
                        value += (double)particle[p].neighb;
                    }
                    value *= 1.0/(double)(nnb + 1);
                    
                    hue = neighbour_color((int)(value + 0.5));
                    hsl_to_rgb(hue, 0.9, 0.5, rgb);
                    break;
                }
                case (BG_CHARGE):
                {
                    avrg_fact = 3;      /* weight of central cell in hashgrid average */
                    value = 0.0;
                    nnb = hashgrid[n].nneighb;
                    for (q=0; q<nnb; q++)
                    {
                        m = hashgrid[n].neighbour[q];
                        for (k=0; k<hashgrid[m].number; k++)
                        {
                            p = hashgrid[m].particles[k];
                            value += particle[p].charge;
                        }
                    }
                    /* hashcell n counts double */
                    for (k=0; k<hashgrid[n].number; k++)
                    {
                        p = hashgrid[n].particles[k];
                        value += (double)(avrg_fact-1)*particle[p].charge;
                    }
                    value *= 1.0/(double)(nnb + avrg_fact);
                    if (CHARGE_OBSTACLES) value += hashgrid[n].charge;
                    hue = (-tanh(BG_CHARGE_SLOPE*value)+1.0)*180.0;
                    hue = p1*oldhue + p2*hue;
                    hsl_to_rgb_twilight(hue, 0.9, 0.5, rgb);
                    break;
                }
                case (BG_EKIN):
                {
                    value = 0.0;
                    nnb = hashgrid[n].nneighb;
                    for (q=0; q<nnb; q++)
                    {
                        m = hashgrid[n].neighbour[q];
                        for (k=0; k<hashgrid[m].number; k++)
                        {
                            p = hashgrid[m].particles[k];
                            value += particle[p].energy;
                        }
                    }
                    /* hashcell n counts double */
                    for (k=0; k<hashgrid[n].number; k++)
                    {
                        p = hashgrid[n].particles[k];
                        value += particle[p].energy;
                    }
                    value *= 1.0/(double)(nnb + 1);
                    hue = ENERGY_HUE_MIN + (ENERGY_HUE_MAX - ENERGY_HUE_MIN)*value/PARTICLE_EMAX;
                    if (hue > ENERGY_HUE_MIN) hue = ENERGY_HUE_MIN;
                    if (hue < ENERGY_HUE_MAX) hue = ENERGY_HUE_MAX;
                    hue = p1*oldhue + p2*hue;
                    hsl_to_rgb_palette(hue, 0.9, 0.5, rgb, COLOR_PALETTE_EKIN);
                    break;
                }
                case (BG_LOG_EKIN):
                {
                    value = 0.0;
                    nnb = hashgrid[n].nneighb;
                    for (q=0; q<nnb; q++)
                    {
                        m = hashgrid[n].neighbour[q];
                        for (k=0; k<hashgrid[m].number; k++)
                        {
                            p = hashgrid[m].particles[k];
                            value += particle[p].energy;
                        }
                    }
                    /* hashcell n counts double */
                    for (k=0; k<hashgrid[n].number; k++)
                    {
                        p = hashgrid[n].particles[k];
                        value += particle[p].energy;
                    }
                    value *= 1.0/(double)(nnb + 1);
                    hue = ENERGY_HUE_MIN + (ENERGY_HUE_MAX - ENERGY_HUE_MIN)*(BG_LOG_EKIN_SHIFT + log(value/PARTICLE_EMAX));
                    if (hue > ENERGY_HUE_MIN) hue = ENERGY_HUE_MIN;
                    if (hue < ENERGY_HUE_MAX) hue = ENERGY_HUE_MAX;
                    hue = p1*oldhue + p2*hue;
                    hsl_to_rgb_palette(hue, 0.9, 0.5, rgb, COLOR_PALETTE_EKIN);
                    break;
                }
                case (BG_EOBSTACLES):
                {
                    value = 0.0;
                    nnb = hashgrid[n].nneighb;
                    for (q=0; q<nnb; q++)
                    {
                        m = hashgrid[n].neighbour[q];
                        if (hashgrid[m].nobs) value += obstacle[hashgrid[m].obstacle].energy;
                    }
                    /* hashcell n counts four times */
                    if (hashgrid[n].nobs) value += 3.0*obstacle[hashgrid[n].obstacle].energy;
                    value *= 1.0/(double)(nnb + 4);
                    
//                     if (hashgrid[n].nobs) value += obstacle[hashgrid[n].obstacle].energy;
                    hue = ENERGY_HUE_MIN + (ENERGY_HUE_MAX - ENERGY_HUE_MIN)*value/OBSTACLE_EMAX;
                    if (hue > ENERGY_HUE_MIN) hue = ENERGY_HUE_MIN;
                    if (hue < ENERGY_HUE_MAX) hue = ENERGY_HUE_MAX;
                    hue = p1*oldhue + p2*hue;
                    hsl_to_rgb_palette(hue, 0.9, 0.5, rgb, COLOR_PALETTE_EKIN);
                    break;
                }
                case (BG_EKIN_OBSTACLES):
                {
                    value = 0.0;
                    nnb = hashgrid[n].nneighb;
                    for (q=0; q<nnb; q++)
                    {
                        m = hashgrid[n].neighbour[q];
                        for (k=0; k<hashgrid[m].number; k++)
                        {
                            p = hashgrid[m].particles[k];
                            value += particle[p].energy;
                        }
                    }
                    /* hashcell n counts double */
                    for (k=0; k<hashgrid[n].number; k++)
                    {
                        p = hashgrid[n].particles[k];
                        value += particle[p].energy;
                    }
                    value *= 1.0/(double)(nnb + 1);
                    if (hashgrid[n].nobs) value += obstacle[hashgrid[n].obstacle].energy;
                    hue = ENERGY_HUE_MIN + (ENERGY_HUE_MAX - ENERGY_HUE_MIN)*value/OBSTACLE_EMAX;
                    if (hue > ENERGY_HUE_MIN) hue = ENERGY_HUE_MIN;
                    if (hue < ENERGY_HUE_MAX) hue = ENERGY_HUE_MAX;
                    hue = p1*oldhue + p2*hue;
                    hsl_to_rgb_palette(hue, 0.9, 0.5, rgb, COLOR_PALETTE_EKIN);
                    break;
                }
                case (BG_DIR_OBSTACLES):
                {
                    valx = 0.0;
                    valy = 0.0;
                    nnb = hashgrid[n].nneighb;
                    for (q=0; q<nnb; q++)
                    {
                        m = hashgrid[n].neighbour[q];
                        if (hashgrid[m].nobs)
                        {
                            valx += obstacle[hashgrid[m].obstacle].vx;
                            valy += obstacle[hashgrid[m].obstacle].vy;
                        }
                    }
                    /* hashcell n counts more */
                    if (hashgrid[n].nobs)
                    {
                        valx += 4.0*obstacle[hashgrid[n].obstacle].vx;
                        valy += 4.0*obstacle[hashgrid[n].obstacle].vy;
                    }
                    valx *= 1.0/(double)(nnb + 4);
                    valy *= 1.0/(double)(nnb + 4);
                    hue = argument(valx, valy);
                    if (hue > DPI) hue -= DPI;
                    if (hue < 0.0) hue += DPI;
                    hue = pp1*oldhue + pp2*hue;
                    lum = module2(valx, valy)/OBSTACLE_VMAX;
                    if (lum > 1.0) lum = 1.0;
                    hsl_to_rgb_palette(hue*360.0/DPI, 0.9, 0.5*lum, rgb, COLOR_PALETTE_DIRECTION);
                    break;
                }
                case (BG_POS_OBSTACLES):
                {
                    valx = 0.0;
                    valy = 0.0;
                    nnb = hashgrid[n].nneighb;
//                     for (q=0; q<nnb; q++)
//                     {
//                         m = hashgrid[n].neighbour[q];
//                         if (hashgrid[m].nobs)
//                         {
//                             valx += obstacle[hashgrid[m].obstacle].xc;
//                             valy += obstacle[hashgrid[m].obstacle].yc;
//                         }
//                     }
                    /* hashcell n counts double */
                    if (hashgrid[n].nobs)
                    {
                        obs = hashgrid[n].obstacle;
                        valx += obstacle[obs].xc - obstacle[obs].xc0;
                        valy += obstacle[obs].yc - obstacle[obs].yc0;
                    }
//                     valx *= 1.0/(double)(nnb + 1);
//                     valy *= 1.0/(double)(nnb + 1);
                    hue = argument(valx, valy);
                    if (hue > DPI) hue -= DPI;
                    if (hue < 0.0) hue += DPI;
                    hue = pp1*oldhue + pp2*hue;
                    lum = module2(valx, valy);
                    if (lum > 1.0) lum = 1.0;
//                     lum = 1.0;
                    hsl_to_rgb_palette(hue*360.0/DPI, 0.9, 0.5*lum, rgb, COLOR_PALETTE_DIRECTION);
                    break;
                }
                case (BG_FORCE):
                {
                    value = 0.0;
                    nnb = hashgrid[n].nneighb;
                    for (q=0; q<nnb; q++)
                    {
                        m = hashgrid[n].neighbour[q];
                        for (k=0; k<hashgrid[m].number; k++)
                        {
                            p = hashgrid[m].particles[k];
                            value += module2(particle[p].fx, particle[p].fy);
                        }
                    }
                    /* hashcell n counts double */
                    for (k=0; k<hashgrid[n].number; k++)
                    {
                        p = hashgrid[n].particles[k];
                        value += module2(particle[p].fx, particle[p].fy);
                    }
                    value *= BG_FORCE_SLOPE/(double)(nnb + 1);
                    hue = (1.0 - tanh(value))*360.0;
                    hue = pp1*oldhue + pp2*hue;
                    hsl_to_rgb_turbo(hue, 0.9, 0.5, rgb);
                    break;
                }
                case (BG_CURRENTX):
                {
                    value = 0.0;
                    nnb = hashgrid[n].nneighb;
                    for (q=0; q<nnb; q++)
                    {
                        m = hashgrid[n].neighbour[q];
                        for (k=0; k<hashgrid[m].number; k++)
                        {
                            p = hashgrid[m].particles[k];
                            value += particle[p].vx*particle[p].charge;
                        }
                    }
                    /* hashcell n counts double */
                    for (k=0; k<hashgrid[n].number; k++)
                    {
                        p = hashgrid[n].particles[k];
                        value += particle[p].vx*particle[p].charge;
                    }
                    value *= 1.0/(double)(nnb + 1);
                    hue = 180.0 + 180.0*value/MAX_CURRENT;
                    if (hue > 360.0) hue = 360.0;
                    if (hue < 0.0) hue = 0.0;
                    hue = p1*oldhue + p2*hue;
                    hsl_to_rgb_palette(hue, 0.9, 0.5, rgb, COLOR_PALETTE_CURRENT);
                    break;
                }
                case (BG_DIRECTION):
                {
                    value = 0.0;
                    nnb = hashgrid[n].nneighb;
                    for (q=0; q<nnb; q++)
                    {
                        m = hashgrid[n].neighbour[q];
                        for (k=0; k<hashgrid[m].number; k++)
                        {
                            p = hashgrid[m].particles[k];
                            value += particle[p].angle*particle[p].spin_freq/DPI;
                            fprintf(lj_log, "value = %.3lg\n", value);
//                             hue -= (double)((int)hue);
                        }
                    }
                    /* hashcell n counts double */
                    for (k=0; k<hashgrid[n].number; k++)
                    {
                        p = hashgrid[n].particles[k];
                        value += particle[p].angle*particle[p].spin_freq/DPI;
                            fprintf(lj_log, "value = %.3lg\n", value);
//                         hue -= (double)((int)hue);
                    }
                    value *= 1.0/(double)(nnb + 1);
                    fprintf(lj_log, "value = %.3lg\n", value);
                    
//                     hue = value*particle[n].spin_freq/DPI;
                    hue = value - (double)((int)value);
                    fprintf(lj_log, "hue = %.3lg\n", hue);

//                     angle = PI - angle;
//                     if (angle < 0.0) angle += DPI;
                    hue = PARTICLE_HUE_MIN + (PARTICLE_HUE_MAX - PARTICLE_HUE_MIN)*(hue);
            
                    hue = p1*oldhue + p2*hue;
                    hsl_to_rgb_palette(hue, 0.9, 0.5, rgb, COLOR_PALETTE_DIRECTION);
                    break;
                }
                
            }
            if (DOUBLE_MOVIE)
            {
                if (counter) hashgrid[n].hue1 = hue;
                else hashgrid[n].hue2 = hue;
            }
            else hashgrid[n].hue1 = hue;
            
            hashgrid[n].r = rgb[0];
            hashgrid[n].g = rgb[1];
            hashgrid[n].b = rgb[2];
        }
    first = 0;
    counter = 1 - counter;
}


/* specific routines for lennardjones on a sphere */

int ij_to_sphere(int i, int j, double r, t_lj_sphere wsphere[NX_SPHERE*NY_SPHERE], double xyz[3], int use_wave_radius)
/* convert spherical to rectangular coordinates */
{
    double pscal, newr;
    static double norm_observer;
    static int first = 1;
    
    if (first)
    {
        norm_observer = sqrt(observer[0]*observer[0] + observer[1]*observer[1] + observer[2]*observer[2]);
        first = 0;
    }
    
    xyz[0] = wsphere[i*NY_SPHERE+j].x;
    xyz[1] = wsphere[i*NY_SPHERE+j].y;
    xyz[2] = wsphere[i*NY_SPHERE+j].z;
    
    pscal = xyz[0]*observer[0] + xyz[1]*observer[1] + xyz[2]*observer[2];
    
    newr = wsphere[i*NY_SPHERE+j].radius;
    xyz[0] *= newr;
    xyz[1] *= newr;
    xyz[2] *= newr;
    
    return(pscal/norm_observer > COS_VISIBLE);
}

double distance_sphere(int i, int k, t_particle* particle, double *ca, double *sa)
/* compute distance and direction on sphere */
{
    double x0, y0, z0, x1, y1, z1, x2, y2, z2, cp1, sp1, ct1, st1, cp2, sp2, ct2, st2;
    double cthat, sphat, cphat, dist, r, reg = 0.001;
    double xx0, yy0, zz0;
    
    cp1 = cos(particle[i].xc);
    sp1 = sin(particle[i].xc);
    ct1 = cos(particle[i].yc);
    st1 = sin(particle[i].yc);
     
    cp2 = cos(particle[k].xc);
    sp2 = sin(particle[k].xc);
    ct2 = cos(particle[k].yc);
    st2 = sin(particle[k].yc);
    
    x0 = cp2*st2;
    y0 = sp2*st2;
    z0 = -ct2;
    
    /* apply first rotation to particle k */
    x1 = cp1*x0 + sp1*y0;
    y1 = -sp1*x0 + cp1*y0;
    
    /* apply second rotation to particle k */
    x2 = -ct1*x1 - st1*z0;
    z2 = st1*x1 - ct1*z0;
    
    /* vector giving direction */
    cthat = -z2;
    r = 1.0/sqrt(y1*y1+x2*x2);
    sphat = y1*r;
    cphat = x2*r;
    
    /* angle of vector */
    *ca = -cthat*sphat;
    *sa = cthat*cphat;

    dist = acos(z2);
    return(dist);
}


double dist_point_to_particle_sphere(int i, double phi, double psi, t_particle* particle, double *ca, double *sa)
/* compute distance and direction between point and particle on sphere */
/* same as distance_sphere but for a general point */
{
    double x0, y0, z0, x1, y1, z1, x2, y2, z2, cp1, sp1, ct1, st1, cp2, sp2, ct2, st2;
    double cthat, sphat, cphat, dist, r, reg = 0.001;
    double xx0, yy0, zz0;
    
    cp1 = cos(particle[i].xc);
    sp1 = sin(particle[i].xc);
    ct1 = cos(particle[i].yc);
    st1 = sin(particle[i].yc);
     
    cp2 = cos(phi);
    sp2 = sin(phi);
    ct2 = cos(psi);
    st2 = sin(psi);
    
    x0 = cp2*st2;
    y0 = sp2*st2;
    z0 = -ct2;
    
    /* apply first rotation to particle k */
    x1 = cp1*x0 + sp1*y0;
    y1 = -sp1*x0 + cp1*y0;
    
    /* apply second rotation to particle k */
    x2 = -ct1*x1 - st1*z0;
    z2 = st1*x1 - ct1*z0;
    
    /* vector giving direction */
    cthat = -z2;
    r = 1.0/sqrt(y1*y1+x2*x2);
    sphat = y1*r;
    cphat = x2*r;
    
    /* angle of vector */
    *ca = -cthat*sphat;
    *sa = cthat*cphat;

    dist = acos(z2);
    return(dist);
}


void init_3d()		/* initialisation of window */
{
    double width, height;
    
    width = 2.0*(WINWIDTH/1760.0);
    height = WINHEIGHT/990.0;
    
    glLineWidth(3);

    glClearColor(0.0, 0.0, 0.0, 1.0);
    glClear(GL_COLOR_BUFFER_BIT);

//     glOrtho(-2.0, 2.0, -1.0, 1.0 , -1.0, 1.0);
    
    glOrtho(-width, width, -height, height, -1.0, 1.0);
}

void xyz_to_xy(double x, double y, double z, double xy_out[2])
{
    int i;
    double s, t, xinter[3];
    static double n2, m2, d, sm2, sn2, v[3], h[2], plane_ratio = 0.5;
    static int first = 1;
    
    if ((first)||(reset_view))
    {
        m2 = observer[0]*observer[0] + observer[1]*observer[1];
        n2 = m2 + observer[2]*observer[2];
        d = plane_ratio*n2;
        sm2 = sqrt(m2);
        sn2 = sqrt(n2);
        h[0] = observer[1]/sm2;
        h[1] = -observer[0]/sm2;
        v[0] = -observer[0]*observer[2]/(sn2*sm2);
        v[1] = -observer[1]*observer[2]/(sn2*sm2);
        v[2] = m2/(sn2*sm2);
        first = 0;
    }
    
    if (z > ZMAX_FACTOR*n2) z = ZMAX_FACTOR*n2;
    z *= Z_SCALING_FACTOR;
    s = observer[0]*x + observer[1]*y + observer[2]*z;
    t = (d - s)/(n2 - s);
    xinter[0] = t*observer[0] + (1.0-t)*x;
    xinter[1] = t*observer[1] + (1.0-t)*y;
    xinter[2] = t*observer[2] + (1.0-t)*z;
            
    xy_out[0] = XSHIFT_3D + FLIPX*XY_SCALING_FACTOR*(xinter[0]*h[0] + xinter[1]*h[1]);
    xy_out[1] = YSHIFT_3D + XY_SCALING_FACTOR*(xinter[0]*v[0] + xinter[1]*v[1] + xinter[2]*v[2]);
}


void draw_vertex_sphere(double xyz[3])
{
    double xy_screen[2];
    
    xyz_to_xy(xyz[0], xyz[1], xyz[2], xy_screen);
    glVertex2d(xy_screen[0], xy_screen[1]);
}

void init_lj_sphere(t_lj_sphere *wsphere)
/* initialize sphere data, taken from sub_sphere.c */
{
    int i, j;
    double dphi, dtheta, theta0, xy[2], phishift, reg_cot;
    
    printf("Initializing wsphere\n");
    
    dphi = DPI/(double)(NX_SPHERE-1);
    dtheta = PI/(double)(NY_SPHERE);
    
    #pragma omp parallel for private(i,j)
    for (i=0; i<NX_SPHERE; i++)
    {                
        for (j=0; j<NY_SPHERE; j++)
        {
            wsphere[i*NY_SPHERE+j].phi = (double)i*dphi;
            wsphere[i*NY_SPHERE+j].theta = (double)j*dtheta;
            
            wsphere[i*NY_SPHERE+j].cphi = cos(wsphere[i*NY_SPHERE+j].phi);
            wsphere[i*NY_SPHERE+j].sphi = sin(wsphere[i*NY_SPHERE+j].phi);
            
            wsphere[i*NY_SPHERE+j].ctheta = cos(wsphere[i*NY_SPHERE+j].theta);
            wsphere[i*NY_SPHERE+j].stheta = sin(wsphere[i*NY_SPHERE+j].theta);
            
            wsphere[i*NY_SPHERE+j].x = wsphere[i*NY_SPHERE+j].cphi*wsphere[i*NY_SPHERE+j].stheta;
            wsphere[i*NY_SPHERE+j].y = wsphere[i*NY_SPHERE+j].sphi*wsphere[i*NY_SPHERE+j].stheta;
            wsphere[i*NY_SPHERE+j].z = -wsphere[i*NY_SPHERE+j].ctheta;
            
            wsphere[i*NY_SPHERE+j].radius = 1.0;
            
            wsphere[i*NY_SPHERE+j].cos_angle_sphere = wsphere[i*NY_SPHERE+j].x*light[0] + wsphere[i*NY_SPHERE+j].y*light[1] + wsphere[i*NY_SPHERE+j].z*light[2];
            
            wsphere[i*NY_SPHERE+j].locked = 0;
        }
    }
}


double test_distance(double phi, double psi, int i, t_particle* particle)
/* determine particle shape */
{
    double d, ca, sa, theta;
    static double c1, c2, c1inv;
    static int first = 1;
    
    if (first)
    {
        switch (INTERACTION) {
            case (I_LJ_DIRECTIONAL):
            {
                c1 = 2.0/PI;
                c2 = 1.0/cos(PI/4.0);
                break;
            }
            case (I_LJ_PENTA):
            {
                c1 = 5.0/DPI;
                c1inv = DPI/5.0;
                c2 = 1.0/cos(PI/5.0);
                break;
            }
        }
        first = 0;
    }
    
    switch (INTERACTION) {
        case (I_LJ_DIRECTIONAL): 
        {
            d = dist_point_to_particle_sphere(i, phi, psi, particle, &ca, &sa);
            theta = (argument(ca, sa) - particle[i].angle + 2.0*DPI)*c1;
            theta -= (double)((int)theta);
            theta *= PID;
            theta -= 0.5*PID;
            return(d*cos(theta)*c2);
        }
        case (I_LJ_PENTA): 
        {
            d = dist_point_to_particle_sphere(i, phi, psi, particle, &ca, &sa);
            theta = (argument(ca, sa) - particle[i].angle + 2.0*DPI)*c1;
            theta -= (double)((int)theta);
            theta *= c1inv;
            theta -= PI/5.0;
            return(d*cos(theta)*c2);
        }
        case (I_LJ_DIPOLE): 
        {
            d = dist_point_to_particle_sphere(i, phi, psi, particle, &ca, &sa);
            theta = (argument(ca, sa) - particle[i].angle);
            return(2.0*d*(1.0 - 0.5*cos(2.0*theta)));
        }
        default:
        {
            if (DRAW_SPIN)
            {
                d = dist_point_to_particle_sphere(i, phi, psi, particle, &ca, &sa);
                theta = (argument(ca, sa) - particle[i].angle);
                return(2.0*d*(1.0 - 0.65*cos(theta)));
            }
            else return(dist_sphere(phi, psi, particle[i].xc, particle[i].yc));
        }
    }
}


void add_particle_to_sphere(int i, int j, int part, t_particle particle[NMAXCIRCLES], t_lj_sphere wsphere[NX_SPHERE*NY_SPHERE])
/* compute the effect at point (i,j) of adding a particle */
{
    int draw_normal_particle = 1, k;
    double x, y, r, dist, ca, sa, x1, y1, x2, y2, d1, d2, r2, omega, h;
    static double dphi, dtheta;
    static int first = 1;
    
    if (first)
    {
        dphi = DPI/(double)NX_SPHERE;
        dtheta = PI/(double)NY_SPHERE;
        first = 0;
    }
    
    x = (double)i*dphi;
    y = (double)j*dtheta;
    r = 1.2*particle[part].radius;
    
    /* special particle shapes for chemical reactions */
    if (REACTION_DIFFUSION)
    {
        switch(RD_REACTION)
        {
            case (CHEM_ABC):
            {
                if (particle[part].type == 3) 
                {
                    dist = dist_point_to_particle_sphere(part, x, y, particle, &ca, &sa);
                    x2 = dist*ca;
                    y2 = dist*sa;
                    
                    x1 = 0.35*r*cos(particle[part].angle);
                    y1 = 0.35*r*sin(particle[part].angle);
                        
                    d1 = (x1-x2)*(x1-x2) + (y1-y2)*(y1-y2);
                    d2 = (x1+x2)*(x1+x2) + (y1+y2)*(y1+y2);
                    
                    if (d2 < d1) d1 = d2;
                    
                    r2 = 0.49*r*r;
                    
                    if ((d1 < r2)&&(!wsphere[i*NY_SPHERE+j].locked))
                    {
                        wsphere[i*NY_SPHERE+j].r = particle[part].rgb[0];
                        wsphere[i*NY_SPHERE+j].g = particle[part].rgb[1];
                        wsphere[i*NY_SPHERE+j].b = particle[part].rgb[2];
                        wsphere[i*NY_SPHERE+j].radius += 1.5*sqrt(r2 - d1);
                    }
                    
                    draw_normal_particle = 0;
                }
                break;
            }
            case (CHEM_AABAA):
            {
                if (particle[part].type == 2) 
                {
                    dist = dist_point_to_particle_sphere(part, x, y, particle, &ca, &sa);
                    x2 = dist*ca;
                    y2 = dist*sa;
                    
                    x1 = 0.35*r*cos(particle[part].angle);
                    y1 = 0.35*r*sin(particle[part].angle);
                        
                    d1 = (x1-x2)*(x1-x2) + (y1-y2)*(y1-y2);
                    d2 = (x1+x2)*(x1+x2) + (y1+y2)*(y1+y2);
                    
                    if (d2 < d1) d1 = d2;
                    
                    r2 = 0.49*r*r;
                    
                    if ((d1 < r2)&&(!wsphere[i*NY_SPHERE+j].locked))
                    {
                        wsphere[i*NY_SPHERE+j].r = particle[part].rgb[0];
                        wsphere[i*NY_SPHERE+j].g = particle[part].rgb[1];
                        wsphere[i*NY_SPHERE+j].b = particle[part].rgb[2];
                        wsphere[i*NY_SPHERE+j].radius += 1.5*sqrt(r2 - d1);
                    }
                    
                    draw_normal_particle = 0;
                }
                break;
            }
            case (CHEM_CATALYTIC_A2D):
            {
                if (particle[part].type == 3) 
                {
                    dist = dist_point_to_particle_sphere(part, x, y, particle, &ca, &sa);
                    x2 = dist*ca;
                    y2 = dist*sa;
                    
                    x1 = 0.5*r*cos(particle[part].angle);
                    y1 = 0.5*r*sin(particle[part].angle);
                        
                    d1 = (x1-x2)*(x1-x2) + (y1-y2)*(y1-y2);
                    d2 = (x1+x2)*(x1+x2) + (y1+y2)*(y1+y2);
                    
                    if (d2 < d1) d1 = d2;
                    
                    r2 = 0.9*r*r;
                    
                    if ((d1 < r2)&&(!wsphere[i*NY_SPHERE+j].locked))
                    {
                        wsphere[i*NY_SPHERE+j].r = particle[part].rgb[0];
                        wsphere[i*NY_SPHERE+j].g = particle[part].rgb[1];
                        wsphere[i*NY_SPHERE+j].b = particle[part].rgb[2];
                        wsphere[i*NY_SPHERE+j].radius += 1.5*sqrt(r2 - d1);
                    }
                    
                    draw_normal_particle = 0;
                }
                else if (particle[part].type == 4) 
                {
                    dist = dist_point_to_particle_sphere(part, x, y, particle, &ca, &sa);
                    x2 = dist*ca;
                    y2 = dist*sa;
                    
                    x1 = 0.4*r*cos(particle[part].angle);
                    y1 = 0.4*r*sin(particle[part].angle);
                        
                    d1 = (x1-x2)*(x1-x2) + (y1-y2)*(y1-y2);
                    d2 = (x1+x2)*(x1+x2) + (y1+y2)*(y1+y2);
                    
                    if (d2 < d1) d1 = d2;
                    
                    r2 = 0.7*r*r;
                    
                    if ((d1 < r2)&&(!wsphere[i*NY_SPHERE+j].locked))
                    {
                        wsphere[i*NY_SPHERE+j].r = particle[part].rgb[0];
                        wsphere[i*NY_SPHERE+j].g = particle[part].rgb[1];
                        wsphere[i*NY_SPHERE+j].b = particle[part].rgb[2];
                        wsphere[i*NY_SPHERE+j].radius += 1.5*sqrt(r2 - d1);
                    }
                    
                    draw_normal_particle = 0;
                }
                break;
            }
            case (CHEM_POLYMER):
            {
                dist = dist_point_to_particle_sphere(part, x, y, particle, &ca, &sa);
                x2 = dist*ca;
                y2 = dist*sa;
                r = 1.2*MU;
                
                if ((dist < r)&&(!wsphere[i*NY_SPHERE+j].locked))
                {
                    wsphere[i*NY_SPHERE+j].r = particle[part].rgb[0];
                    wsphere[i*NY_SPHERE+j].g = particle[part].rgb[1];
                    wsphere[i*NY_SPHERE+j].b = particle[part].rgb[2];
                    wsphere[i*NY_SPHERE+j].radius += 1.5*sqrt(r*r - dist*dist);
                }
                    
                if (particle[part].type > 2)
                {
                    omega = DPI/(double)(particle[part].type - 2);
                
                    for (k=0; k<particle[part].type-2; k++)
                    {
                        x1 = 1.4*r*cos(particle[part].angle + (double)k*omega);
                        y1 = 1.4*r*sin(particle[part].angle + (double)k*omega);
                    
                        d1 = (x1-x2)*(x1-x2) + (y1-y2)*(y1-y2);
                        r2 = 0.9*r*r;
                    
                        if (d1 < r2)
                        {
                            h = 1.0 + 1.5*sqrt(r2 - d1);
                            if ((h > wsphere[i*NY_SPHERE+j].radius)&&(!wsphere[i*NY_SPHERE+j].locked))
                            {
                                wsphere[i*NY_SPHERE+j].r = particle[part].rgb[0];
                                wsphere[i*NY_SPHERE+j].g = particle[part].rgb[1];
                                wsphere[i*NY_SPHERE+j].b = particle[part].rgb[2];
                                wsphere[i*NY_SPHERE+j].radius = h;
                            }
                        }
                    }
                }
                
                draw_normal_particle = 0;
                break;
            }
        }
    }
    
    if (draw_normal_particle)
    {
        dist = test_distance(x, y, part, particle);
                
        if (dist < r)
        {
            h = 1.0 + 1.5*sqrt(r*r - dist*dist);
            if ((h > wsphere[i*NY_SPHERE+j].radius)&&(!wsphere[i*NY_SPHERE+j].locked))
            {
                wsphere[i*NY_SPHERE+j].r = particle[part].rgb[0];
                wsphere[i*NY_SPHERE+j].g = particle[part].rgb[1];
                wsphere[i*NY_SPHERE+j].b = particle[part].rgb[2];
                wsphere[i*NY_SPHERE+j].radius = h;
            }
        }
    }
}

void draw_trajectory_sphere(t_tracer trajectory[TRAJECTORY_LENGTH*N_TRACER_PARTICLES], t_hashgrid hashgrid[HASHX*HASHY], int traj_position, int traj_length, t_particle *particle, t_cluster *cluster, t_lj_sphere wsphere[NX_SPHERE*NY_SPHERE], int *tracer_n, int plot)
/* draw tracer particle trajectory */
{
    int i, j, i0, j0, time, p, q, width, imin, cell, i1, j1;
    double x1, x2, y1, y2, rgb[3], rgbx[3], rgby[3], radius, lum, lum1, rgb_bg[3];
    static double dphi, dtheta;
    static int first = 1;
    
    if (first)
    {
        dphi = DPI/(double)NX_SPHERE;
        dtheta = PI/(double)NY_SPHERE;
        first = 0;
    }
    
    if (traj_length < TRAJECTORY_DRAW_LENGTH) imin = 0;
    else imin = traj_length - TRAJECTORY_DRAW_LENGTH;
    
    if (traj_position < imin) traj_position = imin;
    
    if (traj_length < TRAJECTORY_LENGTH*TRACER_STEPS) 
        for (i=imin; i < traj_length-1; i++)
            for (j=0; j<n_tracers; j++) /*if (particle[tracer_n[j]].active)*/
            {        
//                 printf("Drawing tracer %i\n", j);
                            
                x1 = trajectory[j*TRAJECTORY_LENGTH*TRACER_STEPS + i].xc;
                x2 = trajectory[j*TRAJECTORY_LENGTH*TRACER_STEPS + i+1].xc;
                y1 = trajectory[j*TRAJECTORY_LENGTH*TRACER_STEPS + i].yc;
                y2 = trajectory[j*TRAJECTORY_LENGTH*TRACER_STEPS + i+1].yc;
                
                i0 = (int)(x1/dphi);
                j0 = (int)(y1/dtheta);
                    
                time = traj_length - i;
                lum = 0.9 - TRACER_LUM_FACTOR*(double)time/(double)(TRAJECTORY_LENGTH*TRACER_STEPS);
                if (lum < 0.0) lum = 0.0;
                lum1 = 1.0 - lum;
                
                if (COLOR_BACKGROUND)
                {
                    cell = hash_cell(x1, y1);
                    if (!wsphere[cell].locked)
                    {
                        rgb_bg[0] = hashgrid[cell].r;
                        rgb_bg[1] = hashgrid[cell].g;
                        rgb_bg[2] = hashgrid[cell].b;
                    }
                }
                
                if ((x2 != x1)||(y2 != y1)) for (p=-1; p<2; p++)
                    for (q=-1; q<2; q++)
                    {
                        i1 = i0 + p;
                        if (i1 < 0) i1 = NX_SPHERE-1;
                        if (i1 >= NX_SPHERE) i1 = 0;
                        j1 = j0 + q;
                        if (j1 < 0) j1 = 0;
                        if (j1 >= NY_SPHERE) j1 = NY_SPHERE-1;
                        cell = i1*NY_SPHERE+j1;
                        if (!wsphere[cell].locked)
                        {
                            if (COLOR_BACKGROUND)
                            {
                                wsphere[cell].r = particle[tracer_n[j]].rgb[0]*lum + lum1*rgb_bg[0];
                                wsphere[cell].g = particle[tracer_n[j]].rgb[1]*lum + lum1*rgb_bg[1];
                                wsphere[cell].b = particle[tracer_n[j]].rgb[2]*lum + lum1*rgb_bg[2];
                            }
                            else
                            {
                                wsphere[cell].r = particle[tracer_n[j]].rgb[0]*lum + lum1;
                                wsphere[cell].g = particle[tracer_n[j]].rgb[1]*lum + lum1;
                                wsphere[cell].b = particle[tracer_n[j]].rgb[2]*lum + lum1;
                            }
                        }
                    }
            }
    else 
    {
        for (i = traj_position + 1; i < traj_length-1; i++)
            for (j=0; j<n_tracers; j++)
            {        
                x1 = trajectory[j*TRAJECTORY_LENGTH*TRACER_STEPS + i].xc;
                x2 = trajectory[j*TRAJECTORY_LENGTH*TRACER_STEPS + i+1].xc;
                y1 = trajectory[j*TRAJECTORY_LENGTH*TRACER_STEPS + i].yc;
                y2 = trajectory[j*TRAJECTORY_LENGTH*TRACER_STEPS + i+1].yc;
                
                i0 = (int)(x1/dphi);
                j0 = (int)(y1/dtheta);
                    
                time = traj_position + traj_length - i;
                lum = 0.9 - TRACER_LUM_FACTOR*(double)time/(double)(TRAJECTORY_LENGTH*TRACER_STEPS);
                if (lum < 0.0) lum = 0.0;
                lum1 = 1.0 - lum;
                
                if (COLOR_BACKGROUND)
                {
                    cell = hash_cell(x1, y1);
                    if (!wsphere[cell].locked)
                    {
                        rgb_bg[0] = hashgrid[cell].r;
                        rgb_bg[1] = hashgrid[cell].g;
                        rgb_bg[2] = hashgrid[cell].b;
                    }
                }
                
                if ((x2 != x1)||(y2 != y1)) for (p=-1; p<2; p++)
                    for (q=-1; q<2; q++)
                    {
                        i1 = i0 + p;
                        if (i1 < 0) i1 = NX_SPHERE-1;
                        if (i1 >= NX_SPHERE) i1 = 0;
                        j1 = j0 + q;
                        if (j1 < 0) j1 = 0;
                        if (j1 >= NY_SPHERE) j1 = NY_SPHERE-1;
                        cell = i1*NY_SPHERE+j1;
                        if (!wsphere[cell].locked)
                        {
                            if (COLOR_BACKGROUND)
                            {
                                wsphere[cell].r = particle[tracer_n[j]].rgb[0]*lum + lum1*rgb_bg[0];
                                wsphere[cell].g = particle[tracer_n[j]].rgb[1]*lum + lum1*rgb_bg[1];
                                wsphere[cell].b = particle[tracer_n[j]].rgb[2]*lum + lum1*rgb_bg[2];
                            }
                            else
                            {
                                wsphere[cell].r = particle[tracer_n[j]].rgb[0]*lum + lum1;
                                wsphere[cell].g = particle[tracer_n[j]].rgb[1]*lum + lum1;
                                wsphere[cell].b = particle[tracer_n[j]].rgb[2]*lum + lum1;
                            }
                        }
                    }
            }
        for (i=imin; i < traj_position-1; i++)
            for (j=0; j<n_tracers; j++)
            {        
                x1 = trajectory[j*TRAJECTORY_LENGTH*TRACER_STEPS + i].xc;
                x2 = trajectory[j*TRAJECTORY_LENGTH*TRACER_STEPS + i+1].xc;
                y1 = trajectory[j*TRAJECTORY_LENGTH*TRACER_STEPS + i].yc;
                y2 = trajectory[j*TRAJECTORY_LENGTH*TRACER_STEPS + i+1].yc;
                
                time = traj_position - i;
                lum = 0.9 - TRACER_LUM_FACTOR*(double)time/(double)(TRAJECTORY_LENGTH*TRACER_STEPS);
                if (lum < 0.0) lum = 0.0;
                 lum1 = 1.0 - lum;
                
                if (COLOR_BACKGROUND)
                {
                    cell = hash_cell(x1, y1);
                    if (!wsphere[cell].locked)
                    {
                        rgb_bg[0] = hashgrid[cell].r;
                        rgb_bg[1] = hashgrid[cell].g;
                        rgb_bg[2] = hashgrid[cell].b;
                    }
                }
                
                if ((x2 != x1)||(y2 != y1)) for (p=-1; p<2; p++)
                    for (q=-1; q<2; q++)
                    {
                        i1 = i0 + p;
                        if (i1 < 0) i1 = NX_SPHERE-1;
                        if (i1 >= NX_SPHERE) i1 = 0;
                        j1 = j0 + q;
                        if (j1 < 0) j1 = 0;
                        if (j1 >= NY_SPHERE) j1 = NY_SPHERE-1;
                        cell = i1*NY_SPHERE+j1;
                        if (!wsphere[cell].locked)
                        {
                            if (COLOR_BACKGROUND)
                            {
                                wsphere[cell].r = particle[tracer_n[j]].rgb[0]*lum + lum1*rgb_bg[0];
                                wsphere[cell].g = particle[tracer_n[j]].rgb[1]*lum + lum1*rgb_bg[1];
                                wsphere[cell].b = particle[tracer_n[j]].rgb[2]*lum + lum1*rgb_bg[2];
                            }
                            else
                            {
                                wsphere[cell].r = particle[tracer_n[j]].rgb[0]*lum + lum1;
                                wsphere[cell].g = particle[tracer_n[j]].rgb[1]*lum + lum1;
                                wsphere[cell].b = particle[tracer_n[j]].rgb[2]*lum + lum1;
                            }
                        }
                    }
            }
    }
}


void draw_segments_sphere(t_segment segment[NMAXSEGMENTS], t_lj_sphere wsphere[NX_SPHERE*NY_SPHERE])
/* draw the repelling segments on the sphere */
{
    int s, i, cell, npoints, i0, j0, i1, j1, p, q;
    double x1, y1, x2, y2, x, y, dt, length;
    static double dphi, dtheta;
    static int first = 1;
    
    if (first)
    {
        dphi = DPI/(double)NX_SPHERE;
        dtheta = PI/(double)NY_SPHERE;
        first = 0;
    }
    
    for (s=0; s<nsegments; s++)
    {
        x1 = segment[s].x1;
        y1 = segment[s].y1;
        x2 = segment[s].x2;
        y2 = segment[s].y2;
        
        length = dist_sphere(x1, y1, x2, y2);
        npoints = (int)(length*300.0);
        dt = 1.0/(double)npoints;
        
        for (i=0; i<npoints; i++)
        {
            x = x1 + (x2 - x1)*(double)i*dt;
            y = y1 + (y2 - y1)*(double)i*dt;
            i0 = (int)(x/dphi);
            j0 = (int)(y/dtheta);
            for (p=-1; p<2; p++)
                for (q=-1; q<2; q++)
                {
                    i1 = i0 + p;
                    if (i1 < 0) i1 = NX_SPHERE-1;
                    if (i1 >= NX_SPHERE) i1 = 0;
                    j1 = j0 + q;
                    if (j1 < 0) j1 = 0;
                    if (j1 >= NY_SPHERE) j1 = NY_SPHERE-1;
                    cell = i1*NY_SPHERE+j1;
                    wsphere[cell].r = 0.0;
                    wsphere[cell].g = 0.0;
                    wsphere[cell].b = 0.0;
                    
                    wsphere[cell].radius += 0.005;
                }
        }
    }
}

void draw_absorbers_sphere(t_absorber absorber[NMAX_ABSORBERS], t_lj_sphere wsphere[NX_SPHERE*NY_SPHERE])
/* draw the absorbing discs on the sphere */
{
    int i, j, n, cell;
    double x, y, dist;
    static double dphi, dtheta;
    static int first = 1;
    
    if (first)
    {
        dphi = DPI/(double)NX_SPHERE;
        dtheta = PI/(double)NY_SPHERE;
        first = 0;
    }
    
    for (i=0; i<NX_SPHERE; i++)
        for (j=0; j<NY_SPHERE; j++)
            for (n=0; n<nabsorbers; n++)
            {
                x = XMIN + (double)i*dphi;
                y = YMIN + (double)j*dtheta;
                dist = dist_sphere(x, y, absorber[n].xc, absorber[n].yc);
                if (dist < absorber[n].radius)
                {
                    cell = i*NY_SPHERE+j;
                    wsphere[cell].r = 0.0;
                    wsphere[cell].g = 0.0;
                    wsphere[cell].b = 0.0;
                    
                }
                if (dist < absorber[n].radius*1.02)
                {
                    cell = i*NY_SPHERE+j;
                    wsphere[cell].locked = 1;
                    wsphere[cell].radius = 0.95;
                }
            }
}

void init_sphere_radius(t_particle particle[NMAXCIRCLES], t_hashgrid hashgrid[HASHX*HASHY], t_cluster cluster[NMAXCIRCLES], t_tracer trajectory[TRAJECTORY_LENGTH*N_TRACER_PARTICLES], t_segment segment[NMAXSEGMENTS], int traj_position, int traj_length, t_lj_sphere wsphere[NX_SPHERE*NY_SPHERE], int *tracer_n, int plot, int bg_color, t_absorber absorber[NMAX_ABSORBERS])
/* initialize sphere radius and colors from particles */
{
    int i, j, n, m, i0, j0, i1, j1, p, q, part, width, deltai, imin, imax, deltaj, jmin, jmax, cell;
    double rgb[3], rgbx[3], rgby[3], dist, d2, radius, r, x, y;
    static double dphi, dtheta;
    static int first = 1;
    
    if (first)
    {
        dphi = DPI/(double)NX_SPHERE;
        dtheta = PI/(double)NY_SPHERE;
        first = 0;
    }
    
    /* default radius and color */
    #pragma omp parallel for private(i)
    for (i=0; i<NX_SPHERE*NY_SPHERE; i++) if (!wsphere[i].locked)
    {
        wsphere[i].radius = 1.0;
        
        if ((COLOR_BACKGROUND)&&(bg_color > 0))
        {
            cell = hash_cell(wsphere[i].phi, wsphere[i].theta);
            wsphere[i].r = hashgrid[cell].r;
            wsphere[i].g = hashgrid[cell].g;
            wsphere[i].b = hashgrid[cell].b;
        }
        else
        {
            wsphere[i].r = 1.0;
            wsphere[i].g = 1.0;
            wsphere[i].b = 1.0;
        }
    }
    
    /* add tracer trajectories */
    if (TRACER_PARTICLE) 
        draw_trajectory_sphere(trajectory, hashgrid, traj_position, traj_length, particle, cluster, wsphere, tracer_n, plot);
    
    if (ADD_FIXED_SEGMENTS)
        draw_segments_sphere(segment, wsphere);
    
//     if (ADD_ABSORBERS)
//         draw_absorbers_sphere(absorber, wsphere); 
    
    /* draw zero meridian, for debugging */
//     for (j=0; j<NY_SPHERE; j++)
//         for (i=0; i<20; i++)
//         {
//             wsphere[i*NY_SPHERE + j].r = 0.5;
//             wsphere[i*NY_SPHERE + j].g = 0.5;
//             wsphere[i*NY_SPHERE + j].b = 0.5;
//         }
    
    /* add particles according to hashgrid */
    /* South Pole, q=0 */
    m = mhash(0,0);
    for (n=0; n<hashgrid[m].number; n++)
    {
        part = hashgrid[m].particles[n];
        
        if (particle[part].active)
        {
            deltaj = (int)(2.5*particle[part].radius/dtheta);
            j0 = (int)(particle[part].yc/dtheta);
            jmax = j0 + deltaj;
        
            #pragma omp parallel for private(i,j,x,y,dist,r)
            for (j=0; j<jmax; j++)
                for (i=0; i<NX_SPHERE; i++)
                    add_particle_to_sphere(i, j, part, particle, wsphere);
        }
    }
    
    /* North Pole, q=HASHY-1 */
    m = mhash(0,HASHY-1);
    for (n=0; n<hashgrid[m].number; n++)
    {
        part = hashgrid[m].particles[n];
        
        if (particle[part].active)
        {
            deltaj = (int)(2.5*particle[part].radius/dtheta);
            j0 = (int)(particle[part].yc/dtheta);
            jmin = j0 - deltaj;
        
            #pragma omp parallel for private(i,j,x,y,dist,r)
            for (j=jmin; j<NY_SPHERE; j++)
                for (i=0; i<NX_SPHERE; i++)
                    add_particle_to_sphere(i, j, part, particle, wsphere);
        }
    }
    
    /* bulk */
    for (q=1; q<HASHY-1; q++)
    {
        /* bulk, away from zero meridian */
        for (p=1; p<hashx_sphere[q]-1; p++)
        {
            m = mhash(p,q);
            for (n=0; n<hashgrid[m].number; n++)
            {
                part = hashgrid[m].particles[n];
                
                if (particle[part].active)
                {                    
                    i0 = (int)(particle[part].xc/dphi);
                    j0 = (int)(particle[part].yc/dtheta);
        
                    deltaj = (int)(2.5*particle[part].radius/dtheta);
                    jmin = j0 - deltaj; 
//                     if (jmin < 0) jmin = 0;
                    jmax = j0 + deltaj;
//                     if (jmax > NY_SPHERE-1) jmax = NY_SPHERE-1;
        
                    deltai = (int)(2.5*particle[part].radius/(dphi*sin(particle[part].yc + 0.001)));
                    imin = i0 - deltai; 
//                     if (imin < 0) imin = 0;
                    imax = i0 + deltai;
//                     if (imax > NX_SPHERE-1) imax = NX_SPHERE-1;
                    
                    #pragma omp parallel for private(i,j,x,y,dist,r)
                    for (j=jmin; j<jmax; j++)
                        for (i=imin; i<imax; i++)
                        {
                            add_particle_to_sphere(i, j, part, particle, wsphere);
                        }
                }
            }
            
        }
        /* p=0 */
        m = mhash(0,q);
        for (n=0; n<hashgrid[m].number; n++)
        {
            part = hashgrid[m].particles[n];
                
            if (particle[part].active)
            {
                i0 = (int)(particle[part].xc/dphi);
                j0 = (int)(particle[part].yc/dtheta);
        
                deltaj = (int)(2.5*particle[part].radius/dtheta);
                jmin = j0 - deltaj; 
//                     if (jmin < 0) jmin = 0;
                jmax = j0 + deltaj;
//                     if (jmax > NY_SPHERE-1) jmax = NY_SPHERE-1;
        
                deltai = (int)(2.5*particle[part].radius/(dphi*sin(particle[part].yc + 0.001)));
                imin = NX_SPHERE + i0 - deltai; 
                imax = i0 + deltai;
                    
                #pragma omp parallel for private(i,j,x,y,dist,r)
                for (j=jmin; j<jmax; j++)
                {
                    for (i=0; i<imax; i++)
                        add_particle_to_sphere(i, j, part, particle, wsphere);
                    for (i=imin; i<NX_SPHERE; i++)
                        add_particle_to_sphere(i, j, part, particle, wsphere);
                }
            }
        }
        
        /* p=hashx_sphere[q]-1 */
        m = mhash(hashx_sphere[q]-1,q);
        for (n=0; n<hashgrid[m].number; n++)
        {
            part = hashgrid[m].particles[n];
                
            if (particle[part].active)
            {
                i0 = (int)(particle[part].xc/dphi);
                j0 = (int)(particle[part].yc/dtheta);
        
                deltaj = (int)(2.5*particle[part].radius/dtheta);
                jmin = j0 - deltaj; 
//                     if (jmin < 0) jmin = 0;
                jmax = j0 + deltaj;
//                     if (jmax > NY_SPHERE-1) jmax = NY_SPHERE-1;
        
                deltai = (int)(2.5*particle[part].radius/(dphi*sin(particle[part].yc + 0.001)));
                imin = i0 - deltai; 
                imax = i0 + deltai - NX_SPHERE;
                    
                 #pragma omp parallel for private(i,j,x,y,dist,r)
                for (j=jmin; j<jmax; j++)
                {
                    for (i=imin; i<NX_SPHERE; i++)
                        add_particle_to_sphere(i, j, part, particle, wsphere);
                    for (i=0; i<imax; i++)
                        add_particle_to_sphere(i, j, part, particle, wsphere);
                }
            }
        }
    }
}


void compute_light_angle_sphere(t_lj_sphere wsphere[NX_SPHERE*NY_SPHERE])
/* compute sphere shading, adapted from compute_light_angle_sphere_rde() in sub_rde.c */
{
    int i, j, iplus, b;
    double x, y, z, norm, pscal, deltai[3], deltaj[3], deltar, n[3], r;
    static double dphi, dtheta, vshift;
    static int first = 1;
    
    if (first)
    {
        dphi = DPI/(double)NX_SPHERE;
        dtheta = PI/(double)NY_SPHERE;
        first = 0;
    }
        
    #pragma omp parallel for private(i,j,norm,pscal,deltar,deltai,deltaj,n)
    for (i=0; i<NX_SPHERE-1; i++)
        for (j=0; j<NY_SPHERE-1; j++)
        {
            /* computation of tangent vectors */
            deltar = (wsphere[(i+1)*NY_SPHERE+j].radius - wsphere[i*NY_SPHERE+j].radius)/dphi;
                    
            deltai[0] = -wsphere[i*NY_SPHERE+j].radius*wsphere[i*NY_SPHERE+j].sphi;
            deltai[0] += deltar*wsphere[i*NY_SPHERE+j].cphi;
                    
            deltai[1] = wsphere[i*NY_SPHERE+j].radius*wsphere[i*NY_SPHERE+j].cphi;
            deltai[1] += deltar*wsphere[i*NY_SPHERE+j].sphi;
                    
            deltai[2] = -deltar*wsphere[i*NY_SPHERE+j].cottheta;
            
            deltar = (wsphere[i*NY_SPHERE+j+1].radius - wsphere[i*NY_SPHERE+j].radius)/dtheta;
                    
            deltaj[0] = wsphere[i*NY_SPHERE+j].radius*wsphere[i*NY_SPHERE+j].cphi*wsphere[i*NY_SPHERE+j].ctheta;
            deltaj[0] += deltar*wsphere[i*NY_SPHERE+j].cphi*wsphere[i*NY_SPHERE+j].stheta;
                    
            deltaj[1] = wsphere[i*NY_SPHERE+j].radius*wsphere[i*NY_SPHERE+j].sphi*wsphere[i*NY_SPHERE+j].ctheta;
            deltaj[1] += deltar*wsphere[i*NY_SPHERE+j].sphi*wsphere[i*NY_SPHERE+j].stheta;
                    
            deltaj[2] = wsphere[i*NY_SPHERE+j].radius*wsphere[i*NY_SPHERE+j].stheta;
            deltaj[2] += -deltar*wsphere[i*NY_SPHERE+j].ctheta;
                    
            /* computation of normal vector */
            n[0] = deltai[1]*deltaj[2] - deltai[2]*deltaj[1];
            n[1] = deltai[2]*deltaj[0] - deltai[0]*deltaj[2];
            n[2] = deltai[0]*deltaj[1] - deltai[1]*deltaj[0];
                                        
            norm = sqrt(n[0]*n[0] + n[1]*n[1] + n[2]*n[2]);
                    
            pscal = n[0]*light[0] + n[1]*light[1] + n[2]*light[2];
                
            wsphere[i*NY_SPHERE+j].cos_angle = pscal/norm;
        }
        
    /* i = NX-1 */
    for (j=0; j<NY_SPHERE-1; j++)
    {
        /* computation of tangent vectors */
        deltar = (wsphere[j].radius - wsphere[(NX_SPHERE-1)*NY_SPHERE+j].radius)/dphi;
                    
        deltai[0] = -wsphere[(NX_SPHERE-1)*NY_SPHERE+j].radius*wsphere[(NX_SPHERE-1)*NY_SPHERE+j].sphi;
        deltai[0] += deltar*wsphere[(NX_SPHERE-1)*NY_SPHERE+j].cphi;
                    
        deltai[1] = wsphere[(NX_SPHERE-1)*NY_SPHERE+j].radius*wsphere[(NX_SPHERE-1)*NY_SPHERE+j].cphi;
        deltai[1] += deltar*wsphere[(NX_SPHERE-1)*NY_SPHERE+j].sphi;
                    
        deltai[2] = -deltar*wsphere[(NX_SPHERE-1)*NY_SPHERE+j].cottheta;
                
        deltar = (wsphere[(NX_SPHERE-1)*NY_SPHERE+j+1].radius - wsphere[(NX_SPHERE-1)*NY_SPHERE+j].radius)/dtheta;
                    
        deltaj[0] = wsphere[(NX_SPHERE-1)*NY_SPHERE+j].radius*wsphere[(NX_SPHERE-1)*NY_SPHERE+j].cphi*wsphere[(NX_SPHERE-1)*NY_SPHERE+j].ctheta;
        deltaj[0] += deltar*wsphere[(NX_SPHERE-1)*NY_SPHERE+j].cphi*wsphere[(NX_SPHERE-1)*NY_SPHERE+j].stheta;
                    
        deltaj[1] = wsphere[(NX_SPHERE-1)*NY_SPHERE+j].radius*wsphere[(NX_SPHERE-1)*NY_SPHERE+j].sphi*wsphere[(NX_SPHERE-1)*NY_SPHERE+j].ctheta;
        deltaj[1] += deltar*wsphere[(NX_SPHERE-1)*NY_SPHERE+j].sphi*wsphere[(NX_SPHERE-1)*NY_SPHERE+j].stheta;
                    
        deltaj[2] = wsphere[(NX_SPHERE-1)*NY_SPHERE+j].radius*wsphere[(NX_SPHERE-1)*NY_SPHERE+j].stheta;
        deltaj[2] += -deltar*wsphere[(NX_SPHERE-1)*NY_SPHERE+j].ctheta;
                    
        /* computation of normal vector */
        n[0] = deltai[1]*deltaj[2] - deltai[2]*deltaj[1];
        n[1] = deltai[2]*deltaj[0] - deltai[0]*deltaj[2];
        n[2] = deltai[0]*deltaj[1] - deltai[1]*deltaj[0];
                                        
        norm = sqrt(n[0]*n[0] + n[1]*n[1] + n[2]*n[2]);
                
        pscal = n[0]*light[0] + n[1]*light[1] + n[2]*light[2];
            
        wsphere[(NX_SPHERE-1)*NY_SPHERE+j].cos_angle = pscal/norm;
    }
}

void draw_lj_sphere_ij(int i, int iplus, int j, int jplus, t_lj_sphere wsphere[NX_SPHERE*NY_SPHERE], int fade, double fade_value)
/* draw particles at simulation grid point (i,j) */
{
    int k, l, n, m, s, p, q, i1, j1, prev_cell, cell, iplus1, jplus1;
    short int draw, drawij;
    double xyz[3], ca, rgb[3];
    
    ca = wsphere[i*NY_SPHERE+j].cos_angle;
    ca = (ca + 1.0)*0.4 + 0.2;
    
    glColor3f(wsphere[i*NY_SPHERE+j].r*ca, wsphere[i*NY_SPHERE+j].g*ca, wsphere[i*NY_SPHERE+j].b*ca);
    
    glBegin(GL_TRIANGLE_FAN);
    drawij = ij_to_sphere(i, j, wsphere[i*NY_SPHERE+j].radius, wsphere, xyz, draw);
    if (drawij) draw_vertex_sphere(xyz);
    if (ij_to_sphere(iplus, j, wsphere[iplus*NY_SPHERE+j].radius, wsphere, xyz, draw))
        draw_vertex_sphere(xyz);
    if (ij_to_sphere(iplus, jplus, wsphere[iplus*NY_SPHERE+j+1].radius, wsphere, xyz, draw))
        draw_vertex_sphere(xyz);
    if (ij_to_sphere(i, jplus, wsphere[i*NY_SPHERE+j+1].radius, wsphere, xyz, draw))
        draw_vertex_sphere(xyz);
    glEnd ();
}


void draw_lj_sphere_3d(t_lj_sphere wsphere[NX_SPHERE*NY_SPHERE])
/* draw sphere, adapted from draw_wave_sphere_3d_rde() in sub_rde.c */
{
    int i, j, imax, imin, jmin, jmax, imid, jmid, fade;
    double observer_angle, angle2, observer_latitude, xyz[3], fade_value;
    
    blank();
    
    /* may be defined as parameters later */
    fade = 0;
    fade_value = 1.0;
    
    observer_angle = argument(observer[0], observer[1]);
    if (observer_angle < 0.0) observer_angle += DPI;
    
    angle2 = observer_angle + PI;
    if (angle2 > DPI) angle2 -= DPI;
    
    observer_latitude = asin(observer[2]/module2(observer[0], observer[1]));
    
    imin = (int)(observer_angle*(double)NX_SPHERE/DPI);
    imax = (int)(angle2*(double)NX_SPHERE/DPI);
    if (imin >= NX_SPHERE-1) imin = NX_SPHERE-2;
    if (imax >= NX_SPHERE-1) imax = NX_SPHERE-2;
    
    jmax = NY_SPHERE;
    jmin = 0;
    jmid = (int)((double)NY_SPHERE*(observer_latitude + PID)/PI);
    imid = (imin + imax)/2;
    
    printf("Angle = %.5lg, angle2 = %.5lg, imin = %i, imax = %i\n", observer_angle, angle2, imin, imax);
    
    if (observer[2] > 0.0)
    {
        if (imin < imax)
        {
            for (i=imax; i>imid; i--)
            {
                for (j=jmin; j<=jmax; j++)
                    draw_lj_sphere_ij(i, i+1, j, j+1, wsphere, fade, fade_value);
            }
            for (i=imid; i>imin; i--)
            {
                for (j=jmin; j<=jmid; j++)
                    draw_lj_sphere_ij(i, i+1, j, j+1, wsphere, fade, fade_value);
                for (j=jmax; j>=jmid; j--)
                    draw_lj_sphere_ij(i, i+1, j, j+1, wsphere, fade, fade_value);
            }
            
            for (i=imax+1; i<NX_SPHERE-1; i++)
            {
                for (j=jmin; j<=jmax; j++)
                    draw_lj_sphere_ij(i, i+1, j, j+1, wsphere, fade, fade_value);
            }
        
            for (j=jmin; j<=jmax; j++)
            {
                draw_lj_sphere_ij(NX_SPHERE-1, 0, j, j+1, wsphere, fade, fade_value);
                draw_lj_sphere_ij(0, 1, j, j+1, wsphere, fade, fade_value);
            }
            
            for (i=0; i<=imin; i++)
            {
                for (j=jmin; j<=jmid; j++)
                    draw_lj_sphere_ij(i, i+1, j, j+1, wsphere, fade, fade_value);
                for (j=jmax; j>=jmid; j--)
                    draw_lj_sphere_ij(i, i+1, j, j+1, wsphere, fade, fade_value);
            }
            
            if (imin >= 1)
            {
                for (j=jmax; j>=jmid; j--)
                    draw_lj_sphere_ij(imin-1, imin, j, j+1, wsphere, fade, fade_value);
            }
        }
        else
        {
            for (i=imax; i<imid; i++)
            {
                for (j=jmin; j<=jmax; j++)
                    draw_lj_sphere_ij(i, i+1, j, j+1, wsphere, fade, fade_value);
            }
        
            for (i=imid; i<imin; i++)
            {
                for (j=jmin; j<=jmid; j++)
                    draw_lj_sphere_ij(i, i+1, j, j+1, wsphere, fade, fade_value);
                for (j=jmax; j>=jmid; j--)
                    draw_lj_sphere_ij(i, i+1, j, j+1, wsphere, fade, fade_value);
            }
        
            for (i=imax-1; i>=0; i--)
            {
                for (j=jmin; j<=jmax; j++)
                    draw_lj_sphere_ij(i, i+1, j, j+1, wsphere, fade, fade_value);
            }
            
            for (j=jmin; j<=jmax; j++)
            {
                draw_lj_sphere_ij(NX_SPHERE-1, 0, j, j+1, wsphere, fade, fade_value);
                draw_lj_sphere_ij(0, 1, j, j+1, wsphere, fade, fade_value);
            }
            
            for (i=NX_SPHERE-2; i>=imin; i--)
            {
                for (j=jmin; j<=jmid; j++)
                    draw_lj_sphere_ij(i, i+1, j, j+1, wsphere, fade, fade_value);
                for (j=jmax; j>=jmid; j--)
                    draw_lj_sphere_ij(i, i+1, j, j+1, wsphere, fade, fade_value);
            }

            if (imin >= 1)
            {
                /* experimental */
                for (j=jmid/3; j<=jmid; j++)
                    draw_lj_sphere_ij(imin-1, imin, j, j+1, wsphere, fade, fade_value);
                for (j=jmax; j>=jmid; j--)
                    draw_lj_sphere_ij(imin-1, imin, j, j+1, wsphere, fade, fade_value);
            }

        }
    
        /* North pole */
        for (i=0; i<NX_SPHERE-1; i++) 
            draw_lj_sphere_ij(i, i+1, NY_SPHERE-3, NY_SPHERE-1, wsphere, fade, fade_value);
        
        draw_lj_sphere_ij(NX_SPHERE-1, 0, NY_SPHERE-3, NY_SPHERE-1, wsphere, fade, fade_value);
    }
    else
    {
        if (imin < imax)
        {
            for (i=imax; i>imid; i--)
            {
                for (j=jmax; j>=jmin; j--)
                    draw_lj_sphere_ij(i, i+1, j, j+1, wsphere, fade, fade_value);
            }
            
            for (i=imid; i>imin; i--)
            {
                for (j=jmax; j>=jmid; j--)
                    draw_lj_sphere_ij(i, i+1, j, j+1, wsphere, fade, fade_value);
                for (j=jmin; j<=jmid; j++)
                    draw_lj_sphere_ij(i, i+1, j, j+1, wsphere, fade, fade_value);
            }
            
            for (i=imax+1; i<NX_SPHERE-1; i++)
            {
                for (j=jmax; j>=jmin; j--)
                    draw_lj_sphere_ij(i, i+1, j, j+1, wsphere, fade, fade_value);
            }
            
            for (j=jmax; j>=jmin; j--)
                draw_lj_sphere_ij(NX_SPHERE-1, 0, j, j+1, wsphere, fade, fade_value);
        
            for (i=0; i<=imin; i++)
            {
                for (j=jmax; j>=jmid; j--)
                    draw_lj_sphere_ij(i, i+1, j, j+1, wsphere, fade, fade_value);
                for (j=jmin; j<=jmid; j++)
                    draw_lj_sphere_ij(i, i+1, j, j+1, wsphere, fade, fade_value);
            }
            
            /* TEST */
            for (j=jmin; j<=jmid; j++)
                for (i=-1; i<2; i++) if ((imin+i >= 0)&&(imin+i+1 < NX_SPHERE))
                    draw_lj_sphere_ij(imin+i, imin+i+1, j, j+1, wsphere, fade, fade_value);
            for (j=jmax; j>=jmid; j--)
                for (i=-1; i<2; i++) if ((imin+i >= 0)&&(imin+i+1 < NX_SPHERE))
                    draw_lj_sphere_ij(imin+i, imin+i+1, j, j+1, wsphere, fade, fade_value);
        }
        else
        {
            for (i=imax; i<imid; i++)
            {
                for (j=jmax; j>=jmin; j--)
                    draw_lj_sphere_ij(i, i+1, j, j+1, wsphere, fade, fade_value);
            }
            
            for (i=imid; i<imin-1; i++)
            {
                for (j=jmax; j>=jmid; j--)
                    draw_lj_sphere_ij(i, i+1, j, j+1, wsphere, fade, fade_value);
                for (j=jmin; j<=jmid; j++)
                    draw_lj_sphere_ij(i, i+1, j, j+1, wsphere, fade, fade_value);
            }
            
            for (i=imax-1; i>=0; i--)
            {
                for (j=jmax; j>=jmin; j--)
                    draw_lj_sphere_ij(i, i+1, j, j+1, wsphere, fade, fade_value);
            }
            
            for (j=jmax; j>=jmin; j--)
                draw_lj_sphere_ij(NX_SPHERE-1, 0, j, j+1, wsphere, fade, fade_value);
        
            for (i=NX_SPHERE-2; i>=imin; i--)
            {
                for (j=jmax; j>=jmid; j--)
                    draw_lj_sphere_ij(i, i+1, j, j+1, wsphere, fade, fade_value);
                for (j=jmin; j<=jmid; j++)
                    draw_lj_sphere_ij(i, i+1, j, j+1, wsphere, fade, fade_value);
            }
            
            for (i=0; i<2; i++) if (imin >= i)
            {
                for (j=2*jmax/3; j>=jmid; j--)
                    draw_lj_sphere_ij(imin-i, imin-i+1, j, j+1, wsphere, fade, fade_value);
                for (j=jmin; j<=jmid; j++)
                    draw_lj_sphere_ij(imin-i, imin-i+1, j, j+1, wsphere, fade, fade_value);
            }
            
            /* TEST */
            for (j=jmin; j<=jmid; j++)
                for (i=-1; i<2; i++) if ((imin+i >= 0)&&(imin+i+1 < NX_SPHERE))
                    draw_lj_sphere_ij(imin+i, imin+i+1, j, j+1, wsphere, fade, fade_value);
            for (j=jmax; j>=jmid; j--)
                for (i=-1; i<2; i++) if ((imin+i >= 0)&&(imin+i+1 < NX_SPHERE))
                    draw_lj_sphere_ij(imin+i, imin+i+1, j, j+1, wsphere, fade, fade_value);
        }
    
        /* South pole */
        for (i=0; i<NX_SPHERE-1; i++) for (j=2; j>0; j--)
            draw_lj_sphere_ij(i, i+1, j-1, j, wsphere, fade, fade_value);
        
        for (j=2; j>0; j--)
            draw_lj_sphere_ij(NX_SPHERE-1, 0, j-1, j, wsphere, fade, fade_value);
    }
}

void viewpoint_schedule(int i)
/* change position of observer */
{
    int j;
    double angle, ca, sa, r1, interpolate, rho;
    static double observer_initial[3], r, ratio, rho0, zmax;
    static int first = 1;
    
    if (first)
    {
        for (j=0; j<3; j++) observer_initial[j] = observer[j];
        r1 = observer[0]*observer[0] + observer[1]*observer[1];
        r = sqrt(r1 + observer[2]*observer[2]);
        ratio = r/sqrt(r1);
        rho0 = module2(observer[0], observer[1]);
        if (vabs(rho0) < 0.001) rho0 = 0.001; 
        zmax = r*sin(MAX_LATITUDE*PI/180.0);
        first = 0;
    }
    
    interpolate = (double)i/(double)NSTEPS;
    angle = (ROTATE_ANGLE*DPI/360.0)*interpolate;
//     printf("i = %i, interpolate = %.3lg, angle = %.3lg\n", i, interpolate, angle);
    ca = cos(angle);
    sa = sin(angle);
    switch (VIEWPOINT_TRAJ)
    {
        case (VP_HORIZONTAL):
        {
            observer[0] = ca*observer_initial[0] - sa*observer_initial[1];
            observer[1] = sa*observer_initial[0] + ca*observer_initial[1];
            break;
        }
        case (VP_ORBIT):
        {
            observer[0] = ca*observer_initial[0] - sa*observer_initial[1]*ratio;
            observer[1] = ca*observer_initial[1] + sa*observer_initial[0]*ratio;
            observer[2] = ca*observer_initial[2];
            break;
        }
        case (VP_ORBIT2):
        {
            observer[0] = ca*observer_initial[0] - sa*observer_initial[1]*ratio;
            observer[1] = ca*observer_initial[1] + sa*observer_initial[0]*ratio;
            observer[2] = sa*zmax;
            break;
        }
        case (VP_POLAR):
        {
            rho = -sa*observer_initial[2] + ca*rho0;
            observer[0] = observer_initial[0]*rho/rho0;
            observer[1] = observer_initial[1]*rho/rho0;
            observer[2] = ca*observer_initial[2] + sa*rho0;
            break; 
        }
    }
    
    printf("Angle %.3lg, Observer position (%.3lg, %.3lg, %.3lg)\n", angle, observer[0], observer[1], observer[2]);
}


