/*********************/
/* animation part    */
/*********************/



void init_3d()		/* initialisation of window */
{
    glLineWidth(3);

    glClearColor(0.0, 0.0, 0.0, 1.0);
    glClear(GL_COLOR_BUFFER_BIT);

    glOrtho(XMIN, XMAX, YMIN, YMAX , -1.0, 1.0);
}

void xyz_to_xy(double x, double y, double z, double xy_out[2])
{
    int i;
    double s, t, xinter[3];
    static double n2, m2, d, sm2, sn2, v[3], h[2], plane_ratio = 0.5;
    static int first = 1;
    
    if (((first)&&(REPRESENTATION_3D == REP_PROJ_3D))||(reset_view))
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
        reset_view = 0;
//         printf("h = (%.3lg, %.3lg)\n", h[0], h[1]);
//         printf("v = (%.3lg, %.3lg, %.3lg)\n", v[0], v[1], v[2]);
    }
    
    switch (REPRESENTATION_3D) {
        case (REP_AXO_3D):
        {
            for (i=0; i<2; i++)
                xy_out[i] = x*u_3d[i] + y*v_3d[i] + z*w_3d[i];
            break;
        }
        case (REP_PROJ_3D):
        {
            if (z > ZMAX_FACTOR*n2) z = ZMAX_FACTOR*n2;
            z *= Z_SCALING_FACTOR;
            s = observer[0]*x + observer[1]*y + observer[2]*z;
            t = (d - s)/(n2 - s);
            xinter[0] = t*observer[0] + (1.0-t)*x;
            xinter[1] = t*observer[1] + (1.0-t)*y;
            xinter[2] = t*observer[2] + (1.0-t)*z;
            
            xy_out[0] = XSHIFT_3D + XY_SCALING_FACTOR*(xinter[0]*h[0] + xinter[1]*h[1]);
            xy_out[1] = YSHIFT_3D + XY_SCALING_FACTOR*(xinter[0]*v[0] + xinter[1]*v[1] + xinter[2]*v[2]);
            break;
        }
    }
}


void draw_vertex_ij(int i, int j)
{
    double xy[2];
    
    ij_to_xy(i, j, xy);
//     if (xy[1] > 0.0) printf("(i,j) = (%i,%i), (x,y) = (%.2lg,%.2lg)\n", i, j, xy[0], xy[1]);
    glVertex2d(xy[0], xy[1]);
}


void draw_vertex_xyz(double xy[2], double z)
{
    double xy_screen[2];
    
    xyz_to_xy(xy[0], xy[1], z, xy_screen);
    glVertex2d(xy_screen[0], xy_screen[1]);
}
    
void draw_vertex_xyz_shift(double xy[2], double z, int shiftx, int shifty)
{
    double xy_screen[2];
    
    xyz_to_xy(xy[0] + (double)shiftx*(XMAX - XMIN), xy[1] + (double)shifty*(YMAX - YMIN), z, xy_screen);
    glVertex2d(xy_screen[0], xy_screen[1]);
}
    
void draw_vertex_x_y_z(double x, double y, double z)
{
    double xy_screen[2];
    
    xyz_to_xy(x, y, z, xy_screen);
    glVertex2d(xy_screen[0], xy_screen[1]);
}

void draw_rectangle_3d(double x1, double y1, double x2, double y2)
{
    glBegin(GL_LINE_LOOP);
    draw_vertex_x_y_z(x1, y1, 0.0);
    draw_vertex_x_y_z(x2, y1, 0.0);
    draw_vertex_x_y_z(x2, y2, 0.0);
    draw_vertex_x_y_z(x1, y2, 0.0);    
    glEnd();    
}


void draw_rectangle_noscale(double x1, double y1, double x2, double y2)
{
    glBegin(GL_LINE_LOOP);
    glVertex2d(x1, y1);
    glVertex2d(x2, y1);
    glVertex2d(x2, y2);
    glVertex2d(x1, y2);
    glEnd();    
}

void draw_circle_3d(double x, double y, double r, int nseg)
{
    int i;
    double pos[2], alpha, dalpha;
    
    dalpha = DPI/(double)nseg;
    
    glBegin(GL_LINE_LOOP);
    for (i=0; i<=nseg; i++)
    {
        alpha = (double)i*dalpha;
        draw_vertex_x_y_z(x + r*cos(alpha), y + r*sin(alpha), 0.0);
    }
    glEnd();
}

void tvertex_lineto_3d(t_vertex z)
/* draws boundary segments of isospectral billiard */
{
    draw_vertex_x_y_z(z.x, z.y, 0.0);
}


void draw_billiard_3d(int fade, double fade_value)      /* draws the billiard boundary */
{
    double x0, x, y, x1, y1, dx, dy, phi, r = 0.01, pos[2], pos1[2], alpha, dphi, omega, z, l, width, a, b, c, ymax;
    int i, j, k, k1, k2, mr2;
    static int first = 1, nsides;

    if (BLACK) 
    {
        if (fade) glColor3f(fade_value, fade_value, fade_value);
        else glColor3f(1.0, 1.0, 1.0);
    }
    else 
    {
        if (fade) glColor3f(1.0 - fade_value, 1.0 - fade_value, 1.0 - fade_value);
        else glColor3f(0.0, 0.0, 0.0);
    }
    glLineWidth(BOUNDARY_WIDTH);

    glEnable(GL_LINE_SMOOTH);

    switch (B_DOMAIN) {
        case (D_RECTANGLE):
        {
            glBegin(GL_LINE_LOOP);
            draw_vertex_x_y_z(LAMBDA, -1.0, 0.0);
            draw_vertex_x_y_z(LAMBDA, 1.0, 0.0);
            draw_vertex_x_y_z(-LAMBDA, 1.0, 0.0);
            draw_vertex_x_y_z(-LAMBDA, -1.0, 0.0);
            glEnd();
            break;
        }
        case (D_ELLIPSE):
        {
            glBegin(GL_LINE_LOOP);
            for (i=0; i<=NSEG; i++)
            {
                phi = (double)i*DPI/(double)NSEG;
                x = 1.05*LAMBDA*cos(phi);
                y = 1.05*sin(phi);
                draw_vertex_x_y_z(x, y, 0.0);
            }
            glEnd ();

            /* draw foci */
//             if (FOCI)
//             {
//                 glColor3f(0.3, 0.3, 0.3);
//                 x0 = sqrt(LAMBDA*LAMBDA-1.0);
// 
//                 glLineWidth(2);
//                 glEnable(GL_LINE_SMOOTH);
//                 
//                 draw_circle(x0, 0.0, r, NSEG);
//                 draw_circle(-x0, 0.0, r, NSEG);
//             }
            break;
        }
        case (D_STADIUM):
        {
            glBegin(GL_LINE_LOOP);
            for (i=0; i<=NSEG; i++)
            {
                phi = -PID + (double)i*PI/(double)NSEG;
                x = 0.5*LAMBDA + cos(phi);
                y = sin(phi);
                draw_vertex_x_y_z(x, y, 0.0);
            }
            for (i=0; i<=NSEG; i++)
            {
                phi = PID + (double)i*PI/(double)NSEG;
                x = -0.5*LAMBDA + cos(phi);
                y = sin(phi);
                draw_vertex_x_y_z(x, y, 0.0);
            }
            glEnd();
            break;
        }
        case (D_SINAI):
        {
            draw_circle_3d(0.0, 0.0, LAMBDA, NSEG);
            break;
        }
        case (D_POLYGON):
        {
            omega = DPI/((double)NPOLY);
            glBegin(GL_LINE_LOOP);
            for (i=0; i<=NPOLY; i++)
            {
                x = 1.0075*cos(i*omega + APOLY*PID);
                y = 1.0075*sin(i*omega + APOLY*PID);
                draw_vertex_x_y_z(x, y, 0.0);
            }
            glEnd ();
            break;
        }
        case (D_YOUNG):
        {
            if (FILL_BILLIARD_COMPLEMENT)
            {
                if (fade) glColor3f(0.75*fade_value, 0.75*fade_value, 0.75*fade_value);
                else glColor3f(0.75, 0.75, 0.75);
                
                glBegin(GL_TRIANGLE_FAN);
                draw_vertex_x_y_z(-MU, YMIN, 0.0);
                draw_vertex_x_y_z(-MU, -LAMBDA-MU, 0.0);
                draw_vertex_x_y_z(MU, -LAMBDA-MU, 0.0);
                draw_vertex_x_y_z(MU, YMIN, 0.0);
                glEnd();
            
                glBegin(GL_TRIANGLE_FAN);
                draw_vertex_x_y_z(-MU, YMAX, 0.0);
                draw_vertex_x_y_z(-MU, LAMBDA+MU, 0.0);
                draw_vertex_x_y_z(MU, LAMBDA+MU, 0.0);
                draw_vertex_x_y_z(MU, YMAX, 0.0);
                glEnd();

                glBegin(GL_TRIANGLE_FAN);
                draw_vertex_x_y_z(-MU, -LAMBDA+MU, 0.0);
                draw_vertex_x_y_z(-MU, LAMBDA-MU, 0.0);
                draw_vertex_x_y_z(MU, LAMBDA-MU, 0.0);
                draw_vertex_x_y_z(MU, -LAMBDA+MU, 0.0);
                glEnd();
            }
            else
            {
                glBegin(GL_LINE_STRIP);
                draw_vertex_x_y_z(-MU, YMIN, 0.0);
                draw_vertex_x_y_z(-MU, -LAMBDA-MU, 0.0);
                draw_vertex_x_y_z(MU, -LAMBDA-MU, 0.0);
                draw_vertex_x_y_z(MU, YMIN, 0.0);
                glEnd();
            
                glBegin(GL_LINE_STRIP);
                draw_vertex_x_y_z(-MU, YMAX, 0.0);
                draw_vertex_x_y_z(-MU, LAMBDA+MU, 0.0);
                draw_vertex_x_y_z(MU, LAMBDA+MU, 0.0);
                draw_vertex_x_y_z(MU, YMAX, 0.0);
                glEnd();

                glBegin(GL_LINE_LOOP);
                draw_vertex_x_y_z(-MU, -LAMBDA+MU, 0.0);
                draw_vertex_x_y_z(-MU, LAMBDA-MU, 0.0);
                draw_vertex_x_y_z(MU, LAMBDA-MU, 0.0);
                draw_vertex_x_y_z(MU, -LAMBDA+MU, 0.0);
                glEnd();
            }
        }
        case (D_GRATING):
        {
            k1 = -(int)(-YMIN/LAMBDA);
            k2 = (int)(YMAX/LAMBDA);
            for (i=k1; i<= k2; i++)
            {
                z = (double)i*LAMBDA;
                draw_circle_3d(0.0, z, MU, NSEG);
            }
            break;
        }
        case (D_EHRENFEST):
        {
            alpha = asin(MU/LAMBDA);
            x0 = 1.0 - sqrt(LAMBDA*LAMBDA - MU*MU);
            dphi = 2.0*(PI-alpha)/((double)NSEG); 
            glBegin(GL_LINE_LOOP);
            for (i=0; i<=NSEG; i++)
            {
                phi = -PI + alpha + (double)i*dphi;
                x = 1.0 + (LAMBDA + 0.01)*cos(phi);
                y = (LAMBDA + 0.01)*sin(phi);
                draw_vertex_x_y_z(x, y, 0.0);
            }
            phi = PI - alpha;
            x = 1.0 + (LAMBDA + 0.01)*cos(phi);
            y = 0.01 + (LAMBDA + 0.01)*sin(phi);
            draw_vertex_x_y_z(x, y, 0.0);
            phi = alpha;
            x = -1.0 + (LAMBDA + 0.01)*cos(phi);
            y = 0.01 + (LAMBDA + 0.01)*sin(phi);
            draw_vertex_x_y_z(x, y, 0.0);
            for (i=0; i<=NSEG; i++)
            {
                phi = alpha + (double)i*dphi;
                x = -1.0 + (LAMBDA + 0.01)*cos(phi);
                y = (LAMBDA + 0.01)*sin(phi);
                draw_vertex_x_y_z(x, y, 0.0);
            }
            glEnd ();
            break;
        }
        case (D_TWO_PARABOLAS):
        {
            dy = 3.0*MU/(double)NSEG;
            width = 0.25*MU;
            if (width > 0.2) width = 0.2;
            glBegin(GL_LINE_LOOP);
            for (i = 0; i < NSEG+1; i++) 
            {
                y = -1.5*MU + dy*(double)i;
                x = 0.25*y*y/MU - MU - LAMBDA;
                draw_vertex_x_y_z(x, y, 0.0);
            }
            for (i = 0; i < NSEG+1; i++) 
            {
                y = 1.5*MU - dy*(double)i;
                x = 0.25*y*y/MU - (MU + width) - LAMBDA;
                draw_vertex_x_y_z(x, y, 0.0);
            }
            glEnd ();
            
            glBegin(GL_LINE_LOOP);
            for (i = 0; i < NSEG+1; i++) 
            {
                y = -1.5*MU + dy*(double)i;
                x = LAMBDA + MU - 0.25*y*y/MU;
                draw_vertex_x_y_z(x, y, 0.0);
            }
            for (i = 0; i < NSEG+1; i++) 
            {
                y = 1.5*MU - dy*(double)i;
                x = LAMBDA + (MU + width) - 0.25*y*y/MU;
                draw_vertex_x_y_z(x, y, 0.0);
            }
            glEnd ();
            
//             if (FOCI)
//             {
//                 glColor3f(0.3, 0.3, 0.3);
//                 draw_circle(-LAMBDA, 0.0, r, NSEG);
//                 draw_circle(LAMBDA, 0.0, r, NSEG);
//             }

            break;
        }
        case (D_POLY_PARABOLAS):
        {
            omega = PI/((double)NPOLY);
            a = 0.25/MU;
            b = 1.0/tan(omega);
            c = LAMBDA + MU;
            ymax = (-b + sqrt(b*b + 4.0*a*c))/(2.0*a);
            dy = 2.0*ymax/(double)NSEG; 
            
            glBegin(GL_LINE_LOOP);
            for (k=0; k<NPOLY; k++)  
            {
                alpha = APOLY*PID + (2.0*(double)k+1.0)*omega;
                for (i = 0; i < NSEG+1; i++) 
                {
                    y1 = -ymax + dy*(double)i;
                    x1 = MU + LAMBDA - 0.25*y1*y1/MU;
                    x = x1*cos(alpha) - y1*sin(alpha);
                    y = x1*sin(alpha) + y1*cos(alpha);
                    draw_vertex_x_y_z(x, y, 0.0);
                }
            }
            glEnd ();
            
//             if (FOCI)
//             {
//                 glColor3f(0.3, 0.3, 0.3);
//                 for (k=0; k<NPOLY; k++) 
//                 {
//                     alpha = APOLY*PID + (2.0*(double)k+1.0)*omega;
//                     draw_circle(LAMBDA*cos(alpha), LAMBDA*sin(alpha), r, NSEG);
//                 }
//             }
            
            break;
        }
        case (D_MENGER):
        {
            glLineWidth(3);
//             draw_rectangle(XMIN, -1.0, XMAX, 1.0);
            
            /* level 1 */
            if (MDEPTH > 0)
            {
                glLineWidth(2);
                x = 1.0/((double)MRATIO);
                draw_rectangle_3d(x, x, -x, -x);
            }
            
            /* level 2 */
            if (MDEPTH > 1)
            {
                glLineWidth(1);
                mr2 = MRATIO*MRATIO;
                l = 2.0/((double)mr2);
                
                for (i=0; i<MRATIO; i++)
                    for (j=0; j<MRATIO; j++)
                        if ((i!=MRATIO/2)||(j!=MRATIO/2))
                        {
                            x = -1.0 - 0.5*l + 2.0*((double)i + 0.5)/((double)MRATIO);
                            y = -1.0 - 0.5*l + 2.0*((double)j + 0.5)/((double)MRATIO);
                            draw_rectangle_3d(x, y, x+l, y+l);
                        }
            }
            
            /* level 3 */
            if (MDEPTH > 2)
            {
                glLineWidth(1);
                l = 2.0/((double)(mr2*MRATIO));
                
                for (i=0; i<mr2; i++)
                    for (j=0; j<mr2; j++)
                        if ( (((i%MRATIO!=MRATIO/2))||(j%MRATIO!=MRATIO/2)) && (((i/MRATIO!=MRATIO/2))||(j/MRATIO!=MRATIO/2)) )
                        {
                            x = -1.0 - 0.5*l + 2.0*((double)i + 0.5)/((double)mr2);
                            y = -1.0 - 0.5*l + 2.0*((double)j + 0.5)/((double)mr2);
                            draw_rectangle_3d(x, y, x+l, y+l);
                        }
            }
            
            break;
        }        
        case (D_PENROSE):
        {
            c = sqrt(LAMBDA*LAMBDA - (1.0 - MU)*(1.0 - MU));
            width = 0.1*MU;
            x1 = vabs(x);
            y1 = vabs(y);
            dphi = PI/(double)NSEG;
            
            glBegin(GL_LINE_LOOP);
            /* upper half ellipse */
            for (i=0; i<=NSEG; i++)
            {
                phi = (double)i*dphi;
                x = LAMBDA*cos(phi);
                y = MU + (1.0-MU)*sin(phi);
                draw_vertex_x_y_z(x, y, 0.0);
            }
            
            /* straight parts */
            draw_vertex_x_y_z(-LAMBDA, width, 0.0);
            draw_vertex_x_y_z(-c, width, 0.0);
            draw_vertex_x_y_z(-c, MU, 0.0);
                        
            /* left half ellipse */
            for (i=0; i<=NSEG; i++)
            {
                phi = (double)i*dphi;
                x = -c + 0.5*MU*sin(phi);
                y = MU*cos(phi);
                draw_vertex_x_y_z(x, y, 0.0);
            }
            
            /* straight parts */
            draw_vertex_x_y_z(-c, -width, 0.0);
            draw_vertex_x_y_z(-LAMBDA, -width, 0.0);
            draw_vertex_x_y_z(-LAMBDA, -MU, 0.0);
            
            /* lower half ellipse */
            for (i=0; i<=NSEG; i++)
            {
                phi = (double)i*dphi;
                x = -LAMBDA*cos(phi);
                y = -MU - (1.0-MU)*sin(phi);
                draw_vertex_x_y_z(x, y, 0.0);
            }
            
            /* straight parts */
            draw_vertex_x_y_z(LAMBDA, -width, 0.0);
            draw_vertex_x_y_z(c, -width, 0.0);
            draw_vertex_x_y_z(c, -MU, 0.0);
            
            /* right half ellipse */
            for (i=0; i<=NSEG; i++)
            {
                phi = (double)i*dphi;
                x = c - 0.5*MU*sin(phi);
                y = -MU*cos(phi);
                draw_vertex_x_y_z(x, y, 0.0);
            }
            
            /* straight parts */
            draw_vertex_x_y_z(c, width, 0.0);
            draw_vertex_x_y_z(LAMBDA, width, 0.0);
            draw_vertex_x_y_z(LAMBDA, MU, 0.0);            
            glEnd ();
            break; 
        }
        case (D_TOKA_PRIME):
        {
            glBegin(GL_LINE_LOOP);
            tvertex_lineto_3d(polyline[0]);
            for (i=4; i<43; i++) tvertex_lineto_3d(polyline[i]);
            tvertex_lineto_3d(polyline[3]);
            tvertex_lineto_3d(polyline[2]);
            tvertex_lineto_3d(polyline[1]);

            tvertex_lineto_3d(polyline[44]);
            tvertex_lineto_3d(polyline[45]);
            for (i=84; i>45; i--) tvertex_lineto_3d(polyline[i]);
            glEnd();
            
            /* inner lines */ 
//             glLineWidth(BOUNDARY_WIDTH/2);
            glLineWidth(1);
            glColor3f(0.75, 0.75, 0.75);
            glBegin(GL_LINE_STRIP);
            tvertex_lineto_3d(polyline[0]);
            tvertex_lineto_3d(polyline[1]);
            tvertex_lineto_3d(polyline[2]);
            tvertex_lineto_3d(polyline[0]);
            tvertex_lineto_3d(polyline[3]);
            tvertex_lineto_3d(polyline[4]);
            glEnd();
            
            glBegin(GL_LINE_STRIP);
            tvertex_lineto_3d(polyline[0]);
            tvertex_lineto_3d(polyline[44]);
            tvertex_lineto_3d(polyline[45]);
            tvertex_lineto_3d(polyline[0]);
            tvertex_lineto_3d(polyline[46]);
            tvertex_lineto_3d(polyline[45]);
            glEnd();
            
            for (i=3; i<43; i++)
            {
                glBegin(GL_LINE_STRIP);
                tvertex_lineto_3d(polyline[i]);
                tvertex_lineto_3d(polyline[43]);
                glEnd();
                glBegin(GL_LINE_STRIP);
                tvertex_lineto_3d(polyline[i+42]);
                tvertex_lineto_3d(polyline[85]);
                glEnd();
            }
            
            break;
        }
        case (D_ISOSPECTRAL):
        {
            /* 1st triangle */
            glBegin(GL_LINE_LOOP);
            tvertex_lineto_3d(polyline[0]);
            tvertex_lineto_3d(polyline[4]);
            tvertex_lineto_3d(polyline[7]);
            tvertex_lineto_3d(polyline[1]);
            tvertex_lineto_3d(polyline[5]);
            tvertex_lineto_3d(polyline[8]);
            tvertex_lineto_3d(polyline[2]);
            tvertex_lineto_3d(polyline[3]);
            tvertex_lineto_3d(polyline[6]);
            glEnd();
            
            /* inner lines */
            glBegin(GL_LINE_LOOP);
            tvertex_lineto_3d(polyline[0]);
            tvertex_lineto_3d(polyline[1]);
            tvertex_lineto_3d(polyline[2]);
            tvertex_lineto_3d(polyline[0]);
            tvertex_lineto_3d(polyline[3]);
            tvertex_lineto_3d(polyline[2]);
            tvertex_lineto_3d(polyline[5]);
            tvertex_lineto_3d(polyline[1]);
            tvertex_lineto_3d(polyline[4]);
            glEnd();
            
            /* 2nd triangle */
            glBegin(GL_LINE_LOOP);
            tvertex_lineto_3d( polyline[9]);
            tvertex_lineto_3d(polyline[16]);
            tvertex_lineto_3d(polyline[13]);
            tvertex_lineto_3d(polyline[10]);
            tvertex_lineto_3d(polyline[17]);
            tvertex_lineto_3d(polyline[14]);
            tvertex_lineto_3d(polyline[11]);
            tvertex_lineto_3d(polyline[15]);
            tvertex_lineto_3d(polyline[12]);
            glEnd();
            
            /* inner lines */
            glBegin(GL_LINE_LOOP);
            tvertex_lineto_3d( polyline[9]);
            tvertex_lineto_3d(polyline[10]);
            tvertex_lineto_3d(polyline[11]);
            tvertex_lineto_3d( polyline[9]);
            tvertex_lineto_3d(polyline[13]);
            tvertex_lineto_3d(polyline[10]);
            tvertex_lineto_3d(polyline[14]);
            tvertex_lineto_3d(polyline[11]);
            tvertex_lineto_3d(polyline[12]);
            glEnd();
            break;
        }
        case (D_HOMOPHONIC):
        {
            /* 1st triangle */
            glBegin(GL_LINE_LOOP);
            tvertex_lineto_3d(polyline[1]);
            tvertex_lineto_3d(polyline[3]);
            tvertex_lineto_3d(polyline[4]);
            tvertex_lineto_3d(polyline[5]);
            tvertex_lineto_3d(polyline[6]);
            tvertex_lineto_3d(polyline[8]);
            tvertex_lineto_3d(polyline[9]);
            tvertex_lineto_3d(polyline[10]);
            tvertex_lineto_3d(polyline[12]);
            tvertex_lineto_3d(polyline[13]);
            tvertex_lineto_3d(polyline[15]);
            tvertex_lineto_3d(polyline[16]);
            tvertex_lineto_3d(polyline[17]);
            tvertex_lineto_3d(polyline[18]);
            tvertex_lineto_3d(polyline[20]);
            glEnd();
            
            /* inner lines */
            glLineWidth(BOUNDARY_WIDTH/2);
            glBegin(GL_LINE_STRIP);
            tvertex_lineto_3d(polyline[9]);
            tvertex_lineto_3d(polyline[1]);
            tvertex_lineto_3d(polyline[2]);
            tvertex_lineto_3d(polyline[5]);
            tvertex_lineto_3d(polyline[7]);
            tvertex_lineto_3d(polyline[2]);
            tvertex_lineto_3d(polyline[8]);
            tvertex_lineto_3d(polyline[21]);
            tvertex_lineto_3d(polyline[10]);
            tvertex_lineto_3d(polyline[2]);
            tvertex_lineto_3d(polyline[21]);
            tvertex_lineto_3d(polyline[11]);
            tvertex_lineto_3d(polyline[13]);
            tvertex_lineto_3d(polyline[21]);
            tvertex_lineto_3d(polyline[14]);
            tvertex_lineto_3d(polyline[20]);
            tvertex_lineto_3d(polyline[15]);
            tvertex_lineto_3d(polyline[19]);
            tvertex_lineto_3d(polyline[16]);
            tvertex_lineto_3d(polyline[18]);
            glEnd();
            
            /* 2nd triangle */
            glLineWidth(BOUNDARY_WIDTH);
            glBegin(GL_LINE_LOOP);
            tvertex_lineto_3d(polyline[22+10]);
            tvertex_lineto_3d(polyline[22+16]);
            tvertex_lineto_3d(polyline[22+17]);
            tvertex_lineto_3d(polyline[22+18]);
            tvertex_lineto_3d(polyline[22+12]);
            tvertex_lineto_3d(polyline[22+13]);
            tvertex_lineto_3d(polyline[22+15]);
            tvertex_lineto_3d(polyline[22+19]);
            tvertex_lineto_3d(polyline[22+20]);
            tvertex_lineto_3d(polyline[22+1]);
            tvertex_lineto_3d(polyline[22+4]);
            tvertex_lineto_3d(polyline[22+5]);
            tvertex_lineto_3d(polyline[22+7]);
            tvertex_lineto_3d(polyline[22+8]);
            tvertex_lineto_3d(polyline[22+9]);
            glEnd();
            
            /* inner lines */
            glLineWidth(BOUNDARY_WIDTH/2);
            glBegin(GL_LINE_STRIP);
            tvertex_lineto_3d(polyline[22+2]);
            tvertex_lineto_3d(polyline[22+6]);
            tvertex_lineto_3d(polyline[22+8]);
            tvertex_lineto_3d(polyline[22+2]);
            tvertex_lineto_3d(polyline[22+5]);
            tvertex_lineto_3d(polyline[22+3]);
            tvertex_lineto_3d(polyline[22+2]);
            tvertex_lineto_3d(polyline[22+1]);
            tvertex_lineto_3d(polyline[22+0]);
            tvertex_lineto_3d(polyline[22+21]);
            tvertex_lineto_3d(polyline[22+18]);
            tvertex_lineto_3d(polyline[22+16]);
            tvertex_lineto_3d(polyline[22+13]);
            tvertex_lineto_3d(polyline[22+21]);
            tvertex_lineto_3d(polyline[22+10]);
            tvertex_lineto_3d(polyline[22+12]);
            tvertex_lineto_3d(polyline[22+21]);
            tvertex_lineto_3d(polyline[22+14]);
            tvertex_lineto_3d(polyline[22+20]);
            tvertex_lineto_3d(polyline[22+15]);
            glEnd();
            break;
        }
        default:
        {
            break;
        }   
    }
}

void draw_polyline_visible(int j, int k, double margin)
/* hack to draw the billiard boundary in front of the wave */
/* only parts of the boundary having a small enough angle with the observer vector are drawn */
{
    double x, y, x1, y1, length, length1;
    static int first = 1;
    static double olength;
    
    if (first)
    {
        olength = module2(observer[0], observer[1]);
        first = 0;
    }
    
    x = polyline[j].x;
    y = polyline[j].y;
    x1 = polyline[k].x;
    y1 = polyline[k].y;
    length = module2(x,y);
    length1 = module2(x1,y1);
    if ((x*observer[0] + y*observer[1] > margin*length*olength)&&(x1*observer[0] + y1*observer[1] > margin*length1*olength))
    {
        glBegin(GL_LINE_STRIP);
        tvertex_lineto_3d(polyline[j]);
        tvertex_lineto_3d(polyline[k]);
        glEnd();
    }
}

void draw_vertex_visible(double x, double y, double margin)
/* hack to draw the billiard boundary in front of the wave */
/* only parts of the boundary having a small enough angle with the observer vector are drawn */
{
    static int first = 1;
    static double olength;
    
    if (first)
    {
        olength = module2(observer[0], observer[1]);
        first = 0;
    }
    
    if (x*observer[0] + y*observer[1] > margin*olength*module2(x,y))
        draw_vertex_x_y_z(x, y, 0.0);
    else 
    {
        glEnd();
        glBegin(GL_LINE_STRIP);
    }
}

void draw_billiard_3d_front(int fade, double fade_value)      
/* hack to draw the billiard boundary in front of the wave */
/* only parts of the boundary having a small enough angle with the observer vector are drawn */
{
    double x0, x, y, x1, y1, dx, dy, phi, r = 0.01, pos[2], pos1[2], alpha, dphi, omega, z, l, width, a, b, c, ymax, length, length1;
    int i, j, k, k1, k2, mr2;
    static int first = 1, nsides;
    static double olength;
    
    if (first)
    {
        olength = module2(observer[0], observer[1]);
        first = 0;
    }

    if (BLACK) 
    {
        if (fade) glColor3f(fade_value, fade_value, fade_value);
        else glColor3f(1.0, 1.0, 1.0);
    }
    else 
    {
        if (fade) glColor3f(1.0 - fade_value, 1.0 - fade_value, 1.0 - fade_value);
        else glColor3f(0.0, 0.0, 0.0);
    }
    glLineWidth(BOUNDARY_WIDTH);

    glEnable(GL_LINE_SMOOTH);

    switch (B_DOMAIN) {
        case (D_RECTANGLE):
        {
            glBegin(GL_LINE_STRIP);
            draw_vertex_x_y_z(LAMBDA, -1.0, 0.0);
            draw_vertex_x_y_z(LAMBDA, 1.0, 0.0);
            draw_vertex_x_y_z(-LAMBDA, 1.0, 0.0);
            glEnd();
            break;
        }
        case (D_STADIUM):
        {
            glBegin(GL_LINE_STRIP);
            for (i=0; i<=NSEG; i++)
            {
                phi = -PID + (double)i*PI/(double)NSEG;
                x = 0.5*LAMBDA + cos(phi);
                y = sin(phi);
                draw_vertex_visible(x, y, 0.0);
            }
            for (i=0; i<=NSEG; i++)
            {
                phi = PID + (double)i*PI/(double)NSEG;
                x = -0.5*LAMBDA + cos(phi);
                y = sin(phi);
                draw_vertex_visible(x, y, 0.0);
            }
            glEnd();
            break;
        }
        case (D_POLYGON):
        {
            omega = DPI/((double)NPOLY);
            glBegin(GL_LINE_STRIP);
            for (i=0; i<=NPOLY; i++)
            {
                x = cos(i*omega + APOLY*PID);
                y = sin(i*omega + APOLY*PID);
                draw_vertex_visible(x, y, 0.0);
            }
            glEnd ();
            break;
        }
        case (D_YOUNG):
        {
            glBegin(GL_LINE_STRIP);
            draw_vertex_x_y_z(-MU, YMIN, 0.0);
            draw_vertex_x_y_z(-MU, -LAMBDA-MU, 0.0);
            glEnd();
            
            glBegin(GL_LINE_STRIP);
            draw_vertex_x_y_z(-MU, YMAX, 0.0);
            draw_vertex_x_y_z(-MU, LAMBDA+MU, 0.0);
            glEnd();

            glBegin(GL_LINE_LOOP);
            draw_vertex_x_y_z(-MU, -LAMBDA+MU, 0.0);
            draw_vertex_x_y_z(-MU, LAMBDA-MU, 0.0);
            glEnd();
            break;
        }
        case (D_TOKA_PRIME):
        {
            draw_polyline_visible(0, 4, 0.2);
            for (i=4; i<41; i++) draw_polyline_visible(i, i+1, -0.1);
//             draw_polyline_visible(42, 3, 0.2);
            for (i=84; i>46; i--) draw_polyline_visible(i, i-1, -0.1);
            draw_polyline_visible(46, 0, 0.2);
            break;
        }
        case (D_STAR):
        {
            glLineWidth(BOUNDARY_WIDTH);
            for (i=0; i<npolyline; i++)
                draw_polyline_visible(i%npolyline, (i+1)%npolyline, 0.2);
            break;
        }
        default:
        {
            break;
        }   
    }
}

void compute_energy_field(double phi[NX*NY], double psi[NX*NY], short int xy_in[NX*NY], double energies[NX*NY])
/* computes cosine of angle between normal vector and vector light */
{
    int i, j;
    
    #pragma omp parallel for private(i,j)
    for (i=0; i<NX; i++)
        for (j=0; j<NY; j++)
        {
            if ((TWOSPEEDS)||(xy_in[i*NY+j])) energies[i*NY+j] = PLOT_SCALE_ENERGY*compute_energy_mod(phi, psi, xy_in, i, j);
            else energies[i*NY+j] = 0.0;
        }
}


void compute_phase_field(double phi[NX*NY], double psi[NX*NY], short int xy_in[NX*NY], double phase[NX*NY])
/* computes cosine of angle between normal vector and vector light */
{
    int i, j;
    double angle;
    
    #pragma omp parallel for private(i,j,angle)
    for (i=0; i<NX; i++)
        for (j=0; j<NY; j++)
        {
            if ((TWOSPEEDS)||(xy_in[i*NY+j])) 
            {
                angle = compute_phase(phi, psi, xy_in, i, j)/DPI + PHASE_SHIFT;
                if (angle >= 1.0) angle -= 1.0;
                else if (angle < 0.0) angle += 1.0;
                phase[i*NY+j] = angle;
            }
            else phase[i*NY+j] = 0.0;
        }
}


// void compute_log_energy_field(double *phi[NX], double *psi[NX], short int xy_in[NX*NY], double *energies[NX])
// /* computes cosine of angle between normal vector and vector light */
// {
//     int i, j;
//     
//     for (i=0; i<NX; i++)
//         for (j=0; j<NY; j++)
//         {
//             if ((TWOSPEEDS)||(xy_in[i*NY+j])) energies[i][j] = PLOT_SCALE_ENERGY*compute_energy(phi, psi, xy_in, i, j);
//             else energies[i][j] = 0.0;
//         }
// }


void compute_light_angle(double field[NX*NY], short int xy_in[NX*NY], double cos_angle[NX*NY])
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
                gradx = (field[(i+1)*NY+j] - field[(i-1)*NY+j])/dx;
                grady = (field[i*NY+j+1] - field[i*NY+j-1])/dy;
                norm = sqrt(1.0 + gradx*gradx + grady*grady);
                pscal = -gradx*light[0] - grady*light[1] + 1.0;
                
                cos_angle[i*NY+j] = pscal/norm;
            }
        }
}


void energy_color_scheme(int palette, double energy, double rgb[3])
{
    if (COLOR_PALETTE >= COL_TURBO) color_scheme_asym_palette(COLOR_SCHEME, palette, energy, 1.0, 0, rgb);
    else color_scheme_palette(COLOR_SCHEME, palette, energy, 1.0, 0, rgb);
}

void log_energy_color_scheme(int palette, double energy, double rgb[3])
{
    color_scheme_palette(COLOR_SCHEME, palette, LOG_SHIFT + LOG_SCALE*log(energy), 1.0, 0, rgb);
//     if (energy > 1.0e-8) printf("energy = %.3lg, log energy = %.3lg\n", energy, LOG_SHIFT + LOG_SCALE*log(energy));
}

void phase_color_scheme(int palette, double phase, double rgb[3])
{
//     color_scheme_palette(C_ONEDIM_LINEAR, palette, phase/DPI, 1.0, 0, rgb);
    amp_to_rgb_palette(phase, rgb, palette);
}


void compute_interpolated_colors(int i, int j, double field[NX*NY], double palette, int plot, 
                                 double *z_sw, double *z_se, double *z_nw, double *z_ne, double *z_mid, 
                                 double rgb_e[3], double rgb_w[3], double rgb_n[3], double rgb_s[3])
{
    double zw, ze, zn, zs;
    
    *z_sw = field[i*NY+j];
    *z_se = field[(i+1)*NY+j];
    *z_nw = field[i*NY+j+1];
    *z_ne = field[(i+1)*NY+j+1];
                            
    *z_mid = 0.25*(*z_sw + *z_se + *z_nw + *z_ne);
                            
    zw = (*z_sw + *z_nw + *z_mid)/3.0;
    ze = (*z_se + *z_ne + *z_mid)/3.0;
    zs = (*z_sw + *z_se + *z_mid)/3.0;
    zn = (*z_nw + *z_ne + *z_mid)/3.0;
    
    if (plot == P_3D_ENERGY)
    {
        energy_color_scheme(palette, VSCALE_ENERGY*ze, rgb_e);
        energy_color_scheme(palette, VSCALE_ENERGY*zw, rgb_w);
        energy_color_scheme(palette, VSCALE_ENERGY*zn, rgb_n);
        energy_color_scheme(palette, VSCALE_ENERGY*zs, rgb_s);         
    }
    else if (plot == P_3D_LOG_ENERGY)
    {
        log_energy_color_scheme(palette, ze, rgb_e);
        log_energy_color_scheme(palette, zw, rgb_w);
        log_energy_color_scheme(palette, zn, rgb_n);
        log_energy_color_scheme(palette, zs, rgb_s);         
    }
    else
    {
        color_scheme_palette(COLOR_SCHEME, palette, VSCALE_AMPLITUDE*ze, 1.0, 0, rgb_e);
        color_scheme_palette(COLOR_SCHEME, palette, VSCALE_AMPLITUDE*zw, 1.0, 0, rgb_w);
        color_scheme_palette(COLOR_SCHEME, palette, VSCALE_AMPLITUDE*zn, 1.0, 0, rgb_n);
        color_scheme_palette(COLOR_SCHEME, palette, VSCALE_AMPLITUDE*zs, 1.0, 0, rgb_s); 
    }
}

void compute_interpolated_colors_new(int i, int j, double zfield[NX*NY], double colorfield[NX*NY], double palette, int cplot, 
                                 double *z_sw, double *z_se, double *z_nw, double *z_ne, double *z_mid, 
                                 double rgb_e[3], double rgb_w[3], double rgb_n[3], double rgb_s[3])
{
    double zw, ze, zn, zs, c_sw, c_se, c_nw, c_ne, c_mid;
    
    *z_sw = zfield[i*NY+j];
    *z_se = zfield[(i+1)*NY+j];
    *z_nw = zfield[i*NY+j+1];
    *z_ne = zfield[(i+1)*NY+j+1];
                            
    *z_mid = 0.25*(*z_sw + *z_se + *z_nw + *z_ne);
                            
    c_sw = colorfield[i*NY+j];
    c_se = colorfield[(i+1)*NY+j];
    c_nw = colorfield[i*NY+j+1];
    c_ne = colorfield[(i+1)*NY+j+1];
                            
    c_mid = 0.25*(c_sw + c_se + c_nw + c_ne);

    zw = (c_sw + c_nw + c_mid)/3.0;
    ze = (c_se + c_ne + c_mid)/3.0;
    zs = (c_sw + c_se + c_mid)/3.0;
    zn = (c_nw + c_ne + c_mid)/3.0;

    switch (cplot){
        case (P_3D_ENERGY):
        {
            energy_color_scheme(palette, VSCALE_ENERGY*ze, rgb_e);
            energy_color_scheme(palette, VSCALE_ENERGY*zw, rgb_w);
            energy_color_scheme(palette, VSCALE_ENERGY*zn, rgb_n);
            energy_color_scheme(palette, VSCALE_ENERGY*zs, rgb_s);
            break;
        }
        case (P_3D_LOG_ENERGY):
        {
            log_energy_color_scheme(palette, ze, rgb_e);
            log_energy_color_scheme(palette, zw, rgb_w);
            log_energy_color_scheme(palette, zn, rgb_n);
            log_energy_color_scheme(palette, zs, rgb_s);         
            break;
        }
        case (P_3D_PHASE):
        {
            phase_color_scheme(palette, ze, rgb_e);
            phase_color_scheme(palette, zw, rgb_w);
            phase_color_scheme(palette, zn, rgb_n);
            phase_color_scheme(palette, zs, rgb_s);         
            break;
        }
        default:
        {
//         hsl_to_rgb_palette(VSCALE_AMPLITUDE*ze, 0.9, 0.5, rgb_e, palette);
//         hsl_to_rgb_palette(VSCALE_AMPLITUDE*zw, 0.9, 0.5, rgb_w, palette);
//         hsl_to_rgb_palette(VSCALE_AMPLITUDE*zn, 0.9, 0.5, rgb_n, palette);
//         hsl_to_rgb_palette(VSCALE_AMPLITUDE*zs, 0.9, 0.5, rgb_s, palette);
            color_scheme_palette(COLOR_SCHEME, palette, VSCALE_AMPLITUDE*ze, 1.0, 0, rgb_e);
            color_scheme_palette(COLOR_SCHEME, palette, VSCALE_AMPLITUDE*zw, 1.0, 0, rgb_w);
            color_scheme_palette(COLOR_SCHEME, palette, VSCALE_AMPLITUDE*zn, 1.0, 0, rgb_n);
            color_scheme_palette(COLOR_SCHEME, palette, VSCALE_AMPLITUDE*zs, 1.0, 0, rgb_s); 
        }
    }
}

void compute_interpolated_colors_rde_temp(int i, int j, double zfield[NX*NY], double colorfield[NX*NY], double palette, int cplot, 
                                 double *z_sw, double *z_se, double *z_nw, double *z_ne, double *z_mid, 
                                 double rgb_e[3], double rgb_w[3], double rgb_n[3], double rgb_s[3])
{
    double zw, ze, zn, zs, c_sw, c_se, c_nw, c_ne, c_mid;
    
    *z_sw = zfield[i*NY+j];
    *z_se = zfield[(i+1)*NY+j];
    *z_nw = zfield[i*NY+j+1];
    *z_ne = zfield[(i+1)*NY+j+1];
                            
    *z_mid = 0.25*(*z_sw + *z_se + *z_nw + *z_ne);
                            
    c_sw = colorfield[i*NY+j];
    c_se = colorfield[(i+1)*NY+j];
    c_nw = colorfield[i*NY+j+1];
    c_ne = colorfield[(i+1)*NY+j+1];
                            
    c_mid = 0.25*(c_sw + c_se + c_nw + c_ne);

    zw = (c_sw + c_nw + c_mid)/3.0;
    ze = (c_se + c_ne + c_mid)/3.0;
    zs = (c_sw + c_se + c_mid)/3.0;
    zn = (c_nw + c_ne + c_mid)/3.0;

    switch (cplot){
        case (P_3D_ENERGY):
        {
            energy_color_scheme(palette, VSCALE_ENERGY*ze, rgb_e);
            energy_color_scheme(palette, VSCALE_ENERGY*zw, rgb_w);
            energy_color_scheme(palette, VSCALE_ENERGY*zn, rgb_n);
            energy_color_scheme(palette, VSCALE_ENERGY*zs, rgb_s);
            break;
        }
        case (P_3D_LOG_ENERGY):
        {
            log_energy_color_scheme(palette, ze, rgb_e);
            log_energy_color_scheme(palette, zw, rgb_w);
            log_energy_color_scheme(palette, zn, rgb_n);
            log_energy_color_scheme(palette, zs, rgb_s);         
            break;
        }
        case (P_3D_PHASE):
        {
            phase_color_scheme(palette, ze, rgb_e);
            phase_color_scheme(palette, zw, rgb_w);
            phase_color_scheme(palette, zn, rgb_n);
            phase_color_scheme(palette, zs, rgb_s);         
            break;
        }
        default:
        {
        hsl_to_rgb_palette(VSCALE_AMPLITUDE*ze, 0.9, 0.5, rgb_e, palette);
        hsl_to_rgb_palette(VSCALE_AMPLITUDE*zw, 0.9, 0.5, rgb_w, palette);
        hsl_to_rgb_palette(VSCALE_AMPLITUDE*zn, 0.9, 0.5, rgb_n, palette);
        hsl_to_rgb_palette(VSCALE_AMPLITUDE*zs, 0.9, 0.5, rgb_s, palette);
//             color_scheme_palette(COLOR_SCHEME, palette, VSCALE_AMPLITUDE*ze, 1.0, 0, rgb_e);
//             color_scheme_palette(COLOR_SCHEME, palette, VSCALE_AMPLITUDE*zw, 1.0, 0, rgb_w);
//             color_scheme_palette(COLOR_SCHEME, palette, VSCALE_AMPLITUDE*zn, 1.0, 0, rgb_n);
//             color_scheme_palette(COLOR_SCHEME, palette, VSCALE_AMPLITUDE*zs, 1.0, 0, rgb_s); 
        }
    }
}

void draw_wave_3d(double phi[NX*NY], double psi[NX*NY], short int xy_in[NX*NY], int plot, int cplot, int palette, int fade, double fade_value)
{
    int i, j, k, l, draw = 1;
    double xy[2], xy_screen[2], rgb[3], pos[2], ca, rgb_e[3], rgb_w[3], rgb_n[3], rgb_s[3]; 
    double z_sw, z_se, z_nw, z_ne, z_mid, zw, ze, zn, zs, min = 1000.0, max = 0.0;
    double xy_sw[2], xy_se[2], xy_nw[2], xy_ne[2], xy_mid[2];
    double energy;
    double *cos_angle, *energies, *phases;
    
    blank();
    draw_billiard_3d(fade, fade_value);
    cos_angle = (double *)malloc(NX*NY*sizeof(double));
    energies = (double *)malloc(NX*NY*sizeof(double));
    phases = (double *)malloc(NX*NY*sizeof(double));
    
    
    if ((plot == P_3D_ENERGY)||(plot == P_3D_LOG_ENERGY))
//     if (plot == P_3D_ENERGY)
    {
        compute_energy_field(phi, psi, xy_in, energies);
        compute_light_angle(energies, xy_in, cos_angle);
    }
//     else if (plot == P_3D_LOG_ENERGY)
//     {
//         compute_log_energy_field(phi, psi, xy_in, energies);
//         compute_light_angle(energies, xy_in, cos_angle);
//     }
    else if (plot != P_3D_AMPLITUDE) compute_light_angle(phi, xy_in, cos_angle);
    
    if (cplot == P_3D_PHASE) 
    {
        compute_phase_field(phi, psi, xy_in, phases);
//         compute_energy_field(phi, psi, xy_in, energies);
    }
    
//     if ((plot == P_3D_ANGLE)||(plot == P_3D_AMP_ANGLE)) compute_light_angle(phi, xy_in, cos_angle);
    
//     for (i=0; i<NX-2; i++)
//         for (j=0; j<NY-2; j++)
    for (i=1; i<NX-2; i++)
        for (j=1; j<NY-2; j++)
        {
            if ((TWOSPEEDS)||(xy_in[i*NY+j]))
            {
                switch (plot) {
                    case (P_3D_AMPLITUDE):
                    {
                        if (AMPLITUDE_HIGH_RES)
                        {
                            compute_interpolated_colors(i, j, phi, palette, plot, 
                                 &z_sw, &z_se, &z_nw, &z_ne, &z_mid, rgb_e, rgb_w, rgb_n, rgb_s);
                        }
                        else color_scheme_palette(COLOR_SCHEME, palette, VSCALE_AMPLITUDE*phi[i*NY+j], 1.0, 0, rgb);
                        break;
                    }
                    case (P_3D_ANGLE):
                    {
                        ca = cos_angle[i*NY+j];
                        if (ca < 0.0) ca = 0.0;
                        color_scheme_asym_palette(COLOR_SCHEME, palette, ca, 1.0, 0, rgb);
                        break;
                    }
                    case (P_3D_AMP_ANGLE):
                    {
                        ca = cos_angle[i*NY+j];
                        ca = (ca + 1.0)*0.4 + 0.2;
                        if ((FADE_IN_OBSTACLE)&&(!xy_in[i*NY+j])) ca *= 1.6;
                        if (AMPLITUDE_HIGH_RES)
                        {
                            compute_interpolated_colors_new(i, j, phi, phi, palette, cplot, 
                                 &z_sw, &z_se, &z_nw, &z_ne, &z_mid, rgb_e, rgb_w, rgb_n, rgb_s);
                            for (k=0; k<3; k++) rgb_e[k] *= ca;
                            for (k=0; k<3; k++) rgb_w[k] *= ca;
                            for (k=0; k<3; k++) rgb_n[k] *= ca;
                            for (k=0; k<3; k++) rgb_s[k] *= ca;
                        }
                        else 
                        {
                            color_scheme_palette(COLOR_SCHEME, palette, VSCALE_AMPLITUDE*phi[i*NY+j], 1.0, 0, rgb);
                            for (k=0; k<3; k++) rgb[k] *= ca;
                        }
                        break;
                    }
                    case (P_3D_ENERGY):
                    {
                        ca = cos_angle[i*NY+j];
                        ca = (ca + 1.0)*0.4 + 0.19;
                        if ((FADE_IN_OBSTACLE)&&(!xy_in[i*NY+j])) ca *= 1.6;
                        
                        if (AMPLITUDE_HIGH_RES)
                        {
                            compute_interpolated_colors(i, j, energies, palette, plot, 
                                 &z_sw, &z_se, &z_nw, &z_ne, &z_mid, rgb_e, rgb_w, rgb_n, rgb_s);
                            for (k=0; k<3; k++) 
                            {
                                rgb_e[k] *= ca;
                                rgb_w[k] *= ca;
                                rgb_n[k] *= ca;
                                rgb_s[k] *= ca;
                            }
                        }
                        else 
                        {
                            energy_color_scheme(palette, energies[i*NY+j], rgb);
                            for (k=0; k<3; k++) rgb[k] *= ca;
                        }
                        break;
                    }
                    case (P_3D_LOG_ENERGY):
                    {
                        ca = cos_angle[i*NY+j];
                        ca = (ca + 1.0)*0.4 + 0.19;
                        if ((FADE_IN_OBSTACLE)&&(!xy_in[i*NY+j])) ca *= 1.6;
                        
                        if (AMPLITUDE_HIGH_RES)
                        {
                            compute_interpolated_colors(i, j, energies, palette, plot, 
                                 &z_sw, &z_se, &z_nw, &z_ne, &z_mid, rgb_e, rgb_w, rgb_n, rgb_s);
                            for (k=0; k<3; k++) 
                            {
                                rgb_e[k] *= ca;
                                rgb_w[k] *= ca;
                                rgb_n[k] *= ca;
                                rgb_s[k] *= ca;
                            }
                        }
                        else 
                        {
                            log_energy_color_scheme(palette, energies[i*NY+j], rgb);
                            for (k=0; k<3; k++) rgb[k] *= ca;
                        }
                        break;
                    }
                }
                
                if (fade)
                {
                    for (k=0; k<3; k++)
                    {
                        rgb[k] *= fade_value;
                        rgb_n[k] *= fade_value;
                        rgb_s[k] *= fade_value;
                        rgb_e[k] *= fade_value;
                        rgb_w[k] *= fade_value;
                    }
                    
                }
                
                if ((plot == P_3D_ENERGY)||(plot == P_3D_LOG_ENERGY))
                {                    
                    draw = (xy_in[i*NY+j])&&(xy_in[(i+1)*NY+j])&&(xy_in[i*NY+j+1])&&(xy_in[(i+1)*NY+j+1]);
                }
                
                if ((plot != P_3D_AMPLITUDE)&&(AMPLITUDE_HIGH_RES)&&(draw))
                {
                    ij_to_xy(i, j, xy_sw);
                    ij_to_xy(i+1, j, xy_se);
                    ij_to_xy(i, j+1, xy_nw);
                    ij_to_xy(i+1, j+1, xy_ne);
                    
                    for (k=0; k<2; k++) xy_mid[k] = 0.25*(xy_sw[k] + xy_se[k] + xy_nw[k] + xy_ne[k]);
                    
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
                else
                {
                    glColor3f(rgb[0], rgb[1], rgb[2]);
    
                    glBegin(GL_TRIANGLE_FAN);
                    ij_to_xy(i, j, xy);
                    draw_vertex_xyz(xy, phi[i*NY+j]);
                    ij_to_xy(i+1, j, xy);
                    draw_vertex_xyz(xy, phi[(i+1)*NY+j]);
                    ij_to_xy(i+1, j+1, xy);
                    draw_vertex_xyz(xy, phi[(i+1)*NY+j+1]);
                    ij_to_xy(i, j+1, xy);
                    draw_vertex_xyz(xy, phi[i*NY+j+1]);
                    glEnd ();
                }
            }
            else if (DRAW_OUTSIDE_GRAY)     /* experimental */
            {
                printf("Drawing outside at (%i,%i)\n", i, j);
                glColor3f(0.5, 0.5, 0.5);
                glBegin(GL_TRIANGLE_FAN);
                ij_to_xy(i, j, xy);
                draw_vertex_xyz(xy, 0.0);
                ij_to_xy(i+1, j, xy);
                draw_vertex_xyz(xy, 0.0);
                ij_to_xy(i+1, j+1, xy);
                draw_vertex_xyz(xy, 0.0);
                ij_to_xy(i, j+1, xy);
                draw_vertex_xyz(xy, 0.0);
                glEnd ();
            }
        }

    free(cos_angle);
    free(energies);
    free(phases);
}

void draw_color_scheme_palette_3d(double x1, double y1, double x2, double y2, int plot, double min, double max, int palette, int fade, double fade_value)
{
    int j, k, ij_botleft[2], ij_topright[2], imin, imax, jmin, jmax;
    double y, dy, dy_e, dy_phase, rgb[3], value, lum, amp;
    
    xy_to_ij(x1, y1, ij_botleft);
    xy_to_ij(x2, y2, ij_topright);
    
    rgb[0] = 0.0;   rgb[1] = 0.0;   rgb[2] = 0.0;
    erase_area_rgb(0.5*(x1 + x2), x2 - x1, 0.5*(y1 + y2), y2 - y1, rgb);

    if (ROTATE_COLOR_SCHEME)
    {
        jmin = ij_botleft[0];
        jmax = ij_topright[0];
        imin = ij_botleft[1];
        imax = ij_topright[1];    
    }
    else
    {
        imin = ij_botleft[0];
        imax = ij_topright[0];
        jmin = ij_botleft[1];
        jmax = ij_topright[1];    
    }
        
        
    glBegin(GL_QUADS);
    dy = (max - min)/((double)(jmax - jmin));
    dy_e = max/((double)(jmax - jmin));
    dy_phase = 1.0/((double)(jmax - jmin));
    
    for (j = jmin; j < jmax; j++)
    {
        switch (plot) {
            case (P_3D_AMPLITUDE):
            {
                value = min + 1.0*dy*(double)(j - jmin);
                color_scheme_palette(COLOR_SCHEME, palette, 0.7*value, 1.0, 0, rgb);
                break;
            }
            case (P_3D_ANGLE):
            {
                value = 1.0*dy*(double)(j - jmin);
                color_scheme_asym_palette(COLOR_SCHEME, palette, value, 1.0, 0, rgb);
                break;
            }
            case (P_3D_AMP_ANGLE):
            {
                value = min + 1.0*dy*(double)(j - jmin);
                color_scheme_palette(COLOR_SCHEME, palette, 0.7*value, 1.0, 0, rgb);
                break;
            }
            case (P_3D_ENERGY):
            {
                value = dy_e*(double)(j - jmin)*100.0/E_SCALE;
                if (COLOR_PALETTE >= COL_TURBO) color_scheme_asym_palette(COLOR_SCHEME, palette, value, 1.0, 1, rgb);
                else color_scheme_palette(COLOR_SCHEME, palette, value, 1.0, 1, rgb);
                break;
            }
            case (P_3D_LOG_ENERGY):
            {
//                 value = LOG_SHIFT + LOG_SCALE*dy_e*(double)(j - jmin)*100.0/E_SCALE;
                value = LOG_SCALE*dy_e*(double)(j - jmin)*100.0/E_SCALE;
                color_scheme_asym_palette(COLOR_SCHEME, palette, value, 1.0, 1, rgb);
                break;
            }
            case (P_3D_PHASE):
            {
                value = dy_phase*(double)(j - jmin);
                color_scheme_palette(C_ONEDIM_LINEAR, palette, value, 1.0, 1, rgb);
                break;
            }
            case (Z_MODULE):
            {
                value = min + 1.0*dy*(double)(j - jmin);
                color_scheme_palette(COLOR_SCHEME, palette, 0.7*value, 1.0, 0, rgb);
                break;
            }
            case (Z_ARGUMENT):
            {
                value = dy_phase*(double)(j - jmin);
//                 hsl_to_rgb_palette(value, 0.9, 0.5, rgb, palette);
                color_scheme_palette(C_ONEDIM_LINEAR, palette, value, 1.0, 1, rgb);
                break;
            }
            case (Z_REALPART):
            {
                value = min + 1.0*dy*(double)(j - jmin);
                color_scheme_palette(COLOR_SCHEME, palette, value, 1.0, 1, rgb);
                break;
            }
            case (Z_EULER_VORTICITY):
            {
                value = min + 1.0*dy*(double)(j - jmin);
                color_scheme_palette(COLOR_SCHEME, palette, 0.7*value, 1.0, 0, rgb);
                break;
            }
            case (Z_EULER_LOG_VORTICITY):
            {
                value = min + 1.0*dy*(double)(j - jmin);
                color_scheme_palette(COLOR_SCHEME, palette, 0.7*value, 1.0, 0, rgb);
                break;
            }
            case (Z_EULER_VORTICITY_ASYM):
            {
                value = min + 1.0*dy*(double)(j - jmin);
                color_scheme_palette(COLOR_SCHEME, palette, 0.7*value, 1.0, 0, rgb);
                break;
            }
       }
        if (fade) for (k=0; k<3; k++) rgb[k] *= fade_value;
        glColor3f(rgb[0], rgb[1], rgb[2]);
        if (ROTATE_COLOR_SCHEME)
        {
            draw_vertex_ij(j, imin);
            draw_vertex_ij(j, imax);
            draw_vertex_ij(j+1, imax);
            draw_vertex_ij(j+1, imin);            
        }
        else
        {
            draw_vertex_ij(imin, j);
            draw_vertex_ij(imax, j);
            draw_vertex_ij(imax, j+1);
            draw_vertex_ij(imin, j+1);
        }
    }
    glEnd ();
    
    if (fade) glColor3f(fade_value, fade_value, fade_value);
    else glColor3f(1.0, 1.0, 1.0);
    glLineWidth(BOUNDARY_WIDTH);
    draw_rectangle_noscale(x1, y1, x2, y2);
}

