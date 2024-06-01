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
    
void draw_vertex_x_y_z(double x, double y, double z)
{
    double xy_screen[2];
    
    xyz_to_xy(x, y, z, xy_screen);
    glVertex2d(xy_screen[0], xy_screen[1]);
}

void draw_segment_hsl(double x1, double y1, double x2, double y2, double h, double s, double l)
/* draw line segment (x1,y1)-(x2,y2) in color (h,s,l) */
{
    double rgb[3], pos[2];
    
    glBegin(GL_LINE_STRIP);
    hsl_to_rgb_palette(h, s, l, rgb, COL_JET);
    glColor3f(rgb[0], rgb[1], rgb[2]);
    draw_vertex_x_y_z(x1, y1, 0.0);
    draw_vertex_x_y_z(x2, y2, 0.0);
    glEnd();
}


void draw_segment_rgb(double x1, double y1, double x2, double y2, double r, double g, double b)
/* draw line segment (x1,y1)-(x2,y2) in color (h,s,l) */
{
    double rgb[3], pos[2];
    
    glBegin(GL_LINE_STRIP);
    glColor3f(r, g, b);
    draw_vertex_x_y_z(x1, y1, 0.0);
    draw_vertex_x_y_z(x2, y2, 0.0);
    glEnd();
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

void draw_tpolygon_3d(t_polygon polygon)
{
    int i;
    double pos[2], alpha, dalpha;
    
    dalpha = DPI/(double)polygon.nsides;
    
    glBegin(GL_LINE_LOOP);
    for (i=0; i<=polygon.nsides; i++)
    {
        alpha = PID*polygon.angle + (double)i*dalpha;
        draw_vertex_x_y_z(polygon.xc + polygon.radius*cos(alpha), polygon.yc + polygon.radius*sin(alpha), 0.0);
    }
    glEnd();
}


void draw_billiard_3d(int fade, double fade_value)      /* draws the billiard boundary */
{
    double x0, x, y, x1, y1, dx, dy, phi, r = 0.01, pos[2], pos1[2], alpha, dphi, omega, z, l, width, a, b, c, ymax, padd;
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
                x = LAMBDA*cos(phi);
                y = sin(phi);
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
        case (D_YOUNG):
        {
//             if (FILL_BILLIARD_COMPLEMENT)
//             {
//                 if (fade) glColor3f(0.5*fade_value, 0.5*fade_value, 0.5*fade_value);
//                 else glColor3f(0.5, 0.5, 0.5);
//                 
//                 glBegin(GL_TRIANGLE_FAN);
//                 draw_vertex_x_y_z(-MU, YMIN, 0.0);
//                 draw_vertex_x_y_z(-MU, -LAMBDA-MU, 0.0);
//                 draw_vertex_x_y_z(MU, -LAMBDA-MU, 0.0);
//                 draw_vertex_x_y_z(MU, YMIN, 0.0);
//                 glEnd();
//             
//                 glBegin(GL_TRIANGLE_FAN);
//                 draw_vertex_x_y_z(-MU, YMAX, 0.0);
//                 draw_vertex_x_y_z(-MU, LAMBDA+MU, 0.0);
//                 draw_vertex_x_y_z(MU, LAMBDA+MU, 0.0);
//                 draw_vertex_x_y_z(MU, YMAX, 0.0);
//                 glEnd();
// 
//                 glBegin(GL_TRIANGLE_FAN);
//                 draw_vertex_x_y_z(-MU, -LAMBDA+MU, 0.0);
//                 draw_vertex_x_y_z(-MU, LAMBDA-MU, 0.0);
//                 draw_vertex_x_y_z(MU, LAMBDA-MU, 0.0);
//                 draw_vertex_x_y_z(MU, -LAMBDA+MU, 0.0);
//                 glEnd();
//             }
//             else
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
            break;
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
            if (DRAW_CONSTRUCTION_LINES)
            {
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
            if (DRAW_CONSTRUCTION_LINES)
            {
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
            }
            
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
            if (DRAW_CONSTRUCTION_LINES)
            {
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
            }
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
           if (DRAW_CONSTRUCTION_LINES)
            {
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
            }
            
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
            if (DRAW_CONSTRUCTION_LINES)
            {
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
            }
            break;
        }
        case (D_STAR):
        {
            glLineWidth(BOUNDARY_WIDTH);
            glBegin(GL_LINE_LOOP);
            for (i=0; i<npolyline; i++) tvertex_lineto_3d(polyline[i]);
            glEnd();
            break;
        }
        case (D_NOISEPANEL):
        {
            glLineWidth(BOUNDARY_WIDTH);
            glBegin(GL_LINE_LOOP);
            for (i=0; i<npolyline; i++) tvertex_lineto_3d(polyline[i]);
            glEnd();
            break;
        }
        case (D_POLYGONS):
        {
            glLineWidth(BOUNDARY_WIDTH);
            for (i = 0; i < ncircles; i++) 
                if (polygons[i].active) draw_tpolygon_3d(polygons[i]);
            break;
        }
        case (D_LSHAPE):
        {
            padd = 0.005;
            glLineWidth(BOUNDARY_WIDTH);
            draw_segment_hsl(-LAMBDA - padd, -1.0 - padd, 0.0, -1.0 - padd, 0.0, 1.0, 0.5*fade_value);
            draw_segment_hsl(-LAMBDA - padd,  1.0 + padd, 0.0,  1.0 + padd, 0.0, 1.0, 0.5*fade_value);
            draw_segment_hsl(0.0, -1.0 - padd, LAMBDA + padd, -1.0 - padd, 220.0, 1.0, 0.5*fade_value);
            draw_segment_hsl(0.0, padd, LAMBDA + padd, padd, 220.0, 1.0, 0.5*fade_value);
            draw_segment_hsl(-LAMBDA - padd, -1.0 - padd, -LAMBDA - padd, 0.0, 60.0, 1.0, 0.5*fade_value);
            draw_segment_hsl( LAMBDA + padd, -1.0 - padd,  LAMBDA + padd, 0.0, 60.0, 1.0, 0.5*fade_value);
            draw_segment_hsl(-LAMBDA - padd, 0.0, -LAMBDA - padd, 1.0 + padd, 180.0, 1.0, 0.5*fade_value);
            draw_segment_hsl( padd, padd, padd, LAMBDA + padd, 180.0, 1.0, 0.5*fade_value);
            break;
        }
        case (D_NOTHING):
        {
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
        
void draw_segment_hsl_visible(double x1, double y1, double x2, double y2, double h, double s, double l, double margin)
/* hack to draw the billiard boundary in front of the wave */
/* only parts of the boundary having a small enough angle with the observer vector are drawn */
{
    double length, length1, olength;
    
    olength = module2(observer[0], observer[1]);
    
    length = module2(x1,y1);
    length1 = module2(x2,y2);
    if ((x1*observer[0] + y1*observer[1] > margin*length*olength)&&(x2*observer[0] + y2*observer[1] > margin*length1*olength))
        draw_segment_hsl(x1, y1, x2, y2, h, s, l);
}


void draw_segment_rgb_visible(double x1, double y1, double x2, double y2, double r, double g, double b, double margin)
/* hack to draw the billiard boundary in front of the wave */
/* only parts of the boundary having a small enough angle with the observer vector are drawn */
{
    double length, length1, olength;
    
    olength = module2(observer[0], observer[1]);
    
    length = module2(x1,y1);
    length1 = module2(x2,y2);
    if ((x1*observer[0] + y1*observer[1] > margin*length*olength)&&(x2*observer[0] + y2*observer[1] > margin*length1*olength))
        draw_segment_rgb(x1, y1, x2, y2, r, g, b);
}


void draw_billiard_3d_front(int fade, double fade_value)      
/* hack to draw the billiard boundary in front of the wave */
/* only parts of the boundary having a small enough angle with the observer vector are drawn */
{
    double x0, x, y, x1, y1, dx, dy, phi, r = 0.01, pos[2], pos1[2], alpha, dphi, omega, z, l, width, a, b, c, ymax, length, length1, padd, margin;
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
        case (D_ELLIPSE):
        {
            glBegin(GL_LINE_LOOP);
            dphi = DPI/(double)NSEG;
            margin = 0.0;
            for (i=0; i<=NSEG; i++)
            {
                phi = (double)i*dphi;
                draw_segment_rgb_visible(LAMBDA*cos(phi), sin(phi), LAMBDA*cos(phi+dphi), sin(phi+dphi), fade_value, fade_value, fade_value, margin);
                draw_vertex_x_y_z(LAMBDA*cos(phi), sin(phi), 0.0);
            }
            glEnd ();
            break;
        }
        case (D_YOUNG):
        {
            glBegin(GL_LINE_STRIP);
            draw_vertex_x_y_z(-MU, YMIN, 0.0);
            draw_vertex_x_y_z(-MU, -LAMBDA-MU, 0.0);
//             draw_vertex_x_y_z(MU, -LAMBDA-MU, 0.0);
//             draw_vertex_x_y_z(MU, YMIN, 0.0);
            glEnd();
            
            glBegin(GL_LINE_STRIP);
            draw_vertex_x_y_z(-MU, YMAX, 0.0);
            draw_vertex_x_y_z(-MU, LAMBDA+MU, 0.0);
//             draw_vertex_x_y_z(MU, LAMBDA+MU, 0.0);
//             draw_vertex_x_y_z(MU, YMAX, 0.0);
            glEnd();

            glBegin(GL_LINE_LOOP);
            draw_vertex_x_y_z(-MU, -LAMBDA+MU, 0.0);
            draw_vertex_x_y_z(-MU, LAMBDA-MU, 0.0);
//             draw_vertex_x_y_z(MU, LAMBDA-MU, 0.0);
//             draw_vertex_x_y_z(MU, -LAMBDA+MU, 0.0);
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
        case (D_NOISEPANEL):
        {
            glLineWidth(BOUNDARY_WIDTH);
            for (i=0; i<npolyline; i++)
                draw_polyline_visible(i%npolyline, (i+1)%npolyline, -0.5);
            break;
        }
        case (D_LSHAPE):
        {
            padd = 0.005;
            margin = 0.0;
            glLineWidth(BOUNDARY_WIDTH);
            draw_segment_hsl_visible(-LAMBDA - padd, -1.0 - padd, 0.0, -1.0 - padd, 0.0, 1.0, 0.5*fade_value, margin);
            draw_segment_hsl_visible(-LAMBDA - padd,  1.0 + padd, 0.0,  1.0 + padd, 0.0, 1.0, 0.5*fade_value, margin);
            draw_segment_hsl_visible(0.0, -1.0 - padd, LAMBDA + padd, -1.0 - padd, 220.0, 1.0, 0.5*fade_value, margin);
            draw_segment_hsl_visible(0.0, padd, LAMBDA + padd, padd, 220.0, 1.0, 0.5*fade_value, 0.2);
            draw_segment_hsl_visible(-LAMBDA - padd, -1.0 - padd, -LAMBDA - padd, 0.0, 60.0, 1.0, 0.5*fade_value, margin);
            draw_segment_hsl_visible( LAMBDA + padd, -1.0 - padd,  LAMBDA + padd, 0.0, 60.0, 1.0, 0.5*fade_value, margin);
            draw_segment_hsl_visible(-LAMBDA - padd, 0.0, -LAMBDA - padd, 1.0 + padd, 180.0, 1.0, 0.5*fade_value, margin);
            draw_segment_hsl_visible( padd, padd, padd, LAMBDA + padd, 180.0, 1.0, 0.5*fade_value, 0.6);
            break;
        }
        default:
        {
            break;
        }   
    }
}

void compute_energy_field(double phi[NX*NY], double psi[NX*NY], short int xy_in[NX*NY], t_wave wave[NX*NY])
/* computes cosine of angle between normal vector and vector light */
{
    int i, j, k;
    static int first = 1;
    double energy, logenergy, gx, gy, arg, mod, sum;
    
//     printf("computing energy field\n");
//     printf("COMPUTE_MEAN_ENERGY = %i\n", COMPUTE_MEAN_ENERGY);
    
    #pragma omp parallel for private(i,j)
    for (i=0; i<NX; i++)
        for (j=0; j<NY; j++) 
        {
            if ((TWOSPEEDS)||(xy_in[i*NY+j])) 
            {
                energy = PLOT_SCALE_ENERGY*compute_energy_mod(phi, psi, xy_in, i, j);
                if (energy < 1.0e-100) energy = 0.0;
                wave[i*NY+j].energy = PLOT_SCALE_ENERGY*compute_energy_mod(phi, psi, xy_in, i, j);
                if (COMPUTE_TOTAL_ENERGY)
                {
                    if (first) wave[i*NY+j].total_energy = energy;
                    else wave[i*NY+j].total_energy += energy;
                }
            
                if (COMPUTE_LOG_ENERGY)
                {
                    logenergy = log(wave[i*NY+j].energy);
                    if (logenergy > LOG_ENERGY_FLOOR) wave[i*NY+j].log_energy = LOG_SHIFT + PLOT_SCALE_LOG_ENERGY*logenergy;
                    else wave[i*NY+j].log_energy = LOG_SHIFT + PLOT_SCALE_LOG_ENERGY*LOG_ENERGY_FLOOR;
                }
                
                if (COMPUTE_LOG_TOTAL_ENERGY)
                {
                    logenergy = log(wave[i*NY+j].total_energy);
                    if (logenergy > LOG_ENERGY_FLOOR) wave[i*NY+j].log_total_energy = LOG_SHIFT + PLOT_SCALE_LOG_ENERGY*logenergy;
                    else wave[i*NY+j].log_total_energy = LOG_SHIFT + PLOT_SCALE_LOG_ENERGY*LOG_ENERGY_FLOOR;
                }
                
                if (COMPUTE_MEAN_ENERGY)
                {
                    wave[i*NY+j].mean_energy = wave[i*NY+j].total_energy/((double)(global_time+1));
                }
            
                if (COMPUTE_LOG_MEAN_ENERGY)
                {
                    logenergy = wave[i*NY+j].mean_energy;
                    if (logenergy == 0.0) logenergy = 1.0e-10;
                    logenergy = log(logenergy) + LOG_MEAN_ENERGY_SHIFT;
                    if (logenergy > LOG_ENERGY_FLOOR) wave[i*NY+j].log_mean_energy = LOG_SHIFT + PLOT_SCALE_LOG_ENERGY*logenergy;
                    else wave[i*NY+j].log_mean_energy = LOG_SHIFT + PLOT_SCALE_LOG_ENERGY*LOG_ENERGY_FLOOR;
                }
                
                if (COMPUTE_ENERGY_FLUX)
                {
                    compute_energy_flux_mod(phi, psi, xy_in, i, j, &gx, &gy, &arg, &mod);
                    wave[i*NY+j].flux_direction = arg/DPI;
                    
                    /* compute time-averaged flux intensity */
                    wave[i*NY+j].flux_int_table[wave[i*NY+j].flux_counter] = mod*FLUX_SCALE;
                    sum = 0.0;
                    for (k = 0; k < FLUX_WINDOW; k++) sum += wave[i*NY+j].flux_int_table[k];
                    wave[i*NY+j].flux_intensity = sum/(double)FLUX_WINDOW;
                    wave[i*NY+j].flux_counter++;
                    if (wave[i*NY+j].flux_counter == FLUX_WINDOW) wave[i*NY+j].flux_counter = 0;
                }
            }
            else if (first)
            {
                wave[i*NY+j].energy = 0.0;
                wave[i*NY+j].total_energy = 0.0;
                wave[i*NY+j].log_total_energy = LOG_ENERGY_FLOOR;
                wave[i*NY+j].mean_energy = 0.0;
                wave[i*NY+j].log_mean_energy = LOG_ENERGY_FLOOR;
                wave[i*NY+j].flux_intensity = 0.0;
                wave[i*NY+j].flux_direction = 0.0;
            }
        }
        
    first = 0;
}


void compute_log_energy_field(double phi[NX*NY], double psi[NX*NY], short int xy_in[NX*NY], t_wave wave[NX*NY])
/* computes cosine of angle between normal vector and vector light */
/* TO DO: include into compute_energy_field */
{
    int i, j;
    double value;
    
    #pragma omp parallel for private(i,j,value)
    for (i=0; i<NX; i++)
        for (j=0; j<NY; j++)
        {
            if ((TWOSPEEDS)||(xy_in[i*NY+j])) 
            {
                value = log(compute_energy_mod(phi, psi, xy_in, i, j));
                if (value > LOG_ENERGY_FLOOR) wave[i*NY+j].log_energy = LOG_SHIFT + PLOT_SCALE_LOG_ENERGY*value;
                else wave[i*NY+j].log_energy = LOG_SHIFT + PLOT_SCALE_LOG_ENERGY*LOG_ENERGY_FLOOR;
            }
            else wave[i*NY+j].log_energy = 0.0;
        }
}



void compute_phase_field(double phi[NX*NY], double psi[NX*NY], short int xy_in[NX*NY], t_wave wave[NX*NY])
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
                wave[i*NY+j].phase = angle;
            }
            else wave[i*NY+j].phase = 0.0;
        }
}

void compute_light_angle(short int xy_in[NX*NY], t_wave wave[NX*NY], int movie)
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
                gradx = (*wave[(i+1)*NY+j].p_zfield[movie] - *wave[(i-1)*NY+j].p_zfield[movie])/dx;
                grady = (*wave[i*NY+j+1].p_zfield[movie] - *wave[i*NY+j-1].p_zfield[movie])/dy;
                
                if (ADD_POTENTIAL)
                {
                    gradx -= (*wave[(i+1)*NY+j].potential - *wave[(i-1)*NY+j].potential)*POT_FACT/dx;
                    grady -= (*wave[i*NY+j+1].potential - *wave[i*NY+j-1].potential)*POT_FACT/dy;
                }
                
                norm = sqrt(1.0 + gradx*gradx + grady*grady);
                pscal = -gradx*light[0] - grady*light[1] + 1.0;
                
                wave[i*NY+j].cos_angle = pscal/norm;
            }
        }
}


void compute_field_color(double value, double value2, int cplot, int palette, double rgb[3])
/* compute the color depending on the field value and color palette */
/* value2 is only used for flux representation */
{
    int k;
    
    switch (cplot) {
        case (P_3D_AMP_ANGLE): 
        {
            color_scheme_palette(COLOR_SCHEME, palette, VSCALE_AMPLITUDE*value, 1.0, 0, rgb);
            break;
        }
        case (P_3D_ENERGY):
        {
            if (COLOR_PALETTE >= COL_TURBO) color_scheme_asym_palette(COLOR_SCHEME, palette, VSCALE_ENERGY*value, 1.0, 0, rgb);
            else color_scheme_palette(COLOR_SCHEME, palette, value, 1.0, 0, rgb);
            break;
        }
        case (P_3D_LOG_ENERGY):
        {
            color_scheme_palette(COLOR_SCHEME, palette, LOG_SHIFT + LOG_SCALE*value, 1.0, 0, rgb);
            break;
        }
        case (P_3D_TOTAL_ENERGY):
        {
            if (COLOR_PALETTE >= COL_TURBO) color_scheme_asym_palette(COLOR_SCHEME, palette, VSCALE_ENERGY*value, 1.0, 0, rgb);
            else color_scheme_palette(COLOR_SCHEME, palette, value, 1.0, 0, rgb);
            break;
        }
        case (P_3D_LOG_TOTAL_ENERGY):
        {
            color_scheme_palette(COLOR_SCHEME, palette, LOG_SHIFT + LOG_SCALE*value, 1.0, 0, rgb);
            break;
        }
        case (P_3D_MEAN_ENERGY):
        {
            if (COLOR_PALETTE >= COL_TURBO) color_scheme_asym_palette(COLOR_SCHEME, palette, VSCALE_ENERGY*value, 1.0, 0, rgb);
            else color_scheme_palette(COLOR_SCHEME, palette, value, 1.0, 0, rgb);
            break;
        }
        case (P_3D_LOG_MEAN_ENERGY):
        {
            color_scheme_palette(COLOR_SCHEME, palette, LOG_SHIFT + LOG_SCALE*value, 1.0, 0, rgb);
            break;
        }
        case (P_3D_PHASE):
        {
            amp_to_rgb_palette(value, rgb, palette);
            break;
        }
        case (P_3D_FLUX_INTENSITY):
        {
            color_scheme_asym_palette(COLOR_SCHEME, palette, value*FLUX_SCALE, 1.0, 0, rgb);
            break;
        }
        case (P_3D_FLUX_DIRECTION):
        {
            amp_to_rgb_palette(value, rgb, palette);
            for (k=0; k<3; k++) rgb[k] *= tanh(value2*FLUX_CSCALE); 
            break;
        }
    }
}


double compute_interpolated_colors_wave(int i, int j, short int xy_in[NX*NY], t_wave wave[NX*NY], 
                                        double palette, int cplot, double rgb_e[3], double rgb_w[3], double rgb_n[3], double rgb_s[3], int fade, double fade_value, int movie)
{
    int k;
    double cw, ce, cn, cs, c_sw, c_se, c_nw, c_ne, c_mid, ca, z_mid;
    double cw2, ce2, cn2, cs2, factor;
    double *z_sw, *z_se, *z_nw, *z_ne;
    
    z_sw = wave[i*NY+j].p_zfield[movie];
    z_se = wave[(i+1)*NY+j].p_zfield[movie];
    z_nw = wave[i*NY+j+1].p_zfield[movie];
    z_ne = wave[(i+1)*NY+j+1].p_zfield[movie];

    z_mid = 0.25*(*z_sw + *z_se + *z_nw + *z_ne);
    
    c_sw = *wave[i*NY+j].p_cfield[movie];
    c_se = *wave[(i+1)*NY+j].p_cfield[movie];
    c_nw = *wave[i*NY+j+1].p_cfield[movie];
    c_ne = *wave[(i+1)*NY+j+1].p_cfield[movie];
    
    c_mid = 0.25*(c_sw + c_se + c_nw + c_ne);
                                
    cw = (c_sw + c_nw + c_mid)/3.0;
    ce = (c_se + c_ne + c_mid)/3.0;
    cs = (c_sw + c_se + c_mid)/3.0;
    cn = (c_nw + c_ne + c_mid)/3.0;
    
    /* data for second color parameter */
    if (CHANGE_LUMINOSITY)
    {
        c_sw = *wave[i*NY+j].p_cfield[movie+2];
        c_se = *wave[(i+1)*NY+j].p_cfield[movie+2];
        c_nw = *wave[i*NY+j+1].p_cfield[movie+2];
        c_ne = *wave[(i+1)*NY+j+1].p_cfield[movie+2];
    
        c_mid = 0.25*(c_sw + c_se + c_nw + c_ne);
                                
        cw2 = (c_sw + c_nw + c_mid)/3.0;
        ce2 = (c_se + c_ne + c_mid)/3.0;
        cs2 = (c_sw + c_se + c_mid)/3.0;
        cn2 = (c_nw + c_ne + c_mid)/3.0;
    }
    
    compute_field_color(ce, ce2, cplot, palette, rgb_e);
    compute_field_color(cw, cw2, cplot, palette, rgb_w);
    compute_field_color(cn, cn2, cplot, palette, rgb_n);
    compute_field_color(cs, cs2, cplot, palette, rgb_s);
    
    if (SHADE_3D)
    {
        ca = wave[i*NY+j].cos_angle;
        ca = (ca + 1.0)*0.4 + 0.2;
//         if ((FADE_IN_OBSTACLE)&&(!xy_in[i*NY+j])) ca *= 1.6;
        if ((FADE_IN_OBSTACLE)&&(!xy_in[i*NY+j])) ca = (ca + 0.1)*1.6;
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
    
//     if (ADD_POTENTIAL) 
//     {
//         factor = 0.25*POT_FACT;
//         z_mid += *wave[i*NY+j].potential*factor;
//         z_mid += *wave[(i+1)*NY+j].potential*factor;
//         z_mid += *wave[i*NY+j+1].potential*factor;
//         z_mid += *wave[(i+1)*NY+j+1].potential*factor;
//     }

    
    return(z_mid);
}


void compute_wave_fields(double phi[NX*NY], double psi[NX*NY], short int xy_in[NX*NY], int zplot, int cplot, t_wave wave[NX*NY])
/* compute the necessary auxiliary fields */
{
    int i, j;
    
    if (COMPUTE_ENERGY)
        compute_energy_field(phi, psi, xy_in, wave);
    
    if ((zplot == P_3D_LOG_ENERGY)||(cplot == P_3D_LOG_ENERGY))
        compute_log_energy_field(phi, psi, xy_in, wave);
    
    if ((zplot == P_3D_PHASE)||(cplot == P_3D_PHASE))
        compute_phase_field(phi, psi, xy_in, wave);

}


void init_speed_dissipation(short int xy_in[NX*NY], double tc[NX*NY], double tcc[NX*NY], double tgamma[NX*NY])
/* initialise fields for wave speed and dissipation */
{
    int i, j, k, n, inlens;
    double courant2 = COURANT*COURANT, courantb2 = COURANTB*COURANTB, lambda1, mu1;
    double u, v, u1, x, y, xy[2], norm2, speed, r2, c, salpha, h, ll, ca, sa, x1, y1, dx, dy, sum, sigma, x0, y0, rgb[3];
    double xc[NGRIDX*NGRIDY], yc[NGRIDX*NGRIDY], height[NGRIDX*NGRIDY];
    t_circle circles[NMAXCIRCLES];
    
    if (VARIABLE_IOR)
    {
        switch (IOR) {
            case (IOR_MANDELBROT):
            {
                #pragma omp parallel for private(i,j)
                for (i=0; i<NX; i++){
                    for (j=0; j<NY; j++){
                        ij_to_xy(i, j, xy);
                        x = xy[0];
                        y = xy[1];
                        u = 0.0;
                        v = 0.0;
                        k = 0;
                        while ((k<MANDELLEVEL)&&(u*u+v*v < 1000.0*MANDELLIMIT))
                        {
                            u1 = u*u - v*v + x;
                            v = 2.0*u*v + y;
                            u = u1;
                            k++;
                        }
                        norm2 = u*u + v*v;
                        if (norm2 < MANDELLIMIT)
                        {
                            tc[i*NY+j] = COURANT;
                            tcc[i*NY+j] = courant2;
                            tgamma[i*NY+j] = GAMMA;
                        }
                        else 
                        {
                            speed = 1.0 + MANDEL_IOR_SCALE*log(1.0 + norm2/MANDELLIMIT);
                            if (speed < 0.01) speed = 0.01;
                            tcc[i*NY+j] = courant2*speed;
                            tc[i*NY+j] = COURANT*sqrt(speed);
                            tgamma[i*NY+j] = GAMMA;
                        }
                    }
                }
                break;
            }
            case (IOR_MANDELBROT_LIN):
            {
                #pragma omp parallel for private(i,j)
                for (i=0; i<NX; i++){
                    for (j=0; j<NY; j++){
                        ij_to_xy(i, j, xy);
                        x = xy[0];
                        y = xy[1];
                        u = 0.0;
                        v = 0.0;
                        k = 0;
                        while ((k<MANDELLEVEL)&&(u*u+v*v < 1000.0*MANDELLIMIT))
                        {
                            u1 = u*u - v*v + x;
                            v = 2.0*u*v + y;
                            u = u1;
                            k++;
                        }
                        if (k >= MANDELLEVEL)
                        {
                            tc[i*NY+j] = COURANT;
                            tcc[i*NY+j] = courant2;
                            tgamma[i*NY+j] = GAMMA;
                        }
                        else 
                        {
                            speed = (double)k/(double)MANDELLEVEL;
                            if (speed < 1.0e-10) speed = 1.0e-10;
                            else if (speed > 10.0) speed = 10.0;
                            tcc[i*NY+j] = courant2*speed;
                            tc[i*NY+j] = COURANT*sqrt(speed);
                            tgamma[i*NY+j] = GAMMA;
                        }
                    }
                }
                break;
            }
            case (IOR_MANDELBROT_MOD):
            {
                #pragma omp parallel for private(i,j)
                for (i=0; i<NX; i++){
                    for (j=0; j<NY; j++){
                        ij_to_xy(i, j, xy);
                        x = xy[0];
                        y = xy[1];
                        u = 0.0;
                        v = 0.0;
                        k = 0;
                        while ((k<MANDELLEVEL)&&(u*u+v*v < 1000.0*MANDELLIMIT))
                        {
                            u1 = u*u - v*v + x;
                            v = 2.0*u*v + y;
                            u = u1;
                            k++;
                        }
                        if (k >= MANDELLEVEL)
                        {
                            tc[i*NY+j] = COURANT;
                            tcc[i*NY+j] = courant2;
                            tgamma[i*NY+j] = GAMMA;
                        }
                        else 
                        {
//                             speed = 5.0 - 4.0*pow((double)k/(double)MANDELLEVEL, 0.1);
                            speed = 1.0 + 4.0*log(1.0 - 0.1*log((double)k/(double)MANDELLEVEL));
                            if (speed < 1.0e-10) speed = 1.0e-10;
                            else if (speed > 10.0) speed = 10.0;
                            tcc[i*NY+j] = courantb2*speed;
                            tc[i*NY+j] = COURANTB*sqrt(speed);
                            tgamma[i*NY+j] = GAMMAB;
                        }
                    }
                }
                break;
            }
            case (IOR_EARTH):
            {
                for (i=0; i<NX; i++){
                    for (j=0; j<NY; j++){
                        ij_to_xy(i, j, xy);
                        r2 = xy[0]*xy[0] + xy[1]*xy[1];
                        if (r2 > 1.0) c = 0.0;
                        else if (r2 < 0.25*0.25) c = 0.8*COURANT;
                        else if (r2 < 0.58*0.58) c = COURANT*(0.68 - 0.55*r2);
                        else c = COURANT*(1.3 - 0.9*r2);
                        tc[i*NY+j] = c;
                        tcc[i*NY+j] = c;
                        tgamma[i*NY+j] = GAMMA;
                    }
                }
                break;
            }
            case (IOR_EXPLO_LENSING):
            {
                salpha = DPI/(double)NPOLY;
//                 lambda1 = LAMBDA;
//                 mu1 = LAMBDA;
                lambda1 = 0.5*LAMBDA;
                mu1 = 0.5*LAMBDA;
                h = lambda1*tan(PI/(double)NPOLY);
                if (h < mu1) ll = sqrt(mu1*mu1 - h*h);
                else ll = 0.0;
                
//                 #pragma omp parallel for private(i,j)
                for (i=0; i<NX; i++){
                    for (j=0; j<NY; j++) if (xy_in[i*NY+j]) {
                        ij_to_xy(i, j, xy);
                        x = xy[0];
                        y = xy[1];
                        inlens = 0;
                        for (k=0; k<NPOLY; k++)
                        {
                            ca = cos(((double)k+0.5)*salpha + APOLY*PID);
                            sa = sin(((double)k+0.5)*salpha + APOLY*PID);
                            x1 = x*ca + y*sa;
                            y1 = -x*sa + y*ca;
                            if ((module2(x1 - lambda1 - ll, y1) < mu1)&&(module2(x1 - lambda1 + ll, y1) < mu1)) inlens = 1; 
                        }
                        if (inlens) c = COURANTB;
                        else c = COURANT;
                        tc[i*NY+j] = c;
                        tcc[i*NY+j] = c*c;
                        tgamma[i*NY+j] = GAMMA;
                    }
                    else
                    {
                        tc[i*NY+j] = 0.0;
                        tcc[i*NY+j] = 0.0;
                        tgamma[i*NY+j] = 0.0;
                    }
                }
                break;
            }
            case (IOR_PERIODIC_WELLS):
            {
                dx = (XMAX - XMIN)/(double)NGRIDX;
                dy = (YMAX - YMIN)/(double)NGRIDY;
                sigma = 0.2*dx*dx;
                for (i=0; i<NGRIDX; i++)
                    for (j=0; j<NGRIDY; j++)
                    {
                        
                        n = j*NGRIDX + i;
                        xc[n] = XMIN + dx*((double)i + 0.5);
                        yc[n] = YMIN + dy*((double)j + 0.5);
                        if (j%2 == 1) yc[n] += 0.5*dx;
                    }
                    
                for (i=0; i<NX; i++){
                    for (j=0; j<NY; j++){
                        ij_to_xy(i, j, xy);
                        x = xy[0];
                        y = xy[1];
                        sum = 0.0;
                        for (n = 0; n<NGRIDX*NGRIDY; n++)
                        {
                            r2 = (x - xc[n])*(x - xc[n]) + (y - yc[n])*(y - yc[n]);
                            sum += exp(-r2/(sigma));
                        }
                        tc[i*NY+j] = COURANT*sum + COURANTB*(1.0-sum);
                        tcc[i*NY+j] = COURANT*sum + COURANTB*(1.0-sum);
                        tgamma[i*NY+j] = GAMMA;
                    }
                }
                break;
            }
            case (IOR_RANDOM_WELLS):
            {
                dx = (XMAX - XMIN)/(double)NGRIDX;
                dy = (YMAX - YMIN)/(double)NGRIDY;
                sigma = 0.2*dx*dx;
                for (i=0; i<NGRIDX; i++)
                    for (j=0; j<NGRIDY; j++)
                    {
                        
                        n = j*NGRIDX + i;
                        xc[n] = XMIN + dx*((double)i + 0.5 + 0.1*gaussian());
                        yc[n] = YMIN + dy*((double)j + 0.5 + 0.1*gaussian());
//                         if (j%2 == 1) yc[n] += 0.5*dx;
                        height[n] = 0.5 + 0.5*gaussian();
                        if (height[n] > 1.0) height[n] = 1.0;
                        if (height[n] < 0.0) height[n] = 0.0;
                    }
                    
//                 #pragma omp parallel for private(i,j)
                for (i=0; i<NX; i++){
                    for (j=0; j<NY; j++){
                        ij_to_xy(i, j, xy);
                        x = xy[0];
                        y = xy[1];
                        sum = 0.0;
                        for (n = 0; n<NGRIDX*NGRIDY; n++)
                        {
                            r2 = (x - xc[n])*(x - xc[n]) + (y - yc[n])*(y - yc[n]);
                            sum += exp(-r2/(sigma))*height[n];
                        }
                        sum = tanh(sum);
                        tc[i*NY+j] = COURANT*sum + COURANTB*(1.0-sum);
                        tcc[i*NY+j] = COURANT*sum + COURANTB*(1.0-sum);
                        tgamma[i*NY+j] = GAMMA;
                    }
                }
                break;
            }
            case (IOR_POISSON_WELLS):
            {
                ncircles = init_circle_config_pattern(circles, C_POISSON_DISC);
                for (n = 0; n<ncircles; n++)
                {
                    height[n] = 0.5 + 0.5*gaussian();
                    if (height[n] > 1.0) height[n] = 1.0;
                    if (height[n] < 0.0) height[n] = 0.0;
                }
                
                for (n = 0; n<ncircles; n++) printf("Circle %i at (%.3lg, %.3lg) height %.3lg\n", n, circles[n].xc, circles[n].yc, height[n]);
                
                sigma = 0.2*(XMAX - XMIN)*(YMAX - YMIN)/(double)ncircles;
//                 sigma = MU*MU;
                
                for (i=0; i<NX; i++)
                {
                    if (i%100 == 0) printf("Computing potential for column %i of %i\n", i, NX);
                    for (j=0; j<NY; j++){
                        ij_to_xy(i, j, xy);
                        x = xy[0];
                        y = xy[1];
                        sum = 0.0;
                        for (n = 0; n<ncircles; n++)
                        {
                            r2 = (x - circles[n].xc)*(x - circles[n].xc) + (y - circles[n].yc)*(y - circles[n].yc);
                            sum += exp(-r2/(sigma))*height[n];
                        }
                        sum = tanh(sum);
//                         printf("%.3lg\n", sum);
                        tc[i*NY+j] = COURANT*sum + COURANTB*(1.0-sum);
                        tcc[i*NY+j] = COURANT*sum + COURANTB*(1.0-sum);
                        tgamma[i*NY+j] = GAMMA;
                    }
                }
                
                break;
            }
            case (IOR_PPP_WELLS):
            {
                ncircles = init_circle_config_pattern(circles, C_RAND_POISSON);
                for (n = 0; n<ncircles; n++)
                {
                    height[n] = 0.5 + 0.5*gaussian();
                    if (height[n] > 1.0) height[n] = 1.0;
                    if (height[n] < 0.0) height[n] = 0.0;
                }
                
                for (n = 0; n<ncircles; n++) printf("Circle %i at (%.3lg, %.3lg) height %.3lg\n", n, circles[n].xc, circles[n].yc, height[n]);
                
                sigma = 0.2*(XMAX - XMIN)*(YMAX - YMIN)/(double)ncircles;
//                 sigma = MU*MU;
                
                for (i=0; i<NX; i++)
                {
                    if (i%100 == 0) printf("Computing potential for column %i of %i\n", i, NX);
                    for (j=0; j<NY; j++){
                        ij_to_xy(i, j, xy);
                        x = xy[0];
                        y = xy[1];
                        sum = 0.0;
                        for (n = 0; n<ncircles; n++)
                        {
                            r2 = (x - circles[n].xc)*(x - circles[n].xc) + (y - circles[n].yc)*(y - circles[n].yc);
                            sum += exp(-r2/(sigma))*height[n];
                        }
                        sum = tanh(sum);
//                         printf("%.3lg\n", sum);
                        tc[i*NY+j] = COURANT*sum + COURANTB*(1.0-sum);
                        tcc[i*NY+j] = COURANT*sum + COURANTB*(1.0-sum);
                        tgamma[i*NY+j] = GAMMA;
                    }
                }
                
                break;
            }
            case (IOR_LENS_WALL):
            {
                printf("Case IOR_LENS_WALL in init_speed_dissipation of sub_wave_3d needs to be updated\n");
                exit(1);
                break;
            }
            case (IOR_LENS_CONCAVE):
            {
                printf("Case IOR_LENS_CONCAVE in init_speed_dissipation of sub_wave_3d needs to be updated\n");
                exit(1);
                break;
            }
            case (IOR_LENS_CONVEX_CONCAVE):
            {
                printf("Case IOR_LENS_CONVEX_CONCAVE in init_speed_dissipation of sub_wave_3d needs to be updated\n");
                exit(1);
                break;
            }
            case (IOR_TREE):
            {
                printf("Case IOR_TREE in init_speed_dissipation of sub_wave_3d needs to be updated\n");
                exit(1);
                break;
            }
            case (IOR_WAVE_GUIDES_COUPLED):
            {
                printf("Case IOR_WAVE_GUIDES_COUPLED in init_speed_dissipation of sub_wave_3d needs to be updated\n");
                exit(1);
                break;
            }
            case (IOR_WAVE_GUIDES_COUPLED_B):
            {
                printf("Case IOR_WAVE_GUIDES_COUPLED_B in init_speed_dissipation of sub_wave_3d needs to be updated\n");
                exit(1);
                break;
            }
            case (IOR_WAVE_GUIDE_COATED):
            {
                printf("Case IOR_WAVE_GUIDE_COATED in init_speed_dissipation of sub_wave_3d needs to be updated\n");
                exit(1);
                break;
            }
            case (IOR_MICHELSON):
            {
                printf("Case IOR_MICHELSON in init_speed_dissipation of sub_wave_3d needs to be updated\n");
                exit(1);
                break;
            }
            case (IOR_GRADIENT_INDEX_LENS):
            {
                /* focal distance if f = 1/(4*LAMBDA*n0*MU) */
                /* with n0 = COURANT/COURANTB */
                for (i=0; i<NX; i++){
                    for (j=0; j<NY; j++){
                        ij_to_xy(i, j, xy);
                        x = vabs(xy[0]);
                        y = vabs(xy[1]);
                        tc[i*NY+j] = COURANT;
                        if ((x > LAMBDA)) 
                        {
                            tcc[i*NY+j] = courant2;
                            tgamma[i*NY+j] = GAMMA;
                        }
                        else if (y > 1.0)
                        {
                            tcc[i*NY+j] = 0.0;
                            tgamma[i*NY+j] = 0.0;                            
                        }
                        else
                        {
                            tcc[i*NY+j] = courantb2/(1.0 - MU*y*y);
                            tgamma[i*NY+j] = GAMMAB;
                        }
                    }
                }
                break;
            }
            case (IOR_GRADIENT_INDEX_LENS_B):
            {
                /* focal distance if f = 1/(4*LAMBDA*n0*MU) */
                /* with n0 = COURANT/COURANTB */
                for (i=0; i<NX; i++){
                    for (j=0; j<NY; j++){
                        ij_to_xy(i, j, xy);
                        x = vabs(xy[0]);
                        y = vabs(xy[1]);
                        tc[i*NY+j] = COURANT;
                        if ((x > LAMBDA)) 
                        {
                            tcc[i*NY+j] = courant2;
                            tgamma[i*NY+j] = GAMMA;
                        }
                        else if (y > 1.0)
                        {
                            tcc[i*NY+j] = 0.0;
                            tgamma[i*NY+j] = 0.0;                            
                        }
                        else
                        {
                            speed = 1.0/(1.0 - MU*y*y);
                            speed *= speed;
                            tcc[i*NY+j] = courantb2*speed;
                            tgamma[i*NY+j] = GAMMAB;
                        }
                    }
                }
                break;
            }
            default:
            {
                for (i=0; i<NX; i++){
                    for (j=0; j<NY; j++){
                        tc[i*NY+j] = COURANT;
                        tcc[i*NY+j] = COURANT;
                        tgamma[i*NY+j] = GAMMA;
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
                    tc[i*NY+j] = COURANT;
                    tcc[i*NY+j] = courant2;
                    if (xy_in[i*NY+j] == 1) tgamma[i*NY+j] = GAMMA;
                    else tgamma[i*NY+j] = GAMMAB;
                }
                else if (TWOSPEEDS)
                {
                    tc[i*NY+j] = COURANTB;
                    tcc[i*NY+j] = courantb2;
                    tgamma[i*NY+j] = GAMMAB;
                }
            }
        }
    }
}


void init_zfield(double phi[NX*NY], double psi[NX*NY], short int xy_in[NX*NY], int zplot, t_wave wave[NX*NY], int movie)
/* compute the necessary fields for the z coordinate */
{
    int i, j;
    
    switch(zplot) {
        case (P_3D_AMP_ANGLE):
        {
            #pragma omp parallel for private(i,j)
            for (i=0; i<NX; i++) for (j=0; j<NY; j++) wave[i*NY+j].p_zfield[movie] = &phi[i*NY+j];
            break;
        }
        case (P_3D_ENERGY):
        {
            #pragma omp parallel for private(i,j)
            for (i=0; i<NX; i++) for (j=0; j<NY; j++) wave[i*NY+j].p_zfield[movie] = &wave[i*NY+j].energy;
            break;
        }
        case (P_3D_LOG_ENERGY):
        {
            #pragma omp parallel for private(i,j)
            for (i=0; i<NX; i++) for (j=0; j<NY; j++) wave[i*NY+j].p_zfield[movie] = &wave[i*NY+j].log_energy;
            break;
        }
        case (P_3D_TOTAL_ENERGY):
        {
            #pragma omp parallel for private(i,j)
            for (i=0; i<NX; i++) for (j=0; j<NY; j++) wave[i*NY+j].p_zfield[movie] = &wave[i*NY+j].total_energy;
            break;
        }
        case (P_3D_LOG_TOTAL_ENERGY):
        {
            #pragma omp parallel for private(i,j)
            for (i=0; i<NX; i++) for (j=0; j<NY; j++) wave[i*NY+j].p_zfield[movie] = &wave[i*NY+j].log_total_energy;
            break;
        }
        case (P_3D_MEAN_ENERGY):
        {
            #pragma omp parallel for private(i,j)
            for (i=0; i<NX; i++) for (j=0; j<NY; j++) wave[i*NY+j].p_zfield[movie] = &wave[i*NY+j].mean_energy;
            break;
        }
        case (P_3D_LOG_MEAN_ENERGY):
        {
            #pragma omp parallel for private(i,j)
            for (i=0; i<NX; i++) for (j=0; j<NY; j++) wave[i*NY+j].p_zfield[movie] = &wave[i*NY+j].log_mean_energy;
            break;
        }
        case (P_3D_PHASE):
        {
            #pragma omp parallel for private(i,j)
            for (i=0; i<NX; i++) for (j=0; j<NY; j++) wave[i*NY+j].p_zfield[movie] = &wave[i*NY+j].phase;
            break;
        }
        case (P_3D_FLUX_INTENSITY):
        {
            #pragma omp parallel for private(i,j)
            for (i=0; i<NX; i++) for (j=0; j<NY; j++) wave[i*NY+j].p_zfield[movie] = &wave[i*NY+j].flux_intensity;
            break;
        }
        case (P_3D_FLUX_DIRECTION):
        {
            #pragma omp parallel for private(i,j)
            for (i=0; i<NX; i++) for (j=0; j<NY; j++) wave[i*NY+j].p_zfield[movie] = &wave[i*NY+j].flux_direction;
            break;
        }        
    }
}

void init_cfield(double phi[NX*NY], double psi[NX*NY], short int xy_in[NX*NY], int cplot, t_wave wave[NX*NY], int movie)
/* compute the colors */
{
    int i, j, k;
    
    switch(cplot) {
        case (P_3D_AMP_ANGLE):
        {
            #pragma omp parallel for private(i,j)
            for (i=0; i<NX; i++) for (j=0; j<NY; j++) wave[i*NY+j].p_cfield[movie] = &phi[i*NY+j];
            break;
        }
        case (P_3D_ENERGY):
        {
            #pragma omp parallel for private(i,j)
            for (i=0; i<NX; i++) for (j=0; j<NY; j++) wave[i*NY+j].p_cfield[movie] = &wave[i*NY+j].energy;
            break;
        }
        case (P_3D_LOG_ENERGY):
        {
            #pragma omp parallel for private(i,j)
            for (i=0; i<NX; i++) for (j=0; j<NY; j++) wave[i*NY+j].p_cfield[movie] = &wave[i*NY+j].log_energy;
            break;
        }
        case (P_3D_TOTAL_ENERGY):
        {
            #pragma omp parallel for private(i,j)
            for (i=0; i<NX; i++) for (j=0; j<NY; j++) wave[i*NY+j].p_cfield[movie] = &wave[i*NY+j].total_energy;
            break;
        }
        case (P_3D_LOG_TOTAL_ENERGY):
        {
            #pragma omp parallel for private(i,j)
            for (i=0; i<NX; i++) for (j=0; j<NY; j++) wave[i*NY+j].p_cfield[movie] = &wave[i*NY+j].log_total_energy;
            break;
        }
        case (P_3D_MEAN_ENERGY):
        {
            #pragma omp parallel for private(i,j)
            for (i=0; i<NX; i++) for (j=0; j<NY; j++) wave[i*NY+j].p_cfield[movie] = &wave[i*NY+j].mean_energy;
            break;
        }
        case (P_3D_LOG_MEAN_ENERGY):
        {
            #pragma omp parallel for private(i,j)
            for (i=0; i<NX; i++) for (j=0; j<NY; j++) wave[i*NY+j].p_cfield[movie] = &wave[i*NY+j].log_mean_energy;
            break;
        }
        case (P_3D_PHASE):
        {
            #pragma omp parallel for private(i,j)
            for (i=0; i<NX; i++) for (j=0; j<NY; j++) wave[i*NY+j].p_cfield[movie] = &wave[i*NY+j].phase;
            break;
        }
        case (P_3D_FLUX_INTENSITY):
        {
            #pragma omp parallel for private(i,j)
            for (i=0; i<NX; i++) for (j=0; j<NY; j++) wave[i*NY+j].p_cfield[movie] = &wave[i*NY+j].flux_intensity;
            break;
        }
        case (P_3D_FLUX_DIRECTION):
        {
            #pragma omp parallel for private(i,j)
            for (i=0; i<NX; i++) for (j=0; j<NY; j++) 
            {
                wave[i*NY+j].p_cfield[movie] = &wave[i*NY+j].flux_direction;
                wave[i*NY+j].p_cfield[movie+2] = &wave[i*NY+j].flux_intensity;
                /* information on intensity stored in second pointer to adjust luminosity */
            }
            break;
        }        
    }
    
    if (cplot != P_3D_FLUX_DIRECTION)
    {
        #pragma omp parallel for private(i,j)
            for (i=0; i<NX; i++) for (j=0; j<NY; j++) 
                wave[i*NY+j].p_cfield[movie+2] = wave[i*NY+j].p_cfield[movie];
    }
}


void compute_cfield(short int xy_in[NX*NY], int cplot, int palette, t_wave wave[NX*NY], int fade, double fade_value, int movie)
/* compute the colors */
{
    int i, j, k;
    double ca;
       
    #pragma omp parallel for private(i,j,ca)
    for (i=0; i<NX; i++) for (j=0; j<NY; j++)
    {
        compute_field_color(*wave[i*NY+j].p_cfield[movie], *wave[i*NY+j].p_cfield[movie+2], cplot, palette, wave[i*NY+j].rgb);
        if (SHADE_3D)
        {
            ca = wave[i*NY+j].cos_angle;
            ca = (ca + 1.0)*0.4 + 0.2;
//             if ((FADE_IN_OBSTACLE)&&(!xy_in[i*NY+j])) ca *= 1.6;
            if ((FADE_IN_OBSTACLE)&&(!xy_in[i*NY+j])) ca = (ca + 0.1)*1.6;
            for (k=0; k<3; k++) wave[i*NY+j].rgb[k] *= ca;
        }
        if (fade) for (k=0; k<3; k++) wave[i*NY+j].rgb[k] *= fade_value;
    }
}


void draw_wave_3d_ij(int i, int j, int movie, double phi[NX*NY], double psi[NX*NY], short int xy_in[NX*NY], t_wave wave[NX*NY], int zplot, int cplot, int palette, int fade, double fade_value)
/* draw wave at simulation grid point (i,j) */
{
    int k, l, draw = 1;
    double xy[2], xy2[2], xy_screen[2], rgb[3], pos[2], ca, rgb_e[3], rgb_w[3], rgb_n[3], rgb_s[3]; 
    double z_sw, z_se, z_nw, z_ne, z_mid, zw, ze, zn, zs, min = 1000.0, max = 0.0;
    double xy_sw[2], xy_se[2], xy_nw[2], xy_ne[2], xy_mid[2];
    double energy, zfloor = -10.0;
    
    if (NON_DIRICHLET_BC) 
        draw = (xy_in[i*NY+j])&&(xy_in[(i+1)*NY+j])&&(xy_in[i*NY+j+1])&&(xy_in[(i+1)*NY+j+1]);
    else draw = (TWOSPEEDS)||(xy_in[i*NY+j]);
            
    if (FLOOR_ZCOORD) 
        draw = (draw)&&(*wave[i*NY+j].p_zfield[movie] > zfloor)&&(*wave[(i+1)*NY+j].p_zfield[movie] > zfloor)&&(*wave[i*NY+j+1].p_zfield[movie] > zfloor)&&(*wave[(i+1)*NY+j+1].p_zfield[movie] > zfloor);
            
    if (draw)
    {
        if (AMPLITUDE_HIGH_RES > 0)
        {
            z_mid = compute_interpolated_colors_wave(i, j, xy_in, wave, palette, cplot, 
                                                        rgb_e, rgb_w, rgb_n, rgb_s, fade, fade_value, movie);
            /* TODO: incorporate these in wave structure */
            ij_to_xy(i, j, xy_sw);
            ij_to_xy(i+1, j, xy_se);
            ij_to_xy(i, j+1, xy_nw);
            ij_to_xy(i+1, j+1, xy_ne);
            
//             if (ADD_POTENTIAL) 
//             {
//                 z_mid += *wave[i*NY+j].potential*POT_FACT*0.25;
//                 z_mid += *wave[(i+1)*NY+j].potential*POT_FACT*0.25;
//                 z_mid += *wave[i*NY+j+1].potential*POT_FACT*0.25;
//                 z_mid += *wave[(i+1)*NY+j+1].potential*POT_FACT*0.25;
//             }
                    
            /* TODO: incorporate these in wave structure */
            for (k=0; k<2; k++) xy_mid[k] = 0.25*(xy_sw[k] + xy_se[k] + xy_nw[k] + xy_ne[k]);
                       
            if (AMPLITUDE_HIGH_RES == 1)
            {
                if (ADD_POTENTIAL)
                {
                    glBegin(GL_TRIANGLE_FAN);
                    glColor3f(rgb_w[0], rgb_w[1], rgb_w[2]);
//                     draw_vertex_xyz(xy_mid, z_mid);
                    draw_vertex_xyz(xy_nw, *wave[i*NY+j+1].p_zfield[movie] - *wave[i*NY+j+1].potential*POT_FACT);
                    draw_vertex_xyz(xy_sw, *wave[i*NY+j].p_zfield[movie] - *wave[i*NY+j].potential*POT_FACT);
            
                    glColor3f(rgb_s[0], rgb_s[1], rgb_s[2]);
                    draw_vertex_xyz(xy_se, *wave[(i+1)*NY+j].p_zfield[movie] - *wave[(i+1)*NY+j].potential*POT_FACT);
                
                    glColor3f(rgb_e[0], rgb_e[1], rgb_e[2]);
                    draw_vertex_xyz(xy_ne, *wave[(i+1)*NY+j+1].p_zfield[movie] - *wave[(i+1)*NY+j+1].potential*POT_FACT);
                    
                    glColor3f(rgb_n[0], rgb_n[1], rgb_n[2]);
                    draw_vertex_xyz(xy_nw, *wave[i*NY+j+1].p_zfield[movie] - *wave[i*NY+j].potential*POT_FACT);
                    glEnd ();
                    
                }
                else
                {
                    glBegin(GL_TRIANGLE_FAN);
                    glColor3f(rgb_w[0], rgb_w[1], rgb_w[2]);
                    draw_vertex_xyz(xy_mid, z_mid);
                    draw_vertex_xyz(xy_nw, *wave[i*NY+j+1].p_zfield[movie]);
                    draw_vertex_xyz(xy_sw, *wave[i*NY+j].p_zfield[movie]);
            
                    glColor3f(rgb_s[0], rgb_s[1], rgb_s[2]);
                    draw_vertex_xyz(xy_se, *wave[(i+1)*NY+j].p_zfield[movie]);
                
                    glColor3f(rgb_e[0], rgb_e[1], rgb_e[2]);
                    draw_vertex_xyz(xy_ne, *wave[(i+1)*NY+j+1].p_zfield[movie]);
                    
                    glColor3f(rgb_n[0], rgb_n[1], rgb_n[2]);
                    draw_vertex_xyz(xy_nw, *wave[i*NY+j+1].p_zfield[movie]);
                    glEnd ();
                }
            }
            else /* experimental */
            {
                glColor3f(rgb_w[0], rgb_w[1], rgb_w[2]);
                glBegin(GL_TRIANGLE_STRIP);
                draw_vertex_xyz(xy_mid, z_mid);
                draw_vertex_xyz(xy_nw, *wave[i*NY+j+1].p_zfield[movie]);
                draw_vertex_xyz(xy_sw, *wave[i*NY+j].p_zfield[movie]);
                glEnd ();
                    
                glColor3f(rgb_s[0], rgb_s[1], rgb_s[2]);
                glBegin(GL_TRIANGLE_STRIP);
                draw_vertex_xyz(xy_mid, z_mid);
                draw_vertex_xyz(xy_sw, *wave[i*NY+j].p_zfield[movie]);
                draw_vertex_xyz(xy_se, *wave[(i+1)*NY+j].p_zfield[movie]);
                glEnd ();
                    
                glColor3f(rgb_e[0], rgb_e[1], rgb_e[2]);
                glBegin(GL_TRIANGLE_STRIP);
                draw_vertex_xyz(xy_mid, z_mid);
                draw_vertex_xyz(xy_se, *wave[(i+1)*NY+j].p_zfield[movie]);
                draw_vertex_xyz(xy_ne, *wave[(i+1)*NY+j+1].p_zfield[movie]);
                glEnd ();
                    
                glColor3f(rgb_n[0], rgb_n[1], rgb_n[2]);
                glBegin(GL_TRIANGLE_STRIP);
                draw_vertex_xyz(xy_mid, z_mid);
                draw_vertex_xyz(xy_nw, *wave[i*NY+j+1].p_zfield[movie]);
                draw_vertex_xyz(xy_ne, *wave[(i+1)*NY+j+1].p_zfield[movie]);
                glEnd ();
            }
        }
        else
        {
            glColor3f(wave[i*NY+j].rgb[0], wave[i*NY+j].rgb[1], wave[i*NY+j].rgb[2]);
            glBegin(GL_TRIANGLE_FAN);
            ij_to_xy(i, j, xy);
            draw_vertex_xyz(xy, *wave[i*NY+j].p_zfield[movie]);
            ij_to_xy(i+1, j, xy);
            draw_vertex_xyz(xy, *wave[(i+1)*NY+j].p_zfield[movie]);
            ij_to_xy(i+1, j+1, xy);
            draw_vertex_xyz(xy, *wave[(i+1)*NY+j+1].p_zfield[movie]);
            ij_to_xy(i, j+1, xy);
            draw_vertex_xyz(xy, *wave[i*NY+j+1].p_zfield[movie]);
            glEnd ();
        }
    }
    
    if ((DRAW_OUTSIDE_GRAY)&&((!xy_in[i*NY+j])))
    {
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


void draw_wave_3d(int movie, double phi[NX*NY], double psi[NX*NY], short int xy_in[NX*NY], t_wave wave[NX*NY], int zplot, int cplot, int palette, int fade, double fade_value, int refresh)
{
    int i, j;
    double observer_angle;
    
    blank();
    if (DRAW_BILLIARD) draw_billiard_3d(fade, fade_value);
            
    if (refresh)
    {
        compute_wave_fields(phi, psi, xy_in, zplot, cplot, wave);
        if (SHADE_3D) compute_light_angle(xy_in, wave, movie);
        compute_cfield(xy_in, cplot, palette, wave, fade, fade_value, movie);
    }
    
    if (!ROTATE_VIEW)
    {
        for (i=0; i<NX-2; i++)
            for (j=0; j<NY-2; j++)
                draw_wave_3d_ij(i, j, movie, phi, psi, xy_in, wave, zplot, cplot, palette, fade, fade_value);
    }
    else    /* draw facets in an order depending on the position of the observer */
    {
        observer_angle = argument(observer[0], observer[1]);
        observer_angle -= 0.5*PID;
        if (observer_angle < 0.0) observer_angle += DPI;
        printf("Observer_angle = %.3lg\n", observer_angle*360.0/DPI); 
        
        if ((observer_angle > 0.0)&&(observer_angle < PID))
        {
            for (j=1; j<NY-2; j++)
                for (i=1; i<NX-2; i++)
                    draw_wave_3d_ij(i, j, movie, phi, psi, xy_in, wave, zplot, cplot, palette, fade, fade_value);
        }
        else if (observer_angle < PI)
        {
            for (i=NX-3; i>0; i--)
                for (j=1; j<NY-2; j++)
                    draw_wave_3d_ij(i, j, movie, phi, psi, xy_in, wave, zplot, cplot, palette, fade, fade_value);
        }
        else if (observer_angle < 1.5*PI)
        {
             for (j=NY-3; j>0; j--)
                for (i=1; i<NX-2; i++)
                    draw_wave_3d_ij(i, j, movie, phi, psi, xy_in, wave, zplot, cplot, palette, fade, fade_value);   
        }
        else
        {
            for (i=1; i<NX-2; i++)
                for (j=1; j<NY-2; j++)
                    draw_wave_3d_ij(i, j, movie, phi, psi, xy_in, wave, zplot, cplot, palette, fade, fade_value);
        }
    }
        
    if (DRAW_BILLIARD_FRONT) draw_billiard_3d_front(fade, fade_value);
}

void draw_color_scheme_palette_3d(double x1, double y1, double x2, double y2, int plot, 
                                  double min, double max, int palette, int fade, double fade_value)
{
    int j, k, ij_botleft[2], ij_topright[2], imin, imax, jmin, jmax;
    double y, dy, dy_e, dy_phase, rgb[3], value, lum, amp;
    
//     printf("Drawing color bar\n");
    
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
                value = LOG_SCALE*dy_e*(double)(j - jmin)*100.0/E_SCALE;
                color_scheme_asym_palette(COLOR_SCHEME, palette, value, 1.0, 1, rgb);
                break;
            }
            case (P_3D_TOTAL_ENERGY):
            {
                value = dy_e*(double)(j - jmin)*100.0/E_SCALE;
                if (COLOR_PALETTE >= COL_TURBO) color_scheme_asym_palette(COLOR_SCHEME, palette, value, 1.0, 1, rgb);
                else color_scheme_palette(COLOR_SCHEME, palette, value, 1.0, 1, rgb);
                break;
            }
            case (P_3D_LOG_TOTAL_ENERGY):
            {
                value = LOG_SCALE*dy_e*(double)(j - jmin)*100.0/E_SCALE;
                color_scheme_asym_palette(COLOR_SCHEME, palette, value, 1.0, 1, rgb);
                break;
            }
            case (P_3D_MEAN_ENERGY):
            {
                value = dy_e*(double)(j - jmin)*100.0/E_SCALE;
                if (COLOR_PALETTE >= COL_TURBO) color_scheme_asym_palette(COLOR_SCHEME, palette, value, 1.0, 1, rgb);
                else color_scheme_palette(COLOR_SCHEME, palette, value, 1.0, 1, rgb);
                break;
            }
            case (P_3D_LOG_MEAN_ENERGY):
            {
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
            case (P_3D_FLUX_INTENSITY):
            {
                value = dy_e*(double)(j - jmin)*100.0/E_SCALE;
                if (COLOR_PALETTE >= COL_TURBO) color_scheme_asym_palette(COLOR_SCHEME, palette, value, 1.0, 1, rgb);
                else color_scheme_palette(COLOR_SCHEME, palette, value, 1.0, 1, rgb);
                break;
            }
            case (P_3D_FLUX_DIRECTION):
            {
                value = dy_phase*(double)(j - jmin);
                color_scheme_palette(C_ONEDIM_LINEAR, palette, value, 1.0, 1, rgb);
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
            case (Z_EULER_LPRESSURE):
            {
                value = min + 1.0*dy*(double)(j - jmin);
                color_scheme_palette(COLOR_SCHEME, palette, 0.7*value, 1.0, 0, rgb);
                break;
            }
            case (Z_EULERC_VORTICITY):
             {
                value = min + 1.0*dy*(double)(j - jmin);
                printf("Palette value %.3lg\n", value);
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

void draw_circular_color_scheme_palette_3d(double x1, double y1, double radiusx, double radiusy, int plot, double min, double max, int palette, int fade, double fade_value)
{
    int j, k, ij[2], jmin=0;
    double x, y, dy, dy_e, dy_phase, rgb[3], value, lum, amp, dphi, pos[2], phi, xy[2];
    
    printf("Drawing circular color scheme\n");
    
    glBegin(GL_TRIANGLE_FAN);
//     xy_to_pos(x1, y1, xy);
//     glVertex2d(xy[0], xy[1]);
    xy_to_ij(x1, y1, ij);
    draw_vertex_ij(ij[0], ij[1]);
    dy = (max - min)/360.0;
    dy_e = max/360.0;
    dy_phase = 1.0/360.0;
    dphi = DPI/360.0;
    
    for (j = 0; j <= 360; j++)
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
                value = LOG_SCALE*dy_e*(double)(j - jmin)*100.0/E_SCALE;
                color_scheme_asym_palette(COLOR_SCHEME, palette, value, 1.0, 1, rgb);
                break;
            }
            case (P_3D_TOTAL_ENERGY):
            {
                value = dy_e*(double)(j - jmin)*100.0/E_SCALE;
                if (COLOR_PALETTE >= COL_TURBO) color_scheme_asym_palette(COLOR_SCHEME, palette, value, 1.0, 1, rgb);
                else color_scheme_palette(COLOR_SCHEME, palette, value, 1.0, 1, rgb);
                break;
            }
            case (P_3D_LOG_TOTAL_ENERGY):
            {
                value = LOG_SCALE*dy_e*(double)(j - jmin)*100.0/E_SCALE;
                color_scheme_asym_palette(COLOR_SCHEME, palette, value, 1.0, 1, rgb);
                break;
            }
            case (P_3D_MEAN_ENERGY):
            {
                value = dy_e*(double)(j - jmin)*100.0/E_SCALE;
                if (COLOR_PALETTE >= COL_TURBO) color_scheme_asym_palette(COLOR_SCHEME, palette, value, 1.0, 1, rgb);
                else color_scheme_palette(COLOR_SCHEME, palette, value, 1.0, 1, rgb);
                break;
            }
            case (P_3D_LOG_MEAN_ENERGY):
            {
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
            case (P_3D_FLUX_INTENSITY):
            {
                value = dy_e*(double)(j - jmin)*100.0/E_SCALE;
                if (COLOR_PALETTE >= COL_TURBO) color_scheme_asym_palette(COLOR_SCHEME, palette, value, 1.0, 1, rgb);
                else color_scheme_palette(COLOR_SCHEME, palette, value, 1.0, 1, rgb);
                break;
            }
            case (P_3D_FLUX_DIRECTION):
            {
                value = dy_phase*(double)(j - jmin);
                color_scheme_palette(C_ONEDIM_LINEAR, palette, value, 1.0, 1, rgb);
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
            case (Z_EULER_LPRESSURE):
            {
                value = min + 1.0*dy*(double)(j - jmin);
                color_scheme_palette(COLOR_SCHEME, palette, 0.7*value, 1.0, 0, rgb);
                break;
            }
            case (Z_EULERC_VORTICITY):
             {
                value = min + 1.0*dy*(double)(j - jmin);
                printf("Palette value %.3lg\n", value);
                color_scheme_palette(COLOR_SCHEME, palette, 0.7*value, 1.0, 0, rgb);
                break;
            }
               
        }
        if (fade) for (k=0; k<3; k++) rgb[k] *= fade_value;
        glColor3f(rgb[0], rgb[1], rgb[2]);
        
        xy_to_ij(x1 + radiusx*cos(dphi*(double)j), y1 + radiusy*sin(dphi*(double)j), ij);
        draw_vertex_ij(ij[0], ij[1]);
    }
    
    xy_to_ij(x1 + radiusx*cos(dphi), y1 + radiusy*sin(dphi), ij);
    draw_vertex_ij(ij[0], ij[1]);
    glEnd ();
    
    if (fade) glColor3f(fade_value, fade_value, fade_value);
    else glColor3f(1.0, 1.0, 1.0);
    glLineWidth(BOUNDARY_WIDTH*3/2);
    glEnable(GL_LINE_SMOOTH);
    
    dphi = DPI/(double)NSEG;
    glBegin(GL_LINE_LOOP);
    for (j = 0; j < NSEG; j++)
    {               
        xy_to_ij(x1 + radiusx*cos(dphi*(double)j), y1 + radiusy*sin(dphi*(double)j), ij);
        draw_vertex_ij(ij[0], ij[1]);
    }
    glEnd ();
}


void print_speed_3d(double speed, int fade, double fade_value)
{
    char message[100];
    double y = YMAX - 0.1, pos[2];
    static double xleftbox, xlefttext;
    static int first = 1;
    
    if (first)
    {
        xleftbox = XMIN + 0.3;
        xlefttext = xleftbox - 0.45;
        first = 0;
    }
    
    erase_area_hsl(xleftbox, y + 0.025, 0.22, 0.05, 0.0, 0.9, 0.0);
    if (fade) glColor3f(fade_value, fade_value, fade_value);
    else glColor3f(1.0, 1.0, 1.0);
//     xy_to_pos(xlefttext + 0.28, y, pos);
    sprintf(message, "Mach %.3lg", speed);
    write_text(xlefttext + 0.28, y, message);
}

void init_wave_fields(t_wave wave[NX*NY])
/* initialize some auxiliary fields */
{
    int i, k;
    
    #pragma omp parallel for private(i)
    for (i=0; i<NX*NY; i++)
    {
        wave[i].total_energy = 0.0;
        wave[i].flux_counter = 0;
        for (k=0; k<FLUX_WINDOW; k++)
            wave[i].flux_int_table[k] = 0.0;
    }
    
}
