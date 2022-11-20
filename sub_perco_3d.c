/* routines for 3D representation of percolation clusters, taken from sub_wave_3d_rde.c */

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

// void draw_vertex_ij(int i, int j)
// {
//     double xy[2];
//     
//     ij_to_xy(i, j, xy);
// //     if (xy[1] > 0.0) printf("(i,j) = (%i,%i), (x,y) = (%.2lg,%.2lg)\n", i, j, xy[0], xy[1]);
//     glVertex2d(xy[0], xy[1]);
// }


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

void draw_line_3d(double x1, double y1, double z1, double x2, double y2, double z2)
{
    glBegin(GL_LINE_LOOP);
    draw_vertex_x_y_z(x1, y1, z1);
    draw_vertex_x_y_z(x2, y2, z2);
    glEnd();    
}


double angle_lum(double cosangle)
{
    double mid_lum = 0.6;
    
    if (cosangle > 0.0) return(mid_lum + (1.0 - mid_lum)*cosangle);
    else return (mid_lum + 0.8*(1.0 - mid_lum)*cosangle);
}

void draw_cube(double x, double y, double z, double size, double rgb[3])
/* draw a cube */
{
    double lum_factor;
    
    /* front/back face */
    if (observer[0] > 0.0)
    {
        lum_factor = angle_lum(light[0]);
        glColor3f(rgb[0]*lum_factor, rgb[1]*lum_factor, rgb[2]*lum_factor);
        glBegin(GL_TRIANGLE_FAN);
        draw_vertex_x_y_z(x + size, y, z);
        draw_vertex_x_y_z(x + size, y + size, z);
        draw_vertex_x_y_z(x + size, y + size, z + size);
        draw_vertex_x_y_z(x + size, y, z + size);
        glEnd();
    }
    else
    {
//         cosangle = -light[0];
        lum_factor = angle_lum(-light[0]);
        glColor3f(rgb[0]*lum_factor, rgb[1]*lum_factor, rgb[2]*lum_factor);
        glBegin(GL_TRIANGLE_FAN);
        draw_vertex_x_y_z(x, y, z);
        draw_vertex_x_y_z(x, y + size, z);
        draw_vertex_x_y_z(x, y + size, z + size);
        draw_vertex_x_y_z(x, y, z + size);
        glEnd();
    }
        
    
    /* right/left face */
    if (observer[1] > 0.0)
    {
//         cosangle = light[1];
        lum_factor = angle_lum(light[1]);
        glColor3f(rgb[0]*lum_factor, rgb[1]*lum_factor, rgb[2]*lum_factor);
        glBegin(GL_TRIANGLE_FAN);
        draw_vertex_x_y_z(x + size, y + size, z);
        draw_vertex_x_y_z(x, y + size, z);
        draw_vertex_x_y_z(x, y + size, z + size);
        draw_vertex_x_y_z(x + size, y + size, z + size);
        glEnd();
    }
    else
    {
//         cosangle = -light[1];
        lum_factor = angle_lum(-light[1]);
        glColor3f(rgb[0]*lum_factor, rgb[1]*lum_factor, rgb[2]*lum_factor);
        glBegin(GL_TRIANGLE_FAN);
        draw_vertex_x_y_z(x + size, y, z);
        draw_vertex_x_y_z(x, y, z);
        draw_vertex_x_y_z(x, y, z + size);
        draw_vertex_x_y_z(x + size, y, z + size);
        glEnd();
    }
        
    
    /* top face */
//     cosangle = light[2];
    lum_factor = angle_lum(light[2]);
    glColor3f(rgb[0]*lum_factor, rgb[1]*lum_factor, rgb[2]*lum_factor);
    glBegin(GL_TRIANGLE_FAN);
    draw_vertex_x_y_z(x + size, y, z + size);
    draw_vertex_x_y_z(x + size, y + size, z + size);
    draw_vertex_x_y_z(x, y + size, z + size);
    draw_vertex_x_y_z(x, y, z + size);
    glEnd();
}

void viewpoint_schedule(int i)
/* change position of observer */
{
    int j;
    double angle, ca, sa;
    static double observer_initial[3];
    static int first = 1;
    
    if (first)
    {
        for (j=0; j<3; j++) observer_initial[j] = observer[j];
        first = 0;
    }
    
    angle = (ROTATE_ANGLE*DPI/360.0)*(double)i/(double)NSTEPS;
    ca = cos(angle);
    sa = sin(angle);
    observer[0] = ca*observer_initial[0] - sa*observer_initial[1];
    observer[1] = sa*observer_initial[0] + ca*observer_initial[1];
    printf("Angle %.3lg, Observer position (%.3lg, %.3lg, %.3lg)\n", angle, observer[0], observer[1], observer[2]);
}
