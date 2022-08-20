/*********************************************************************************/
/*                                                                               */
/*  Simulation of percolation in 2D                                              */
/*                                                                               */
/*  N. Berglund, July 2022                                                       */
/*                                                                               */
/*  Feel free to reuse, but if doing so it would be nice to drop a               */
/*  line to nils.berglund@univ-orleans.fr - Thanks!                              */
/*                                                                               */
/*  compile with                                                                 */
/*  gcc -o percolation percolation.c                                             */
/* -L/usr/X11R6/lib -ltiff -lm -lGL -lGLU -lX11 -lXmu -lglut -O3 -fopenmp        */
/*                                                                               */
/*  OMP acceleration may be more effective after executing                       */
/*  export OMP_NUM_THREADS=2 in the shell before running the program             */
/*                                                                               */
/*  To make a video, set MOVIE to 1 and create subfolder tif_perc                */
/*  It may be possible to increase parameter PAUSE                               */
/*                                                                               */
/*  create movie using                                                           */
/*  ffmpeg -i perc.%05d.tif -vcodec libx264 perc.mp4                             */
/*                                                                               */
/*********************************************************************************/

/*********************************************************************************/
/*                                                                               */
/* NB: The algorithm used to simulate the wave equation is highly paralellizable */
/* One could make it much faster by using a GPU                                  */
/*                                                                               */
/*********************************************************************************/

#include <math.h>
#include <string.h>
#include <GL/glut.h>
#include <GL/glu.h>
#include <unistd.h>
#include <sys/types.h>
#include <tiffio.h>     /* Sam Leffler's libtiff library. */
#include <omp.h>
#include <time.h>

#define MOVIE 1         /* set to 1 to generate movie */

/* General geometrical parameters */

// #define WINWIDTH 	1920  /* window width */
// #define WINHEIGHT 	1000  /* window height */
// #define NX 1920          /* number of grid points on x axis */
// #define NY 992           /* number of grid points on y axis */
// 
// #define XMIN -2.0
// #define XMAX 2.0	/* x interval  */
// #define YMIN -1.041666667
// #define YMAX 1.041666667	/* y interval for 9/16 aspect ratio */

#define HIGHRES 0        /* set to 1 if resolution of grid is double that of displayed image */

#define WINWIDTH 	1280  /* window width */
#define WINHEIGHT 	720   /* window height */

#define NX 1280          /* number of grid points on x axis */
#define NY 720          /* number of grid points on y axis */

#define XMIN -2.0
#define XMAX 2.0	/* x interval  */
#define YMIN -1.125
#define YMAX 1.125	/* y interval for 9/16 aspect ratio */

/* Boundary conditions, see list in global_pdes.c  */

#define LATTICE 2

#define FLOOD_LEFT_BOUNDARY 0   /* set to 1 to flood cells on left boundary */
#define FIND_ALL_CLUSTERS 1     /* set to 1 to find all open clusters */

#define PLOT_CLUSTER_SIZE 0     /* set to 1 to add a plot for the size of the percolation cluster */
#define PLOT_CLUSTER_NUMBER 0   /* set to 1 to add a graph of the number of clusters */
#define PLOT_CLUSTER_HISTOGRAM 1    /* set to 1 to add a histogram of the number of clusters */
#define PRINT_LARGEST_CLUSTER_SIZE 1    /* set to 1 to print size of largest cluster */

#define MAX_CLUSTER_NUMBER 6    /* vertical scale of the cluster number plot */
#define HISTO_BINS 30           /* number of bins in histogram */

#define NSTEPS 100        /* number of frames of movie */
// #define NSTEPS 700        /* number of frames of movie */
// #define NSTEPS 830        /* number of frames of movie */

#define PAUSE 200       /* number of frames after which to pause */
#define PSLEEP 2         /* sleep time during pause */
#define SLEEP1  1        /* initial sleeping time */
#define SLEEP2  1        /* final sleeping time */
#define MID_FRAMES 20    /* number of still frames between parts of two-part movie */
#define END_FRAMES 100    /* number of still frames at end of movie */
#define FADE 1           /* set to 1 to fade at end of movie */

/* Color schemes */

#define COLOR_PALETTE 15     /* Color palette, see list in global_pdes.c  */

#define BLACK 1          /* background */

#define COLOR_CLUSTERS_BY_SIZE 1    /* set to 1 to link cluster color to their size */

#define SCALE 0          /* set to 1 to adjust color scheme to variance of field */
#define SLOPE 1.0        /* sensitivity of color on wave amplitude */
#define ATTENUATION 0.0  /* exponential attenuation coefficient of contrast with time */

#define HUE_CLOSED 350.0   /* color hue of closed cells */
#define HUE_OPEN 200.0     /* color hue of open (dry) cells */
#define HUE_FLOODED 45.0   /* color hue of open flooded cells */
#define HUE_GRAPH 250.0    /* color hue in graph of cluster size */

#define CLUSTER_HUEMIN 60.0     /* minimal color hue of clusters */
#define CLUSTER_HUEMAX 300.0    /* maximal color hue of clusters */
#define N_CLUSTER_COLORS 20     /* number of different colors of clusters */

#define COLORHUE 260     /* initial hue of water color for scheme C_LUM */
#define COLORDRIFT 0.0   /* how much the color hue drifts during the whole simulation */
#define LUMMEAN 0.5      /* amplitude of luminosity variation for scheme C_LUM */
#define LUMAMP 0.3       /* amplitude of luminosity variation for scheme C_LUM */
#define HUEMEAN 180.0    /* mean value of hue for color scheme C_HUE */
#define HUEAMP -180.0      /* amplitude of variation of hue for color scheme C_HUE */

/* debugging options */
#define VERBOSE 0       /* set to 1 to print more messages in shell */
#define DEBUG 0         /* set to 1 for some debugging features */
#define DEBUG_SLEEP_TIME 1  /* sleep time between frames when debugging */
#define TEST_GRAPH 0   /* set to 1 to test graph connectivity matrix */
#define TEST_START 2210    /* start position of connectivity test */

#define ADD_PLOT ((PLOT_CLUSTER_SIZE)||(PLOT_CLUSTER_NUMBER)||(PLOT_CLUSTER_HISTOGRAM))
#define FIND_CLUSTER_SIZES ((COLOR_CLUSTERS_BY_SIZE)||(PLOT_CLUSTER_HISTOGRAM))

#include "global_perc.c"        /* constants and global variables */
#include "sub_perco.c"          



void animation(int size)
{
    int i, j, k, s, nx, ny, nmaxcells, maxsize, nopen, nflooded, nstack, nclusters, maxclustersize = 0, maxclusterlabel;
    int *plot_cluster_number, *cluster_sizes;
    int ncells;
    double p, *plot_cluster_size;
    t_perco *cell;
    t_perco **pstack;
    
    compute_nxny(size, &nx, &ny);
    
    nmaxcells = cell_number(NX, NY);
    
    cell = (t_perco *)malloc(nmaxcells*sizeof(t_perco));
    if (PLOT_CLUSTER_SIZE) plot_cluster_size = (double *)malloc(NSTEPS*sizeof(double));
    if (PLOT_CLUSTER_NUMBER) plot_cluster_number = (int *)malloc(NSTEPS*sizeof(double));
//     if (FIND_CLUSTER_SIZES) 
        cluster_sizes = (int *)malloc(2*nmaxcells*sizeof(int));
    
    ncells = init_cell_lattice(cell, nx, ny);
    printf("nx = %i, ny = %i, ncells = %i, maxcells = %i\n", nx, ny, ncells, nmaxcells);
    
    pstack = (t_perco* *)malloc(ncells*sizeof(struct t_perco *));
    
    if (TEST_GRAPH) test_neighbours(TEST_START, cell, nx, ny, size, ncells);
    
    init_cell_probabilities(cell, ncells);
    
    for (i=0; i<NSTEPS; i++)
    {
        p = p_schedule(i);
        printf("\ni = %i, p = %.4lg\n", i, p);
        
        init_cell_state(cell, p, ncells, (i == 0));
        
        if (FLOOD_LEFT_BOUNDARY) nstack = init_flooded_cells(cell, ncells, nx, ny, pstack);
        nopen = count_open_cells(cell, ncells);
        
        printf("Flooded cells, %i open cells, nstack = %i\n", nopen, nstack);
        
        if (FLOOD_LEFT_BOUNDARY) 
        {
            nflooded = find_percolation_cluster(cell, ncells, pstack, nstack);
            printf("Found percolation cluster with %i flooded cells\n", nflooded);
        }

        if (FIND_ALL_CLUSTERS) 
        {
            nclusters = find_all_clusters(cell, ncells, (i == 0), &maxclusterlabel);
            printf("Found %i clusters\n", nclusters);
        }
        
        if (FIND_CLUSTER_SIZES) 
        {
            maxclustersize = find_cluster_sizes(cell, ncells, cluster_sizes, &maxclusterlabel); 
            printf("Max cluster size %i, max cluster label %i, ncells %i\n", maxclustersize, maxclusterlabel, ncells);
        }
        
//         print_cluster_sizes(cell, ncells, cluster_sizes);
        
        draw_configuration(cell, cluster_sizes, ncells, nx, ny, size, ncells);
        print_p(p);
        if (PRINT_LARGEST_CLUSTER_SIZE) print_largest_cluster_size(maxclustersize);
        
        printf("%i open cells, %i flooded cells out of %i cells\n", nopen, nflooded, ncells);
        
        if (PLOT_CLUSTER_SIZE) 
        {
            plot_cluster_size[i] = (double)(nflooded)/(double)(nopen);
            draw_size_plot(plot_cluster_size, i,  pcritical(LATTICE));
        }
        
        if (PLOT_CLUSTER_NUMBER) 
        {
            plot_cluster_number[i] = nclusters;
            draw_cluster_number_plot(plot_cluster_number, ncells/MAX_CLUSTER_NUMBER, i,  pcritical(LATTICE));
            print_nclusters(nclusters);
        }
        
        if (PLOT_CLUSTER_HISTOGRAM) draw_cluster_histogram(ncells, cluster_sizes, maxclustersize, maxclusterlabel);
                
        glutSwapBuffers();
        
        if (DEBUG) 
        {
            printf("\n\n");
            sleep(DEBUG_SLEEP_TIME);
        }
        
        if (MOVIE) save_frame_perc();
    }
    
    if (MOVIE) 
    {
        for (i=0; i<MID_FRAMES; i++) save_frame_perc();
        s = system("mv perc*.tif tif_perc/");
    }
    
    free(cell);
    free(pstack);
    
    if (PLOT_CLUSTER_SIZE) free(plot_cluster_size);
    if (PLOT_CLUSTER_NUMBER) free(plot_cluster_number);
//     if (FIND_CLUSTER_SIZES) 
        free(cluster_sizes);
}




void display(void)
{
    time_t rawtime;
    struct tm * timeinfo;

    time(&rawtime);
    timeinfo = localtime(&rawtime);
    
    glPushMatrix();

    blank();
    glutSwapBuffers();
    blank();
    glutSwapBuffers();

//     animation(128);
//     animation(64);
//     animation(32);
//     animation(16);
//     animation(8);
    animation(4);
//     animation(2);
//     animation(1);
    
    
    sleep(SLEEP2);

    glPopMatrix();

    glutDestroyWindow(glutGetWindow());
    
    printf("Start local time and date: %s", asctime(timeinfo));
    time(&rawtime);
    timeinfo = localtime(&rawtime);
    printf("Current local time and date: %s", asctime(timeinfo));
}


int main(int argc, char** argv)
{
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH);
    glutInitWindowSize(WINWIDTH,WINHEIGHT);
    glutCreateWindow("Percolation in a planar domain");

    init();

    glutDisplayFunc(display);

    glutMainLoop();

    return 0;
}

