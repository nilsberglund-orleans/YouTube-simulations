/* The function init_maze has been improved and should return a maze with a solution  */
/* The current algorithm uses a self-avoiding random walk. A better option may be */
/* to give random weights to the dual graph, and finite a maximal spanning tree   */

/* Change constant RAND_SHIFT to change the maze */

#define MAZE_TYPE_SQUARE 0      /* maze with square cells */
#define MAZE_TYPE_CIRCLE 1      /* circular maze */
#define MAZE_TYPE_HEX 2         /* honeycomb maze */
#define MAZE_TYPE_OCT 3         /* maze with octagonal and square cells */

typedef struct
{
    short int nneighb;                      /* number of neighbours */
    int neighb[MAZE_MAX_NGBH];              /* neighbour cells */
    short int directions[MAZE_MAX_NGBH];    /* direction of neighbours */
    short int north, east, south, west;     /* closed walls */
    short int northeast, northwest, southeast, southwest;  /* closed walls */
    short int active;                       /* takes value 1 if currently active in RW path */
    short int tested;                       /* takes value 1 if tested */
    short int connected;                    /* takes value 1 if connected to exit */
    short int closed;                       /* takes value 1 if no untested neighbours */
} t_maze;


int nmaze(int i, int j)
{
    return(NXMAZE*j + i);
}

void init_maze_graph(t_maze maze[NXMAZE*NYMAZE])
{
    int i, j, k, n;
    
    printf("Initializing maze\n");

    /* initialize neighbours */
    /* in the bulk */
    for (i=1; i<NXMAZE-1; i++)
        for (j=1; j<NYMAZE-1; j++)
        {
            n = nmaze(i, j);
            maze[n].nneighb = 4;
            maze[n].neighb[0] = nmaze(i, j+1);
            maze[n].neighb[1] = nmaze(i+1, j);
            maze[n].neighb[2] = nmaze(i, j-1);
            maze[n].neighb[3] = nmaze(i-1, j);
            for (k=0; k<4; k++) maze[n].directions[k] = k;
        }
    
    /* left side */
    for (j=1; j<NYMAZE-1; j++)
    {
        n = nmaze(0, j);
        maze[n].nneighb = 3;
        maze[n].neighb[0] = nmaze(0, j+1);
        maze[n].neighb[1] = nmaze(1, j);
        maze[n].neighb[2] = nmaze(0, j-1);
        for (k=0; k<3; k++) maze[n].directions[k] = k;
    }
    /* right side */
    for (j=1; j<NYMAZE-1; j++)
    {
        n = nmaze(NXMAZE-1, j);
        maze[n].nneighb = 3;
        maze[n].neighb[0] = nmaze(NXMAZE-1, j+1);
        maze[n].neighb[1] = nmaze(NXMAZE-2, j);
        maze[n].neighb[2] = nmaze(NXMAZE-1, j-1);
        maze[n].directions[0] = 0;
        maze[n].directions[1] = 3;
        maze[n].directions[2] = 2;
    }
    /* bottom side */
    for (i=1; i<NXMAZE-1; i++)
    {
        n = nmaze(i, 0);
        maze[n].nneighb = 3;
        maze[n].neighb[0] = nmaze(i, 1);
        maze[n].neighb[1] = nmaze(i+1, 0);
        maze[n].neighb[2] = nmaze(i-1, 0);
        maze[n].directions[0] = 0;
        maze[n].directions[1] = 1;
        maze[n].directions[2] = 3;
    }
    /* top side */
    for (i=1; i<NXMAZE-1; i++)
    {
        n = nmaze(i, NYMAZE-1);
        maze[n].nneighb = 3;
        maze[n].neighb[0] = nmaze(i, NYMAZE-2);
        maze[n].neighb[1] = nmaze(i+1, NYMAZE-1);
        maze[n].neighb[2] = nmaze(i-1, NYMAZE-1);
        maze[n].directions[0] = 2;
        maze[n].directions[1] = 1;
        maze[n].directions[2] = 3;
    }
    /* corners */
    n = nmaze(0,0);
    maze[n].nneighb = 2;
    maze[n].neighb[0] = nmaze(1,0);
    maze[n].neighb[1] = nmaze(0,1);
    maze[n].directions[0] = 1;
    maze[n].directions[1] = 0;
    
    n = nmaze(NXMAZE-1,0);
    maze[n].nneighb = 2;
    maze[n].neighb[0] = nmaze(NXMAZE-2,0);
    maze[n].neighb[1] = nmaze(NXMAZE-1,1);
    maze[n].directions[0] = 3;
    maze[n].directions[1] = 0;

    n = nmaze(0,NYMAZE-1);
    maze[n].nneighb = 2;
    maze[n].neighb[0] = nmaze(1,NYMAZE-1);
    maze[n].neighb[1] = nmaze(0,NYMAZE-2);
    maze[n].directions[0] = 1;
    maze[n].directions[1] = 2;
    
    n = nmaze(NXMAZE-1,NYMAZE-1);
    maze[n].nneighb = 2;
    maze[n].neighb[0] = nmaze(NXMAZE-2,NYMAZE-1);
    maze[n].neighb[1] = nmaze(NXMAZE-1,NYMAZE-2);
    maze[n].directions[0] = 3;
    maze[n].directions[1] = 2;

    /* initialize other parameters */
    for (i=0; i<NXMAZE; i++)
        for (j=0; j<NYMAZE; j++)
        {
            n = nmaze(i, j);
            maze[n].active = 0;
            maze[n].tested = 0;
            maze[n].connected = 0;
            maze[n].closed = 0;
            maze[n].north = 1;
            maze[n].east = 1;
            maze[n].south = 1;
            maze[n].west = 1;
        }
}


void init_circular_maze_graph(t_maze maze[NXMAZE*NYMAZE])
/* initialise graph of circular maze */
/* NXMAZE should be a power of 2, and NYMAZE a multiple of NXMAZE */
/* row number i represents the radial coordinate, and column number j the angular one */
/* the maze is split into square blocks of size NXMAZE times NXMAZE */
/* in each block, the first column has 1 cell, the second has 2 cells */ 
/* then there are 2 columns of 4 cells, 4 columns of 8 cells, etc */
{
    int i, j, k, n, p, q, block, nblocks, width, nextblock, prevblock;
    
    printf("Initializing maze\n");
    if (MAZE_MAX_NGBH < 5)
    {
        printf("Error: MAZE_MAX_NGBH should be at least 5 for circular maze\n");
        exit(0);
    }
    if (NYMAZE%NXMAZE != 0) printf("Warning: NYMAZE should be a multiple of NXMAZE\n");
    
    /* set dummy variables for potentially unused cells */
    for (i=0; i<NXMAZE*NYMAZE; i++) maze[i].nneighb = 0;

    nblocks = NYMAZE/NXMAZE;
    
    /* initialize neighbours */
    for (block=0; block<nblocks; block++)
    {
        nextblock = (block+1)%nblocks;
        prevblock = block-1;    if (prevblock < 0) prevblock = nblocks-1;
        
        /* first column */
        j = block*NXMAZE;
        n = nmaze(0,j);
        maze[n].nneighb = 4;
        maze[n].neighb[0] = nmaze(1, j);
        maze[n].neighb[1] = nmaze(1, j+1);
        maze[n].neighb[2] = nmaze(0, nextblock*NXMAZE);
        maze[n].neighb[3] = nmaze(0, prevblock*NXMAZE);
        maze[n].directions[0] = 1;
        maze[n].directions[1] = 4;
        maze[n].directions[2] = 0;
        maze[n].directions[3] = 2;
        
        /* second column */
        for (q=0; q<2; q++)
        {
            j = block*NXMAZE + q;
            n = nmaze(1,j);
            maze[n].nneighb = 5;
            maze[n].neighb[0] = nmaze(2, j + q);
            maze[n].neighb[1] = nmaze(2, j + q + 1);
            if (q == 1)  maze[n].neighb[2] = nmaze(1, nextblock*NXMAZE);
            else maze[n].neighb[2] = nmaze(1, j+1);
            if (q == 0) maze[n].neighb[3] = nmaze(1, prevblock*NXMAZE + 1);
            else maze[n].neighb[3] = nmaze(1, j-1);
            maze[n].neighb[4] = nmaze(0, block*NXMAZE);
            
            maze[n].directions[0] = 1;
            maze[n].directions[1] = 4;
            maze[n].directions[2] = 0;
            maze[n].directions[3] = 2;
            maze[n].directions[4] = 3;
        }
        
        /* other columns */
        width = 2;
        i = 2;
        while (width < NXMAZE)
        {
            /* left column of block */
            for (q = 0; q < 2*width; q++)
            {
                j = block*NXMAZE + q;
                n = nmaze(i,j);
                maze[n].nneighb = 4;
                maze[n].neighb[0] = nmaze(i+1, j);
                if (q == 2*width-1) maze[n].neighb[1] = nmaze(i, nextblock*NXMAZE);
                else maze[n].neighb[1] = nmaze(i, j+1);
                if (q == 0) maze[n].neighb[2] = nmaze(i, prevblock*NXMAZE + 2*width - 1);
                else maze[n].neighb[2] = nmaze(i, j-1);
                maze[n].neighb[3] = nmaze(i-1, block*NXMAZE + q/2);
            
                maze[n].directions[0] = 1;
                maze[n].directions[1] = 0;
                maze[n].directions[2] = 2;
                maze[n].directions[3] = 3;
            }
            
            /* middle columns of block */
            for (p = 1; p < width-1; p++)
            {
                i++;
                for (q = 0; q < 2*width; q++)
                {
                    j = block*NXMAZE + q;
                    n = nmaze(i,j);
                    maze[n].nneighb = 4;
                    maze[n].neighb[0] = nmaze(i+1, j);
                    if (q == 2*width-1) maze[n].neighb[1] = nmaze(i, nextblock*NXMAZE);
                    else maze[n].neighb[1] = nmaze(i, j+1);
                    if (q == 0) maze[n].neighb[2] = nmaze(i, prevblock*NXMAZE + 2*width - 1);
                    else maze[n].neighb[2] = nmaze(i, j-1);
                    maze[n].neighb[3] = nmaze(i-1, j);
            
                    maze[n].directions[0] = 1;
                    maze[n].directions[1] = 0;
                    maze[n].directions[2] = 2;
                    maze[n].directions[3] = 3;
                }
            }
            
            /* right column of block */
            i++;
            for (q = 0; q < 2*width; q++)
            {
                j = block*NXMAZE + q;
                n = nmaze(i,j);
                if (i<NXMAZE-1) maze[n].nneighb = 5;
                else maze[n].nneighb = 3;
                
                maze[n].neighb[0] = nmaze(i-1, j); 
                if (q == 2*width-1) maze[n].neighb[1] = nmaze(i, nextblock*NXMAZE);
                else maze[n].neighb[1] = nmaze(i, j+1);
                if (q == 0) maze[n].neighb[2] = nmaze(i, prevblock*NXMAZE + 2*width - 1);
                else maze[n].neighb[2] = nmaze(i, j-1);
                    
                maze[n].directions[0] = 3;
                maze[n].directions[1] = 0;
                maze[n].directions[2] = 2;
                
                if (i<NXMAZE-1)
                {
                    printf("i = %i, q = %i, j+ = %i\n", i, q, block*NXMAZE + 2*q);
                    maze[n].neighb[3] = nmaze(i+1, block*NXMAZE + 2*q);
                    printf("i = %i, q = %i, j+ = %i\n", i, q, block*NXMAZE + 2*q + 1);
                    maze[n].neighb[4] = nmaze(i+1, block*NXMAZE + 2*q + 1);
                    
                    maze[n].directions[3] = 1;
                    maze[n].directions[4] = 4;
                }
            }
            i++;
            width *= 2;
        }
    }
    
    /* initialize other parameters */
    for (i=0; i<NXMAZE; i++)
        for (j=0; j<NYMAZE; j++)
        {
            n = nmaze(i, j);
            maze[n].active = 0;
            if (maze[n].nneighb == 0) maze[n].tested = 1;
            else maze[n].tested = 0;
            maze[n].connected = 0;
            maze[n].closed = 0;
            maze[n].north = 1;
            maze[n].east = 1;
            maze[n].south = 1;
            maze[n].west = 1;
            maze[n].northeast = 1;
        }
        
    /* for debugging */
//     for (i=0; i<NXMAZE; i++)
//         for (j=0; j<NYMAZE; j++)
//         {
//             n = nmaze(i, j);
//             q = maze[n].nneighb;
//             if (q > 0) printf("Cell (%i, %i)\n", i, j);
//             for (k = 0; k <q; k++)
//             {
//                 p = maze[n].neighb[k];
//                 printf("Neighbour %i at (%i, %i)\t", k, p%NXMAZE, p/NXMAZE);
//                 printf("Direction %i\n", maze[n].directions[k]);
//             }
//         }
//     sleep(5);
}

void init_hex_maze_graph(t_maze maze[NXMAZE*NYMAZE])
/* initialise graph of maze with honeycomb cells */
/* NXMAZE and NYMAZE are assumed to be even */
/* directions are: 0 - north, 1 - NE, 2 - SE, 3 - south, 4 - SW, 5 - NW */
{
    int i, j, k, n;
    
    printf("Initializing maze\n");
    if (MAZE_MAX_NGBH < 6)
    {
        printf("Error: MAZE_MAX_NGBH should be at least 5 for circular maze\n");
        exit(0);
    }
    
    /* initialize neighbours */
    /* in the bulk */
    for (i=1; i<NXMAZE-1; i++)
        for (j=1; j<NYMAZE-1; j++)
        {
            n = nmaze(i, j);
            maze[n].nneighb = 6;
            maze[n].neighb[0] = nmaze(i, j+1);
            maze[n].neighb[3] = nmaze(i, j-1);
            if (i%2 == 0)
            {
                maze[n].neighb[1] = nmaze(i+1, j+1);
                maze[n].neighb[2] = nmaze(i+1, j);
                maze[n].neighb[4] = nmaze(i-1, j);
                maze[n].neighb[5] = nmaze(i-1, j+1);
            }
            else
            {
                maze[n].neighb[1] = nmaze(i+1, j);
                maze[n].neighb[2] = nmaze(i+1, j-1);
                maze[n].neighb[4] = nmaze(i-1, j-1);
                maze[n].neighb[5] = nmaze(i-1, j);                
            }
            for (k=0; k<6; k++) maze[n].directions[k] = k;
        }
    
    /* left side */
    for (j=1; j<NYMAZE-1; j++)
    {
        n = nmaze(0, j);
        maze[n].nneighb = 4;
        maze[n].neighb[0] = nmaze(0, j+1);
        maze[n].neighb[1] = nmaze(1, j+1);
        maze[n].neighb[2] = nmaze(1, j);
        maze[n].neighb[3] = nmaze(0, j-1);
        for (k=0; k<4; k++) maze[n].directions[k] = k;
    }
    /* right side */
    for (j=1; j<NYMAZE-1; j++)
    {
        n = nmaze(NXMAZE-1, j);
        maze[n].nneighb = 4;
        maze[n].neighb[0] = nmaze(NXMAZE-1, j+1);
        maze[n].neighb[1] = nmaze(NXMAZE-1, j-1);
        maze[n].neighb[2] = nmaze(NXMAZE-2, j-1);
        maze[n].neighb[3] = nmaze(NXMAZE-2, j);
        maze[n].directions[0] = 0;
        maze[n].directions[1] = 3;
        maze[n].directions[2] = 4;
        maze[n].directions[3] = 5;
    }
    /* bottom side */
    for (i=1; i<NXMAZE-1; i++)
    {
        n = nmaze(i, 0);
        if (i%2 == 0)
        {
            maze[n].nneighb = 5;
            maze[n].neighb[0] = nmaze(i, 1);
            maze[n].neighb[1] = nmaze(i+1, 1);
            maze[n].neighb[2] = nmaze(i+1, 0);
            maze[n].neighb[3] = nmaze(i-1, 0);
            maze[n].neighb[4] = nmaze(i-1, 1);
            maze[n].directions[0] = 0;
            maze[n].directions[1] = 1;
            maze[n].directions[2] = 2;
            maze[n].directions[3] = 4;
            maze[n].directions[4] = 5;
        }
        else
        {
            maze[n].nneighb = 3;
            maze[n].neighb[0] = nmaze(i, 1);
            maze[n].neighb[1] = nmaze(i+1, 0);
            maze[n].neighb[2] = nmaze(i-1, 0);
            maze[n].directions[0] = 0;
            maze[n].directions[1] = 1;
            maze[n].directions[2] = 5;
        }
    }
    /* top side */
    for (i=1; i<NXMAZE-1; i++)
    {
        n = nmaze(i, NYMAZE-1);
        if (i%2 == 0)
        {
            maze[n].nneighb = 3;
            maze[n].neighb[0] = nmaze(i+1, NYMAZE-1);
            maze[n].neighb[1] = nmaze(i, NYMAZE-2);
            maze[n].neighb[2] = nmaze(i-1, NYMAZE-1);
            maze[n].directions[0] = 2;
            maze[n].directions[1] = 3;
            maze[n].directions[2] = 4;
        }
        else
        {
            maze[n].nneighb = 5;
            maze[n].neighb[0] = nmaze(i+1, NYMAZE-1);
            maze[n].neighb[1] = nmaze(i+1, NYMAZE-2);
            maze[n].neighb[2] = nmaze(i, NYMAZE-2);
            maze[n].neighb[3] = nmaze(i-1, NYMAZE-2);
            maze[n].neighb[4] = nmaze(i-1, NYMAZE-1);
            maze[n].directions[0] = 1;
            maze[n].directions[1] = 2;
            maze[n].directions[2] = 3;
            maze[n].directions[3] = 4;
            maze[n].directions[4] = 5;
        }
    }
    /* corners */
    n = nmaze(0,0);
    maze[n].nneighb = 3;
    maze[n].neighb[0] = nmaze(0,1);
    maze[n].neighb[1] = nmaze(1,1);
    maze[n].neighb[2] = nmaze(1,0);
    maze[n].directions[0] = 0;
    maze[n].directions[1] = 1;
    maze[n].directions[2] = 2;
    
    n = nmaze(NXMAZE-1,0);
    maze[n].nneighb = 2;
    maze[n].neighb[0] = nmaze(NXMAZE-1,1);
    maze[n].neighb[1] = nmaze(NXMAZE-2,0);
    maze[n].directions[0] = 0;
    maze[n].directions[1] = 5;

    n = nmaze(0,NYMAZE-1);
    maze[n].nneighb = 2;
    maze[n].neighb[0] = nmaze(1,NYMAZE-1);
    maze[n].neighb[1] = nmaze(0,NYMAZE-2);
    maze[n].directions[0] = 2;
    maze[n].directions[1] = 3;
    
    n = nmaze(NXMAZE-1,NYMAZE-1);
    maze[n].nneighb = 3;
    maze[n].neighb[0] = nmaze(NXMAZE-1,NYMAZE-2);
    maze[n].neighb[1] = nmaze(NXMAZE-2,NYMAZE-2);
    maze[n].neighb[2] = nmaze(NXMAZE-2,NYMAZE-1);
    maze[n].directions[0] = 3;
    maze[n].directions[1] = 4;
    maze[n].directions[2] = 5;

    /* initialize other parameters */
    for (i=0; i<NXMAZE; i++)
        for (j=0; j<NYMAZE; j++)
        {
            n = nmaze(i, j);
            maze[n].active = 0;
            maze[n].tested = 0;
            maze[n].connected = 0;
            maze[n].closed = 0;
            maze[n].north = 1;
            maze[n].northeast = 1;
            maze[n].southeast = 1;
            maze[n].south = 1;
            maze[n].southwest = 1;
            maze[n].northwest = 1;
        }
}

void init_oct_maze_graph(t_maze maze[NXMAZE*NYMAZE])
/* initialise graph of maze made of octagons and squares */
{
    int i, j, k, n, p, q;
    
    printf("Initializing maze\n");
    if (MAZE_MAX_NGBH < 8)
    {
        printf("Error: MAZE_MAX_NGBH should be at least 8 for circular maze\n");
        exit(0);
    }

    /* initialize neighbours */
    /* in the bulk */
    for (i=1; i<NXMAZE-1; i++)
        for (j=1; j<NYMAZE-1; j++)
        {
            n = nmaze(i, j);
            maze[n].nneighb = 4;
            maze[n].neighb[0] = nmaze(i, j+1);
            maze[n].neighb[1] = nmaze(i+1, j);
            maze[n].neighb[2] = nmaze(i, j-1);
            maze[n].neighb[3] = nmaze(i-1, j);
            for (k=0; k<4; k++) maze[n].directions[k] = k;
            
            if ((i+j)%2 == 0)
            {
                maze[n].nneighb = 8;
                maze[n].neighb[4] = nmaze(i+1, j+1);
                maze[n].neighb[5] = nmaze(i+1, j-1);
                maze[n].neighb[6] = nmaze(i-1, j-1);
                maze[n].neighb[7] = nmaze(i-1, j+1);
                for (k=4; k<8; k++) maze[n].directions[k] = k;
            }
        }
    
    /* left side */
    for (j=1; j<NYMAZE-1; j++)
    {
        n = nmaze(0, j);
        maze[n].nneighb = 3;
        maze[n].neighb[0] = nmaze(0, j+1);
        maze[n].neighb[1] = nmaze(1, j);
        maze[n].neighb[2] = nmaze(0, j-1);
        for (k=0; k<3; k++) maze[n].directions[k] = k;
        
        if (j%2 == 0)
        {
            maze[n].nneighb = 5;
            maze[n].neighb[3] = nmaze(1, j+1);
            maze[n].neighb[4] = nmaze(1, j-1);
            for (k=3; k<5; k++) maze[n].directions[k] = k+1;
        }
    }
    /* right side */
    for (j=1; j<NYMAZE-1; j++)
    {
        n = nmaze(NXMAZE-1, j);
        maze[n].nneighb = 3;
        maze[n].neighb[0] = nmaze(NXMAZE-1, j+1);
        maze[n].neighb[1] = nmaze(NXMAZE-2, j);
        maze[n].neighb[2] = nmaze(NXMAZE-1, j-1);
        maze[n].directions[0] = 0;
        maze[n].directions[1] = 3;
        maze[n].directions[2] = 2;
        
        if ((NXMAZE-1+j)%2 == 0)
        {
            maze[n].nneighb = 5;
            maze[n].neighb[3] = nmaze(NXMAZE-2, j-1);
            maze[n].neighb[4] = nmaze(NXMAZE-2, j+1);
            for (k=3; k<5; k++) maze[n].directions[k] = k+3;
        }
    }
    /* bottom side */
    for (i=1; i<NXMAZE-1; i++)
    {
        n = nmaze(i, 0);
        maze[n].nneighb = 3;
        maze[n].neighb[0] = nmaze(i, 1);
        maze[n].neighb[1] = nmaze(i+1, 0);
        maze[n].neighb[2] = nmaze(i-1, 0);
        maze[n].directions[0] = 0;
        maze[n].directions[1] = 1;
        maze[n].directions[2] = 3;
        
        if (i%2 == 0)
        {
            maze[n].nneighb = 5;
            maze[n].neighb[3] = nmaze(i+1, 1);
            maze[n].neighb[4] = nmaze(i-1, 1);
            maze[n].directions[3] = 4;
            maze[n].directions[4] = 7;
        }
    }
    /* top side */
    for (i=1; i<NXMAZE-1; i++)
    {
        n = nmaze(i, NYMAZE-1);
        maze[n].nneighb = 3;
        maze[n].neighb[0] = nmaze(i, NYMAZE-2);
        maze[n].neighb[1] = nmaze(i+1, NYMAZE-1);
        maze[n].neighb[2] = nmaze(i-1, NYMAZE-1);
        maze[n].directions[0] = 2;
        maze[n].directions[1] = 1;
        maze[n].directions[2] = 3;
        
        if ((i+NXMAZE-1)%2 == 0)
        {
            maze[n].nneighb = 5;
            maze[n].neighb[3] = nmaze(i+1, NYMAZE-2);
            maze[n].neighb[4] = nmaze(i-1, NYMAZE-2);
            maze[n].directions[3] = 5;
            maze[n].directions[4] = 6;
        }
    }
    /* corners */
    n = nmaze(0,0);
    maze[n].nneighb = 3;
    maze[n].neighb[0] = nmaze(1,0);
    maze[n].neighb[1] = nmaze(0,1);
    maze[n].neighb[2] = nmaze(1,1);
    maze[n].directions[0] = 1;
    maze[n].directions[1] = 0;
    maze[n].directions[2] = 4;
    
    n = nmaze(NXMAZE-1,0);
    maze[n].nneighb = 2;
    maze[n].neighb[0] = nmaze(NXMAZE-2,0);
    maze[n].neighb[1] = nmaze(NXMAZE-1,1);
    maze[n].directions[0] = 3;
    maze[n].directions[1] = 0;
    if ((NXMAZE-1)%2 == 0)
    {
        maze[n].nneighb = 3;
        maze[n].neighb[2] = nmaze(NXMAZE-2,1);
        maze[n].directions[2] = 7;
    }

    n = nmaze(0,NYMAZE-1);
    maze[n].nneighb = 2;
    maze[n].neighb[0] = nmaze(1,NYMAZE-1);
    maze[n].neighb[1] = nmaze(0,NYMAZE-2);
    maze[n].directions[0] = 1;
    maze[n].directions[1] = 2;
    if ((NYMAZE-1)%2 == 0)
    {
        maze[n].nneighb = 3;
        maze[n].neighb[2] = nmaze(1,NYMAZE-2);
        maze[n].directions[2] = 5;
    }
    
    n = nmaze(NXMAZE-1,NYMAZE-1);
    maze[n].nneighb = 2;
    maze[n].neighb[0] = nmaze(NXMAZE-2,NYMAZE-1);
    maze[n].neighb[1] = nmaze(NXMAZE-1,NYMAZE-2);
    maze[n].directions[0] = 3;
    maze[n].directions[1] = 2;
    if ((NXMAZE+NYMAZE)%2 == 0)
    {
        maze[n].nneighb = 3;
        maze[n].neighb[2] = nmaze(NXMAZE-2,NYMAZE-2);
        maze[n].directions[2] = 6;
    }

    /* initialize other parameters */
    for (i=0; i<NXMAZE; i++)
        for (j=0; j<NYMAZE; j++)
        {
            n = nmaze(i, j);
            maze[n].active = 0;
            maze[n].tested = 0;
            maze[n].connected = 0;
            maze[n].closed = 0;
            maze[n].north = 1;
            maze[n].east = 1;
            maze[n].south = 1;
            maze[n].west = 1;
            maze[n].northeast = 1;
            maze[n].southeast = 1;
            maze[n].southwest = 1;
            maze[n].northwest = 1;
        }
        /* for debugging */
//     for (i=0; i<NXMAZE; i++)
//         for (j=0; j<NYMAZE; j++)
//         {
//             n = nmaze(i, j);
//             q = maze[n].nneighb;
//             if (q > 0) printf("Cell (%i, %i)\n", i, j);
//             for (k = 0; k <q; k++)
//             {
//                 p = maze[n].neighb[k];
//                 printf("Neighbour %i at (%i, %i)\t", k, p%NXMAZE, p/NXMAZE);
//                 printf("Direction %i\n", maze[n].directions[k]);
//             }
//         }
//     sleep(5);

}


int find_maze_path(t_maze maze[NXMAZE*NYMAZE], int n0, int *path, int *pathlength, int mazetype)
/* find a random walk path in the maze */
/* returns 0 or 1 depending on whether path reaches a tested cell or a deadend */
{
    int active_counter = 0, i, n = n0, npaths, inext, nextcell, trial, nnext, deadend = 1, length = 0;
    int next_table[MAZE_MAX_NGBH];
    
    /* contruct random walk */
    npaths = maze[n].nneighb;
    path[0] = n0;
    
//     while ((npaths > 0)&&(!maze[n].tested))
    while ((npaths > 0))
    {
        maze[n].active = 1;
        
        printf("Cell (%i, %i) ", n%NXMAZE, n/NXMAZE);
        
        nnext = 0;
        for (i=0; i<npaths; i++)
        {
            nextcell = maze[n].neighb[i];
            if ((!maze[nextcell].active)&&((maze[nextcell].connected)||(!maze[nextcell].tested)))
            {
                next_table[nnext] = i;
                nnext++;
            }
        }
        
        if (nnext == 0) 
        {
            deadend = 1;
            printf("Ended path\n");
//             sleep(5);
            npaths = 0;
            maze[n].closed = 1;
        }
        else
        {
            deadend = 0;
            inext = next_table[rand()%nnext];
            nextcell = maze[n].neighb[inext];
            /* square and circular maze */
            if (mazetype < MAZE_TYPE_HEX) switch(maze[n].directions[inext]){
                case(0): 
                {
                    printf("Moving north\n");
                    maze[n].north = 0;
                    maze[nextcell].south = 0;
                    break;
                }
                case(1): 
                {
                    printf("Moving east\n");
                    maze[n].east = 0;
                    maze[nextcell].west = 0;
                    break;
                }
                case(2): 
                {
                    printf("Moving south\n");
                    maze[n].south = 0;
                    maze[nextcell].north = 0;
                    break;
                }
                case(3): 
                {
                    printf("Moving west\n");
                    maze[n].west = 0;
                    /* TODO find which is which */
                    maze[nextcell].east = 0;
                    maze[nextcell].northeast = 0;
                    break;
                }
                case(4): /* for circular maze */
                {
                    printf("Moving north-east\n");
                    maze[n].northeast = 0;
                    maze[nextcell].west = 0;
                    break;
                }
            }
            /* case of hexagonal maze */
            else if (mazetype == MAZE_TYPE_HEX) switch(maze[n].directions[inext]){
                case(0): 
                {
                    printf("Moving north\n");
                    maze[n].north = 0;
                    maze[nextcell].south = 0;
                    break;
                }
                case(1): 
                {
                    printf("Moving north-east\n");
                    maze[n].northeast = 0;
                    maze[nextcell].southwest = 0;
                    break;
                }
                case(2): 
                {
                    printf("Moving south-east\n");
                    maze[n].southeast = 0;
                    maze[nextcell].northwest = 0;
                    break;
                }
                case(3): 
                {
                    printf("Moving south\n");
                    maze[n].south = 0;
                    maze[nextcell].north = 0;
                    break;
                }
                case(4): 
                {
                    printf("Moving south-west\n");
                    maze[n].southwest = 0;
                    maze[nextcell].northeast = 0;
                    break;
                }
                case(5): 
                {
                    printf("Moving north-west\n");
                    maze[n].northwest = 0;
                    maze[nextcell].southeast = 0;
                    break;
                }
            }
            /* case of octagonal maze */
            else if (mazetype == MAZE_TYPE_OCT) switch(maze[n].directions[inext]){
                case(0): 
                {
                    printf("Moving north\n");
                    maze[n].north = 0;
                    maze[nextcell].south = 0;
                    break;
                }
                case(1): 
                {
                    printf("Moving east\n");
                    maze[n].east = 0;
                    maze[nextcell].west = 0;
                    break;
                }
                case(2): 
                {
                    printf("Moving south\n");
                    maze[n].south = 0;
                    maze[nextcell].north = 0;
                    break;
                }
                case(3): 
                {
                    printf("Moving west\n");
                    maze[n].west = 0;
                    maze[nextcell].east = 0;
                    break;
                }
                case(4): 
                {
                    printf("Moving north-east\n");
                    maze[n].northeast = 0;
                    maze[nextcell].southwest = 0;
                    break;
                }
                case(5): 
                {
                    printf("Moving south-east\n");
                    maze[n].southeast = 0;
                    maze[nextcell].northwest = 0;
                    break;
                }
                case(6): 
                {
                    printf("Moving south-west\n");
                    maze[n].southwest = 0;
                    maze[nextcell].northeast = 0;
                    break;
                }
                case(7): 
                {
                    printf("Moving north-west\n");
                    maze[n].northwest = 0;
                    maze[nextcell].southeast = 0;
                    break;
                }
            }
            n = nextcell;
            if (maze[n].tested) npaths = 0;
            else npaths = maze[n].nneighb;
            active_counter++;
            
            if (length < NXMAZE*NYMAZE) 
            {
                length++;
                path[length] = n;
            }
            deadend = 0;
        }
    }
    printf("Reached tested cell (%i, %i)\n", n%NXMAZE, n/NXMAZE);
    
    if (!maze[n].connected) deadend = 1;
    
    /* update cell status */
    for (n=0; n<NXMAZE*NYMAZE; n++) if (maze[n].active)
    {
        maze[n].active = 0;
        maze[n].tested = 1;
    }
    
    printf("Ended path\n");
    if (deadend) printf("Deadend\n"); 
    *pathlength = length;
    printf("Path length %i \n", length);
    
    return(deadend);
//     return(active_counter);
}

void init_maze_old(t_maze maze[NXMAZE*NYMAZE])
/* init a maze */
{
    int i, pathlength, *path;
    
    init_maze_graph(maze);
    
    for (i=0; i<RAND_SHIFT; i++) rand();
    
    for (i=0; i<NXMAZE*NYMAZE; i++) if (!maze[i].tested) find_maze_path(maze, i, path, &pathlength, 0);
    
}

void init_maze_oftype(t_maze maze[NXMAZE*NYMAZE], int type)
/* init a maze of given type */
{
    int i, j, n, deadend, pathlength, newpathlength;
    int *path, *newpath;
    
    path = (int *)malloc(2*NXMAZE*NYMAZE*sizeof(short int));
    newpath = (int *)malloc(2*NXMAZE*NYMAZE*sizeof(short int));
    
    switch (type) {
        case (MAZE_TYPE_SQUARE):
        {
            init_maze_graph(maze);
            break;
        }
        case (MAZE_TYPE_CIRCLE):
        {
            init_circular_maze_graph(maze);
            break;
        }
        case (MAZE_TYPE_HEX):
        {
            init_hex_maze_graph(maze);
            break;
        }
        case (MAZE_TYPE_OCT):
        {
            init_oct_maze_graph(maze);
            break;
        }
    }
    
    for (i=0; i<RAND_SHIFT; i++) rand();
    
    find_maze_path(maze, 0, path, &pathlength, type);
    
    for (n=0; n<pathlength; n++) maze[path[n]].connected = 1;
    
    for (i=0; i<NXMAZE*NYMAZE; i++) if ((!maze[i].tested)&&(!maze[i].connected))
    {
        deadend = find_maze_path(maze, i, path, &pathlength, type);
        if (!deadend) for (n=0; n<pathlength; n++) maze[path[n]].connected = 1;
        j = 0;
        printf("deadend = %i, pathlength = %i\n", deadend, pathlength);
        
//         while ((deadend)&&(j < pathlength))
        while (deadend)
        {
            j++;
            if (j > pathlength) j = 0; 
            printf("j = %i\n", j);
//             while (deadend) 
            if (!maze[path[j]].connected) deadend = find_maze_path(maze, path[j], newpath, &newpathlength, type);
            if (!deadend) for (n=0; n<newpathlength; n++) maze[newpath[n]].connected = 1;
        }
        
//         for (n=0; n<newpathlength; n++) maze[newpath[n]].connected = 1;
        for (j=0; j<NXMAZE*NYMAZE; j++) if (maze[j].tested) maze[j].connected = 1;
    }
    
    free(path);
    free(newpath);
}

void init_maze(t_maze maze[NXMAZE*NYMAZE])
/* init a maze */
{
    init_maze_oftype(maze, MAZE_TYPE_SQUARE);
}

void init_circular_maze(t_maze maze[NXMAZE*NYMAZE])
/* init a circular maze */
{
//     int i, j, n, q; 
    
    init_maze_oftype(maze, MAZE_TYPE_CIRCLE);
    
        /* for debugging */
//     for (i=0; i<NXMAZE; i++)
//         for (j=0; j<NYMAZE; j++)
//         {
//             n = nmaze(i, j);
//             q = maze[n].nneighb;
//             if (q > 0) 
//             {
//                 printf("Cell (%i, %i)\t", i, j);
//                 if (maze[n].north) printf("North ");
//                 if (maze[n].northeast) printf("N-E ");
//                 if (maze[n].east) printf("East ");
//                 if (maze[n].south) printf("South ");
//                 if (maze[n].west) printf("West ");
//                 printf("\n");
//             }
//         }
//     sleep(5);
}

void init_hex_maze(t_maze maze[NXMAZE*NYMAZE])
/* init a maze with hexagonal cells */
{
    init_maze_oftype(maze, MAZE_TYPE_HEX);
}

void init_oct_maze(t_maze maze[NXMAZE*NYMAZE])
/* init a maze with hexagonal cells */
{
    init_maze_oftype(maze, MAZE_TYPE_OCT);
}

void init_maze_exit(int nx, int ny, t_maze maze[NXMAZE*NYMAZE])
/* init a maze with exit at (nx, ny) */
{
    int i, j, n, deadend, pathlength, newpathlength;
    int *path, *newpath;
    
    path = (int *)malloc(2*NXMAZE*NYMAZE*sizeof(short int));
    newpath = (int *)malloc(2*NXMAZE*NYMAZE*sizeof(short int));
    
    init_maze_graph(maze);
    
    for (i=0; i<RAND_SHIFT; i++) rand();
    
    find_maze_path(maze, nmaze(nx, ny), path, &pathlength, 0);
    for (n=0; n<pathlength; n++) maze[path[n]].connected = 1;
    
    for (i=0; i<NXMAZE*NYMAZE; i++) if ((!maze[i].tested)&&(!maze[i].connected))
    {
        deadend = find_maze_path(maze, i, path, &pathlength, 0);
        if (!deadend) for (n=0; n<pathlength; n++) maze[path[n]].connected = 1;
        j = 0;
        printf("deadend = %i, pathlength = %i\n", deadend, pathlength);
        
//         while ((deadend)&&(j < pathlength))
        while (deadend)
        {
            j++;
            if (j > pathlength) j = 0; 
            printf("j = %i\n", j);
//             while (deadend) 
            if (!maze[path[j]].connected) deadend = find_maze_path(maze, path[j], newpath, &newpathlength, 0);
            if (!deadend) for (n=0; n<newpathlength; n++) maze[newpath[n]].connected = 1;
        }
        
//         for (n=0; n<newpathlength; n++) maze[newpath[n]].connected = 1;
        for (j=0; j<NXMAZE*NYMAZE; j++) if (maze[j].tested) maze[j].connected = 1;
    }
    
    free(path);
    free(newpath);
}

