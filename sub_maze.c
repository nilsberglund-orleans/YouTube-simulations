/* Warning: the function init_maze does not always return a maze with a solution  */
/* The current algorithm uses a self-avoiding random walk. A better option may be */
/* to give random weights to the dual graph, and finite a maximal spanning tree   */

/* Change constant RAND_SHIFT to change the maze */

typedef struct
{
    short int nneighb;                      /* number of neighbours */
    int neighb[MAZE_MAX_NGBH];              /* neighbour cells */
    short int directions[MAZE_MAX_NGBH];    /* direction of neighbours */
    short int north, east, south, west;     /* closed walls */
    short int active;                       /* takes value 1 if currently active in RW path */
    short int tested;                       /* takes value 1 if tested */
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
            maze[n].north = 1;
            maze[n].east = 1;
            maze[n].south = 1;
            maze[n].west = 1;
        }
}

int find_maze_path(t_maze maze[NXMAZE*NYMAZE], int n0)
/* find a random walk path in the maze */
{
    int active_counter = 0, i, n = n0, npaths, inext, nextcell, trial, nnext;
    int next_table[4];
    
    /* contruct random walk */
    npaths = maze[n].nneighb;
//     while ((npaths > 0)&&(!maze[n].tested))
    while ((npaths > 0))
    {
        maze[n].active = 1;
        
        printf("Cell (%i, %i) ", n%NXMAZE, n/NXMAZE);
        
        nnext = 0;
        for (i=0; i<npaths; i++)
        {
            nextcell = maze[n].neighb[i];
            if (!maze[nextcell].active)
            {
                next_table[nnext] = i;
                nnext++;
            }
        }
        
        if (nnext == 0) 
        {
            printf("Ended path\n");
//             sleep(5);
            npaths = 0;
        }
        else
        {
            inext = next_table[rand()%nnext];
            nextcell = maze[n].neighb[inext];
            switch(maze[n].directions[inext]){
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
            }
        
            n = nextcell;
            if (maze[n].tested) npaths = 0;
            else npaths = maze[n].nneighb;
            active_counter++;
        }
    }
    
    /* update cell status */
    for (n=0; n<NXMAZE*NYMAZE; n++) if (maze[n].active)
    {
        maze[n].active = 0;
        maze[n].tested = 1;
    }
    
    printf("Ended path\n");
    
    return(active_counter);
}

void init_maze(t_maze maze[NXMAZE*NYMAZE])
/* init a maze */
{
    int i;
    
    init_maze_graph(maze);
    
    for (i=0; i<RAND_SHIFT; i++) rand();
    
    for (i=0; i<NXMAZE*NYMAZE; i++) if (!maze[i].tested) find_maze_path(maze, i);
    
}


