### Tool to create videos of particles or waves in different 2D domains.

Created by **Nils Berglund** and optimized by **Marco Mancini**

C code for videos on YouTube Channel https://www.youtube.com/c/NilsBerglund

Parameter values used in specific simulations will be gradually added to file `Parameters.md`, `Parameters_June21.md` and so on.

There are four groups of 6 files, 19 files, 5 files and 4 files. 
In addition the following files handling color schemes have been included:

1. `hsluv.c`and `hsluv.h` from https://github.com/adammaj1/hsluv-color-gradient 
2. `turbo_colormap.c` from https://gist.github.com/mikhailov-work/6a308c20e494d9e0ccc29036b28faa7a
3. `colormaps.c` containing look-up tables from https://github.com/yuki-koyama/tinycolormap

The following file (beta version) provides support for creating mazes:

4. `sub_maze.c`

The file 

5. `Earth_Map_Blue_Marble_2002_large.ppm.gz` 

is required by `wave_sphere.c` and should be unzipped before compiling. 

### Simulations of classical particles in billiards.

1. *particle_billiard.c*:   simulation of a collection of non-interacting particles in a billiard
2. *drop_billiard.c*:       simulation of an expanding front of particles
3. *particle_pinball.c*:    variant of `particle_billiard` with some extra statistics plots 
4. *global_particles.c*:    global variables and parameters
5. *sub_part_billiard.c*:   drawing/computation routines common to `particle_billiard` and `drop_billiard`
6. *sub_part_pinball.c*:    additional drawing/computation routines for `particle_pinball`


- Create subfolders `tif_part`, `tif_drop`
- Customize constants at beginning of .c file
- Compile with 

`gcc -o particle_billiard particle_billiard.c-O3 -L/usr/X11R6/lib -ltiff -lm -lGL -lGLU -lX11 -lXmu -lglut`

`gcc -o drop_billiard drop_billiard.c-O3 -L/usr/X11R6/lib -ltiff -lm -lGL -lGLU -lX11 -lXmu -lglut`

- Many laptops claim to have 4 cores, but two of those are virtual. OMP acceleration may be more effective after executing           

`export OMP_NUM_THREADS=2` 

in the shell before running the program

- Generate movie with 

`ffmpeg -i part.%05d.tif -vcodec libx264 part.mp4`

### Simulations of wave equation and reaction-diffusion equations, including the Schrodinger equation.

1. *wave_billiard.c*:    simulation of the (linear) wave equation
2. *wave_3d.c*:          3d rendering of wave equation
3. *wave_sphere.c*:      wave equation on a sphere (3D and 2D render)
4. *wave_comparison.c*:  comparison of the wave equation in two different domains
5. *wave_energy.c*:      a version of `wave_billiard` plotting the energy profile of the wave
6. *mangrove.c*:         a version of `wave_billiard` with additional features to animate mangroves
7. *heat.c*:             simulation of the heat equation, with optional drawing of gradient field lines
8. *rde.c*:              simulation of reaction-diffusion equations, plots in 2D and 3D (including Schrödinger equation, 
                         Euler equation, and shallow water equation)
9. *schrodinger.c*:      simulation of the Schrodinger equation in 2D (old version)
10. *global_pdes.c*:      global variables and parameters
11. *global_3d.c*:        additional global variables for 3d version
12. *sub_wave.c*:         drawing/computation routines common to `wave_billiard`, `heat` and `schrodinger`
13. *sub_wave_comp.c*:    some modified functions needed by `wave_comparison`
14. *sub_wave_3d.c*:      additional functions for 3d version
15. *common_wave.c*:      common functions of `wave_billiard` and `wave_comparison`
16. *colors_waves.c*:     colormaps used by wave simulations
17. *sub_rde.c*:          additional routines for rde.c
18. *sub_wave_rde_3d.c*:  additional 3d drawing routines for rde.c
19. *sub_sphere.c*:       additional routines for wave_sphere.c

- Create subfolders `tif_wave`, `tif_heat`, `tif_bz`, `tif_schrod`
- Customize constants at beginning of .c file
- Compile with 

`gcc -o wave_billiard wave_billiard.c -L/usr/X11R6/lib -ltiff -lm -lGL -lGLU -lX11 -lXmu -lglut -O3 -fopenmp`

`gcc -o wave_comparison wave_comparison.c -L/usr/X11R6/lib -ltiff -lm -lGL -lGLU -lX11 -lXmu -lglut -O3 -fopenmp`

`gcc -o heat heat.c -L/usr/X11R6/lib -ltiff -lm -lGL -lGLU -lX11 -lXmu -lglut -O3 -fopenmp`

`gcc -o schrodinger schrodinger.c -L/usr/X11R6/lib -ltiff -lm -lGL -lGLU -lX11 -lXmu -lglut -O3 -fopenmp`

- Many laptops claim to have 4 cores, but two of those are virtual. OMP acceleration may be more effective after executing           

`export OMP_NUM_THREADS=2` 

in the shell before running the program

- Generate movie with 

`ffmpeg -i wave.%05d.tif -vcodec libx264 wave.mp4`

### Molecular dynamics simulations.

1. *lennardjones.c*:      simulation of molecular dynamics
2. *global_ljones.c*:     global variables and parameters
3. *sub_lj.c*:            drawing and initialization routines
4. *sub_hashgrid.c*:      hashgrid manipulation routines
5. *lj_movie.c*:          render movie with precomputed particle positions 
                          (requires files lj_time_series.dat and lj_final_positions.dat generated by lennardjones)

- Create subfolder `tif_ljones`
- Customize constants at beginning of .c file
- Compile with 

`gcc -o lennardjones lennardjones.c -L/usr/X11R6/lib -ltiff -lm -lGL -lGLU -lX11 -lXmu -lglut -O3 -fopenmp`

- Generate movie with 

`ffmpeg -i lennardjones.%05d.tif -vcodec libx264 lennardjones.mp4`

### Percolation simulations.

1. *percolation.c*:     simulation of Bernoulli percolation 
2. *global_perc.c*:     global variables and parameters
3. *sub_perco.c*:       drawing and cluster finding routines
4. *sub_perco_3d.c*:    3D drawing routines

- Create subfolder `tif_perc`
- Customize constants at beginning of .c file
- Compile with 

`gcc -o percolation percolation.c -L/usr/X11R6/lib -ltiff -lm -lGL -lGLU -lX11 -lXmu -lglut -O3 -fopenmp`

- Generate movie with 

`ffmpeg -i percolation.%05d.tif -vcodec libx264 percolation.mp4`

#### Some references ####

- Discretizing the wave equation: https://hplgit.github.io/fdm-book/doc/pub/wave/pdf/wave-4print.pdf
- Absorbing boundary conditions: https://hal.archives-ouvertes.fr/hal-01374183
- Cloaking device: https://www.sciencedirect.com/science/article/pii/S0165212514001759
- Poisson disc sampling: https://bl.ocks.org/mbostock/dbb02448b0f93e4c82c3
- Thermostat algorithm: https://doi.org/10.1007/s10955-009-9734-0
or http://www.maths.warwick.ac.uk/~theil/HL12-3-2009.pdf

