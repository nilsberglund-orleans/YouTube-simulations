### Tool to create videos of particles or waves in different 2D domains.

Created by **Nils Berglund** and optimized by **Marco Mancini**

C code for videos on YouTube Channel https://www.youtube.com/c/NilsBerglund

Parameter values used in specific simulations will be gradually added to file `Parameters.md`, `Parameters_June21.md` and so on.

There are two groups of 4 files and 8 files:

### Simulations of classical particles in billiards.

1. *global_particles.c*:    global variables and parameters
2. *sub_part_billiard.c*:   drawing/computation routines common to `particle_billiard` and `drop_billiard`
3. *particle_billiard.c*:   simulation of a collection of non-interacting particles in a billiard
4. *drop_billiard.c*:       simulation of an expanding front of particles

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

### Simulations of wave equation, heat equation and Schrodinger equation.

1. *global_pdes.c*:     global variables and parameters
2. *sub_wave.c*:        drawing/computation routines common to `wave_billiard`, `heat` and `schrodinger`
3. *sub_wave_comp.c*:   some modified functions needed by `wave_comparison`
4. *common_wave.c*:     common functions of `wave_billiard` and `wave_comparison`
5. *wave_billiard.c*:   simulation of the (linear) wave equation
6. *wave_comparison.c*: comparison of the wave equation in two different domains
7. *heat.c*:            simulation of the heat equation, with optional drawing of gradient field lines
8. *schrodinger.c*:     simulation of the Schrodinger equation

- Create subfolders `tif_wave`, `tif_heat`, `tif_schrod`
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

#### Some references ####

- Discretizing the wave equation: https://hplgit.github.io/fdm-book/doc/pub/wave/pdf/wave-4print.pdf
- Absorbing boundary conditions: https://hal.archives-ouvertes.fr/hal-01374183
- Cloaking device: https://www.sciencedirect.com/science/article/pii/S0165212514001759

