### Tool to create videos of particles or waves in different 2D domains.

Created by **Nils Berglund** and optimized by **Marco Mancini**

C code for videos on YouTube Channel https://www.youtube.com/c/NilsBerglund

I plan to add a list of parameter values later on. 

There are two groups of 4 files:

### Simulations of classical particles in billiards.

1. global_particles.c
2. sub_part_billiard.c
3. particle_billiard.c
4. drop_billiard.c

- Create subfolders 'tif_part', 'tif_drop'
- Customize constants at beginning of .c file
- Compile with 

gcc -o particle_billiard particle_billiard.c-O3 -L/usr/X11R6/lib -ltiff -lm -lGL -lGLU -lX11 -lXmu -lglut

gcc -o drop_billiard drop_billiard.c-O3 -L/usr/X11R6/lib -ltiff -lm -lGL -lGLU -lX11 -lXmu -lglut

- Generate movie with 

ffmpeg -i part.%05d.tif -vcodec libx264 part.mp4

### Simulations of wave equation, heat equation and Schrodinger equation.

1. global_pdes.c
2. sub_wave.c
3. wave_billiard.c
4. heat.c
5. schrodinger.c

- Create subfolders 'tif_wave', 'tif_heat', 'tif_schrod'
- Customize constants at beginning of .c file
- Compile with 

gcc -o wave_billiard wave_billiard.c -L/usr/X11R6/lib -ltiff -lm -lGL -lGLU -lX11 -lXmu -lglut -O3 -fopenmp

gcc -o heat heat.c -L/usr/X11R6/lib -ltiff -lm -lGL -lGLU -lX11 -lXmu -lglut -O3 -fopenmp

gcc -o schrodinger schrodinger.c -L/usr/X11R6/lib -ltiff -lm -lGL -lGLU -lX11 -lXmu -lglut -O3 -fopenmp

- Generate movie with 

ffmpeg -i wave.%05d.tif -vcodec libx264 wave.mp4
