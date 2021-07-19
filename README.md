### Tool to create videos of particles of waves in different 2D domains.

Created by **Nils Berglund** and optimized by **Marco Mancini**

C code for videos on YouTube Channel https://www.youtube.com/c/NilsBerglund

I plan to add a list of parameter values later on. 

There are two groups of 4 files:

### Simulations of classical particles in billiards.

global_particles.c

sub_part_billiard.c

particle_billiard.c

drop_billiard.c

Create subfolders 'tif_part', 'tif_drop'

Customize constants at beginning of .c file

Compile with 

gcc -o particle_billiard particle_billiard.c-O3 -L/usr/X11R6/lib -ltiff -lm -lGL -lGLU -lX11 -lXmu -lglut

gcc -o drop_billiard drop_billiard.c-O3 -L/usr/X11R6/lib -ltiff -lm -lGL -lGLU -lX11 -lXmu -lglut

Generate movie with 

ffmpeg -i part.%05d.tif -vcodec libx264 part.mp4

### Simulations of wave equation, heat equation and Schrodinger equation.

global_pdes.c

sub_wave.c

wave_billiard.c

heat.c

schrodinger.c

Create subfolders 'tif_wave', 'tif_heat', 'tif_schrod'

Customize constants at beginning of .c file

Compile with 

gcc -o wave_billiard wave_billiard.c -L/usr/X11R6/lib -ltiff -lm -lGL -lGLU -lX11 -lXmu -lglut -O3 -fopenmp

gcc -o heat heat.c -L/usr/X11R6/lib -ltiff -lm -lGL -lGLU -lX11 -lXmu -lglut -O3 -fopenmp

gcc -o schrodinger schrodinger.c -L/usr/X11R6/lib -ltiff -lm -lGL -lGLU -lX11 -lXmu -lglut -O3 -fopenmp

Generate movie with 

ffmpeg -i wave.%05d.tif -vcodec libx264 part.mp4
