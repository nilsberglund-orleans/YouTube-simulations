LIBS = -L/usr/X11R6/lib -ltiff -lm -lGL -lGLU -lX11 -lglut -fopenmp

.c:
	gcc -o $@ $< -O3 $(LIBS)
