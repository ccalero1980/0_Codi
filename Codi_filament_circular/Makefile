.SUFFIXES:
.SUFFIXES: .c .o
.c.o:
	gcc -c ona_filament_circ_albeniz.c

filament_circ: ona_filament_circ_albeniz.o
	gcc -o filament_circ ona_filament_circ_albeniz.o -lm

clean:
	rm *.o
