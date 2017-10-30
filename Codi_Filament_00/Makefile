.SUFFIXES:
.SUFFIXES: .c .o
.c.o:
	gcc -c ona_filament.c

filament: ona_filament.o
	gcc -o filament ona_filament.o -lm

clean:
	rm -f *.o
