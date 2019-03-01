CFLAGS= -O2 -Ofast -Wall -std=gnu11

default: simpar

all: simpar


simpar: simpar.c physics.c init_program.c
	gcc physics.c init_program.c simpar.c -o simpar $(CFLAGS)

clean:
	rm -f *.o simpar 
