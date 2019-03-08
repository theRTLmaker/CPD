CC=gcc
CFLAGS= -O2 -I. -lm -Wall -std=gnu11
DEPS = debug.h init_program.h linkedList.h physics.h
OBJ =  debug.c init_program.c linkedList.c physics.c simpar.c

default: simpar

all: simpar

%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

simpar: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS)

clean:
	rm -f *.o simpar 