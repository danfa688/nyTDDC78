CC=mpicc

CFLAGS = -Nmpi -g

LFLAGS= -Nmpi -g -lpthread -lrt -lm

all: main 

clean:
	-$(RM) *.o main extra physics

main: extra.o physics.o main.o
	$(CC) -o $@ extra.o physics.o main.o $(LFLAGS)

extra:
	$(CC) -o $@ extra.o $(LFLAGS)

physics: 
	$(CC) -o $@ physics.o $(LFLAGS)
