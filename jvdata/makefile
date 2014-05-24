CC=g++
CFLAGS=-O2 -ansi -pedantic -Wall
GLFLAGS=`./glflags.bash`

EXE=binboxing
all : $(EXE)

binboxing.o: binboxing.cc binboxing.hh
	$(CC) -c binboxing.cc $(CFLAGS)

clean:
	rm *.o fit.log $(EXE) >&-
