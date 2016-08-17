CC=g++
CFLAGS=-w -Wall -c -O2
all: BROWNIES

BROWNIES: BROWNIES.o classes.o file_functions.o  icc.o  input_output.o  \
			ions_functions.o  PACO.o  physics_functions.o  retrace.o  \
			sim_domain.o  utils.o
		g++ -g  -w -Wall -O2 -Wl,--allow-multiple-definition  \
			BROWNIES.o classes.o file_functions.o  icc.o  input_output.o  \
			ions_functions.o  PACO.o  physics_functions.o  retrace.o  \
			sim_domain.o  utils.o -o BROWNIES \
			-lm -llapack -lblas -lgomp

BROWNIES.o: BROWNIES.cc
		$(CC) $(CFLAGS) BROWNIES.cc
classes.o: classes.cc classes.h
		$(CC) $(CFLAGS) classes.cc
file_functions.o: file_functions.cc file_functions.h
		$(CC) $(CFLAGS) file_functions.cc
icc.o: icc.cc icc.h
		$(CC) $(CFLAGS) icc.cc
input_output.o: input_output.cc input_output.h
		$(CC) $(CFLAGS) input_output.cc
ions_functions.o: ions_functions.cc ions_functions.h
		$(CC) $(CFLAGS) ions_functions.cc
PACO.o: PACO.cc PACO.h
		$(CC) $(CFLAGS) PACO.cc
physics_functions.o: physics_functions.cc physics_functions.h
		$(CC) $(CFLAGS) physics_functions.cc
retrace.o: retrace.cc retrace.h
		$(CC) $(CFLAGS) retrace.cc
sim_domain.o: sim_domain.cc sim_domain.h
		$(CC) $(CFLAGS) sim_domain.cc
utils.o: utils.cc utils.h
		$(CC) $(CFLAGS) utils.cc

clean:
		rm -rf *o BROWNIES
