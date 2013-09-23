# C compiler to run, change if necessary.
CC = /usr/bin/gcc

# compiler flags, assuming gcc compiler is used.
CFLAGS = -Wall -pedantic -O3

# If the C compiler does not find the FFTW header and library files, 
# replace the line above with something like this:
# CFLAGS = -Wall -pedantic -O3 -I. -I<FFTW directory>/api -L<FFTW directory>/.libs
#
# The names of the appropriate directories seem to have changed between different
# FFTW versions. 
#
# The -I flag indicates directories that are searched for header files to be included 
# in the compilation. It should be followed by the path to a directory containing 
# a file called fftw3.h
#
# If the wrong include directory is given, output from the compiler will contain
# the error message
# error: fftw3.h: No such file or directory
#
# The -L flag indicates directories where the linker looks for libraries to link.
# It should be followed by the path to a directory containing a file called libfftw3.a
# 
# If the -I flag is correct, but the -L flag is wrong, output from the compiler will
# contain the error message
# ld: library not found for -lfftw3

all: TEM-simulator
TEM-simulator: TEM-simulator.c array.o detector.o electronbeam.o functions.o geometry.o input.o log.o matrix.o misc.o mrc.o optics.o particle.o particleset.o pdb.o random.o sample.o simulation.o volume.o wavefunction.o
	$(CC) $(CFLAGS) -o $@ $^ -lfftw3 -lm
clean:
	rm *.o TEM-simulator
array.o: array.c array.h
	$(CC) -c $(CFLAGS) -o $@ $<
detector.o: detector.c detector.h array.h functions.h geometry.h input.h log.h macros.h misc.h optics.h random.h simulation.h structs.h wavefunction.h
	$(CC) -c $(CFLAGS) -o $@ $<
electronbeam.o: electronbeam.c electronbeam.h geometry.h input.h log.h macros.h matrix.h misc.h random.h simulation.h structs.h
	$(CC) -c $(CFLAGS) -o $@ $<
functions.o: functions.c functions.h array.h matrix.h misc.h
	$(CC) -c $(CFLAGS) -o $@ $<
geometry.o: geometry.c geometry.h input.h log.h macros.h matrix.h misc.h particle.h particleset.h random.h sample.h simulation.h structs.h
	$(CC) -c $(CFLAGS) -o $@ $<
input.o: input.c input.h log.h macros.h misc.h
	$(CC) -c $(CFLAGS) -o $@ $<
log.o: log.c log.h macros.h misc.h
	$(CC) -c $(CFLAGS) -o $@ $<
matrix.o: matrix.c matrix.h log.h macros.h misc.h
	$(CC) -c $(CFLAGS) -o $@ $<
misc.o: misc.c misc.h
	$(CC) -c $(CFLAGS) -o $@ $<
mrc.o: mrc.c mrc.h array.h macros.h log.h misc.h
	$(CC) -c $(CFLAGS) -o $@ $<
optics.o: optics.c optics.h geometry.h input.h log.h macros.h matrix.h misc.h random.h simulation.h 
	$(CC) -c $(CFLAGS) -o $@ $<
particle.o: particle.c particle.h array.h electronbeam.h functions.h geometry.h input.h log.h macros.h matrix.h misc.h mrc.h pdb.h random.h sample.h simulation.h structs.h wavefunction.h
	$(CC) -c $(CFLAGS) -o $@ $<
particleset.o: particleset.c particleset.h array.h input.h log.h macros.h matrix.h misc.h particle.h random.h sample.h simulation.h structs.h
	$(CC) -c $(CFLAGS) -o $@ $<
pdb.o: pdb.c pdb.h electronbeam.h input.h log.h macros.h misc.h sample.h
	$(CC) -c $(CFLAGS) -o $@ $<
random.o: random.c random.h macros.h
	$(CC) -c $(CFLAGS) -o $@ $<
sample.o: sample.c sample.h array.h electronbeam.h functions.h geometry.h input.h log.h macros.h matrix.h misc.h particle.h particleset.h pdb.h simulation.h structs.h wavefunction.h
	$(CC) -c $(CFLAGS) -o $@ $<
simulation.o: simulation.c simulation.h detector.h electronbeam.h geometry.h input.h log.h macros.h misc.h optics.h particle.h sample.h structs.h volume.h wavefunction.h
	$(CC) -c $(CFLAGS) -o $@ $<
volume.o: volume.c volume.h array.h electronbeam.h geometry.h input.h log.h macros.h matrix.h misc.h mrc.h particle.h particleset.h sample.h simulation.h
	$(CC) -c $(CFLAGS) -o $@ $<
wavefunction.o: wavefunction.c wavefunction.h array.h electronbeam.h functions.h geometry.h log.h misc.h optics.h particle.h particleset.h sample.h simulation.h structs.h
	$(CC) -c $(CFLAGS) -o $@ $<