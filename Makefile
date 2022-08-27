CC = g++
flags = -O3 -fopenmp

BUILDDIR = build

objects = $(addprefix $(BUILDDIR)\,testing.o transform.o substrate.o read_myocytes.o process_signal.o polygon.o walkers.o)

run_sim: $(objects) libmatio-11.dll
	$(CC) $(flags) $(objects) run_sim.cpp libmatio-11.dll -o run_sim.exe

# Montecarlo
$(addprefix $(BUILDDIR)\,walkers.o): montecarlo\walkers.h montecarlo\walkers.cpp
	$(CC) $(flags) -c montecarlo\walkers.cpp -o $(addprefix $(BUILDDIR)\,walkers.o)

# Geometry
$(addprefix $(BUILDDIR)\,polygon.o): geometry\polygon.h geometry\polygon.cpp
	$(CC) $(flags) -c geometry\polygon.cpp -o $(addprefix $(BUILDDIR)\,polygon.o)

# MRI
$(addprefix $(BUILDDIR)\,process_signal.o): MRI\process_signal.h MRI\process_signal.cpp
	$(CC) $(flags) -c MRI\process_signal.cpp -o $(addprefix $(BUILDDIR)\,process_signal.o)

# Substrate
$(addprefix $(BUILDDIR)\,read_myocytes.o): substrate\read_myocytes.h substrate\read_myocytes.cpp
	$(CC) $(flags) -c substrate\read_myocytes.cpp -o $(addprefix $(BUILDDIR)\,read_myocytes.o)

$(addprefix $(BUILDDIR)\,substrate.o): substrate\substrate.h substrate\substrate.cpp
	$(CC) $(flags) -c substrate\substrate.cpp -o $(addprefix $(BUILDDIR)\,substrate.o)

$(addprefix $(BUILDDIR)\,transform.o): substrate\transform.h substrate\transform.cpp
	$(CC) $(flags) -c substrate\transform.cpp -o $(addprefix $(BUILDDIR)\,transform.o)

# Testing
$(addprefix $(BUILDDIR)\,testing.o): testing\testing.h testing\testing.cpp
	$(CC) $(flags) -c testing\testing.cpp -o $(addprefix $(BUILDDIR)\,testing.o)


.PHONY : clean
clean:
	del $(objects) run_sim.exe