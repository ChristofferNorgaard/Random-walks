# Fortran Compiler
FC = gfortran

# Compiler flags
# -Wall: Enable all warnings
# -Wextra: Enable extra warnings (often useful for Fortran)
# -O3: Optimization level 3 (for performance)
# -g: Include debugging information (can be removed for release builds)
FFLAGS = -Wall -Wextra -O3 -g -fopenmp

# Directory where LAPACK and BLAS libraries are located
# IMPORTANT: Change this to your actual library directory
LIBDIR = ../../LAPACK

# Libraries to link
# -L$(LIBDIR): Specifies the directory to search for libraries
# -llapack: Links the LAPACK library (e.g., liblapack.a or liblapack.so)
# -lblas: Links the BLAS library (e.g., libblas.a or libblas.so)
# The order is often important: LAPACK depends on BLAS, so BLAS should come after LAPACK.
LDFLAGS = -L$(LIBDIR) -L/opt/homebrew/opt/libomp/lib -llapack -lblas -fopenmp

# Name of the executable
TARGET = gamma

# List of all Fortran source files (.f90, .f, .for)
SRCS = gamma-par.f90 

# Automatically generate object file names from source files
OBJS = $(SRCS:.f90=.o)

.PHONY: all clean

all: $(TARGET)

$(TARGET): $(OBJS)
	$(FC) $(FFLAGS) $^ -o $@ $(LDFLAGS)

# Rule for compiling Fortran source files into object files
# $<: The name of the prerequisite (the .f90 file)
# $@: The name of the target (the .o file)
%.o: %.f90
	$(FC) $(FFLAGS) -c $< -o $@

clean:
	rm -f $(OBJS) $(TARGET) *.mod
