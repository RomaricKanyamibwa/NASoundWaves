CXX=gfortran -O3 -fimplicit-none -fopenmp -Werror -Wall -ffast-math -ftree-vectorize -march=native -funroll-loops

DEVOPT=-Wextra -fPIC -fmax-errors=5 -fcheck=all -fbacktrace
# -Wimplicit-interface
FLAGS=-Iinc -g

LIB=liblapack.a librefblas.a liblapacke.a 
#LDFLAGS=-lboost_program_options -lm
#$(CXX) -o $@ $^ $(LDFLAGS)

FILES=finite_diff finite_diff_diffusion inv_matrix plot_finite_diff finite_diff_acoustic finite_diff_continuity integrate_trapezoid
OBJECTS=$(addsuffix .o, $(FILES))
BIN=bin/SoundWaves

all : $(BIN)

$(BIN) : $(addprefix obj/, $(OBJECTS))
	$(CXX) -o $@ $^ $(LIB)

obj/%.o : src/%.f90
	$(CXX) -c -o $@ $^ $(FLAGS) $(DEVOPT)


clean :
	rm -f bin/* obj/*
