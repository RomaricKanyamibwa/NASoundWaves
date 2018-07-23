CXX=gfortran -O3 -fimplicit-none -fopenmp -Werror -Wall -O3  
#-ffast-math -ftree-vectorize -march=native -funroll-loops

FLAGS=-Iinc -g

LIB=liblapack.a librefblas.a liblapacke.a 
#LDFLAGS=-lboost_program_options -lm
#$(CXX) -o $@ $^ $(LDFLAGS)

FILES=finite_diff
OBJECTS=$(addsuffix .o, $(FILES))
BIN=bin/SoundWaves

all : $(BIN)

$(BIN) : $(addprefix obj/, $(OBJECTS))
	$(CXX) -o $@ $^ $(LIB)

obj/%.o : src/%.f90
	$(CXX) -c -o $@ $^ $(FLAGS)


clean :
	rm -f bin/* obj/*
