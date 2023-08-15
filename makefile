PETSC_DIR = /Users/kuanrenqian/petsc
EIGEN_DIR = /Users/kuanrenqian/eigen-3.4.0

include ${PETSC_DIR}/lib/petsc/conf/variables
 
# Compiler=/Users/kuanrenqian/petsc/arch-darwin-c-debug/bin/mpicc
Compiler=/opt/homebrew/Cellar/gcc/13.1.0/bin/aarch64-apple-darwin22-g++-13
CFLAGS= ${PETSC_CC_INCLUDES} ${CXX_FLAGS} ${CXXFLAGS} ${CPPFLAGS}  ${PSOURCECXX}  -std=c++17 -Wall -fopenmp -I ${EIGEN_DIR}

all: 2DNG

2DNG: ng_func.o main.o
	$(Compiler) ng_func.o main.o -o 2DNG ${PETSC_LIB} $(CFLAGS)
	
main.o: main.cpp
	$(Compiler) -c main.cpp $(CFLAGS)

ng_func.o: ng_func.cpp
	$(Compiler) -c ng_func.cpp $(CFLAGS)	

clean:
	clear
	rm *o 2DNG
