###################################################################################################	
UNAME_S := $(shell uname -s)
HOSTNAME := $(shell hostname)
ifeq ($(UNAME_S),Linux)
	ifeq ($(HOSTNAME),Kuanren-G14) # WSL
		PETSC_DIR = /home/kuanrenqian/Research/petsc
		Compiler = /home/kuanrenqian/Research/petsc/arch-linux-c-opt/bin/mpic++
	else
		PETSC_DIR = /ocean/projects/eng170006p/ussqww/petsc
		Compiler = /ocean/projects/eng170006p/ussqww/petsc/arch-linux-c-opt/bin/mpic++
	endif
endif
ifeq ($(UNAME_S),Darwin)
	PETSC_DIR = /Users/kuanrenqian/petsc
	Compiler=/opt/homebrew/Cellar/gcc/14.2.0_1/bin/aarch64-apple-darwin24-g++-14
endif

EIGEN_DIR = ../eigen-3.4.0
CFLAGS= ${PETSC_CC_INCLUDES} ${CXX_FLAGS} ${CXXFLAGS} ${CPPFLAGS}  ${PSOURCECXX} -Ofast -ffast-math -flto -std=c++11 -fopenmp -Wall -I ${EIGEN_DIR}

all: THS3D

THS3D: BasicDataStructure.o BSplineBasis.o kernel.o KnotInsertion.o Laplace.o LeastSquare.o main.o SingularEval.o TTSP_3D.o TTSP_2D.o
	$(Compiler) BasicDataStructure.o BSplineBasis.o kernel.o KnotInsertion.o Laplace.o LeastSquare.o main.o SingularEval.o TTSP_3D.o TTSP_2D.o -o THS3D ${PETSC_LIB} $(CFLAGS)
	
BasicDataStructure.o: BasicDataStructure.cpp
	$(Compiler)	-c BasicDataStructure.cpp $(CFLAGS)
	
BSplineBasis.o: BSplineBasis.cpp
	$(Compiler)	-c BSplineBasis.cpp $(CFLAGS)
	
Elasiticity_3D.o: Elasiticity_3D.cpp
	$(Compiler)	-c Elasiticity_3D.cpp $(CFLAGS)
	
kernel.o: kernel.cpp
	$(Compiler)	-c kernel.cpp $(CFLAGS)
	
KnotInsertion.o: KnotInsertion.cpp
	$(Compiler)	-c KnotInsertion.cpp $(CFLAGS)

Laplace.o: Laplace.cpp
	$(Compiler)	-c Laplace.cpp $(CFLAGS)

LeastSquare.o: LeastSquare.cpp
	$(Compiler)	-c LeastSquare.cpp $(CFLAGS)
	
main.o: main.cpp
	$(Compiler)	-c main.cpp $(CFLAGS)
	
SingularEval.o: SingularEval.cpp
	$(Compiler)	-c SingularEval.cpp $(CFLAGS)

TTSP_2D.o: TTSP_2D.cpp
	$(Compiler)	-c TTSP_2D.cpp $(CFLAGS)

TTSP_3D.o: TTSP_3D.cpp
	$(Compiler)	-c TTSP_3D.cpp $(CFLAGS)

clean:
	rm *o THS3D

