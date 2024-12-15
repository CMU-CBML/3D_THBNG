# I am a comment, and I want to say that the variable CC will be
# the compiler to use.
# Hey!, I am comment number 2. I want to say that CFLAGS will be the
# options I'll pass to the compiler.

###################################################################################################	
UNAME_S := $(shell uname -s)
HOSTNAME := $(shell hostname)
ifeq ($(UNAME_S),Linux)
	ifeq ($(HOSTNAME),Kuanren-G14) # WSL
		Compiler = /home/kuanrenqian/Research/petsc/arch-linux-c-opt/bin/mpic++
	else
		Compiler = /ocean/projects/eng170006p/ussqww/petsc/arch-linux-c-opt/bin/mpic++
	endif
endif
ifeq ($(UNAME_S),Darwin)
	Compiler=/opt/homebrew/Cellar/gcc/14.2.0_1/bin/aarch64-apple-darwin24-g++-14
endif

EIGEN_DIR = ../eigen-3.4.0

CFLAGS= -Ofast -ffast-math -flto -c -std=c++17 -I ${EIGEN_DIR} -Wall 
#-DMKL_ILP64
FLAGS_LINK = -fopenmp -lpthread -lm -ldl


all: spline

spline: BasicDataStructure.o kernel.o KnotInsertion.o  main.o 
	$(Compiler) BasicDataStructure.o kernel.o KnotInsertion.o  main.o -o spline $(FLAGS_LINK)

BasicDataStructure.o: BasicDataStructure.cpp
	$(Compiler)	$(CFLAGS) BasicDataStructure.cpp $(FLAGS_LINK)
	
kernel.o: kernel.cpp
	$(Compiler)	$(CFLAGS) kernel.cpp $(FLAGS_LINK)

KnotInsertion.o: KnotInsertion.cpp
	$(Compiler)	$(CFLAGS) KnotInsertion.cpp $(FLAGS_LINK)

main.o: main.cpp
	$(Compiler)	$(CFLAGS) main.cpp $(FLAGS_LINK)

clean:
	rm *o *Zone.Identifier

