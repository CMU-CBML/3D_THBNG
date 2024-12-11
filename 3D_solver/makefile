###################################################################################################	
UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Linux)
	PETSC_DIR = /ocean/projects/eng170006p/ussqww/petsc
	Compiler = /ocean/projects/eng170006p/ussqww/petsc/arch-linux-c-opt/bin/mpic++
endif
ifeq ($(UNAME_S),Darwin)
	PETSC_DIR = /Users/kuanrenqian/petsc
	Compiler=/opt/homebrew/Cellar/gcc/14.2.0_1/bin/aarch64-apple-darwin24-g++-14
endif

EIGEN_DIR = ../eigen-3.4.0

include ${PETSC_DIR}/lib/petsc/conf/variables
CFLAGS= ${PETSC_CC_INCLUDES} ${CXX_FLAGS} ${CXXFLAGS} ${CPPFLAGS} ${PSOURCECXX} -Ofast -march=native -march=native -std=c++17 -fopenmp -flto -I ${EIGEN_DIR}

all: 3DNG

3DNG: BasicDataStructure.o NeuronGrowth.o utils.o main.o 
	$(Compiler) BasicDataStructure.o NeuronGrowth.o utils.o main.o -o 3DNG ${PETSC_LIB} $(CFLAGS)
	
BasicDataStructure.o: BasicDataStructure.cpp
	$(Compiler)	-c BasicDataStructure.cpp $(CFLAGS)
	
# UserSetting.o: UserSetting.cpp
# 	$(Compiler)	-c UserSetting.cpp $(CFLAGS)
	
NeuronGrowth.o: NeuronGrowth.cpp
	$(Compiler)	-c NeuronGrowth.cpp $(CFLAGS)

# ngfunc.o: ngfunc.cpp
# 	$(Compiler)	-c ngfunc.cpp $(CFLAGS)

# NeuronGrowth.o: NeuronGrowth.cpp
# 	$(Compiler)	-c NeuronGrowth.cpp $(CFLAGS)

utils.o: utils.cpp
	$(Compiler)	-c utils.cpp $(CFLAGS)

main.o: main.cpp
	$(Compiler)	-c main.cpp $(CFLAGS)

clean:
	rm *o 3DNG

