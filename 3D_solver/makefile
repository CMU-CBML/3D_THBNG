###################################################################################################	
UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Linux)
	PETSC_DIR = /home/kuanrenqian/Research/petsc
	EIGEN_DIR = /home/kuanrenqian/Research/eigen-3.4.0
	Compiler=/home/kuanrenqian/Research/petsc/arch-linux-c-opt/bin/mpicxx
endif
ifeq ($(UNAME_S),Darwin)
	PETSC_DIR = /Users/kuanrenqian/petsc
	EIGEN_DIR = /Users/kuanrenqian/eigen-3.4.0
	Compiler=/opt/homebrew/Cellar/gcc/13.2.0/bin/aarch64-apple-darwin23-g++-13
endif

include ${PETSC_DIR}/lib/petsc/conf/variables
CFLAGS= ${PETSC_CC_INCLUDES} ${CXX_FLAGS} ${CXXFLAGS} ${CPPFLAGS} ${PSOURCECXX} -O2 -march=native -march=native -std=c++17 -fopenmp -I ${EIGEN_DIR}
# -ffast-math -O3 

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

