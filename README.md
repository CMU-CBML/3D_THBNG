# 3D Phase Field Neuron Growth Simulation (IGA/THB-Splines)

## Overview

This C++ code simulates 3D neuron growth and neurodevelopmental disorder (NDD) deterioration using an Isogeometric Analysis (IGA) based phase field model. It employs Truncated Hierarchical B-splines (THB-splines) for adaptive local mesh refinement and couples phase field evolution with tubulin and neurotrophin dynamics. The framework supports modeling NDDs by adjusting simulation parameters.

## Key Features

* 3D Coupled Phase Field Model (Phase Field, Tubulin, Neurotrophin)
* Isogeometric Analysis (IGA) using THB-Splines
* Adaptive Multi-Level Local Refinement (Phase-field driven)
* Dynamic Domain Expansion during simulation
* 3D Neuron Identification and Tip Detection algorithms
* NDD Deterioration Modeling capability
* Parallel Execution using MPI and PETSc (SNES/KSP solvers)
* Checkpoint/Restart Functionality

## Example Simulations

*(Note: Store GIF files in a `media/` folder within the repository for these links to work).*

### Single Neuron Growth

Examples showcasing the growth and branching of a single neuron over time.

| Case 1                       | Case 2                       |
| :--------------------------: | :--------------------------: |
| ![Single Neuron 1](media/ref0.gif) | ![Single Neuron 2](media/ref1.gif) |
| **Case 3** | **Case 4** |
| ![Single Neuron 3](media/ref2.gif) | ![Single Neuron 4](media/ref9.gif) |

### Multi-Neuron Interaction

Examples showing the interaction and connection formation between two initially separate neurons.

| Case 1                       | Case 2                       |
| :--------------------------: | :--------------------------: |
| ![Multi Neuron 1](media/2n22.gif) | ![Multi Neuron 2](media/2n23.gif) |
| **Case 3** | **Case 4** |
| ![Multi Neuron 3](media/2n24.gif) | ![Multi Neuron 4](media/2n25.gif) |

*(Tip: Optimize your GIFs for size to ensure the README loads reasonably quickly).*

## Dependencies

* **C++ Compiler:** C++11 or later (GCC, Clang, Intel).
* **MPI:** MPI library (MPICH, OpenMPI).
* **PETSc:** Scientific computation library (\url{https://petsc.org/}), built with MPI support. Requires `PETSC_DIR` and `PETSC_ARCH` environment variables.
* **Nanoflann:** Header-only KD-tree library (included via relative path `../nanoflann/`). (\url{https://github.com/jlblancoc/nanoflann})
* **METIS (Recommended):** For graph partitioning (`mpmetis` command needed). (\url{http://glaros.dtc.umn.edu/gkhome/metis/metis/overview})
* **External Spline Tools (Required):** Pre-compiled executables needed in relative paths or system PATH:
    * `../spline3D_src/spline`: For Bezier extraction / spline processing.
    * `../THS3D/THS3D`: For THB-spline local refinement.
    *(Build these external tools separately).*

## Building the Code

*(Requires a PETSc-compatible build environment. Makefile not included).*

1.  Ensure all dependencies are installed and accessible.
2.  Set `PETSC_DIR` and `PETSC_ARCH` environment variables.
3.  Compile using `mpicxx` and link against PETSc libraries.

    **Example Makefile Snippet:**
    ```makefile
    include ${PETSC_DIR}/lib/petsc/conf/variables
    include ${PETSC_DIR}/lib/petsc/conf/rules

    TARGET = 3DNG
    SRCS = main.cpp NeuronGrowth.cpp utils.cpp BasicDataStructure.cpp
    OBJS = $(SRCS:.cpp=.o)
    INCLUDES = -I${PETSC_DIR}/${PETSC_ARCH}/include -I../nanoflann/1.5.5/include # Adjust Nanoflann path

    $(TARGET): $(OBJS)
        $(CXX) $(CFLAGS) -o $(TARGET) $(OBJS) ${PETSC_LIB}

    %.o: %.cpp
        $(CXX) $(CFLAGS) $(INCLUDES) -c $< -o $@

    clean:
        rm -f $(OBJS) $(TARGET)
    ```
    Run `make`.

## Input Files

Requires an input directory (`--path_in`) containing:

1.  **`simulation_parameters.txt`**: Key-value pairs defining model parameters (e.g., `dt`, `M_phi`). See example file.
2.  **Conditional Files:**
    * `phi.txt`: Needed if local refinement is active (specifies refinement flags).
    * `domain_size.txt`: Needed for restarting (contains grid dimensions/origin).
    * External tools might generate intermediate files here (`bzmeshinfo.txt`, `cmat.txt`, etc.).

## Running the Simulation

Execute via `mpirun`:

```bash
mpirun -np <N> ./3DNG \
       --numNeuron=<int> \
       --end_iter=<int> \
       --solver=<snes|ksp> \
       --path_in=<dir> \
       [--restart=<yes|no>]