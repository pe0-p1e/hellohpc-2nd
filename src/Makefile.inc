CXX := $(if $(filter default,$(origin CXX)),mpicxx,$(CXX))
FC := $(if $(filter default,$(origin FC)),mpifort,$(FC))
F77 := $(if $(filter default,$(origin F77)),mpif77,$(F77))

CPPFLAGS ?= -Dfortran3 -Dnewc
CXXFLAGS ?= -O3 -march=native -mtune=native -fno-math-errno -fno-trapping-math -DNDEBUG
FCFLAGS ?= -O3 -march=native -cpp
LDLIBS ?= -lstdc++

Cu = nvcc
CUDA_LIB_PATH = -L/usr/lib/cuda/lib64 -I/usr/include -I/usr/lib/cuda/include
CUDA_APP_FLAGS = -c -g -O3 --ptxas-options=-v -Dfortran3 -Dnewc
