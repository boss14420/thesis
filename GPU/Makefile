all: hydrodynamic-flow-cuda

CXX = g++
#CXX = g++-4.7

#CXXFLAGS_BASIC = -std=c++11 -pipe
#CXXFLAGS_DEBUG = $(CXXFLAGS_BASIC) -g
#CXXFLAGS_RELEASE = $(CXXFLAGS_BASIC) -O3 -march=native
#CXXFLAGS = $(CXXFLAGS_RELEASE)
CXXFLAGS ?=

LIBS ?=
INCLUDES ?=
NVCC ?= nvcc
INCLUDES ?=

ALL_CXXFLAGS =
#CCFLAGS = -ccbin $(GCC) -m64
CCFLAGS_DEBUG = -ccbin $(CXX) -m64 --device-debug -g
CCFLAGS_RELEASE = -ccbin $(CXX) -m64 -O3 --use_fast_math
ALL_CXXFLAGS += $(CCFLAGS_RELEASE)
ifneq (,$(CXXFLAGS))
    ALL_CXXFLAGS += -Xcompiler $(CXXFLAGS)
endif

#GENCODE_SM20    := -gencode arch=compute_20,code=sm_20
# GENCODE_SM30    := -gencode arch=compute_30,code=sm_30 -gencode arch=compute_35,code=\"sm_35,compute_35\" -dc
#GENCODE_SM30    := -gencode arch=compute_30,code=sm_30
GENCODE_SM30    := -gencode arch=compute_35,code=\"sm_35,compute_35\" -dc
GENCODE_FLAGS   := $(GENCODE_SM20) $(GENCODE_SM30)

#LDFLAGS = -locelot
LDFLAGS ?= -arch=sm_35 -rdc=true -lcudadevrt

hydrodynamic-flow-cuda: hydrodynamic-flow-cuda.o
	$(NVCC) $(ALL_CXXFLAGS) $(LDFLAGS) $(LIBS) $(INCLUDES) -o $@ $<

hydrodynamic-flow-cuda.o: hydrodynamic-flow-cuda.cu
	$(NVCC) $(ALL_CXXFLAGS) $(INCLUDES) $(GENCODE_FLAGS) -o $@ -c $<

.PHONY: clean

clean:
	rm -rf *.o
