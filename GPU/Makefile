all: hydrodynamic-flow-cuda

CXX = g++
GCC = g++-4.7
#CXXFLAGS_BASIC = -std=c++11 -pipe
#CXXFLAGS_DEBUG = $(CXXFLAGS_BASIC) -g
#CXXFLAGS_RELEASE = $(CXXFLAGS_BASIC) -O3 -march=native
#CXXFLAGS = $(CXXFLAGS_RELEASE)
LIBS =
LDFLAGS =
INCLUDES =
NVCC = nvcc
INCLUDES ?=

ALL_CXXFLAGS =
CCFLAGS = -ccbin $(GCC) -m64
ALL_CXXFLAGS += $(CCFLAGS)
ifneq ($(CXXFLAGS),)
    ALL_CXXFLAGS = -Xcompiler $(CXXFLAGS)
endif

GENCODE_SM20    := -gencode arch=compute_20,code=sm_20
GENCODE_SM30    := -gencode arch=compute_30,code=sm_30 -gencode arch=compute_35,code=\"sm_35,compute_35\"
GENCODE_FLAGS   := $(GENCODE_SM20) $(GENCODE_SM30)

LDFLAGS = -locelot

hydrodynamic-flow-cuda: hydrodynamic-flow-cuda.o
	$(NVCC) $(ALL_CXXFLAGS) $(LDFLAGS) $(LIBS) $(INCLUDES) -o $@ $<

hydrodynamic-flow-cuda.o: hydrodynamic-flow-cuda.cu
	$(NVCC) $(ALL_CXXFLAGS) $(INCLUDES) $(GENCODE_FLAGS) -o $@ -c $<

.PHONY: clean

clean:
	rm -rf *.o
