all: hydrodynamic-flow

CXX = g++
CXXFLAGS_BASIC = -std=c++11 -pipe
CXXFLAGS_DEBUG = $(CXXFLAGS_BASIC) -g
CXXFLAGS_RELEASE = $(CXXFLAGS_BASIC) -O3 -march=native
CXXFLAGS = $(CXXFLAGS_RELEASE)
LIBS =
LDFLAGS =
INCLUDES =

hydrodynamic-flow: hydrodynamic-flow.o
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $(LIBS) $(INCLUDES) -o $@ $<

hydrodynamic-flow.o: hydrodynamic-flow.cc
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c -o $@ $<

.PHONY: clean

clean:
	rm -rf *.o
