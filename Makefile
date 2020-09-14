CXX := g++
CPPFLAGS := -march=native -fopenmp -O3 --std=c++0x -o scalescan

all: 
	make scalescan; make csr

scalescan: main.cpp graph.cpp
	$(CXX) $(CPPFLAGS) main.cpp graph.cpp

csr:
	$(CXX) -O3 -std=c++11 csr_uint.cpp -o csr_gen

clean:
	$(RM) scalescan csr_gen
