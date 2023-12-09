#scl enable devtoolset-9 bash

CC=g++ -std=c++17 -Ofast -lrt -DNDEBUG  -DHAVE_CXX0X -march=native -fpic -w -ftree-vectorize -ftree-vectorizer-verbose=0

CCOMP=g++ -std=c++17 -mavx512f -Ofast -lrt -DNDEBUG  -DHAVE_CXX0X -openmp -march=native -fpic -w -fopenmp -ftree-vectorize -ftree-vectorizer-verbose=0
SRCS=$(wildcard ./src/*.cpp)

TARGET=sos
RESDIR=./results
OUTPUT=./output
.PHONY: $(TARGET)

sos:$(SRCS)
	rm -rf $(TARGET)
	@test -d $(RESDIR) | mkdir -p $(RESDIR)
	@test -d $(OUTPUT) | mkdir -p $(OUTPUT)
	$(CCOMP) $(SRCS) -o $(TARGET)


clean:
	rm -rf $(TARGET)
	rm -rf $(wildcard *.txt)
	rm -rf ./datasets
	rm -rf ./output
	rm -rf ./results

cleancodes:
	rm -rf $(wildcard *.cpp) $(wildcard *.h) $(wildcard ./src/*.cpp) $(wildcard ./src/*.h)
