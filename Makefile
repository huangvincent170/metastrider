CXX = g++
CFLAGS  = -g -Wall

all: spgemm

spgemm: spgemm.cpp
	$(CXX) $(CFLAGS) -o spgemm spgemm.cpp matrix.hpp metastrider.hpp

clean:
	rm spgemm