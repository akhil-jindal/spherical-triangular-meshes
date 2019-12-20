CXX = g++
CXXFLAGS = -std=c++11 -Wall -g -rdynamic -fno-omit-frame-pointer
VPATH = src

SRC = manip.cpp
OBJ = mainp.o


main: main.o manip.o
	$(CXX) $(CXXFLAGS) $^ -o $@

main.o: main.cpp manip.o
	$(CXX) -c $(CXXFLAGS) $< -o $@

${OBJ}: $(SRC) mainp.hpp conn.hpp co.hpp
	$(CXX) -c $(CXXFLAGS) $< -o $@

test: main
	@echo
	@./$<

clean:
	rm -f *.o main

debug: main
	@echo
	valgrind -v ./$<
