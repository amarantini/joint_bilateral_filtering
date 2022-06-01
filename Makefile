# Makefile

CXX=clang++

CXXFLAGS = -std=c++17 

BINARIES=main

all: ${BINARIES}

main:main.o pfm.o Vector.o 
	${CXX} $^ -o $@

clean:
	/bin/rm -f ${BINARIES} *.o