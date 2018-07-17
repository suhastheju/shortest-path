##
 # This file is part of an experimental software implementation for solving the
 # single-source shortest path problem in edge-weighted graphs. The algorithm
 # considered for the implementation is Dijkstra's algorithm for the shortest
 # path problem and it runs in linear time with respect to the size of the host
 # graph.
 #
 # The source code is subject to the following license.
 #
 # The MIT License (MIT)
 #
 # Copyright (c) 2018 Suhas Thejaswi
 # Copyright (c) 2014 A. Björklund, P. Kaski, Ł. Kowalik, J. Lauri
 #
 # Permission is hereby granted, free of charge, to any person obtaining a copy
 # of this software and associated documentation files (the "Software"), to deal
 # in the Software without restriction, including without limitation the rights
 # to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 # copies of the Software, and to permit persons to whom the Software is
 # furnished to do so, subject to the following conditions:
 #
 # The above copyright notice and this permission notice shall be included in all
 # copies or substantial portions of the Software.
 #
 # THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 # IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 # FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 # AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 # LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 # OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 # SOFTWARE.
 #
 ##

MAKE = make
CC = gcc
CFLAGS = -O3 -Wall -march=native -std=c99 -fopenmp

all: graph-gen

EXE=graph-gen

graph-gen: graph-gen.c ffprng.h
	$(CC) $(CFLAGS) -o $(EXE) graph-gen.c -lm

.PHONY: $(EXE)
clean:
	rm -f *~ $(EXE)
