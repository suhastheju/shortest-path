# Shortest-path

## Overview

This software repository contains an implementation of the Dijkstra's algorithm
for solving the single-source shortest path problem in edge-weighted graphs. The
software is written in C programming language. The accompanying software also
contain a graph generator which is capable of generating regular, powerlaw and
clique graphs.

The source code is configured for a gcc build for Intel microarchitectures.
Other builds are possible but require manual configuration of the 'Makefile'.

The source code is subject to MIT license, see 'LICENSE' for details.

## Build

Use GNU make to build the software.

Our implementation makes use of preprocessor directives to enable conditional
compilation to enable resource tracking (memory bandwidth). Use
'TRACK_RESOURCES' compilation flag to enable resource tracking.

Check 'Makefile' for building the software.

### Usage
```
./shortest-path -in input-file -src source -dst destination

arguments:
    -in  : input file in DIMACS format
           reads from standard input if the argument is not specified
    -src : source vertex
    -dst : destination vertex
           -1, to report shortest path to all vertices

Example output:
./shortest-path -in input.gr -src 1 -dst 5
[source: 1] [destination: 5] [cost: 24]
E 5 2
E 2 10
E 10 9
E 9 1
```

## Input graphs
The implementation accepts input graphs in DIMACS format. An example graph input
is specified below.

```
section graph
nodes 10
edges 10
e 2 5 9
e 3 6 3
e 4 7 9
e 9 10 8
e 3 7 2
e 1 8 9
e 2 10 1
e 4 8 4
e 5 6 7
e 1 9 6
end
eof
```

## Performance
The implementation can find a shortest path in a 20-regular one-million vertex
edge-weighted graph (ten million edges) in less than one second using a single
core of an Intel(R) Core(TM) i7-7660U CPU @ 2.50GHz processor and 16GiB of main
memory (An Apple macbook pro laptop). For dense graphs the performance is even 
better, for example it takes less than 70 milliseconds to find a shortest path 
in a two-thousand-regular ten-thousand vertex edge-weighted graph (ten million 
edges). The graph preprocessing time (graph read and initialisation) is less 
than ten seconds.

### Sparse graphs
```
$ ../graph-gen/graph-gen regular 1000000 20 100 1234 | ./shortest-path-perf -src 1 -dst 100

invoked as: ../graph-gen/graph-gen regular 1000000 20 100 1234
invoked as: ./shortest-path-perf -src 1 -dst 100
Input file not specified, redirecting to standard input stream
gen-unique [regular]: n = 1000000, m = 10000000, seed = 1234
input: n = 1000000, m = 10000000 [8126.07ms] {peak: 0.22GiB} {curr: 0.15GiB}
root build: [zero: 0.15ms] [pos: 58.55ms] [adj: 1984.37ms] done. [2043.12ms] {peak: 0.31GiB} {curr: 0.31GiB}
dijkstra: [init: 11.61ms] [visit: 962.51ms] [950.89ms]

tracepath: [source: 1] [destination: 100] [cost: 69]
E 100 704842
E 704842 113871
E 113871 353568
E 353568 977365
E 977365 589807
E 589807 486424
E 486424 188775
E 188775 565251
E 565251 1

dijkstra-query: [query: 962.82ms 0.74GiB/s] [trace: 0.03ms] [962.86ms] {peak: 0.18GiB} {curr: 0.18GiB}
grand total [11145.96ms] {peak: 0.31GiB}
host: maagha
build: single thread, binary heap
compiler: gcc 7.3.0
```

### Dense graphs
```
$ ../graph-gen/graph-gen regular 10000 2000 100 1234 | ./shortest-path-perf -src 1 -dst 100
invoked as: ../graph-gen/graph-gen regular 10000 2000 100 1234
invoked as: ./shortest-path-perf -src 1 -dst 100
Input file not specified, redirecting to standard input stream
gen-unique [regular]: n = 10000, m = 10000000, seed = 1234
input: n = 10000, m = 10000000 [7275.85ms] {peak: 0.22GiB} {curr: 0.15GiB}
root build: [zero: 0.01ms] [pos: 12.18ms] [adj: 620.70ms] done. [632.92ms] {peak: 0.30GiB} {curr: 0.30GiB}
dijkstra: [init: 0.11ms] [visit: 67.21ms] [67.09ms]

tracepath: [source: 1] [destination: 100] [cost: 3]
E 100 4846
E 4846 3608
E 3608 1

dijkstra-query: [query: 67.22ms 3.37GiB/s] [trace: 0.03ms] [67.26ms] {peak: 0.15GiB} {curr: 0.15GiB}
grand total [7988.76ms] {peak: 0.30GiB}
host: maagha
build: single thread, binary heap
compiler: gcc 7.3.0
```

## Bug-report
For bug fixes and support, mail me at suhas.thejaswi@aalto.fi

## Citation
We encourage you to cite our work if you find this implementation useful. You
can use the following bibtex:

```
@mastersthesis{dijkstra-implementation,
title={Scalable Parameterised Algorithms for two Steiner Problems},
author={Thejaswi, Suhas},
year={2017},
pages={104 + 8},
url={http://urn.fi/URN:NBN:fi:aalto-201709046839},
}
```
