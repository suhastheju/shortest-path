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
           reads input from standard input if the argument is not specified
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
edge-weighted graph (ten million edges) in less than two seconds using a single
core of an Intel(R) Core(TM) i7-7660U CPU @ 2.50GHz processor and 16GiB of main
memory (An Apple macbook pro laptop). The graph preprocessing time (graph read
and initialisation) is less than ten seconds.

```
invoked as: ../graph-gen/graph-gen regular 1000000 20 100 1234
invoked as: ./shortest-path-perf -src 1 -dst 100
Input file not specified, redirecting to standard input stream
gen-unique [regular]: n = 1000000, m = 10000000, seed = 1234
input: n = 1000000, m = 10000000 [8305.68ms] {peak: 0.22GiB} {curr: 0.15GiB}
root build: [zero: 0.16ms] [pos: 71.77ms] [adj: 2081.41ms] done. [2153.41ms] {peak: 0.31GiB} {curr: 0.31GiB}
dijkstra: [init: 11.48ms] [visit: 1265.98ms] [1254.49ms]

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

dijkstra-query: [query: 1266.33ms 0.56GiB/s] [trace: 0.03ms] [1266.39ms] {peak: 0.18GiB} {curr: 0.18GiB}
grand total [11740.87ms] {peak: 0.31GiB}
host: maagha
build: single thread, binary heap
compiler: gcc 7.3.0
```

## Bug-fixes and support
Mail me at suhas.muniyappa@aalto.fi
