# Finding-Shortest-path
A parallel implementation of Shortest Path problem in C, in order to Run the C code, we need to have a gcc compiler and Open MPI in our Mac System. we can go to the following link: https://wiki.helsinki.fi/display/HUGG/Open+MPI+install+on+Mac+OS+X to know, how to install Open MPI in Mac.To install gcc compiler on Mac OS X, we need to download and install “Command Line Tools for Xcode”, which is available in Apple’s developer page. After the installation of Xcode, we can compile and run the C program as follows:

1. To compile the parallel code run the command “mpicc -o m shortestpath.c” on Terminal
2. mpirun -np 1 ./m
