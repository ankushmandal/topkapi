# Topkapi
This is a C++ implementaion of "Topkapi" algorithm from the following work:

Ankush Mandal, Cary Jiang, Anshumali Shrivastava, and Vivek Sarkar. **Topkapi: Parallel and Fast Sketches for Finding Top-K Frequent Elements**. In *Neural Information Processing Systems*(*NIPS*), Montreal, Canada, 2018.

This implementation finds Top-K frequent words from text data. It supports multi-threaded execution using OpenMP and distributed computation using MPI.

## Prerequisites
- GNU Compiler Collection (GCC)
- OpenMPI
<a/>
The code was developed and tested using GCC 6.2.0 and OpenMPI 1.10.3 on Red Hat Enterprise Linux Server release 6.5 (Santiago).

## Quick Start on Linux
1. Download the code from Github: `git clone https://github.com/ankushmandal/topkapi.git`
2. Go to `src` directory and create the executable using `make` command
3. For details on usage, type `./topkapi -h` or `./topkapi --help`

## Data
Instructions on preprocessing any text data is given in `utils` directory. The experiments on the paper were carried out using two data sets:
- Gutenberg dataset from [Project Gutenberg](https://www.gutenberg.org/). Useful instructions on downloading the data set can be found at [Nico's Blog](http://blog.ditullio.fr/2015/10/31/mini-cluster-part-iv-word-count-benchmark/).
- Puma datasets under "Wikipedia" section from [here](https://engineering.purdue.edu/~puma/datasets.htm).
