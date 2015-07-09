# pypraxis #
## Introduction ##
This code provides a Python interface to Richard P. Brent's [Principal Axis algorithm](http://wwwmaths.anu.edu.au/~brent/pub/pub011.html) for minimization without derivatives. It consists of the following parts:
  1. John Burkardt's Fortran90 implementation of the algorithm, retrieved from [here](http://people.sc.fsu.edu/~burkardt/f_src/praxis/praxis.html).
  1. Optionally, the original Fortran77 implementation from [netlib](http://www.netlib.org/opt/) can be used.
  1. An f2py signature file (`praxis.pyf`).
  1. A simple script to compile a shared library that can be imported to Python.
  1. A function wrapper and some test cases in `test.py` to illustrate usage.
Brent's Principal Axis algorithm typically requires less function calls than other gradient- and derivative-free multidimensional minimizing algorithms, such as the Nelder-Mead Simplex or Powell's Direction Set method.

## Building ##
Building requires [f2py](http://cens.ioc.ee/projects/f2py2e/) (part of [NumPy](http://numpy.scipy.org)), a Fortran compiler (I've tested gfortran and g95) and the Python development files.
To build, edit the file `compile.sh` to reflect your Fortran compiler, and then run `./compile.sh`.
Known building issues:
  1. Using g95 on a 64bit system requires passing the `-fpic` flag to the compiler.
  1. Optimization has to be switched off when using gfortran on both 32bit and 64bit systems.

## Usage ##
The test case provided in test.py is probably the best place to start. It contains a function wrapper (`fmin_praxis`) that has a similar signature as the other functions in [scipy.optimize](http://docs.scipy.org/doc/scipy/reference/tutorial/optimize.html).

## Reference ##
R. P. Brent, Algorithms for Minimization without Derivatives, Prentice-Hall, Englewood Cliffs, New Jersey, 1973, 195 pp. ISBN 0-13-022335-2. CR 15#26544; MR 49#4251, 51#7283. Reprinted in paperback by Dover Publications, Mineola, New York, January 2002. ISBN 0-486-41998-3.