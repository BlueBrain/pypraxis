#! /bin/bash

# f2py --fcompiler=g95 --f90flags=-fpic --f77flags=-fpic -c praxis.pyf praxis.f
f2py --fcompiler=gfortran --noopt -c praxis.pyf praxis.f

