#! /bin/bash

# This works on 32bit systems:
# f2py --fcompiler=g95 -c praxis.pyf praxis.f

# This works on 32bit and 64bit systems:
f2py --fcompiler=gnu95 --f77flags=-fomit-frame-pointer --noopt -c praxis.pyf praxis.f

# This works on 64bit systems:
# f2py --fcompiler=g95 --f90flags=-fpic --f77flags=-fpic -c praxis.pyf praxis.f

# This works on 32 and 64bit systems:
# f2py --fcompiler=gfortran --noopt -c praxis.pyf praxis.f

