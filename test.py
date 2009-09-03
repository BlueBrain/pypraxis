import numpy as np
from scipy.optimize import *
from praxis import praxis

def fmin_praxis(func, x0, args=(), ftol=1e-5, maxstep=1.0, disp=1):
    """Minimize a function using Brent's principal axis method.

    :Parameters:

      func : callable f(x,*args)
          Objective function to be minimized.
      x0 : ndarray
          Initial guess.
      args : tuple
          Extra arguments passed to func.

    :Returns: (xopt, fopt)

        xopt : ndarray
            Parameter which minimizes `func`.
        fopt : number
            Value of function at minimum: ``fopt = func(xopt)``.

    *Other Parameters*:

      ftol : float
          fmin_praxis attempts to return fopt such that if x0 is the
          true local minimum near xopt, then 
          norm(x-x0) < ftol + sqrt(macheps)*norm(xopt)
      maxstep : int
          maximum step size. maxstep should be set to about the
          maximum distance from the initial guess to the minimum. If
          maxstep is set too large or too small, the initial rate of
          convergence may be slow.
      disp : int
          controls the printing of intermediate results. Id disp=0,
          nothing is printed. Will be more and more verbose up to
          prin=4.

    """

    eps = np.finfo(np.double).eps

    fopt, xopt = praxis( t0=ftol, h0=maxstep, x=x0, f=func,
                         f_extra_args = args, prin=disp, machep=eps)

    return xopt, fopt

def main():
    import time

    times = []
    algor = []
    x0 = [0.8,1.2,0.7]
    print "Nelder-Mead Simplex"
    print "==================="
    start = time.time()
    x = fmin(rosen,x0)
    print x
    times.append(time.time() - start)
    algor.append('Nelder-Mead Simplex\t')

    print
    print "Powell Direction Set Method"
    print "==========================="
    start = time.time()
    x = fmin_powell(rosen,x0)
    print x
    times.append(time.time() - start)
    algor.append('Powell Direction Set Method.')

    print
    print "Brent Principal Axis Method"
    print "==========================="
    start = time.time()
    x = fmin_praxis(rosen,x0)
    print x
    times.append(time.time() - start)
    algor.append('Brent Principal Axis Method.')

    print
    print "Nonlinear CG"
    print "============"
    start = time.time()
    x = fmin_cg(rosen, x0, fprime=rosen_der, maxiter=200)
    print x
    times.append(time.time() - start)
    algor.append('Nonlinear CG     \t')

    print
    print "BFGS Quasi-Newton"
    print "================="
    start = time.time()
    x = fmin_bfgs(rosen, x0, fprime=rosen_der, maxiter=80)
    print x
    times.append(time.time() - start)
    algor.append('BFGS Quasi-Newton\t')

    print
    print "BFGS approximate gradient"
    print "========================="
    start = time.time()
    x = fmin_bfgs(rosen, x0, gtol=1e-4, maxiter=100)
    print x
    times.append(time.time() - start)
    algor.append('BFGS without gradient\t')


    print
    print "Newton-CG with Hessian product"
    print "=============================="
    start = time.time()
    x = fmin_ncg(rosen, x0, rosen_der, fhess_p=rosen_hess_prod, maxiter=80)
    print x
    times.append(time.time() - start)
    algor.append('Newton-CG with hessian product')


    print
    print "Newton-CG with full Hessian"
    print "==========================="
    start = time.time()
    x = fmin_ncg(rosen, x0, rosen_der, fhess=rosen_hess, maxiter=80)
    print x
    times.append(time.time() - start)
    algor.append('Newton-CG with full hessian')

    print
    print "\nMinimizing the Rosenbrock function of order 3\n"
    print " Algorithm \t\t\t       Seconds"
    print "===========\t\t\t      ========="
    for k in range(len(algor)):
        print algor[k], "\t -- ", times[k]

if __name__ == "__main__":
    main()
