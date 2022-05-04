CL1NORM: Simplex L1 Solver Function
-----------------------------------

Based on the extremely well done and valuable work of Barrodale and Roberts:

    I. Barrodale and F. D. K. Roberts. 1980.
    Algorithm 552:
    Solution of the Constrained L1 Linear Approximation Problem [F4].
    ACM Trans. Math. Softw. 6, 2 (June 1980), 231-235.

This module, written in C (with a C++ wrapper), is compiled to a MATLAB
.mex module `cl1norm` to solve linear programming problems in the following
form:

    A*x = b

**A** is a known design matrix, **b** is a known vector (e.g. from measurements).
The vector **x** is the resulting L1 solution calculated by the `cl1norm` function:

    x = cl1norm(A, b);

![gif](doc/cl1example.gif?raw=1)

The sum of the absolute values of the residuals is minimized (L1 norm):

    min. Σ( |b - A*x| )

Optionally linear constraints and linear inequality constraints can be specified:

    C*x = d

and

    E*x <= f

The vector **x** is then calculated by the `cl1norm` function:

    x = cl1norm(A, b, C, d, E, f);

Re-build MATLAB Mex file:
-------------------------

Pre-compiled binaries for Windows and Linux are available in the build/ folder.

To re-build from source, type in the MATLAB console:

    mex -O ./src/cl1norm.cpp

![gif](doc/mexcompile.gif?raw=1)

On Linux the following script might help to compile the .mex file:

    ./mexcompile_linux.sh

Details on the cl1norm(...) MATLAB command
------------------------------------------

    Usage: [x, res, info] = cl1norm(A, B, C, D, E, F, tol, maxiter);

    Input: A, B, C (optional), D (optional), E (optional), F (optional),
    tolerance (optional), maxiter (optional).
    A*x = B, C*x = D, E*x <= F
    C and D can be empty matrices [] and/or E and F can be empty matrices [].
    tolerance: A small positive tolerance. Default: 1e-9.
    maxiter: Maximum number of iterations for the algorithm.
    Output: x, residuals (optional), simplexinfo (optional)
    simplexinfo(1): 0 - optimal solution found. >=1 no solution found.
    simplexinfo(2): Minimum sum of absolute values of the residuals.
    simplexinfo(3): Number of simplex iterations.

C++ environment
---------------

The function can be used in a C++ program with the Eigen math library
(link: https://eigen.tuxfamily.org/).

For `Eigen::Matrix<double>` the following function can be used:

    int cl1_double(const Matrix<double, Dynamic, Dynamic>& A,
                   const Matrix<double, Dynamic, 1>& B,
                   Matrix<double, Dynamic, 1>& X,
                   const Matrix<double, Dynamic, Dynamic>* C,
                   const Matrix<double, Dynamic, 1>* D,
                   const Matrix<double, Dynamic, Dynamic>* E,
                   const Matrix<double, Dynamic, 1>* F,
                   ... )

and for `Eigen::Matrix<float>`:

    int cl1_float( ... )

Note: optional input arguments (C,D,E,F) are handed over as pointers and can be set to NULL.


C environment
-------------

The basic functionality is in the `cl1` function that can be extracted from cl1norm.cpp and
embedded into a pure C project.

    void cl1(const int *k, const int *l, const int *m,
             const int *n, const int *klmd, const int *klm2d,
             const int *nklmd, const int *n2d, real *q,
             int *kode, const real *toler, int *iter, real *x,
             real *res, real *error, real *cu, int *iu, int *s)
    {
        ...
    }

Note that a `real` typedef is used here that can be setup according to project needs, e.g.:

    typedef double real;

The **A**, **b**, **C**, **d**, **E** and **f** matrices are collected in a single **q** matrix:

        A b
    q = C d
        E f

