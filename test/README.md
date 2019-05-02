
7/4/2018

This directory contains several julia scripts for running
tests when developing STARRY/limbdark.

These still need to be made into a set of tests for the
code.

7/11/2018

Comparison with Pal's code, ntiq-fortran.f, is accomplished
as follows:

1).  Compile code from the prompt with:

tests$ gfortran -ffixed-line-length-none ntiq-fortran.f -o ntiq-fortran

2).  Run the program:

tests$ ./ntiq-fortran > ntiq_test.out

3).  Plot the results.
