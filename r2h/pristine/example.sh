gcc -c -Df2cFortran -I/cern/pro/include/cfortran chbook-example.c
g77 chbook-example.o `cernlib packlib,mathlib` -o chbook-example
./chbook-example
