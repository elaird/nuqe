NAME="chbook-example_ted"
gcc -c -Df2cFortran $NAME.c
gfortran $NAME.o `cernlib packlib,mathlib` -o $NAME
./$NAME
