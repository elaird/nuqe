NAME="hbook-ntuple_ted"
FLAGS="-O2 -Wall -Wextra -pedantic -Wno-long-long -ggdb"
LDFLAGS="$FLAGS"
gcc $FLAGS -c -Df2cFortran $NAME.c \
&& gfortran $LDFLAGS $NAME.o `cernlib packlib,mathlib` -o $NAME \
&& ./$NAME
