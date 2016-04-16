lorenz:
	gfortran -c lorenz.f03
	gfortran flow.f03 lorenz.o -o flow
	gfortran verify.f03 lorenz.o -o verify
	gfortran finite-diff.f03 lorenz.o -o finite-diff
	gfortran adj.f03 lorenz.o -o adj
	mv adj ../adj/

vanderpol:
	gfortran -c vanderpol.f03
	gfortran flow.f03 vanderpol.o -o flow
	gfortran verify.f03 vanderpol.o -o verify
	gfortran finite-diff.f03 vanderpol.o -o finite-diff
	gfortran adj.f03 vanderpol.o -o adj
	mv adj ../adj/
