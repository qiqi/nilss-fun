lorenz:
	gfortran -c lorenz.f03
	gfortran flow.f03 lorenz.o -o flow
	gfortran verify.f03 lorenz.o -o verify
	gfortran adj.f03 lorenz.o -o adj
	mv adj ../adj/
